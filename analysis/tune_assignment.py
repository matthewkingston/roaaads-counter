"""
Powell's-method tuning of gravity model parameters against combined chi-squared
objective across all traffic count observations.

K (global scale) and per-slot hourly fractions {f_s} are jointly calibrated at each
evaluation via alternating minimisation (not in the optimizer).  Walking observations
are compared to the model in count space (numerator: Jeffreys count n_eff = n + 0.5;
denominator: n_eff, i.e. Poisson weight).  Official AADT sites remain in AADT space.
A Gaussian prior on each f_s anchors it to the NI-average hourly fraction; the prior
mean and std are derived from hourly_fractions.csv via the law of total variance,
grouping days into weekday / Saturday / Sunday.  Replaces the former Woodbury rank-1
correction; both approaches preserve one df per time slot (N_eff = N - N_slots).

All optimizer parameters are stored in log-space to enforce positivity.

Stage 1 (--gravity, default): tune W_BIZ, P, ALPHA (3 params; 4 with THETA if k=3 cache).
Stage 2 (--full): also tune city-level pop/wp and sub-1 dampings (26 params)
  with L2 regularisation relative to simulation/tuner_config.json references.

Node 180 (pop=50, wp=0) is excluded from stage 2 — too small to tune.

Results:
  simulation/tuned_params.json   best params from this run (read by build_assignment.py)
  simulation/tuning_history.jsonl  appended record of every run
  simulation/gravity_model_curve.png  kernel shape plot

Usage:
  python3 analysis/tune_assignment.py
  python3 analysis/tune_assignment.py --full
  python3 analysis/tune_assignment.py --note "added-june-counts"
"""

import csv, json, math, os, secrets, subprocess, sys, time, xml.etree.ElementTree as ET
from datetime import datetime, timezone

import numpy as np
import scipy.optimize

sys.path.insert(0, "simulation")
from model import (COUNT_SITES, EXCLUDE_LINKS, PATHS_CACHE, WEIGHTS_FILE,
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS,
                   print_chi2_table)

# ── Paths ─────────────────────────────────────────────────────────────────────

CONS_GRAPH        = "simulation/newtownards_consolidated.graphml"
HISTORY_FILE      = "simulation/tuning_history.jsonl"
CURVE_PNG         = "simulation/gravity_model_curve.png"
HOURLY_FRACS_FILE = "analysis/hourly_fractions.csv"

# ── CLI args ──────────────────────────────────────────────────────────────────

stage = "gravity"
note  = None
argv  = sys.argv[1:]
i = 0
while i < len(argv):
    if argv[i] == "--full":
        stage = "full"
    elif argv[i] == "--note" and i + 1 < len(argv):
        i += 1
        note = argv[i]
    i += 1

run_id = secrets.token_hex(4)

try:
    git_hash = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
    ).decode().strip()
except Exception:
    git_hash = "unknown"

print(f"Run ID: {run_id}  stage: {stage}  git: {git_hash}" +
      (f"  note: {note}" if note else ""))

# ── Load previous best from history ──────────────────────────────────────────

prev_chi2_per_n = None
prev_id         = None
if os.path.exists(HISTORY_FILE):
    with open(HISTORY_FILE) as _hf:
        _lines = [ln for ln in _hf if ln.strip()]
    if _lines:
        _last = json.loads(_lines[-1])
        prev_chi2_per_n = _last.get("chi2_per_n")
        prev_id         = _last.get("id", "")
        print(f"Previous run:  χ²/N={prev_chi2_per_n:.4f}  (id={prev_id})")

# ── Load paths cache ──────────────────────────────────────────────────────────

if not os.path.exists(PATHS_CACHE):
    print(f"ERROR: {PATHS_CACHE} not found — run build_paths.py first")
    raise SystemExit(1)

print("Loading paths cache …")
cache        = np.load(PATHS_CACHE)
node_ids_arr = cache["node_ids"]
od_src       = cache["od_src"]
od_dst       = cache["od_dst"]
od_dist      = cache["od_dist"].astype(np.float64)
pair_idx     = cache["pair_idx"]
link_idx_arr = cache["link_idx"]
link_u       = cache["link_u"]
link_v       = cache["link_v"]
N_links      = len(link_u)
node_ids     = [int(nid) for nid in node_ids_arr]
N_nodes      = len(node_ids)

link_index = {(int(link_u[k]), int(link_v[k])): k for k in range(N_links)}

_has_stoch = "pair_idx_2" in cache
if _has_stoch:
    _od_dist_2  = cache["od_dist_2"].astype(np.float64)
    _pair_idx_2 = cache["pair_idx_2"]
    _link_idx_2 = cache["link_idx_2"]
    _od_dist_3  = cache["od_dist_3"].astype(np.float64)
    _pair_idx_3 = cache["pair_idx_3"]
    _link_idx_3 = cache["link_idx_3"]
    # Per-OD-pair float32 distances and precomputed log(d) for the gravity kernel.
    # u^(ALPHA+1) = exp((ALPHA+1)*(log_d - log_P)); log_d fixed, log_P cheap scalar.
    _od_dist_f32   = od_dist.astype(np.float32)
    _od_dist_2_f32 = _od_dist_2.astype(np.float32)
    _od_dist_3_f32 = _od_dist_3.astype(np.float32)
    _log_od_dist   = np.log(od_dist).astype(np.float32)
    _log_od_dist_2 = np.log(_od_dist_2).astype(np.float32)
    _log_od_dist_3 = np.log(_od_dist_3).astype(np.float32)
    # Precompute CSR sparse link-pair matrices (one per path).
    # SpMV `A_r @ w_r` (N_links × N_OD) × (N_OD,) replaces the
    # gather+bincount scatter, cutting per-eval scatter cost ~9×.
    from scipy.sparse import csr_matrix as _csr_build
    print("  Building sparse link-pair matrices …", end=" ", flush=True)
    _t_csr = time.time()
    _N_OD = len(od_src)
    _A1 = _csr_build((np.ones(len(pair_idx),    np.float32), (link_idx_arr, pair_idx)),    shape=(N_links, _N_OD))
    _A2 = _csr_build((np.ones(len(_pair_idx_2), np.float32), (_link_idx_2,  _pair_idx_2)), shape=(N_links, _N_OD))
    _A3 = _csr_build((np.ones(len(_pair_idx_3), np.float32), (_link_idx_3,  _pair_idx_3)), shape=(N_links, _N_OD))
    print(f"done ({time.time()-_t_csr:.1f}s)")
else:
    _od_dist_2 = _pair_idx_2 = _link_idx_2 = None
    _od_dist_3 = _pair_idx_3 = _link_idx_3 = None

print(f"  {N_nodes} nodes  {N_links} links  {len(od_src):,} OD pairs"
      + ("  stochastic k=3 paths loaded" if _has_stoch else "  (no stochastic paths — run build_paths.py)"))

# ── Load node weights ─────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    wdata = json.load(f)

node_pop_full = {int(k): v for k, v in wdata["node_population"].items()}
node_biz_full = {int(k): v for k, v in wdata["node_business_demand"].items()}

# Precomputed base weight arrays (used in stage 1 and as fallback in stage 2)
base_w_pop = np.array([node_pop_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
base_w_biz = np.array([node_biz_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)

# ── Load street names from consolidated GraphML (sequential node IDs) ─────────

link_name = {}
if os.path.exists(CONS_GRAPH):
    _ns  = "http://graphml.graphstruct.org/graphml"
    _ns2 = "http://graphml.graphstruct.org/graphml"
    try:
        _tree = ET.parse(CONS_GRAPH)
        _root = _tree.getroot()
        _nsmap = {"g": _root.tag.split("}")[0].lstrip("{")} if "}" in _root.tag else {"g": ""}
        _pfx = "{" + _nsmap["g"] + "}" if _nsmap["g"] else ""
        for _edge in _root.iter(f"{_pfx}edge"):
            _u = int(_edge.get("source"))
            _v = int(_edge.get("target"))
            for _data in _edge:
                if _data.get("key") == "d14" and _data.text:
                    link_name[(_u, _v)] = _data.text
                    break
    except Exception as _e:
        print(f"Warning: could not load street names ({_e})")


def _link_label(u, v):
    name = link_name.get((u, v), "")
    return f"{u}→{v}  {name}" if name else f"{u}→{v}"

# ── Load tuner config ─────────────────────────────────────────────────────────

with open(TUNER_CONFIG) as f:
    config = json.load(f)

lam       = config["lambda"]
city_list = list(config["cities"].items())

grav_ref = config.get("gravity_ref", {})
grav_lam = config.get("gravity_lambda", 0.0)
_grav_ref_vals = [
    math.log(max(grav_ref.get("W_BIZ", 1.0),   1e-4)),
    math.log(max(grav_ref.get("P",     300.0),  1e-4)),
    math.log(max(grav_ref.get("ALPHA", 2.0),    1e-4)),
]
if _has_stoch:
    _grav_ref_vals.append(math.log(max(grav_ref.get("THETA", 1.0), 1e-4)))
log_grav_ref = np.array(_grav_ref_vals)

# External nodes covered by tuner_config (node 180 excluded)
external_nodes = set()
for _, city_cfg in city_list:
    external_nodes.update(city_cfg["nodes"])

# Node array indices for external nodes (for fast overriding in stage 2)
ext_indices = [(i, nid) for i, nid in enumerate(node_ids) if nid in external_nodes]

# Tunable dampings: those with value < 1.0, in deterministic order
tunable_dampings = []  # [(city_name, node_id, ref_damping), ...]
for city_name, city_cfg in city_list:
    for node_str in sorted(city_cfg["dampings"], key=int):
        damp = city_cfg["dampings"][node_str]
        if damp < 1.0:
            tunable_dampings.append((city_name, int(node_str), damp))

n_gravity = 4 if _has_stoch else 3
n_city    = len(city_list) * 2
n_damp    = len(tunable_dampings)
n_ext     = n_city + n_damp  # external params in stage 2

print(f"  {len(city_list)} cities  {len(tunable_dampings)} tunable dampings")

# Override external node weights from tuner_config so gravity-only tuning
# always uses the same city values as the --full stage, regardless of what
# node_weights.json contains.
_ext_pop_cfg = {}
_ext_biz_cfg = {}
for _city_name, _city_cfg in city_list:
    for _nid in _city_cfg["nodes"]:
        _damp = _city_cfg["dampings"][str(_nid)]
        _ext_pop_cfg[_nid] = _city_cfg["ref_pop"] * _damp
        _ext_biz_cfg[_nid] = _city_cfg["ref_wp"]  * _damp
for _arr_i, _nid in ext_indices:
    base_w_pop[_arr_i] = _ext_pop_cfg[_nid]
    base_w_biz[_arr_i] = _ext_biz_cfg[_nid]

# ── Precompute distance-bin link matrices (all-or-nothing mode only) ──────────
# The 20M-entry scatter (pair_idx gather + bincount) dominates each evaluation.
# Precompute link-bin accumulation matrices so each eval is a small matmul instead.
#
# flow[l] = Σ_k M[l,k] · w_prod[k] · f(d_k; P,ALPHA)
#         ≈ Σ_b f(d_b; P,ALPHA) · Σ_{k in bin b} M[l,k] · w_prod[k]
#         = (link_bin_pp + W·link_bin_pb + W²·link_bin_bb) @ f_b
#
# Stage 1: all OD-pair weights fixed → use all-pair matrices.
# Stage 2: external-node weights vary → precompute internal-only matrices and
#   compute the ~23K external-involved OD pairs exactly per eval.
#
# Skipped when _has_stoch=True: the stochastic path uses exact per-pair scatter
# for every eval and never calls the bin-matrix path.

if not _has_stoch:
    print("Precomputing link-bin matrices …")
    _t_pc = time.time()

    from scipy.sparse import coo_matrix as _coo

    N_BINS      = 300
    _d_lo       = float(od_dist.min())
    _d_hi       = float(od_dist.max()) * 1.001
    _bin_edges  = np.logspace(np.log10(_d_lo), np.log10(_d_hi), N_BINS + 1)
    _bin_centers = np.sqrt(_bin_edges[:-1] * _bin_edges[1:])  # geometric mean per bin

    # Bin index for each OD pair (0 … N_BINS-1)
    _pair_bin  = np.clip(np.digitize(od_dist, _bin_edges) - 1, 0, N_BINS - 1).astype(np.int32)
    _col_all   = _pair_bin[pair_idx]   # bin index for every link-pair entry

    # Weight products per OD pair under base weights (ext nodes already at ref values)
    _pp_od = base_w_pop[od_src] * base_w_pop[od_dst]
    _pb_od = (base_w_pop[od_src] * base_w_biz[od_dst]
              + base_w_biz[od_src] * base_w_pop[od_dst])
    _bb_od = base_w_biz[od_src] * base_w_biz[od_dst]

    # Mask: True for OD pairs that involve at least one external node
    _is_ext_node = np.zeros(N_nodes, dtype=bool)
    for _ai, _ in ext_indices:
        _is_ext_node[_ai] = True
    _ext_pair_mask  = _is_ext_node[od_src] | _is_ext_node[od_dst]  # (N_OD,)
    _entry_is_ext   = _ext_pair_mask[pair_idx]                       # (N_entries,)

    _int_sel = np.where(~_entry_is_ext)[0]   # link-pair entries for internal-only pairs
    _ext_sel = np.where(_entry_is_ext)[0]    # link-pair entries for external-involved pairs


    def _link_bin(dat_od, sel=None):
        """Build a dense (N_links, N_BINS) accumulation matrix via COO→dense."""
        if sel is None:
            r, c, d = link_idx_arr, _col_all, dat_od[pair_idx]
        else:
            pi = pair_idx[sel]
            r, c, d = link_idx_arr[sel], _col_all[sel], dat_od[pi]
        return _coo((d, (r, c)), shape=(N_links, N_BINS)).toarray()


    # Stage-1 matrices: all OD pairs (external nodes at reference weights, fixed)
    # Stored as float32: f32 × f32 BLAS SGEMV is ~30× faster than f64 DGEMV.
    all_bin_pp = _link_bin(_pp_od).astype(np.float32)
    all_bin_pb = _link_bin(_pb_od).astype(np.float32)
    all_bin_bb = _link_bin(_bb_od).astype(np.float32)

    # Stage-2 matrices: internal-internal pairs only
    int_bin_pp = _link_bin(_pp_od, _int_sel).astype(np.float32)
    int_bin_pb = _link_bin(_pb_od, _int_sel).astype(np.float32)
    int_bin_bb = _link_bin(_bb_od, _int_sel).astype(np.float32)

    # External pair arrays for exact per-eval scatter in stage 2
    _ext_p     = np.unique(pair_idx[_ext_sel])                # unique global pair indices
    _ext_dist  = od_dist[_ext_p]                              # distances for ext OD pairs
    _ext_src   = od_src[_ext_p]                               # node-array src index
    _ext_dst   = od_dst[_ext_p]                               # node-array dst index
    _ext_local = np.empty(len(od_src), dtype=np.int32)        # global→local pair map
    _ext_local[_ext_p] = np.arange(len(_ext_p), dtype=np.int32)
    _ext_link  = link_idx_arr[_ext_sel]                       # link idx per ext entry
    _ext_lp    = _ext_local[pair_idx[_ext_sel]]               # local pair idx per ext entry

    del _is_ext_node, _ext_pair_mask, _entry_is_ext, _int_sel, _ext_sel, _ext_local
    del _pp_od, _pb_od, _bb_od, _col_all

    print(f"  {N_BINS} bins  {len(_ext_p):,} ext pairs  {len(_ext_link):,} ext entries"
          f"  ({time.time()-_t_pc:.1f}s)")

    def _kern_b(P, ALPHA):
        u = _bin_centers / P
        return ((ALPHA + 1) * u / (ALPHA + u ** (ALPHA + 1))).astype(np.float32)

# ── Per-slot hourly fraction priors ───────────────────────────────────────────
# Group days into weekday (0), Saturday (1), Sunday (2).
# Prior mean and std for each (day_type, hour) slot are derived from
# hourly_fractions.csv via the law of total variance:
#   total_var = between_day_var(means) + mean(within_day_var)
# This equals the pooled variance computed directly from the raw NI count data.

_DOW_TO_TYPE = {d: (0 if d < 5 else (1 if d == 5 else 2)) for d in range(7)}
_DT_DOWS     = {0: list(range(5)), 1: [5], 2: [6]}

_raw_fracs = {}  # {(dow, hour): (mean_f, std_f)}
with open(HOURLY_FRACS_FILE, newline="") as _fh:
    for _row in csv.DictReader(_fh):
        _dow  = int(_row["day_of_week"])
        _hour = int(_row["hour"].split(":")[0])
        _raw_fracs[(_dow, _hour)] = (float(_row["mean_fraction"]), float(_row["std_fraction"]))

slot_prior = {}  # {(day_type, hour): (mean_f, std_f)}
for _dt, _dows in _DT_DOWS.items():
    for _h in range(24):
        _means = [_raw_fracs[(_d, _h)][0] for _d in _dows if (_d, _h) in _raw_fracs]
        _stds  = [_raw_fracs[(_d, _h)][1] for _d in _dows if (_d, _h) in _raw_fracs]
        if not _means:
            continue
        _mf      = sum(_means) / len(_means)
        _between = sum((_m - _mf) ** 2 for _m in _means) / len(_means)  # population var
        _within  = sum(_s ** 2 for _s in _stds) / len(_stds)
        slot_prior[(_dt, _h)] = (_mf, math.sqrt(_between + _within))

# ── Build observation list ────────────────────────────────────────────────────

observations  = []
obs_slot_keys = []   # (day_type, hour) per obs, or None for official/unslotted
obs_n_eff_lst = []   # Jeffreys count (n + 0.5) per obs, None for official
obs_dur_lst   = []   # duration_s per obs, None for official

for s in COUNT_SITES:
    observations.append(("official", s["node"], s["links"],
                         float(s["observed"]), 0.10 * s["observed"]))
    obs_slot_keys.append(None)
    obs_n_eff_lst.append(None)
    obs_dur_lst.append(None)

if os.path.exists(LINK_AADT):
    with open(LINK_AADT) as f:
        link_aadt_data = json.load(f)["links"]
    for key, entry in sorted(link_aadt_data.items()):
        u, v = map(int, key.split(","))
        if (u, v) in EXCLUDE_LINKS:
            continue
        for sess_obs in entry["observations"]:
            observations.append(("walking", (u, v), None,
                float(sess_obs["aadt"]), float(sess_obs["aadt_uncertainty"])))
            ts   = sess_obs.get("time_slot")
            neff = sess_obs.get("n_eff")
            dur  = sess_obs.get("duration_s")
            if ts is not None:
                sk = (_DOW_TO_TYPE[ts[0]], ts[1])
                obs_slot_keys.append(sk if sk in slot_prior else None)
            else:
                obs_slot_keys.append(None)
            obs_n_eff_lst.append(float(neff) if neff is not None else None)
            obs_dur_lst.append(float(dur)  if dur  is not None else None)

n_obs    = len(observations)
obs_arr  = np.array([o[3] for o in observations], dtype=np.float64)
sig_arr  = np.array([o[4] for o in observations], dtype=np.float64)

# Count-space arrays for walking obs (T in fractional hours)
obs_n_eff = np.zeros(n_obs)
obs_Th    = np.zeros(n_obs)   # T / 3600
for i in range(n_obs):
    if obs_n_eff_lst[i] is not None:
        obs_n_eff[i] = obs_n_eff_lst[i]
    if obs_dur_lst[i] is not None:
        obs_Th[i] = obs_dur_lst[i] / 3600.0

# Group slotted observations by (day_type, hour)
slot_groups = {}
for i, sk in enumerate(obs_slot_keys):
    if sk is not None:
        slot_groups.setdefault(sk, []).append(i)
slot_list = list(slot_groups.items())

unslotted_idxs = [i for i, sk in enumerate(obs_slot_keys) if sk is None]

n_slots   = len(slot_list)
n_slotted = sum(len(idxs) for _, idxs in slot_list)
n_eff     = n_obs - n_slots
n_official = len(COUNT_SITES)
print(f"  {n_obs} observations ({n_official} official, {n_obs - n_official} walking"
      f" in {n_slots} time slot(s))")
if n_slotted < n_obs - n_official:
    print(f"  Warning: {n_obs - n_official - n_slotted} walking obs have no time-slot"
          f" data — treated as independent (AADT space, original uncertainty)")
for sk, idxs in sorted(slot_list):
    if len(idxs) > 2:
        print(f"  Slot {sk}: {len(idxs)} observations")

# Precompute link index sets per observation for fast model flow extraction
obs_link_idxs = []
for kind, target, links, *_ in observations:
    if kind == "official":
        if links:
            idxs = [link_index[lnk] for lnk in links if lnk in link_index]
        else:
            idxs = [k for k in range(N_links)
                    if int(link_u[k]) == target or int(link_v[k]) == target]
    else:
        k = link_index.get(target, -1)
        idxs = [k] if k >= 0 else []
    obs_link_idxs.append(idxs)

# ── Precomputed arrays for vectorised model_obs / calibrate_K_and_fracs / chi2 ─

# Walking obs each have exactly 1 link; store as a numpy index array
_walk_link_raw = np.array(
    [obs_link_idxs[i][0] if obs_link_idxs[i] else -1
     for i in range(n_official, n_obs)],
    dtype=np.int64)
_walk_link_safe = np.where(_walk_link_raw >= 0, _walk_link_raw, 0)
_walk_valid     = _walk_link_raw >= 0

# Unslotted obs (official sites) as numpy arrays
_uns_arr    = np.array(unslotted_idxs, dtype=np.int64)
_uns_sig_sq = sig_arr[_uns_arr] ** 2
_uns_obs    = obs_arr[_uns_arr]

# Per-slot count-space arrays for calibrate_K_and_fracs and chi²
# Each entry: (slot_key, idxs, n_effs_arr, Ths_arr, mean_f, inv_var_f)
# inv_var_f = 1/std_f² is the Gaussian prior precision on f_s
_slot_data = [
    (sk,
     np.array(idxs, dtype=np.int64),
     obs_n_eff[np.array(idxs)],
     obs_Th[np.array(idxs)],
     slot_prior[sk][0],
     1.0 / slot_prior[sk][1] ** 2)
    for sk, idxs in slot_list
]

# ── Assignment and chi-squared helpers ───────────────────────────────────────

def run_assignment(W_BIZ, P, ALPHA, w_pop, w_biz, THETA=None):
    """Gravity assignment with optional stochastic logit routing.

    THETA=None  → fast binned all-or-nothing:
      Stage 1: ~0.12 ms (3 f32 matmuls).  Stage 2: ~5.5 ms (+ exact ext scatter).
    THETA given → CSR SpMV scatter over k=3 paths with logit weights (~150 ms).
      share(r) ∝ exp(−THETA · d_r / P).  Works identically for stage 1 and 2.
    """
    if THETA is None or not _has_stoch:
        # Fast bin-matrix path (all-or-nothing)
        f_b = _kern_b(P, ALPHA)
        W   = np.float32(W_BIZ)
        W2  = np.float32(W_BIZ ** 2)
        if stage == "full":
            flow = (int_bin_pp @ f_b
                    + W  * (int_bin_pb @ f_b)
                    + W2 * (int_bin_bb @ f_b)).astype(np.float64)
            w_src = w_pop[_ext_src] + W_BIZ * w_biz[_ext_src]
            w_dst = w_pop[_ext_dst] + W_BIZ * w_biz[_ext_dst]
            u_e   = _ext_dist / P
            t_e   = w_src * w_dst * (ALPHA + 1) * u_e / (ALPHA + u_e ** (ALPHA + 1))
            flow += np.bincount(_ext_link, weights=t_e[_ext_lp], minlength=N_links)
        else:
            flow = (all_bin_pp @ f_b
                    + W  * (all_bin_pb @ f_b)
                    + W2 * (all_bin_bb @ f_b)).astype(np.float64)
        return flow

    # Stochastic logit: exact scatter over 3 paths.
    # All per-OD-pair (820K) arrays use float32.
    # Scatter via prebuilt CSR SpMV (~38ms × 3) rather than gather+bincount
    # (~102+350ms × 3), for ~9× scatter speedup.
    w_vec = (w_pop + W_BIZ * w_biz).astype(np.float32)
    t_ij  = w_vec[od_src] * w_vec[od_dst]  # float32, N_OD

    # Logit shares (float32, stable via row-wise max subtraction)
    d_mat  = np.stack([_od_dist_f32, _od_dist_2_f32, _od_dist_3_f32], axis=1)
    log_w  = np.float32(-THETA / P) * d_mat
    log_w -= log_w.max(axis=1, keepdims=True)
    shares = np.exp(log_w)
    shares /= shares.sum(axis=1, keepdims=True)  # float32, (N_OD, 3)

    # Kernel via precomputed log(d): avoids log() call per eval
    log_P   = np.float32(math.log(P))
    alpha1  = np.float32(ALPHA + 1)
    alpha_f = np.float32(ALPHA)

    flow = np.zeros(N_links, dtype=np.float64)
    for r, (A_r, d_f32, log_d) in enumerate([
        (_A1, _od_dist_f32,   _log_od_dist),
        (_A2, _od_dist_2_f32, _log_od_dist_2),
        (_A3, _od_dist_3_f32, _log_od_dist_3),
    ]):
        log_u = log_d - log_P
        u_pow = np.exp(alpha1 * log_u)
        u_r   = d_f32 * np.float32(1.0 / P)
        f_r   = alpha1 * u_r / (alpha_f + u_pow)   # float32, N_OD
        w_r   = t_ij * shares[:, r] * f_r           # float32, N_OD
        flow += A_r @ w_r                            # SpMV: (N_links, N_OD) × (N_OD,)
    return flow


def model_obs(raw_flow):
    m = np.empty(n_obs)
    for i, idxs in enumerate(obs_link_idxs[:n_official]):
        m[i] = raw_flow[idxs].sum() if idxs else 0.0
    m[n_official:] = np.where(_walk_valid, raw_flow[_walk_link_safe], 0.0)
    return m


def calibrate_K_and_fracs(m_arr):
    """Jointly calibrate K and per-slot hourly fractions via alternating minimisation.

    Walking obs are in count space: chi²_walk = Σ (K·m·T·f_s/3600 − n_eff)²/n_eff.
    Official obs remain in AADT space: chi²_off = Σ (K·m − y)²/σ².
    K and each f_s appear as K·f_s products, so coordinate descent alternates:
      K-step  (quadratic in K for fixed {f_s})
      f_s-step (quadratic in f_s for fixed K, one per slot)
    Converges in 3–5 iterations; 10 iterations is ample.
    Returns (K, slot_fracs) where slot_fracs = {(day_type, hour): f_s}.
    """
    slot_fracs = {sk: mean_f for sk, _, _, _, mean_f, _ in _slot_data}
    K = 1.0
    for _ in range(10):
        # K-step: K = B / A
        # A = Σ_off m²/σ²  + Σ_walk (m·T/3600·f_s)²/n_eff
        # B = Σ_off m·y/σ² + Σ_walk m·T/3600·f_s        (n_eff/n_eff = 1)
        A = float(np.sum(m_arr[_uns_arr] ** 2 / _uns_sig_sq))
        B = float(np.sum(m_arr[_uns_arr] * _uns_obs / _uns_sig_sq))
        for sk, ia, n_effs, Ths, mean_f, inv_var_f in _slot_data:
            B_i = m_arr[ia] * Ths * slot_fracs[sk]  # m · (T/3600) · f_s
            A  += float(np.sum(B_i ** 2 / n_effs))
            B  += float(np.sum(B_i))
        K = B / A if A > 0 else 1.0

        # f_s-step: f_s = (Σ C_i + mean_f/std_f²) / (Σ C_i²/n_eff_i + 1/std_f²)
        # where C_i = K · m_i · T_i/3600
        for sk, ia, n_effs, Ths, mean_f, inv_var_f in _slot_data:
            C_i = K * m_arr[ia] * Ths
            num = float(np.sum(C_i)) + mean_f * inv_var_f
            den = float(np.sum(C_i ** 2 / n_effs)) + inv_var_f
            slot_fracs[sk] = num / den if den > 0 else mean_f

    return K, slot_fracs


# ── Objective function ────────────────────────────────────────────────────────

eval_count = [0]
best       = {"chi2": float("inf"), "log_params": None, "K": 1.0}
t0         = time.time()


def objective(log_params, log_ref=None):
    log_params = np.clip(log_params, -100, 100)
    W_BIZ = math.exp(log_params[0])
    P     = math.exp(log_params[1])
    ALPHA = math.exp(log_params[2])
    THETA = math.exp(log_params[3]) if _has_stoch else None

    if stage == "full":
        # Build per-node weight arrays from city-level params and dampings
        w_pop = base_w_pop.copy()
        w_biz = base_w_biz.copy()

        idx = n_gravity
        city_pops = {}
        city_wps  = {}
        for city_name, _ in city_list:
            city_pops[city_name] = math.exp(log_params[idx])
            city_wps[city_name]  = math.exp(log_params[idx + 1])
            idx += 2

        # Start from reference dampings, then override tunable ones
        curr_dampings = {}
        for city_name, city_cfg in city_list:
            for node_str, damp in city_cfg["dampings"].items():
                curr_dampings[int(node_str)] = damp
        for i_td, (_, node_id, _) in enumerate(tunable_dampings):
            curr_dampings[node_id] = math.exp(log_params[n_gravity + n_city + i_td])

        # Apply to node arrays
        for city_name, city_cfg in city_list:
            for node_id in city_cfg["nodes"]:
                damp = curr_dampings[node_id]
                # Find this node's index in node_ids
                pass  # done below via ext_indices

        # Use ext_indices for efficient override
        ext_pop_map = {}
        ext_biz_map = {}
        for city_name, city_cfg in city_list:
            for node_id in city_cfg["nodes"]:
                damp = curr_dampings[node_id]
                ext_pop_map[node_id] = city_pops[city_name] * damp
                ext_biz_map[node_id] = city_wps[city_name]  * damp
        for arr_i, nid in ext_indices:
            w_pop[arr_i] = ext_pop_map.get(nid, base_w_pop[arr_i])
            w_biz[arr_i] = ext_biz_map.get(nid, base_w_biz[arr_i])
    else:
        w_pop = base_w_pop
        w_biz = base_w_biz

    raw_flow      = run_assignment(W_BIZ, P, ALPHA, w_pop, w_biz, THETA)
    m_arr         = model_obs(raw_flow)
    K, slot_fracs = calibrate_K_and_fracs(m_arr)
    # Official obs: AADT space
    r_u  = K * m_arr[_uns_arr] - _uns_obs
    chi2 = float(np.sum(r_u ** 2 / _uns_sig_sq))
    # Walking obs: count space + Gaussian prior on each slot fraction
    for sk, ia, n_effs, Ths, mean_f, inv_var_f in _slot_data:
        pred  = K * m_arr[ia] * Ths * slot_fracs[sk]  # expected counts
        chi2 += float(np.sum((pred - n_effs) ** 2 / n_effs))
        chi2 += (slot_fracs[sk] - mean_f) ** 2 * inv_var_f

    if stage == "full" and log_ref is not None:
        chi2 += lam * float(np.sum((log_params[n_gravity:] - log_ref[n_gravity:]) ** 2))
    if grav_lam > 0:
        chi2 += grav_lam * float(np.sum((log_params[:n_gravity] - log_grav_ref) ** 2))

    eval_count[0] += 1
    if chi2 < best["chi2"]:
        best["chi2"]       = chi2
        best["log_params"] = log_params.copy()
        best["K"]          = K

    if eval_count[0] % 20 == 0:
        elapsed = time.time() - t0
        print(f"  {eval_count[0]:4d}  χ²={chi2:.2f}  χ²/N={chi2/n_obs:.3f}"
              f"  best={best['chi2']/n_obs:.3f}  ({elapsed:.0f}s)")

    return chi2


# ── Build initial parameter vector ────────────────────────────────────────────

# Gravity start: from tuned_params.json if available, else hardcoded defaults
grav_start = {"W_BIZ": 1.0, "P": 300.0, "ALPHA": 2.0, "THETA": 1.0}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        tp = json.load(f)
    for k in ("W_BIZ", "P", "ALPHA", "THETA"):
        if k in tp:
            grav_start[k] = tp[k]
    print(f"Starting gravity params from {TUNED_PARAMS}")

# Clamp to a safe minimum before log-transform (guards against degenerate prior runs)
_LOG_MIN = 1e-4
_log_p0_vals = [
    math.log(max(grav_start["W_BIZ"], _LOG_MIN)),
    math.log(max(grav_start["P"],     _LOG_MIN)),
    math.log(max(grav_start["ALPHA"], _LOG_MIN)),
]
if _has_stoch:
    _log_p0_vals.append(math.log(max(grav_start["THETA"], _LOG_MIN)))
log_p0 = np.array(_log_p0_vals, dtype=np.float64)

log_ref = None

if stage == "full":
    # City pop/wp: reference values from config
    for city_name, city_cfg in city_list:
        log_p0 = np.append(log_p0, [
            math.log(city_cfg["ref_pop"]),
            math.log(max(city_cfg["ref_wp"], 1)),  # guard against 0 (not present in config)
        ])
    # Tunable dampings: reference values from config
    for _, _, ref_damp in tunable_dampings:
        log_p0 = np.append(log_p0, math.log(ref_damp))

    log_ref = log_p0.copy()
    print(f"Full stage: {len(log_p0)} params "
          f"({n_gravity} gravity + {n_city} city pop/wp + {n_damp} dampings)")
else:
    _theta_note = "  [W_BIZ, P, ALPHA, THETA]" if _has_stoch else "  [W_BIZ, P, ALPHA]"
    print(f"Gravity stage: {len(log_p0)} params{_theta_note}")

# ── Run optimization ──────────────────────────────────────────────────────────

print(f"\nRunning Powell's method (λ={lam}) …")
print(f"  {'eval':>4s}  χ²/N(curr)  χ²/N(best)  elapsed")

# Evaluate initial point
objective(log_p0, log_ref)

result = scipy.optimize.minimize(
    lambda p: objective(p, log_ref),
    log_p0,
    method="powell",
    options={"maxiter": 5000, "ftol": 1e-5 if stage == "gravity" else 1e-4,
             "xtol":  1e-5 if stage == "gravity" else 1e-4},
)

# Use best params seen (Powell may backtrack at convergence)
log_best = best["log_params"]

# ── Unpack best params ────────────────────────────────────────────────────────

W_BIZ = math.exp(log_best[0])
P     = math.exp(log_best[1])
ALPHA = math.exp(log_best[2])
THETA = math.exp(log_best[3]) if _has_stoch else None

ext_pop_map  = {}
ext_biz_map  = {}
city_pops_out = {}
city_wps_out  = {}
dampings_out  = {}

if stage == "full":
    idx = n_gravity
    for city_name, city_cfg in city_list:
        city_pops_out[city_name] = math.exp(log_best[idx])
        city_wps_out[city_name]  = math.exp(log_best[idx + 1])
        idx += 2

    curr_dampings = {}
    for city_name, city_cfg in city_list:
        for node_str, damp in city_cfg["dampings"].items():
            curr_dampings[int(node_str)] = damp
    for i_td, (_, node_id, _) in enumerate(tunable_dampings):
        d = math.exp(log_best[n_gravity + n_city + i_td])
        curr_dampings[node_id] = d
        dampings_out[str(node_id)] = d

    for city_name, city_cfg in city_list:
        for node_id in city_cfg["nodes"]:
            damp = curr_dampings[node_id]
            ext_pop_map[node_id] = city_pops_out[city_name] * damp
            ext_biz_map[node_id] = city_wps_out[city_name]  * damp

# Final evaluation for clean chi2, K, and slot_fracs (no L2 term)
w_pop_f = base_w_pop.copy() if ext_pop_map else base_w_pop
w_biz_f = base_w_biz.copy() if ext_biz_map else base_w_biz
for arr_i, nid in ext_indices:
    if nid in ext_pop_map:
        w_pop_f[arr_i] = ext_pop_map[nid]
        w_biz_f[arr_i] = ext_biz_map[nid]

raw_flow      = run_assignment(W_BIZ, P, ALPHA, w_pop_f, w_biz_f, THETA)
m_arr         = model_obs(raw_flow)
K, slot_fracs = calibrate_K_and_fracs(m_arr)
# Official obs: AADT space
r_u   = K * m_arr[_uns_arr] - _uns_obs
chi2  = float(np.sum(r_u ** 2 / _uns_sig_sq))
# Walking obs: count space + Gaussian prior on each slot fraction
for sk, ia, n_effs, Ths, mean_f, inv_var_f in _slot_data:
    pred  = K * m_arr[ia] * Ths * slot_fracs[sk]
    chi2 += float(np.sum((pred - n_effs) ** 2 / n_effs))
    chi2 += (slot_fracs[sk] - mean_f) ** 2 * inv_var_f
chi2_per_n = chi2 / n_obs

# Effective AADT and sigma for each obs (using inferred slot fractions for walking obs).
# z-score is identical in count and AADT space (numerator and denominator scale the same way).
obs_eff = obs_arr.copy()
sig_eff = sig_arr.copy()
for i in range(n_official, n_obs):
    sk = obs_slot_keys[i]
    if sk is not None and sk in slot_fracs and obs_Th[i] > 0 and obs_n_eff[i] > 0:
        f_s = slot_fracs[sk]
        th  = obs_Th[i]
        obs_eff[i] = obs_n_eff[i] / (th * f_s)              # effective AADT
        sig_eff[i] = math.sqrt(obs_n_eff[i]) / (th * f_s)   # Poisson sigma in AADT space
resid = (K * m_arr - obs_eff) / sig_eff

elapsed = time.time() - t0
print(f"\nResult  ({eval_count[0]} evals, {elapsed:.0f}s)")
_theta_str = f"  THETA={THETA:.4f}" if THETA is not None else ""
print(f"  K={K:.4e}  W_BIZ={W_BIZ:.4f}  P={P:.2f}s  ALPHA={ALPHA:.4f}{_theta_str}")
print(f"  χ²={chi2:.2f}  χ²/N={chi2_per_n:.4f}  χ²/N_eff={chi2/n_eff:.3f}  (N={n_obs}, N_eff={n_eff})")
if prev_chi2_per_n is not None:
    delta = chi2_per_n - prev_chi2_per_n
    direction = "improvement" if delta < 0 else "regression"
    print(f"  vs previous ({prev_id}):  Δχ²/N={delta:+.4f}  ({direction})")

# ── City parameter delta table (full stage) ───────────────────────────────────

if stage == "full":
    print(f"\n  {'City':<12}  {'ref_pop':>10}  {'tuned_pop':>10}  {'Δpop%':>7}  "
          f"{'ref_wp':>8}  {'tuned_wp':>8}  {'Δwp%':>7}")
    for city_name, city_cfg in city_list:
        rp = city_cfg["ref_pop"]
        tp = city_pops_out[city_name]
        rw = city_cfg["ref_wp"]
        tw = city_wps_out[city_name]
        dp = 100.0 * (tp - rp) / rp
        dw = 100.0 * (tw - rw) / rw if rw > 0 else float("nan")
        print(f"  {city_name:<12}  {rp:>10,.0f}  {tp:>10,.0f}  {dp:>+7.1f}%  "
              f"{rw:>8,.0f}  {tw:>8,.0f}  {dw:>+7.1f}%")

# ── Per-slot fraction table ───────────────────────────────────────────────────

_DT_NAMES = {0: "Wkday", 1: "Sat", 2: "Sun"}
print(f"\n  Per-slot hourly fractions (inferred vs NI prior):")
print(f"  {'Type':<5}  {'Hr':>2}  {'Prior':>9}  {'Inferred':>9}  {'Δ%':>6}  {'|Δ|/σ':>6}  N")
for sk in sorted(slot_fracs):
    dt, h         = sk
    mean_f, std_f = slot_prior[sk]
    f_s           = slot_fracs[sk]
    n_in_slot     = len(slot_groups.get(sk, []))
    delta_pct     = 100.0 * (f_s - mean_f) / mean_f if mean_f > 0 else 0.0
    pull          = abs(f_s - mean_f) / std_f if std_f > 0 else 0.0
    print(f"  {_DT_NAMES[dt]:<5}  {h:>2d}  {mean_f:>9.6f}  {f_s:>9.6f}"
          f"  {delta_pct:>+5.1f}%  {pull:>6.2f}σ  {n_in_slot}")

# ── Goodness-of-fit table (sorted by |z|) ────────────────────────────────────

fit_rows = []
for i_obs, (kind, target, links, _, _) in enumerate(observations):
    if kind == "official":
        lbl = next(s["label"] for s in COUNT_SITES if s["node"] == target)
    else:
        lbl = _link_label(target[0], target[1])
    mod = K * m_arr[i_obs]
    z   = resid[i_obs]
    fit_rows.append((kind, lbl, obs_eff[i_obs], sig_eff[i_obs], mod, z))

fit_rows.sort(key=lambda r: abs(r[5]), reverse=True)
print_chi2_table(fit_rows, chi2, len(fit_rows), n_eff=n_eff)

# ── Save tuned_params.json ────────────────────────────────────────────────────

tuned = {
    "kernel": "rational",
    "K":      round(K, 6),
    "W_BIZ":  round(W_BIZ, 6),
    "P":      round(P, 4),
    "ALPHA":  round(ALPHA, 6),
    **( {"THETA": round(THETA, 6)} if THETA is not None else {} ),
    "slot_fracs": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs.items()},
    "external_node_pop": {str(k): round(v) for k, v in ext_pop_map.items()},
    "external_node_biz": {str(k): round(v) for k, v in ext_biz_map.items()},
    "chi2":       round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "n_obs":      n_obs,
    "n_slots":    n_slots,
    "n_eff":      n_eff,
    "stage":      stage,
}
if stage == "full":
    tuned["external_city_pop"] = {k: round(v) for k, v in city_pops_out.items()}
    tuned["external_city_wp"]  = {k: round(v) for k, v in city_wps_out.items()}
    tuned["external_dampings"] = {k: round(v, 4) for k, v in dampings_out.items()}

with open(TUNED_PARAMS, "w") as f:
    json.dump(tuned, f, indent=2)
print(f"\nSaved → {TUNED_PARAMS}")

# ── Gravity model kernel plot ─────────────────────────────────────────────────

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d_sec  = np.logspace(np.log10(30), np.log10(7200), 500)  # 30s – 120 min
    d_min  = d_sec / 60.0
    u      = d_sec / P
    kernel = K * (ALPHA + 1) * u / (ALPHA + u ** (ALPHA + 1))

    d_peak_sec = P          # peak is exactly at d = P
    k_peak     = K          # f(P) = 1, so kernel peak = K

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(d_min, kernel, linewidth=1.8)
    ax.axvline(d_peak_sec / 60.0, color="r", linestyle="--", alpha=0.7,
               label=f"peak  {d_peak_sec/60.0:.1f} min  (K={k_peak:.3e})")
    ax.set_xlabel("Travel time (minutes)")
    ax.set_ylabel("K · (ALPHA+1) · u / (ALPHA + u^(ALPHA+1))  where u = d/P")
    ax.set_title(
        f"Gravity kernel (rational)  "
        f"P={P:.1f}s  ALPHA={ALPHA:.3f}  K={K:.3e}\n"
        f"stage={stage}  χ²/N={chi2_per_n:.4f}  id={run_id}  git={git_hash}"
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(CURVE_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved → {CURVE_PNG}")
except Exception as _e:
    print(f"Warning: could not save gravity curve plot ({_e})")

# ── Append to tuning history ──────────────────────────────────────────────────

params = {
    "kernel": "rational",
    "K":      round(K, 6),
    "W_BIZ":  round(W_BIZ, 6),
    "P":      round(P, 4),
    "ALPHA":  round(ALPHA, 6),
    **( {"THETA": round(THETA, 6)} if THETA is not None else {} ),
    "slot_fracs": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs.items()},
}
if stage == "full":
    params["external_node_pop"] = {str(k): round(v) for k, v in ext_pop_map.items()}
    params["external_node_biz"] = {str(k): round(v) for k, v in ext_biz_map.items()}
    params["external_city_pop"] = {k: round(v) for k, v in city_pops_out.items()}
    params["external_city_wp"]  = {k: round(v) for k, v in city_wps_out.items()}
    params["external_dampings"] = {k: round(v, 4) for k, v in dampings_out.items()}

history_entry = {
    "id":        run_id,
    "git_hash":  git_hash,
    "timestamp": datetime.now(timezone.utc).isoformat(),
    "stage":     stage,
    "n_evals":   eval_count[0],
    "n_obs":     n_obs,
    "n_slots":   n_slots,
    "n_eff":     n_eff,
    "n_params":  len(log_p0),
    "chi2":      round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "params":    params,
    "observations": [
        {
            "kind":     observations[i_obs][0],
            "target":   (list(observations[i_obs][1])
                         if isinstance(observations[i_obs][1], tuple)
                         else observations[i_obs][1]),
            "observed": round(float(obs_eff[i_obs]), 1),
            "sigma":    round(float(sig_eff[i_obs]), 1),
            "model":    round(float(K * m_arr[i_obs]), 1),
            "z":        round(float(resid[i_obs]), 3),
        }
        for i_obs in range(n_obs)
    ],
}
if note:
    history_entry["note"] = note

with open(HISTORY_FILE, "a") as f:
    f.write(json.dumps(history_entry) + "\n")
print(f"Appended → {HISTORY_FILE}")
