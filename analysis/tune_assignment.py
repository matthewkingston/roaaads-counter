"""
Powell's-method tuning of gravity model parameters against combined chi-squared
objective across all traffic count observations.

Two-component model: residential (pop→pop) and business-adjacent (pb+bb) flows each
carry their own temporal profile (f_s_res, f_s_biz) and global scale (K_res, K_biz).
All four are jointly calibrated at each evaluation via 4-block alternating minimisation:
  K-step   (1D solve for total K; phi-step splits into K_res, K_biz)
  f_res-step  (per-slot analytical update, anchored by NTS residential prior)
  f_biz-step  (per-slot analytical update, anchored by NTS business prior)
  + aggregate coupling γ·(f_res + f_biz − 2·f_agg)² per slot to prevent collective drift.

Observations:
  Official sites: hourly count obs from data/official_hourly.json (24 h × 3 day-types
    × 3 sites = 216 obs), with Gaussian error (sigma from between-day variance).
  Walking obs: per-session count obs from data/link_aadt.json, Poisson error (n_eff).
Both types are in count space, unified in _slot_data with per-obs weights and rhs.

All optimizer parameters are stored in log-space to enforce positivity.

Stage 1 (--gravity, default): tune W_BIZ, P, ALPHA (3 params; 4 with THETA if k=3 cache).
Stage 2 (--full): also tune city-level pop/wp and sub-1 dampings (26 params)
  with L2 regularisation relative to simulation/tuner_config.json references.

Node 180 (pop=50, wp=0) is excluded from stage 2 — too small to tune.

Results:
  simulation/tuned_params.json   best params from this run (read by build_assignment.py)
  simulation/tuning_history.jsonl  appended record of every run
  reports/gravity_model_curve.png  kernel shape plot

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
from model import (EXCLUDE_LINKS, PATHS_CACHE, WEIGHTS_FILE,
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS,
                   print_chi2_table)

# ── Paths ─────────────────────────────────────────────────────────────────────

CONS_GRAPH        = "simulation/newtownards_consolidated.graphml"
HISTORY_FILE      = "simulation/tuning_history.jsonl"
CURVE_PNG         = "reports/gravity_model_curve.png"
HOURLY_FRACS_FILE = "analysis/hourly_fractions.csv"
OFFICIAL_HOURLY   = "data/official_hourly.json"

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

lam                  = config["lambda"]
gamma_coupling_scale = config.get("gamma_coupling_scale",
                                   config.get("gamma_coupling", 1.0))
phi_prior            = config.get("phi_prior", 0.35)
phi_std              = config.get("phi_std",   0.15)
city_list    = list(config["cities"].items())

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

_raw_fracs = {}  # {(dow, hour): (mean_f, std_f, mean_f_res, mean_f_biz)}
with open(HOURLY_FRACS_FILE, newline="") as _fh:
    for _row in csv.DictReader(_fh):
        _dow  = int(_row["day_of_week"])
        _hour = int(_row["hour"].split(":")[0])
        _mfr  = float(_row["mean_fraction_res"]) if "mean_fraction_res" in _row else None
        _mfb  = float(_row["mean_fraction_biz"]) if "mean_fraction_biz" in _row else None
        _raw_fracs[(_dow, _hour)] = (float(_row["mean_fraction"]), float(_row["std_fraction"]),
                                     _mfr, _mfb)

# slot_prior[key] = (mean_f_agg, std_f_agg, mean_f_res, mean_f_biz)
# std_f is the same for both components (aggregate day-to-day variability as proxy)
slot_prior = {}
for _dt, _dows in _DT_DOWS.items():
    for _h in range(24):
        _entries = [_raw_fracs[(_d, _h)] for _d in _dows if (_d, _h) in _raw_fracs]
        if not _entries:
            continue
        _means = [e[0] for e in _entries]
        _stds  = [e[1] for e in _entries]
        _mf      = sum(_means) / len(_means)
        _between = sum((_m - _mf) ** 2 for _m in _means) / len(_means)
        _within  = sum(_s ** 2 for _s in _stds) / len(_stds)
        _std     = math.sqrt(_between + _within)
        # component means: average over day-type's days
        _mfr = (sum(e[2] for e in _entries) / len(_entries)
                if _entries[0][2] is not None else _mf)
        _mfb = (sum(e[3] for e in _entries) / len(_entries)
                if _entries[0][3] is not None else _mf)
        slot_prior[(_dt, _h)] = (_mf, _std, _mfr, _mfb)

# ── Build observation list ────────────────────────────────────────────────────
# All observations are slotted (day_type, hour) in count space:
#   official_hourly: Gaussian error, weight = 1/sigma², T_s = 3600 s
#   walking:         Poisson error,  weight = 1/n_eff,  T_s = duration_s

observations  = []   # (kind, link_key, link_idxs_placeholder, rhs, sigma, Ts_s)
obs_slot_keys = []   # (day_type, hour) per obs
obs_weights   = []   # 1/sigma² (official) or 1/n_eff (walking)
obs_rhs       = []   # target value in count space
obs_Th_lst    = []   # T_s / 3600

# Official hourly obs from ODS-derived JSON (replace single AADT constraints)
n_official_hourly = 0
if os.path.exists(OFFICIAL_HOURLY):
    with open(OFFICIAL_HOURLY) as _f:
        _oh = json.load(_f)
    for _site_id, _site in _oh.items():
        _node  = _site["node"]
        _links = [tuple(lnk) for lnk in _site["links"]] if _site["links"] else None
        for _obs in _site["observations"]:
            _ts    = tuple(_obs["time_slot"])
            _sk    = (_ts[0], _ts[1])   # time_slot already encodes (day_type, hour)
            _count = float(_obs["count"])
            _sigma = float(_obs["sigma"])
            observations.append(("official_hourly", _node, _links, _count, _sigma, 3600.0))
            obs_slot_keys.append(_sk if _sk in slot_prior else None)
            obs_weights.append(1.0 / (_sigma ** 2))
            obs_rhs.append(_count)
            obs_Th_lst.append(1.0)  # T=3600 s → T/3600 = 1
            n_official_hourly += 1
else:
    print(f"  Warning: {OFFICIAL_HOURLY} not found — no official hourly observations")

# Walking obs from link_aadt.json
if os.path.exists(LINK_AADT):
    with open(LINK_AADT) as _f:
        link_aadt_data = json.load(_f)["links"]
    for _key, _entry in sorted(link_aadt_data.items()):
        _u, _v = map(int, _key.split(","))
        if (_u, _v) in EXCLUDE_LINKS:
            continue
        for _sess in _entry["observations"]:
            _neff = float(_sess.get("n_eff", 0.5))
            _dur  = float(_sess.get("duration_s", 0.0))
            _ts   = _sess.get("time_slot")
            _aadt = float(_sess["aadt"])
            _aadt_unc = float(_sess["aadt_uncertainty"])
            observations.append(("walking", (_u, _v), None, _aadt, _aadt_unc, _dur))
            if _ts is not None:
                _sk = (_DOW_TO_TYPE[_ts[0]], _ts[1])
                obs_slot_keys.append(_sk if _sk in slot_prior else None)
            else:
                obs_slot_keys.append(None)
            obs_weights.append(1.0 / _neff)    # Poisson: weight = 1/n_eff
            obs_rhs.append(_neff)              # Poisson: compare to n_eff
            obs_Th_lst.append(_dur / 3600.0)

n_obs      = len(observations)
obs_Th     = np.array(obs_Th_lst,    dtype=np.float64)
obs_w_arr  = np.array(obs_weights,   dtype=np.float64)
obs_rhs_arr= np.array(obs_rhs,       dtype=np.float64)

# Precompute link index sets per observation for fast model flow extraction
obs_link_idxs = []
for kind, target, links, *_ in observations:
    if kind == "official_hourly":
        if links:
            idxs = [link_index[lnk] for lnk in links if lnk in link_index]
        else:
            idxs = [k for k in range(N_links)
                    if int(link_u[k]) == target or int(link_v[k]) == target]
    else:  # walking
        k = link_index.get(target, -1)
        idxs = [k] if k >= 0 else []
    obs_link_idxs.append(idxs)

# Walking obs: vectorized link index array
_walk_link_raw = np.array(
    [obs_link_idxs[i][0] if obs_link_idxs[i] else -1
     for i in range(n_official_hourly, n_obs)],
    dtype=np.int64)
_walk_link_safe = np.where(_walk_link_raw >= 0, _walk_link_raw, 0)
_walk_valid     = _walk_link_raw >= 0

# Group slotted observations by (day_type, hour)
slot_groups = {}
for i, sk in enumerate(obs_slot_keys):
    if sk is not None:
        slot_groups.setdefault(sk, []).append(i)
slot_list = list(slot_groups.items())

n_slots   = len(slot_list)
n_walking = n_obs - n_official_hourly
n_slotted = sum(len(idxs) for _, idxs in slot_list)
# Two df per slot (f_s_res and f_s_biz each consume one df)
n_eff = n_obs - 2 * n_slots
print(f"  {n_obs} observations ({n_official_hourly} official hourly, {n_walking} walking"
      f" in {n_slots} time slot(s))  N_eff={n_eff}")
for sk, idxs in sorted(slot_list):
    if len(idxs) > 3:
        print(f"  Slot {sk}: {len(idxs)} observations")

# Per-slot data for calibrate_Ks_and_fracs and chi²
# Each entry: (slot_key, ia, weights, rhs, Ths,
#              mean_f_res, inv_var_res, mean_f_biz, inv_var_biz,
#              mean_f_agg, gamma)
_slot_data = []
for sk, idxs in slot_list:
    ia       = np.array(idxs, dtype=np.int64)
    mfa, std_f, mfr, mfb = slot_prior[sk]
    inv_var  = 1.0 / (std_f ** 2)
    _slot_data.append((
        sk, ia,
        obs_w_arr[ia],     # weights: 1/sigma² or 1/n_eff per obs
        obs_rhs_arr[ia],   # rhs: count or n_eff per obs
        obs_Th[ia],        # T/3600 per obs
        mfr, inv_var,      # residential prior
        mfb, inv_var,      # business prior (same std as aggregate)
        mfa, gamma_coupling_scale * inv_var,  # aggregate coupling: scale/std_f² per slot
    ))

# ── Assignment and chi-squared helpers ───────────────────────────────────────

def run_assignment(W_BIZ, P, ALPHA, w_pop, w_biz, THETA=None):
    """Gravity assignment → (flow_res, flow_biz).

    flow_res = pop→pop component (purely residential).
    flow_biz = W_BIZ·pb + W_BIZ²·bb component (home↔work/retail).

    THETA=None  → fast binned all-or-nothing (bin-matrix matmul path).
    THETA given → CSR SpMV scatter over k=3 paths with logit weights.
    """
    if THETA is None or not _has_stoch:
        f_b = _kern_b(P, ALPHA)
        W   = W_BIZ
        W2  = W_BIZ ** 2
        if stage == "full":
            flow_res = (int_bin_pp @ f_b).astype(np.float64)
            flow_biz = (W * (int_bin_pb @ f_b) + W2 * (int_bin_bb @ f_b)).astype(np.float64)
            # Exact scatter for external-involved OD pairs
            w_src = w_pop[_ext_src]
            w_dst = w_pop[_ext_dst]
            b_src = w_biz[_ext_src]
            b_dst = w_biz[_ext_dst]
            u_e   = _ext_dist / P
            f_e   = (ALPHA + 1) * u_e / (ALPHA + u_e ** (ALPHA + 1))
            t_pp  = w_src * w_dst * f_e
            t_pb  = (w_src * b_dst + b_src * w_dst) * f_e
            t_bb  = b_src * b_dst * f_e
            flow_res += np.bincount(_ext_link, weights=t_pp[_ext_lp], minlength=N_links)
            flow_biz += np.bincount(_ext_link,
                                    weights=(W_BIZ * t_pb + W_BIZ**2 * t_bb)[_ext_lp],
                                    minlength=N_links)
        else:
            flow_res = (all_bin_pp @ f_b).astype(np.float64)
            flow_biz = (W * (all_bin_pb @ f_b) + W2 * (all_bin_bb @ f_b)).astype(np.float64)
        return flow_res, flow_biz

    # Stochastic logit: exact scatter over 3 paths, split by component.
    # Per-OD weight products in float64 — float32 overflows when the optimizer
    # explores extreme W_BIZ or wp values during Powell's method.
    # Distance arrays stay float32 (large, always finite); upcasting happens
    # naturally when multiplied against float64 weight products.
    pp_od = w_pop[od_src] * w_pop[od_dst]
    pb_od = w_pop[od_src] * w_biz[od_dst] + w_biz[od_src] * w_pop[od_dst]
    bb_od = w_biz[od_src] * w_biz[od_dst]

    # Logit shares (float32 distances sufficient here)
    d_mat  = np.stack([_od_dist_f32, _od_dist_2_f32, _od_dist_3_f32], axis=1)
    log_w  = (-THETA / P) * d_mat
    log_w -= log_w.max(axis=1, keepdims=True)
    shares = np.exp(log_w)
    shares /= shares.sum(axis=1, keepdims=True)

    log_P   = math.log(P)
    alpha1  = ALPHA + 1
    alpha_f = ALPHA

    flow_res = np.zeros(N_links, dtype=np.float64)
    flow_biz = np.zeros(N_links, dtype=np.float64)
    for r, (A_r, d_f32, log_d) in enumerate([
        (_A1, _od_dist_f32,   _log_od_dist),
        (_A2, _od_dist_2_f32, _log_od_dist_2),
        (_A3, _od_dist_3_f32, _log_od_dist_3),
    ]):
        log_u  = log_d - log_P
        u_pow  = np.exp(alpha1 * log_u)
        u_r    = d_f32 * (1.0 / P)
        f_r    = alpha1 * u_r / (alpha_f + u_pow)      # (N_OD,) kernel
        sr     = shares[:, r]
        flow_res += A_r @ (pp_od * sr * f_r)
        flow_biz += A_r @ ((W_BIZ * pb_od + W_BIZ**2 * bb_od) * sr * f_r)
    return flow_res, flow_biz


def model_obs_2c(flow_res, flow_biz):
    """Extract per-observation modelled flows for both components."""
    m_r = np.empty(n_obs)
    m_b = np.empty(n_obs)
    for i, idxs in enumerate(obs_link_idxs[:n_official_hourly]):
        m_r[i] = flow_res[idxs].sum() if idxs else 0.0
        m_b[i] = flow_biz[idxs].sum() if idxs else 0.0
    m_r[n_official_hourly:] = np.where(_walk_valid, flow_res[_walk_link_safe], 0.0)
    m_b[n_official_hourly:] = np.where(_walk_valid, flow_biz[_walk_link_safe], 0.0)
    return m_r, m_b


def calibrate_Ks_and_fracs(m_res, m_biz):
    """4-block alternating minimisation: (K, phi, f_s_res, f_s_biz).

    Reparameterised as K_total = K_res + K_biz and phi = K_biz / K_total.
    This breaks the K_biz / W_BIZ degeneracy that otherwise collapses K_biz → 0.

    K-step:   1D solve (linear in K, same as single-component case).
    phi-step: 1D solve with Gaussian prior phi ~ N(PHI_PRIOR, PHI_STD²).
    f_res-step: per-slot analytical, anchored by NTS residential prior.
    f_biz-step: symmetric.
    Coupling: γ·(f_res + f_biz − 2·f_agg)² per slot (2 lines each f-step).
    Converges in 3–5 iterations; 10 iterations is ample.
    Returns (K_res, K_biz, slot_fracs_res, slot_fracs_biz).
    """
    PHI_PRIOR = phi_prior
    PHI_STD   = phi_std
    INV_PHI_V = 1.0 / (PHI_STD ** 2)

    slot_fracs_res = {sk: mfr for sk, _, _, _, _, mfr, _, _, _, _, _ in _slot_data}
    slot_fracs_biz = {sk: mfb for sk, _, _, _, _, _, _, mfb, _, _, _ in _slot_data}
    K   = 1.0
    phi = PHI_PRIOR

    for _ in range(10):
        K_res = K * (1 - phi)
        K_biz = K * phi

        # ── K-step: 1D solve (linear in K) ───────────────────────────────────
        A = B = 0.0
        for sk, ia, w, rhs, Ths, mfr, ivr, mfb, ivb, mfa, gam in _slot_data:
            f_r = slot_fracs_res[sk]
            f_b = slot_fracs_biz[sk]
            # combined coefficient per obs: (1-phi)*c_r*f_r + phi*c_b*f_b
            coeff = ((1 - phi) * m_res[ia] * Ths * f_r
                     + phi     * m_biz[ia] * Ths * f_b)
            A += float(np.dot(w * coeff, coeff))
            B += float(np.dot(w * coeff, rhs))
        K = max(B / A, 1e-30) if A > 0 else K

        # ── phi-step: 1D solve (linear in phi) ───────────────────────────────
        A_p = B_p = 0.0
        for sk, ia, w, rhs, Ths, mfr, ivr, mfb, ivb, mfa, gam in _slot_data:
            f_r  = slot_fracs_res[sk]
            f_b  = slot_fracs_biz[sk]
            base = m_res[ia] * Ths * f_r             # c_r·f_r
            delt = m_biz[ia] * Ths * f_b - base      # c_b·f_b − c_r·f_r
            # pred = K*(base + phi*delt); minimise Σ w*(K*(base+phi*delt)−rhs)²
            A_p += float(np.dot(w * delt, delt)) * K ** 2
            B_p += float(np.dot(w * delt, rhs - K * base)) * K
        A_p += INV_PHI_V
        B_p += PHI_PRIOR * INV_PHI_V
        phi = max(0.01, min(0.99, B_p / A_p)) if A_p > 0 else phi
        K_res = K * (1 - phi)
        K_biz = K * phi

        # ── f_res-step (per slot, holding K_res, K_biz, f_biz fixed) ─────────
        for sk, ia, w, rhs, Ths, mfr, ivr, mfb_prior, ivb, mfa, gam in _slot_data:
            f_b = slot_fracs_biz[sk]
            cr  = m_res[ia] * Ths
            cb  = m_biz[ia] * Ths
            h_i = rhs - K_biz * cb * f_b
            num = K_res * float(np.dot(w * cr, h_i)) + mfr * ivr + gam * (2*mfa - f_b)
            den = K_res**2 * float(np.dot(w * cr, cr)) + ivr + gam
            slot_fracs_res[sk] = max(num / den, 1e-12) if den > 0 else mfr

        # ── f_biz-step (symmetric) ────────────────────────────────────────────
        for sk, ia, w, rhs, Ths, mfr_prior, ivr, mfb, ivb, mfa, gam in _slot_data:
            f_r = slot_fracs_res[sk]
            cr  = m_res[ia] * Ths
            cb  = m_biz[ia] * Ths
            h_i = rhs - K_res * cr * f_r
            num = K_biz * float(np.dot(w * cb, h_i)) + mfb * ivb + gam * (2*mfa - f_r)
            den = K_biz**2 * float(np.dot(w * cb, cb)) + ivb + gam
            slot_fracs_biz[sk] = max(num / den, 1e-12) if den > 0 else mfb

    return K_res, K_biz, slot_fracs_res, slot_fracs_biz


# ── Objective function ────────────────────────────────────────────────────────

eval_count = [0]
best       = {"chi2": float("inf"), "log_params": None, "K_res": 1.0, "K_biz": 1.0}
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

    flow_res, flow_biz                          = run_assignment(W_BIZ, P, ALPHA, w_pop, w_biz, THETA)
    m_res, m_biz                               = model_obs_2c(flow_res, flow_biz)
    K_res, K_biz, slot_fracs_res, slot_fracs_biz = calibrate_Ks_and_fracs(m_res, m_biz)
    chi2 = 0.0
    for sk, ia, w, rhs, Ths, mfr, ivr, mfb, ivb, mfa, gam in _slot_data:
        f_r  = slot_fracs_res[sk]
        f_b  = slot_fracs_biz[sk]
        pred = K_res * m_res[ia] * Ths * f_r + K_biz * m_biz[ia] * Ths * f_b
        chi2 += float(np.dot(w * (pred - rhs), pred - rhs))
        chi2 += (f_r - mfr) ** 2 * ivr
        chi2 += (f_b - mfb) ** 2 * ivb
        chi2 += gam * (f_r + f_b - 2 * mfa) ** 2

    if stage == "full" and log_ref is not None:
        chi2 += lam * float(np.sum((log_params[n_gravity:] - log_ref[n_gravity:]) ** 2))
    if grav_lam > 0:
        chi2 += grav_lam * float(np.sum((log_params[:n_gravity] - log_grav_ref) ** 2))

    eval_count[0] += 1
    if chi2 < best["chi2"]:
        best["chi2"]       = chi2
        best["log_params"] = log_params.copy()
        best["K_res"]      = K_res
        best["K_biz"]      = K_biz

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

# Capture starting gravity params for history (before any optimization)
initial_gravity = {k: grav_start[k] for k in ("W_BIZ", "P", "ALPHA") if k in grav_start}
if _has_stoch and "THETA" in grav_start:
    initial_gravity["THETA"] = grav_start["THETA"]

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

# Final evaluation for clean chi2, K_res, K_biz, and slot_fracs (no L2 term)
w_pop_f = base_w_pop.copy() if ext_pop_map else base_w_pop
w_biz_f = base_w_biz.copy() if ext_biz_map else base_w_biz
for arr_i, nid in ext_indices:
    if nid in ext_pop_map:
        w_pop_f[arr_i] = ext_pop_map[nid]
        w_biz_f[arr_i] = ext_biz_map[nid]

flow_res, flow_biz                             = run_assignment(W_BIZ, P, ALPHA, w_pop_f, w_biz_f, THETA)
m_res, m_biz                                   = model_obs_2c(flow_res, flow_biz)
K_res, K_biz, slot_fracs_res, slot_fracs_biz  = calibrate_Ks_and_fracs(m_res, m_biz)
K = K_res + K_biz   # total scale (for display and backward compat)
chi2 = 0.0
for sk, ia, w, rhs, Ths, mfr, ivr, mfb, ivb, mfa, gam in _slot_data:
    f_r  = slot_fracs_res[sk]
    f_b  = slot_fracs_biz[sk]
    pred = K_res * m_res[ia] * Ths * f_r + K_biz * m_biz[ia] * Ths * f_b
    chi2 += float(np.dot(w * (pred - rhs), pred - rhs))
    chi2 += (f_r - mfr) ** 2 * ivr
    chi2 += (f_b - mfb) ** 2 * ivb
    chi2 += gam * (f_r + f_b - 2 * mfa) ** 2
chi2_per_n = chi2 / n_obs

# Build per-obs residuals for the fit table.
# Convert count-space obs to an effective hourly-average for display:
#   for walking obs: n_eff / (T/3600) = expected hourly count (un-slots the fraction)
#   for official hourly obs: count is already vehicles/hour
# The total modelled hourly count is K_res*m_res*f_res + K_biz*m_biz*f_biz.
obs_eff = np.empty(n_obs)
sig_eff = np.empty(n_obs)
mod_eff = np.empty(n_obs)
for i, (kind, target, links, _obs, _sig, Ts) in enumerate(observations):
    sk = obs_slot_keys[i]
    Th = obs_Th[i]
    if kind == "official_hourly":
        obs_eff[i] = obs_rhs_arr[i]                   # vehicles/hour
        sig_eff[i] = _sig
        if sk is not None:
            f_r = slot_fracs_res.get(sk, slot_prior[sk][2])
            f_b = slot_fracs_biz.get(sk, slot_prior[sk][3])
            mod_eff[i] = K_res * m_res[i] * Th * f_r + K_biz * m_biz[i] * Th * f_b
        else:
            mod_eff[i] = (K_res * m_res[i] + K_biz * m_biz[i]) * Th
    else:  # walking: show as vehicles/hour using slot fraction
        n_eff_i = obs_rhs_arr[i]
        if sk is not None and Th > 0 and n_eff_i > 0:
            f_r = slot_fracs_res.get(sk, slot_prior[sk][2])
            f_b = slot_fracs_biz.get(sk, slot_prior[sk][3])
            # Combined f_s for this obs (weighted by component flow at this link)
            m_r_i, m_b_i = m_res[i], m_biz[i]
            m_tot = m_r_i * f_r + m_b_i * f_b
            f_eff = m_tot / (m_r_i + m_b_i) if (m_r_i + m_b_i) > 0 else (f_r + f_b) / 2
            obs_eff[i] = n_eff_i / (Th * f_eff)
            sig_eff[i] = math.sqrt(n_eff_i) / (Th * f_eff)
            mod_eff[i] = (K_res * m_r_i * f_r + K_biz * m_b_i * f_b) / f_eff
        else:
            obs_eff[i] = n_eff_i / max(Th, 1e-9)
            sig_eff[i] = math.sqrt(n_eff_i) / max(Th, 1e-9)
            mod_eff[i] = (K_res * m_res[i] + K_biz * m_biz[i])
resid = np.where(sig_eff > 0, (mod_eff - obs_eff) / sig_eff, 0.0)

elapsed = time.time() - t0
print(f"\nResult  ({eval_count[0]} evals, {elapsed:.0f}s)")
_theta_str = f"  THETA={THETA:.4f}" if THETA is not None else ""
print(f"  K_res={K_res:.4e}  K_biz={K_biz:.4e}  (K={K:.4e})")
print(f"  W_BIZ={W_BIZ:.4f}  P={P:.2f}s  ALPHA={ALPHA:.4f}{_theta_str}")
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
print(f"\n  Per-slot hourly fractions (res / biz vs NTS priors):")
print(f"  {'Type':<5}  {'Hr':>2}  {'PriorAgg':>9}  {'f_res':>9}  {'f_biz':>9}  {'Δres%':>6}  {'Δbiz%':>6}  N")
for sk in sorted(slot_fracs_res):
    dt, h              = sk
    mfa, std_f, mfr, mfb = slot_prior[sk]
    f_r                = slot_fracs_res[sk]
    f_b                = slot_fracs_biz[sk]
    n_in_slot          = len(slot_groups.get(sk, []))
    dr  = 100.0 * (f_r - mfr) / mfr if mfr > 0 else 0.0
    db  = 100.0 * (f_b - mfb) / mfb if mfb > 0 else 0.0
    print(f"  {_DT_NAMES[dt]:<5}  {h:>2d}  {mfa:>9.6f}  {f_r:>9.6f}  {f_b:>9.6f}"
          f"  {dr:>+5.1f}%  {db:>+5.1f}%  {n_in_slot}")

# ── Goodness-of-fit table (sorted by |z|) ────────────────────────────────────

fit_rows = []
for i_obs, (kind, target, links, _, _, _) in enumerate(observations):
    if kind == "official_hourly":
        sk  = obs_slot_keys[i_obs]
        lbl = f"{target} h{sk[1]:02d}{'WD' if sk[0]==0 else ('Sa' if sk[0]==1 else 'Su')}"
    else:
        lbl = _link_label(target[0], target[1])
    fit_rows.append((kind, lbl, obs_eff[i_obs], sig_eff[i_obs], mod_eff[i_obs], resid[i_obs]))

fit_rows.sort(key=lambda r: abs(r[5]), reverse=True)
print_chi2_table(fit_rows, chi2, len(fit_rows), n_eff=n_eff)

# ── Save tuned_params.json ────────────────────────────────────────────────────

tuned = {
    "kernel":  "rational",
    "K":       round(K, 6),      # K_res + K_biz (backward compat for build_assignment.py)
    "K_res":   round(K_res, 6),
    "K_biz":   round(K_biz, 6),
    "W_BIZ":   round(W_BIZ, 6),
    "P":       round(P, 4),
    "ALPHA":   round(ALPHA, 6),
    **( {"THETA": round(THETA, 6)} if THETA is not None else {} ),
    "slot_fracs_res": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_res.items()},
    "slot_fracs_biz": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_biz.items()},
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
    os.makedirs("reports", exist_ok=True)
    fig.savefig(CURVE_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved → {CURVE_PNG}")
except Exception as _e:
    print(f"Warning: could not save gravity curve plot ({_e})")

# ── Append to tuning history ──────────────────────────────────────────────────

params = {
    "kernel": "rational",
    "K":      round(K, 6),
    "K_res":  round(K_res, 6),
    "K_biz":  round(K_biz, 6),
    "W_BIZ":  round(W_BIZ, 6),
    "P":      round(P, 4),
    "ALPHA":  round(ALPHA, 6),
    **( {"THETA": round(THETA, 6)} if THETA is not None else {} ),
    "slot_fracs_res": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_res.items()},
    "slot_fracs_biz": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_biz.items()},
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
    "tuner_hyperparams": {
        "phi_prior":            phi_prior,
        "phi_std":              phi_std,
        "gamma_coupling_scale": gamma_coupling_scale,
        "gravity_lambda":       grav_lam,
        "lambda":               lam,
    },
    "initial_gravity": {k: round(v, 6) for k, v in initial_gravity.items()},
    "slot_prior": {
        f"{dt},{h}": [round(mfa, 8), round(std_f, 8), round(mfr, 8), round(mfb, 8)]
        for (dt, h), (mfa, std_f, mfr, mfb) in slot_prior.items()
    },
    "observations": [
        {
            "kind":     observations[i_obs][0],
            "label":    (f"{observations[i_obs][1]} h{obs_slot_keys[i_obs][1]:02d}"
                         if observations[i_obs][0] == "official_hourly"
                         else _link_label(observations[i_obs][1][0],
                                          observations[i_obs][1][1])),
            "observed": round(float(obs_eff[i_obs]), 1),
            "sigma":    round(float(sig_eff[i_obs]), 1),
            "model":    round(float(mod_eff[i_obs]), 1),
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
