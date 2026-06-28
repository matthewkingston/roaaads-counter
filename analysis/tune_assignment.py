"""
Powell's-method tuning of gravity model parameters against combined chi-squared
objective across all traffic count observations.

Three-component model: residential (pop→pop), business-adjacent (pb+bb), and school
(pop→school) flows each carry their own temporal profile and global scale.
Jointly calibrated at each evaluation via 5-block alternating minimisation:
  K-step     (1D solve for total K)
  phi_biz-step  (1D solve for business fraction, with Gaussian prior)
  phi_sch-step  (1D solve for school fraction, with Gaussian prior)
  f_res-step    (per-slot analytical, anchored by NTS residential prior)
  f_biz-step / f_school-step  (symmetric)
  + aggregate coupling γ·(f_res + f_biz + f_school − f_agg)² per slot.

Observations:
  Official sites: hourly count obs from data/official_hourly.json (24 h × 3 day-types
    × 3 sites = 216 obs), with Gaussian error (sigma from between-day variance).
  Walking obs: per-session count obs from data/link_aadt.json, Poisson error (n_eff).
Both types are in count space, unified in _slot_data with per-obs weights and rhs.

All optimizer parameters are stored in log-space to enforce positivity.

Tunes W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, P_school, ALPHA_school (8 params with school).
Production-constrained per component: each origin's trip production is fixed by its
producing weight, independent of accessibility (school magnitude is K_sch; W_SCHOOL removed
as redundant under the constraint).
External zone values are fixed from census data and are not tuned.

Results:
  simulation/tuned_params.json   best params from this run (read by build_assignment.py)
  simulation/tuning_history.jsonl  appended record of every run
  reports/gravity_model_curve.png  kernel shape plot

Usage:
  python3 analysis/tune_assignment.py
  python3 analysis/tune_assignment.py --note "added-june-counts"
  python3 analysis/tune_assignment.py --fast          # looser tolerances, ~2× faster
"""

import csv, json, math, os, secrets, subprocess, sys, time, xml.etree.ElementTree as ET
from datetime import datetime, timezone

import numpy as np
import scipy.optimize

sys.path.insert(0, "simulation")
from model import (EXCLUDE_LINKS, PATHS_CACHE, WEIGHTS_FILE,
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS,
                   constrained_od_flows, scatter_od_to_links, load_self_terms,
                   print_chi2_table, assert_paths_cache_fresh,
                   format_slot_time, nice_official)

# ── Paths ─────────────────────────────────────────────────────────────────────

CONS_GRAPH        = "simulation/newtownards_reduced.graphml"  # dead-end-reduced (street names)
HISTORY_FILE      = "simulation/tuning_history.jsonl"
CURVE_PNG         = "reports/gravity_model_curve.png"
HOURLY_FRACS_FILE = "analysis/hourly_fractions.csv"
OFFICIAL_HOURLY   = "data/official_hourly.json"

# ── CLI args ──────────────────────────────────────────────────────────────────

note     = None
fast     = False
argv  = sys.argv[1:]
i = 0
while i < len(argv):
    if argv[i] == "--fast":
        fast = True
    elif argv[i] == "--f-frozen":
        # Deprecated no-op: temporal fractions are now ALWAYS pinned at the NTS
        # profile (f is never tuned), so this flag is the default behaviour.
        print("note: --f-frozen is deprecated and now the default "
              "(temporal fractions are always pinned at the NTS profile); ignoring.")
    elif argv[i] == "--note" and i + 1 < len(argv):
        i += 1
        note = argv[i]
    i += 1

# Temporal fractions f_res/f_biz/f_school are pinned at the NTS profile
# (mean_fraction_*) and never tuned — every residual is purely spatial. The inner
# calibration therefore solves only the three component scales (K_res, K_biz, K_sch)
# directly via a convex Newton solve (see solve_scales); no φ reparam, no f-steps,
# no alternation, no best-iterate bookkeeping.

run_id = secrets.token_hex(4)

try:
    git_hash = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
    ).decode().strip()
except Exception:
    git_hash = "unknown"

print(f"Run ID: {run_id}  git: {git_hash}" +
      ("  FAST" if fast else "") +
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
cache        = np.load(PATHS_CACHE, allow_pickle=True)
assert_paths_cache_fresh(cache)
node_ids_arr = cache["node_ids"]
od_src       = cache["od_src"]
od_dst       = cache["od_dst"]
od_dist      = cache["od_dist"].astype(np.float64)
pair_idx     = cache["pair_idx"]
link_idx_arr = cache["link_idx"]
link_u       = cache["link_u"]
link_v       = cache["link_v"]
N_links      = len(link_u)
node_ids     = list(node_ids_arr)
N_nodes      = len(node_ids)

# External intra-zonal self-term (denominator-only; from build_intra_times.py).
# Constant across evals — only the kernel F(self_dist) recomputes inside run_assignment.
_self_src, _self_dist, _self_w = load_self_terms(node_ids)

link_index = {(int(link_u[k]), int(link_v[k])): k for k in range(N_links)}

# Probit fractional link weights (new cache format); None for legacy binary caches.
_link_weight = cache["link_weight"].astype(np.float32) if "link_weight" in cache else None

_has_stoch = "pair_idx_2" in cache
if _has_stoch:
    _od_dist_2  = cache["od_dist_2"].astype(np.float64)
    _pair_idx_2 = cache["pair_idx_2"]
    _link_idx_2 = cache["link_idx_2"]
    _od_dist_3  = cache["od_dist_3"].astype(np.float64)
    _pair_idx_3 = cache["pair_idx_3"]
    _link_idx_3 = cache["link_idx_3"]
    # Per-OD-pair float32 distances and precomputed log(d) for the gravity kernel.
    # u^(ALPHA+BETA) = exp((ALPHA+BETA)*(log_d - log_P)); log_d fixed, log_P cheap scalar.
    _od_dist_f32   = od_dist.astype(np.float32)
    _od_dist_2_f32 = _od_dist_2.astype(np.float32)
    _od_dist_3_f32 = _od_dist_3.astype(np.float32)
    _log_od_dist   = np.log(od_dist).astype(np.float32)
    _log_od_dist_2 = np.log(_od_dist_2).astype(np.float32)
    _log_od_dist_3 = np.log(_od_dist_3).astype(np.float32)
    # Stacked constants for batch logit-share and kernel computation (allocated once).
    _d_mat_f32   = np.column_stack([_od_dist_f32, _od_dist_2_f32, _od_dist_3_f32])
    _log_d_stack = np.column_stack([_log_od_dist, _log_od_dist_2, _log_od_dist_3])
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

if _link_weight is not None:
    _n_passes = int(cache["probit_n_passes"])
    _cv       = float(cache["probit_cv"])
    print(f"  {N_nodes} nodes  {N_links} links  {len(od_src):,} OD pairs"
          f"  probit loading: {_n_passes} passes  CV={_cv:.2f}")
else:
    print(f"  {N_nodes} nodes  {N_links} links  {len(od_src):,} OD pairs"
          + ("  stochastic k=3 paths loaded" if _has_stoch else "  (no stochastic paths — run build_paths.py)"))

# ── Load node weights ─────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    wdata = json.load(f)

_pnid = lambda k: (int(k) if k.lstrip("-").isdigit() else k)
node_pop_full    = {_pnid(k): v for k, v in wdata["node_population"].items()}
node_biz_full    = {_pnid(k): v for k, v in wdata["node_business_demand"].items()}
node_school_full = {_pnid(k): v for k, v in wdata.get("node_school_demand", {}).items()}
node_school_prod_full = {_pnid(k): v for k, v in wdata.get("node_school_producers", {}).items()}

# Precomputed base weight arrays (from census + OSM demand; external zones fixed)
base_w_pop    = np.array([node_pop_full.get(nid, 0.0)    for nid in node_ids], dtype=np.float64)
base_w_biz    = np.array([node_biz_full.get(nid, 0.0)    for nid in node_ids], dtype=np.float64)
base_w_school = np.array([node_school_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
base_w_school_prod = np.array([node_school_prod_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
if base_w_school_prod.sum() == 0:
    base_w_school_prod = None   # fall back to population producer (legacy weights)

_has_school = base_w_school.sum() > 0
if not _has_school:
    print("  Warning: no node_school_demand in weights — school component disabled")

# ── Load street names from consolidated GraphML (sequential node IDs) ─────────

link_name = {}
if os.path.exists(CONS_GRAPH):
    try:
        _tree = ET.parse(CONS_GRAPH)
        _root = _tree.getroot()
        _nsmap = {"g": _root.tag.split("}")[0].lstrip("{")} if "}" in _root.tag else {"g": ""}
        _pfx = "{" + _nsmap["g"] + "}" if _nsmap["g"] else ""
        # Resolve the data-key id for the edge "name" attribute.  GraphML dN key
        # ids are assigned in attribute-appearance order and are NOT stable across
        # network regenerations, so look the id up from the <key> header rather
        # than hardcoding it (a hardcoded "d14" silently became "oneway").
        _name_key = next(
            (_k.get("id") for _k in _root.iter(f"{_pfx}key")
             if _k.get("for") == "edge" and _k.get("attr.name") == "name"),
            None)
        if _name_key is None:
            print("Warning: no edge 'name' attribute in GraphML — labels omit street names")
        else:
            for _edge in _root.iter(f"{_pfx}edge"):
                _u = int(_edge.get("source"))
                _v = int(_edge.get("target"))
                for _data in _edge:
                    if _data.get("key") == _name_key and _data.text:
                        link_name[(_u, _v)] = _data.text
                        break
    except Exception as _e:
        print(f"Warning: could not load street names ({_e})")


def _link_label(u, v):
    """Human-readable street name for a directed link ('' if unnamed).

    Used as the model.compute_chi2 label_fn and to populate the report Label column;
    the precise u→v reference now lives in the table's separate Link column.
    """
    return link_name.get((u, v), "")

# ── Load tuner config ─────────────────────────────────────────────────────────

with open(TUNER_CONFIG) as f:
    config = json.load(f)

lam                  = config["lambda"]
gamma_coupling_scale = config.get("gamma_coupling_scale",
                                   config.get("gamma_coupling", 1.0))
phi_prior            = config.get("phi_biz_prior", config.get("phi_prior",  0.35))
phi_std              = config.get("phi_biz_std",   config.get("phi_std",    0.15))
phi_school_prior     = config.get("phi_school_prior",  0.10)
phi_school_std       = config.get("phi_school_std",    0.08)

grav_ref = config.get("gravity_ref", {})
grav_lam_raw = config.get("gravity_lambda", 0.0)
_grav_param_names = ["W_BIZ", "P", "ALPHA", "BETA", "P_biz", "ALPHA_biz"]
if _has_stoch:
    _grav_param_names.append("THETA")
if _has_school:
    # W_SCHOOL removed: under the production constraint the school component is
    # K_sch·(per-origin-normalised pop↔school split), so a separate W_SCHOOL scale is
    # exactly redundant with K_sch. Only the school KERNEL SHAPE (P_school, ALPHA_school)
    # is tuned; its magnitude comes from K_sch (analytical).
    _grav_param_names += ["P_school", "ALPHA_school"]
_grav_ref_vals = [
    math.log(max(grav_ref.get("W_BIZ",     1.0),  1e-4)),
    math.log(max(grav_ref.get("P",         600.0), 1e-4)),
    math.log(max(grav_ref.get("ALPHA",     2.0),   1e-4)),
    math.log(max(grav_ref.get("BETA",      1.0),   1e-4)),
    math.log(max(grav_ref.get("P_biz",     600.0), 1e-4)),
    math.log(max(grav_ref.get("ALPHA_biz", 2.0),   1e-4)),
]
if _has_stoch:
    _grav_ref_vals.append(math.log(max(grav_ref.get("THETA", 1.0), 1e-4)))
if _has_school:
    _grav_ref_vals += [
        math.log(max(grav_ref.get("P_school",     600.0), 1e-4)),
        math.log(max(grav_ref.get("ALPHA_school", 2.0),   1e-4)),
    ]
log_grav_ref = np.array(_grav_ref_vals)
if isinstance(grav_lam_raw, dict):
    log_grav_lam = np.array([grav_lam_raw.get(k, 0.0) for k in _grav_param_names])
else:
    log_grav_lam = np.full(len(_grav_param_names), float(grav_lam_raw))

n_gravity = (7 if _has_stoch else 6) + (2 if _has_school else 0)

# External node weights come from node_weights.json (census data, fixed — not tuned).
# Stochastic-path weight products (constant across all evaluations).
if _has_stoch:
    _pp_od_s1 = base_w_pop[od_src] * base_w_pop[od_dst]
    _pb_od_s1 = (base_w_pop[od_src] * base_w_biz[od_dst]
                 + base_w_biz[od_src] * base_w_pop[od_dst])
    _bb_od_s1 = base_w_biz[od_src] * base_w_biz[od_dst]

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

# NOTE: The legacy distance-bin link matrices are gone. The production-constrained
# assignment normalises each origin's flow by its own (per-eval, kernel-dependent)
# denominator D_i, which cannot be pre-summed into static distance bins the way the
# unconstrained w_i·w_j product could. run_assignment now calls model.constrained_od_flows
# (per-pair flows + cheap per-origin bincount denominators) and scatters via the probit
# routing incidence — see the production-constrained gravity plan / project memory note.

N_nodes = len(node_ids)   # for the per-origin denominator bincounts

# ── Per-slot hourly fraction priors ───────────────────────────────────────────
# Group days into weekday (0), Saturday (1), Sunday (2).
# Prior mean and std for each (day_type, hour) slot are derived from
# hourly_fractions.csv via the law of total variance:
#   total_var = between_day_var(means) + mean(within_day_var)
# This equals the pooled variance computed directly from the raw NI count data.

_DOW_TO_TYPE = {d: (0 if d < 5 else (1 if d == 5 else 2)) for d in range(7)}
_DT_DOWS     = {0: list(range(5)), 1: [5], 2: [6]}

_raw_fracs = {}  # {(dow, hour): (mean_f, std_f, mean_f_res, mean_f_biz, mean_f_school)}
with open(HOURLY_FRACS_FILE, newline="") as _fh:
    for _row in csv.DictReader(_fh):
        _dow  = int(_row["day_of_week"])
        _hour = int(_row["hour"].split(":")[0])
        _mfr  = float(_row["mean_fraction_res"])    if "mean_fraction_res"    in _row else None
        _mfb  = float(_row["mean_fraction_biz"])    if "mean_fraction_biz"    in _row else None
        _mfs  = float(_row["mean_fraction_school"]) if "mean_fraction_school" in _row else None
        _raw_fracs[(_dow, _hour)] = (float(_row["mean_fraction"]), float(_row["std_fraction"]),
                                     _mfr, _mfb, _mfs)

# slot_prior[key] = (mean_f_agg, std_f_agg, mean_f_res, mean_f_biz, mean_f_school)
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
        _mfr = (sum(e[2] for e in _entries) / len(_entries)
                if _entries[0][2] is not None else _mf * 0.5)
        _mfb = (sum(e[3] for e in _entries) / len(_entries)
                if _entries[0][3] is not None else _mf * 0.35)
        _mfs = (sum(e[4] for e in _entries) / len(_entries)
                if _entries[0][4] is not None else 0.0)
        slot_prior[(_dt, _h)] = (_mf, _std, _mfr, _mfb, _mfs)

# ── Build observation list ────────────────────────────────────────────────────
# All observations are slotted (day_type, hour) in count space:
#   official_hourly: Gaussian error, weight = 1/sigma², T_s = 3600 s
#   walking:         Poisson error,  weight = 1/n_eff,  T_s = duration_s

observations  = []   # (kind, link_key, link_idxs_placeholder, rhs, sigma, Ts_s)
obs_slot_keys = []   # (day_type, hour) per obs
obs_weights   = []   # 1/sigma² (official) or 1/n_eff (walking)
obs_rhs       = []   # target value in count space
obs_Th_lst    = []   # T_s / 3600
obs_meta      = []   # {"label": display name, "link": precise u→v / node ref} per obs

# Official hourly obs from ODS-derived JSON (replace single AADT constraints)
n_official_hourly = 0
if os.path.exists(OFFICIAL_HOURLY):
    with open(OFFICIAL_HOURLY) as _f:
        _oh = json.load(_f)
    for _site_id, _site in _oh.items():
        _node  = _site["node"]
        _links = [tuple(lnk) for lnk in _site["links"]] if _site["links"] else None
        _site_label, _site_link = nice_official(_site)
        for _obs in _site["observations"]:
            _ts    = tuple(_obs["time_slot"])
            _sk    = (_ts[0], _ts[1])   # time_slot already encodes (day_type, hour)
            _count = float(_obs["count"])
            _sigma = float(_obs["sigma"])
            observations.append(("official_hourly", _node, _links, _count, _sigma, 3600.0))
            obs_meta.append({"label": _site_label, "link": _site_link})
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
            obs_meta.append({"label": link_name.get((_u, _v), "") or "(unnamed)",
                             "link":  f"{_u}→{_v}"})
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

# ── Observed-link scatter restriction (tuner-only speedup) ────────────────────
# model_obs_3c reads modelled flow on only the OBSERVED links (the official-site +
# walking links in obs_link_idxs / _walk_link_safe). run_assignment, however, used
# to scatter the full ~62M probit-incidence entries onto ALL N_links every eval,
# then discard flow on the ~1400 links no observation ever touches. The scatter
# dominates each eval (~3×2.8 s), so we precompute the subset of incidence entries
# landing on observed links and scatter only those into a COMPACT observed-link
# space (link id → 0..n_obs_links-1). Identical result on observed links; build
# time only (no per-eval cost). build_assignment.py keeps the full scatter for the
# map/flow outputs — this restriction is local to the tuner's objective.
_obs_link_ids = np.array(sorted({k for idxs in obs_link_idxs for k in idxs}),
                         dtype=np.int64)
_N_obs_links  = len(_obs_link_ids)
_link_remap   = np.full(N_links, -1, dtype=np.int64)   # full link id → compact id (or -1)
_link_remap[_obs_link_ids] = np.arange(_N_obs_links)

_compact_link_idx = _link_remap[link_idx_arr]          # -1 where the entry's link is unobserved
_scatter_keep     = _compact_link_idx >= 0
_pair_idx_obs     = np.ascontiguousarray(pair_idx[_scatter_keep])
_link_idx_obs     = np.ascontiguousarray(_compact_link_idx[_scatter_keep])
_link_weight_obs  = (np.ascontiguousarray(_link_weight[_scatter_keep])
                     if _link_weight is not None else None)

# Readback indices remapped into the compact observed-link space.
_obs_link_idxs_compact = [[int(_link_remap[k]) for k in idxs] for idxs in obs_link_idxs]
_walk_link_safe_compact = np.where(_walk_valid, _link_remap[_walk_link_safe], 0)
print(f"  Observed-link scatter: {_N_obs_links} links carry "
      f"{_scatter_keep.sum():,}/{len(link_idx_arr):,} incidence entries "
      f"({100*_scatter_keep.sum()/len(link_idx_arr):.0f}%) — scatter restricted to these")

# Group slotted observations by (day_type, hour)
slot_groups = {}
for i, sk in enumerate(obs_slot_keys):
    if sk is not None:
        slot_groups.setdefault(sk, []).append(i)
slot_list = list(slot_groups.items())

n_slots   = len(slot_list)
n_walking = n_obs - n_official_hourly
n_slotted = sum(len(idxs) for _, idxs in slot_list)
# Temporal fractions are pinned at the NTS profile (never fitted), so NO per-slot
# temporal degrees of freedom are consumed — N_eff counts all observations. (The
# few global df — gravity shape params + 3 scales — are not subtracted, consistent
# with the prior convention that only ever counted per-slot temporal df.) This
# corrects the old `n_obs − 3·n_slots`, which subtracted df for fractions the
# --f-frozen path never actually fit; the χ²/N basis therefore changes and is NOT
# comparable to pre-2026-06-27 history.
n_eff = n_obs
print(f"  {n_obs} observations ({n_official_hourly} official hourly, {n_walking} walking"
      f" in {n_slots} time slot(s))  N_eff={n_eff}")
print("  Objective: Gaussian chi² for official hourly; Poisson deviance for walking obs."
      "  χ²/N is a mixed criterion — not comparable to pre-Poisson runs.")
for sk, idxs in sorted(slot_list):
    if len(idxs) > 3:
        print(f"  Slot {sk}: {len(idxs)} observations")

# Per-slot data: builds the constant NTS f vectors (below) and the chi² weight arrays
# Each entry: (slot_key, ia, weights, rhs, Ths,
#              mean_f_res, inv_var_res, mean_f_biz, inv_var_biz,
#              mean_f_school, inv_var_school, mean_f_agg, gamma)
_slot_data = []
for sk, idxs in slot_list:
    ia       = np.array(idxs, dtype=np.int64)
    mfa, std_f, mfr, mfb, mfs = slot_prior[sk]
    inv_var  = 1.0 / (std_f ** 2)
    _slot_data.append((
        sk, ia,
        obs_w_arr[ia],     # weights: 1/sigma² or 1/n_eff per obs
        obs_rhs_arr[ia],   # rhs: count or n_eff per obs
        obs_Th[ia],        # T/3600 per obs
        mfr, inv_var,      # residential prior
        mfb, inv_var,      # business prior (same std as aggregate)
        mfs, inv_var,      # school prior (same std as aggregate)
        mfa, gamma_coupling_scale * inv_var,  # aggregate coupling: scale/std_f² per slot
    ))

# Per-obs weight array with unslotted obs zeroed out (used in vectorised chi² in objective).
_slot_w_arr = np.zeros(n_obs, dtype=np.float64)
for _, _ia, _w, *_ in _slot_data:
    _slot_w_arr[_ia] = _w

# ── Poisson deviance setup for walking obs ────────────────────────────────────
# Walking obs occupy positions [n_official_hourly:].
# Actual integer count: n = n_eff − 0.5 (inverse of the Jeffreys n_eff = n + 0.5).
# Official hourly obs are Gaussian; walking obs use Poisson deviance 2(pred − n·log pred).
_walk_n_arr = np.zeros(n_obs, dtype=np.float64)
_walk_n_arr[n_official_hourly:] = obs_rhs_arr[n_official_hourly:] - 0.5

# Slotted walking obs mask (unslotted obs have _slot_w_arr = 0; mirror that here).
_walk_slotted = np.zeros(n_obs, dtype=bool)
_walk_slotted[n_official_hourly:] = _slot_w_arr[n_official_hourly:] > 0

# Gaussian-only weight array: official hourly weights unchanged, walking zeroed.
# Used in K-step and f_s-steps to isolate the Gaussian chi-squared contribution.
_gauss_w_arr = _slot_w_arr.copy()
_gauss_w_arr[n_official_hourly:] = 0.0

# Precomputed per-obs arrays (slot id, weights) used by solve_scales / the objective.
# _obs_slot_id[i] = slot index in slot_list (0..n_slots-1) if slotted, n_slots if not.
_obs_slot_id = np.full(n_obs, n_slots, dtype=np.int32)
for _si, (_, _sidxs) in enumerate(slot_list):
    _obs_slot_id[np.array(_sidxs, dtype=np.int64)] = _si

# Slot IDs and counts for slotted walking obs (constant; used in per-slot bincounts).
_walk_sl_sid = _obs_slot_id[_walk_slotted]
_walk_sl_n   = _walk_n_arr[_walk_slotted]

_slot_mfr = np.array([e[5]  for e in _slot_data], dtype=np.float64)
_slot_ivr = np.array([e[6]  for e in _slot_data], dtype=np.float64)
_slot_mfb = np.array([e[7]  for e in _slot_data], dtype=np.float64)
_slot_ivb = np.array([e[8]  for e in _slot_data], dtype=np.float64)
_slot_mfs = np.array([e[9]  for e in _slot_data], dtype=np.float64)
_slot_ivs = np.array([e[10] for e in _slot_data], dtype=np.float64)
_slot_mfa = np.array([e[11] for e in _slot_data], dtype=np.float64)
_slot_gam = np.array([e[12] for e in _slot_data], dtype=np.float64)

# ── Frozen temporal fractions (NTS profile) ───────────────────────────────────
# f is pinned at the NTS profile (mean_fraction_res/biz/school) and never tuned.
# Build the constant per-obs fraction vectors ONCE (the per-eval f-steps are gone).
# Unslotted obs map to the sentinel slot index n_slots → fraction 0; they carry
# zero objective weight anyway (_gauss_w_arr / _walk_slotted are 0 there).
_f_r_by_slot = np.append(_slot_mfr, 0.0)
_f_b_by_slot = np.append(_slot_mfb, 0.0)
_f_s_by_slot = np.append(_slot_mfs, 0.0)
_obs_f_r = _f_r_by_slot[_obs_slot_id]
_obs_f_b = _f_b_by_slot[_obs_slot_id]
_obs_f_s = _f_s_by_slot[_obs_slot_id]

# Constant NTS slot-fraction dicts, returned by solve_scales for persistence to
# tuned_params.json / tuning_history.jsonl (downstream consumers unchanged).
slot_fracs_res_nts    = {sk: float(_slot_mfr[si]) for si, (sk, _) in enumerate(slot_list)}
slot_fracs_biz_nts    = {sk: float(_slot_mfb[si]) for si, (sk, _) in enumerate(slot_list)}
slot_fracs_school_nts = {sk: float(_slot_mfs[si]) for si, (sk, _) in enumerate(slot_list)}

# ── Assignment and chi-squared helpers ───────────────────────────────────────

def run_assignment(W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz,
                   P_school, ALPHA_school,
                   w_pop, w_biz, THETA=None):
    """Production-constrained gravity assignment → (flow_res, flow_biz, flow_school).

    Each component is singly (production) constrained: T^c_ij = K_c·p^c_i·a^c_j·F_c/D^c_i
    (see model.constrained_od_flows). This returns the PRE-K per-link flows; the
    K_res/K_biz/K_sch scaling is applied analytically by solve_scales (D_i has no K,
    so flow stays linear in K — with f pinned at NTS the scale solve is convex).

    flow_res    = pop→pop,              kernel (P, ALPHA, BETA).
    flow_biz    = symmetric pop↔biz split (per-origin normalised) + W_BIZ·(biz×biz),
                  kernel (P_biz, ALPHA_biz, BETA).
    flow_school = symmetric pop↔school split (per-origin normalised), kernel
                  (P_school, ALPHA_school, BETA). Magnitude is K_sch (no separate
                  W_SCHOOL scale — it would be redundant with K_sch under the constraint).

    THETA is accepted but ignored (probit cache; the legacy k=3 logit scatter is retired).
    """
    t_res, t_biz, t_sch = constrained_od_flows(
        od_src, od_dst, od_dist, N_nodes, w_pop, w_biz, base_w_school,
        W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz,
        P_school=P_school, ALPHA_school=ALPHA_school, with_school=_has_school,
        self_src=_self_src, self_dist=_self_dist, self_w=_self_w,
        w_school_prod=base_w_school_prod)
    # Scatter only onto observed links (compact space) — see the observed-link
    # scatter-restriction block above. flow_* are indexed by COMPACT link id
    # (0.._N_obs_links-1); model_obs_3c reads them via the compact remap.
    flow_res = scatter_od_to_links(t_res, _pair_idx_obs, _link_idx_obs, _link_weight_obs, _N_obs_links)
    flow_biz = scatter_od_to_links(t_biz, _pair_idx_obs, _link_idx_obs, _link_weight_obs, _N_obs_links)
    if _has_school:
        flow_school = scatter_od_to_links(t_sch, _pair_idx_obs, _link_idx_obs, _link_weight_obs, _N_obs_links)
    else:
        flow_school = np.zeros(_N_obs_links, dtype=np.float64)
    return flow_res, flow_biz, flow_school


def model_obs_3c(flow_res, flow_biz, flow_school):
    """Extract per-observation modelled flows for all three components."""
    m_r = np.empty(n_obs)
    m_b = np.empty(n_obs)
    m_s = np.empty(n_obs)
    # flow_* are indexed by compact observed-link id (see run_assignment).
    for i, idxs in enumerate(_obs_link_idxs_compact[:n_official_hourly]):
        m_r[i] = flow_res[idxs].sum()    if idxs else 0.0
        m_b[i] = flow_biz[idxs].sum()    if idxs else 0.0
        m_s[i] = flow_school[idxs].sum() if idxs else 0.0
    m_r[n_official_hourly:] = np.where(_walk_valid, flow_res[_walk_link_safe_compact],    0.0)
    m_b[n_official_hourly:] = np.where(_walk_valid, flow_biz[_walk_link_safe_compact],    0.0)
    m_s[n_official_hourly:] = np.where(_walk_valid, flow_school[_walk_link_safe_compact], 0.0)
    return m_r, m_b, m_s


def solve_scales(m_res, m_biz, m_school):
    """Direct convex solve for the three component scales (K_res, K_biz, K_sch).

    Temporal fractions f are pinned at the NTS profile (never tuned), so each
    observation prediction is LINEAR in the scales:

        pred_i = K_res·a_i + K_biz·b_i + K_sch·d_i
        a_i = m_res_i·Th_i·f_res[s_i]   (b, d analogously) — f_* constant (_obs_f_*).

    The objective is therefore CONVEX over K ≥ 0:
      • Gaussian WLS over the official-hourly obs (convex quadratic),
      • Poisson identity-link deviance 2·Σ(n·log(n/pred)+pred−n) over the slotted
        walking obs (convex in the mean pred>0),
      • a scale-share prior on φ_biz = K_biz/ΣK (and φ_sch) carried over from the
        old φ-steps — it breaks the K_biz×W_BIZ degeneracy by anchoring the
        business/school share (it regularises the inner K-solve only and is NOT
        part of the reported χ², exactly as the old φ-prior was not).

    Solved by damped (Levenberg) Newton with a backtracking line search on the
    full objective — monotone by construction, so there is no K-collapse and no
    best-iterate bookkeeping. Replaces the old 5-block alternating minimisation.

    Returns (K_res, K_biz, K_sch, slot_fracs_res, slot_fracs_biz, slot_fracs_school)
    — the slot_fracs are the constant NTS profile dicts (for persistence).
    """
    a = m_res    * obs_Th * _obs_f_r
    b = m_biz    * obs_Th * _obs_f_b
    d = m_school * obs_Th * _obs_f_s

    C  = np.column_stack((a, b, d))          # (n_obs, 3): pred = C @ K
    y  = obs_rhs_arr
    wg = _gauss_w_arr                        # Gaussian weights (walking zeroed)
    wm = _walk_slotted                       # slotted-walking (Poisson) mask
    Cw = C[wm]
    nw = _walk_sl_n
    nw_pos = np.maximum(nw, 1e-300)

    mu_b, iv_b = phi_prior, 1.0 / (phi_std ** 2)
    mu_s, iv_s = ((phi_school_prior, 1.0 / (phi_school_std ** 2))
                  if _has_school else (0.0, 0.0))
    active = np.array([True, True, _has_school])

    def objval(K):
        pred = C @ K
        r  = pred - y
        Lg = float((wg * r) @ r)
        pw = np.maximum(pred[wm], 1e-30)
        Lp = float(np.sum(2.0 * np.where(
            nw > 0, nw * np.log(nw_pos / pw) + (pw - nw), pw)))
        S  = K[0] + K[1] + K[2]
        Lpr = 0.0
        if S > 1e-300:
            Lpr = iv_b * (K[1] / S - mu_b) ** 2
            if _has_school:
                Lpr += iv_s * (K[2] / S - mu_s) ** 2
        return Lg + Lp + Lpr

    def grad_hess(K):
        pred = C @ K
        # Gaussian
        g = 2.0 * (C.T @ (wg * (pred - y)))
        H = 2.0 * ((C.T * wg) @ C)
        # Poisson (identity link) over slotted walking obs
        pw = np.maximum(pred[wm], 1e-30)
        g += Cw.T @ (2.0 * (1.0 - nw / pw))
        H += (Cw.T * (2.0 * nw / (pw * pw))) @ Cw
        # Scale-share prior (Gauss-Newton PSD surrogate for its Hessian)
        S = K[0] + K[1] + K[2]
        if S > 1e-300:
            S2 = S * S
            jb = np.array([-K[1] / S2, (S - K[1]) / S2, -K[1] / S2])
            g += 2.0 * iv_b * (K[1] / S - mu_b) * jb
            H += 2.0 * iv_b * np.outer(jb, jb)
            if _has_school:
                js = np.array([-K[2] / S2, -K[2] / S2, (S - K[2]) / S2])
                g += 2.0 * iv_s * (K[2] / S - mu_s) * js
                H += 2.0 * iv_s * np.outer(js, js)
        return g, H

    # ── Init: single global scale from a moment fit, split by the prior shares.
    # K spans many orders of magnitude as the kernel params change, so a data-driven
    # magnitude is essential to land in the right basin before Newton refines.
    c     = a + b + d
    den_g = float((wg * c) @ c)
    num_g = float((wg * c) @ y)
    if den_g > 0 and num_g > 0:
        s0 = num_g / den_g                   # Gaussian (official) least-squares scale
    else:
        cw_sum = float(c[wm].sum())
        s0 = (float(nw.sum()) / cw_sum) if cw_sum > 0 else 1.0   # Poisson moment
    s0   = max(s0, 1e-30)
    sh_s = mu_s if _has_school else 0.0
    K = np.array([s0 * max(1.0 - mu_b - sh_s, 1e-6), s0 * mu_b, s0 * sh_s])
    K = np.maximum(K, 1e-30)
    if not _has_school:
        K[2] = 0.0

    f_cur = objval(K)
    for _ in range(60):
        g, H = grad_hess(K)
        gi = g[active]
        Hi = H[np.ix_(active, active)]
        lam_lm = 1e-9 * (np.trace(Hi) / max(gi.size, 1) + 1.0)   # Levenberg damping
        try:
            step_i = np.linalg.solve(Hi + lam_lm * np.eye(gi.size), -gi)
        except np.linalg.LinAlgError:
            step_i = -gi
        step = np.zeros(3)
        step[active] = step_i
        # Backtracking line search on the FULL objective → guaranteed descent.
        t, improved = 1.0, False
        for _ls in range(40):
            Kn = np.maximum(K + t * step, 1e-30)
            if not _has_school:
                Kn[2] = 0.0
            fn = objval(Kn)
            if fn < f_cur - 1e-12 * abs(f_cur):
                K, f_cur, improved = Kn, fn, True
                break
            t *= 0.5
        if not improved:
            break
        if np.max(np.abs(t * step_i)) < 1e-10 * (np.max(np.abs(K[active])) + 1e-30):
            break

    return (float(K[0]), float(K[1]), float(K[2]),
            slot_fracs_res_nts, slot_fracs_biz_nts, slot_fracs_school_nts)


# ── Objective function ────────────────────────────────────────────────────────

eval_count = [0]
best       = {"chi2": float("inf"), "log_params": None,
              "K_res": 1.0, "K_biz": 1.0, "K_sch": 0.0}
t0         = time.time()


def objective(log_params, log_ref=None):
    log_params = np.clip(log_params, -100, 100)
    W_BIZ     = math.exp(log_params[0])
    P         = math.exp(log_params[1])
    ALPHA     = math.exp(log_params[2])
    BETA      = math.exp(log_params[3])
    P_biz     = math.exp(log_params[4])
    ALPHA_biz = math.exp(log_params[5])
    THETA     = math.exp(log_params[6]) if _has_stoch else None
    _sch_off  = 7 if _has_stoch else 6
    if _has_school:
        P_school     = math.exp(log_params[_sch_off])
        ALPHA_school = math.exp(log_params[_sch_off + 1])
    else:
        P_school = ALPHA_school = 1.0

    w_pop = base_w_pop
    w_biz = base_w_biz

    flow_res, flow_biz, flow_school = run_assignment(
        W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz,
        P_school, ALPHA_school,
        w_pop, w_biz, THETA)
    m_res, m_biz, m_school = model_obs_3c(flow_res, flow_biz, flow_school)
    K_res, K_biz, K_sch = solve_scales(m_res, m_biz, m_school)[:3]
    # f is pinned at the NTS profile, so the per-obs fraction vectors _obs_f_* are
    # module-level constants and the old f-prior/coupling penalty is identically
    # zero (f ≡ NTS prior and res+biz+sch = agg per slot); it is dropped.
    _pred = (K_res * m_res * obs_Th * _obs_f_r
             + K_biz * m_biz * obs_Th * _obs_f_b
             + K_sch * m_school * obs_Th * _obs_f_s)
    # Gaussian chi-squared for official hourly obs (_gauss_w_arr zeroes walking obs).
    _resid    = _pred - obs_rhs_arr
    chi2_data = float((_gauss_w_arr * _resid) @ _resid)
    # Poisson deviance 2*(n*log(n/pred) + pred - n) for slotted walking obs.
    # Always ≥ 0; minimum 0 at pred=n.  The "raw" form 2*(pred-n*log(pred)) omits
    # the saturated term 2*(n*log(n)-n) and goes negative for n > e ≈ 2.718.
    _pred_w   = np.maximum(_pred[_walk_slotted], 1e-30)
    _pos_n    = np.maximum(_walk_sl_n, 1e-300)
    _pois_dev = 2.0 * np.where(
        _walk_sl_n > 0,
        _walk_sl_n * np.log(_pos_n / _pred_w) + (_pred_w - _walk_sl_n),
        _pred_w)
    chi2 = chi2_data + float(_pois_dev.sum())

    if np.any(log_grav_lam > 0):
        chi2 += float(np.dot(log_grav_lam, (log_params[:n_gravity] - log_grav_ref) ** 2))

    eval_count[0] += 1
    if chi2 < best["chi2"]:
        best["chi2"]       = chi2
        best["log_params"] = log_params.copy()
        best["K_res"]      = K_res
        best["K_biz"]      = K_biz
        best["K_sch"]      = K_sch

    if eval_count[0] % 20 == 0:
        elapsed = time.time() - t0
        print(f"  {eval_count[0]:4d}  χ²={chi2:.2f}  χ²/N={chi2/n_obs:.3f}"
              f"  best={best['chi2']/n_obs:.3f}  ({elapsed:.0f}s)")

    return chi2


# ── Build initial parameter vector ────────────────────────────────────────────

# Gravity start: from tuned_params.json if available, else hardcoded defaults
grav_start = {"W_BIZ": 1.0, "P": 300.0, "ALPHA": 2.0, "BETA": 1.0,
              "P_biz": 600.0, "ALPHA_biz": 2.0, "THETA": 1.0,
              "P_school": 300.0, "ALPHA_school": 2.0}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        tp = json.load(f)
    for k in ("W_BIZ", "P", "ALPHA", "BETA", "P_biz", "ALPHA_biz", "THETA",
              "P_school", "ALPHA_school"):
        if k in tp:
            grav_start[k] = tp[k]
    print(f"Starting gravity params from {TUNED_PARAMS}")

# Clamp to a safe minimum before log-transform (guards against degenerate prior runs)
_LOG_MIN = 1e-4
_log_p0_vals = [
    math.log(max(grav_start["W_BIZ"],     _LOG_MIN)),
    math.log(max(grav_start["P"],         _LOG_MIN)),
    math.log(max(grav_start["ALPHA"],     _LOG_MIN)),
    math.log(max(grav_start["BETA"],      _LOG_MIN)),
    math.log(max(grav_start["P_biz"],     _LOG_MIN)),
    math.log(max(grav_start["ALPHA_biz"], _LOG_MIN)),
]
if _has_stoch:
    _log_p0_vals.append(math.log(max(grav_start["THETA"], _LOG_MIN)))
if _has_school:
    _log_p0_vals += [
        math.log(max(grav_start["P_school"],     _LOG_MIN)),
        math.log(max(grav_start["ALPHA_school"], _LOG_MIN)),
    ]
log_p0 = np.array(_log_p0_vals, dtype=np.float64)

# Guard against param-vector / index drift (loud failure, not a silent mismatch).
assert len(log_p0) == n_gravity == len(_grav_param_names) == len(log_grav_ref), (
    f"param vector length mismatch: log_p0={len(log_p0)} n_gravity={n_gravity} "
    f"names={len(_grav_param_names)} ref={len(log_grav_ref)}")

# Capture starting gravity params for history (before any optimization)
initial_gravity = {k: grav_start[k]
                   for k in ("W_BIZ", "P", "ALPHA", "BETA", "P_biz", "ALPHA_biz")
                   if k in grav_start}
if _has_stoch and "THETA" in grav_start:
    initial_gravity["THETA"] = grav_start["THETA"]
if _has_school:
    for k in ("P_school", "ALPHA_school"):
        initial_gravity[k] = grav_start[k]

log_ref = None

if _has_stoch:
    _grav_note = "  [W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, THETA]"
elif _has_school:
    _grav_note = "  [W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, P_school, ALPHA_school]"
else:
    _grav_note = "  [W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz]"
print(f"Gravity: {len(log_p0)} params{_grav_note}")

# ── Calibration convergence probe (CALIBRATE_PROBE=1) ─────────────────────────
# Diagnostic only: at the start params, run the direct-K solve and report the
# residual global scale λ that would further reduce the data χ² (λ≈1 ⇒ K is at its
# optimum — the convex solver has converged). No optimization, no writes.
if os.environ.get("CALIBRATE_PROBE"):
    print("\n=== CALIBRATE PROBE (no optimization, no writes) ===")
    _W=math.exp(log_p0[0]);_P=math.exp(log_p0[1]);_A=math.exp(log_p0[2]);_B=math.exp(log_p0[3])
    _Pb=math.exp(log_p0[4]);_Ab=math.exp(log_p0[5])
    _so=7 if _has_stoch else 6
    _Ps=math.exp(log_p0[_so])   if _has_school else 1.0
    _As=math.exp(log_p0[_so+1]) if _has_school else 1.0
    _fr,_fb,_fs=run_assignment(_W,_P,_A,_B,_Pb,_Ab,_Ps,_As,base_w_pop,base_w_biz,None)
    _mr,_mb,_ms=model_obs_3c(_fr,_fb,_fs)
    def _probe_eval(Kr,Kb,Ks,sfr,sfb,sfs):
        fr=np.array([sfr.get(sk,1/24) if sk else 1/24 for sk in obs_slot_keys])
        fb=np.array([sfb.get(sk,1/24) if sk else 1/24 for sk in obs_slot_keys])
        fs=np.array([sfs.get(sk,0.0) if sk else 0.0 for sk in obs_slot_keys])
        pred=Kr*_mr*obs_Th*fr+Kb*_mb*obs_Th*fb+Ks*_ms*obs_Th*fs
        def ch(lam):
            p=np.maximum(lam*pred,1e-30)
            c=float((_gauss_w_arr*(p-obs_rhs_arr))@(p-obs_rhs_arr))
            pw=p[_walk_slotted]; n=_walk_sl_n; pos=np.maximum(n,1e-300)
            c+=float(np.sum(2*np.where(n>0,n*np.log(pos/pw)+(pw-n),pw)))
            return c
        ls=np.logspace(-1,2,2000); cs=[ch(l) for l in ls]; lo=ls[int(np.argmin(cs))]
        return ch(1.0),lo,ch(lo)
    print(f"  params: W_BIZ={_W:.3f} P={_P:.1f} ALPHA={_A:.3f} BETA={_B:.3f} P_biz={_Pb:.1f} ALPHA_biz={_Ab:.3f}")
    Kr,Kb,Ks,sfr,sfb,sfs=solve_scales(_mr,_mb,_ms)
    c1,lo,clo=_probe_eval(Kr,Kb,Ks,sfr,sfb,sfs)
    print(f"  direct solve  K_res={Kr:.4e} K_biz={Kb:.4e} K_sch={Ks:.4e}"
          f"  data χ²/N={c1/n_obs:7.3f}  opt_λ={lo:6.3f}  χ²/N@optλ={clo/n_obs:7.3f}")
    sys.exit(0)

# ── Run optimization ──────────────────────────────────────────────────────────

_tol = 5e-5 if fast else 1e-5
print(f"\nRunning Powell's method (λ={lam}"
      + (f"  fast: ftol/xtol={_tol:.0e}" if fast else "") + ") …")
print(f"  {'eval':>4s}  χ²/N(curr)  χ²/N(best)  elapsed")

# Evaluate initial point
objective(log_p0, log_ref)

result = scipy.optimize.minimize(
    lambda p: objective(p, log_ref),
    log_p0,
    method="powell",
    options={"maxiter": 5000, "ftol": _tol, "xtol": _tol},
)

# Use best params seen (Powell may backtrack at convergence)
log_best = best["log_params"]

# ── Unpack best params ────────────────────────────────────────────────────────

W_BIZ     = math.exp(log_best[0])
P         = math.exp(log_best[1])
ALPHA     = math.exp(log_best[2])
BETA      = math.exp(log_best[3])
P_biz     = math.exp(log_best[4])
ALPHA_biz = math.exp(log_best[5])
THETA     = math.exp(log_best[6]) if _has_stoch else None
_sch_off  = 7 if _has_stoch else 6
if _has_school:
    P_school     = math.exp(log_best[_sch_off])
    ALPHA_school = math.exp(log_best[_sch_off + 1])
else:
    P_school = ALPHA_school = None

# Final evaluation for clean chi2, K_res, K_biz, K_sch, and slot_fracs (no L2 term)
w_pop_f = base_w_pop
w_biz_f = base_w_biz

_ps_final = 1.0 if not _has_school else P_school
_as_final = 1.0 if not _has_school else ALPHA_school
flow_res, flow_biz, flow_school = run_assignment(
    W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz,
    _ps_final, _as_final,
    w_pop_f, w_biz_f, THETA)
m_res, m_biz, m_school = model_obs_3c(flow_res, flow_biz, flow_school)
K_res, K_biz, K_sch, slot_fracs_res, slot_fracs_biz, slot_fracs_school = \
    solve_scales(m_res, m_biz, m_school)
K = K_res + K_biz + K_sch   # total scale (for display and backward compat)
# f pinned at the NTS profile → constant fraction vectors; the f-prior/coupling
# penalty is identically zero and is dropped (see objective()).
_pred = (K_res * m_res * obs_Th * _obs_f_r
         + K_biz * m_biz * obs_Th * _obs_f_b
         + K_sch * m_school * obs_Th * _obs_f_s)
_resid    = _pred - obs_rhs_arr
chi2_data = float((_gauss_w_arr * _resid) @ _resid)
_pred_w   = np.maximum(_pred[_walk_slotted], 1e-30)
_pos_n    = np.maximum(_walk_sl_n, 1e-300)
_pois_dev = 2.0 * np.where(
    _walk_sl_n > 0,
    _walk_sl_n * np.log(_pos_n / _pred_w) + (_pred_w - _walk_sl_n),
    _pred_w)
chi2 = chi2_data + float(_pois_dev.sum())
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
            f_s = slot_fracs_school.get(sk, slot_prior[sk][4])
            mod_eff[i] = (K_res * m_res[i] * Th * f_r
                          + K_biz * m_biz[i] * Th * f_b
                          + K_sch * m_school[i] * Th * f_s)
        else:
            mod_eff[i] = (K_res * m_res[i] + K_biz * m_biz[i] + K_sch * m_school[i]) * Th
    else:  # walking: show as vehicles/hour using slot fraction
        n_eff_i = obs_rhs_arr[i]
        if sk is not None and Th > 0 and n_eff_i > 0:
            f_r = slot_fracs_res.get(sk, slot_prior[sk][2])
            f_b = slot_fracs_biz.get(sk, slot_prior[sk][3])
            f_s = slot_fracs_school.get(sk, slot_prior[sk][4])
            m_r_i, m_b_i, m_s_i = m_res[i], m_biz[i], m_school[i]
            raw_wtd = m_r_i * f_r + m_b_i * f_b + m_s_i * f_s
            denom   = m_r_i + m_b_i + m_s_i
            f_eff   = raw_wtd / denom if denom > 0 else (f_r + f_b + f_s) / 3
            obs_eff[i] = n_eff_i / (Th * f_eff)
            sig_eff[i] = math.sqrt(n_eff_i) / (Th * f_eff)
            mod_eff[i] = K_res * m_r_i + K_biz * m_b_i + K_sch * m_s_i   # combined AADT
        else:
            obs_eff[i] = n_eff_i / max(Th, 1e-9)
            sig_eff[i] = math.sqrt(n_eff_i) / max(Th, 1e-9)
            mod_eff[i] = (K_res * m_res[i] + K_biz * m_biz[i] + K_sch * m_school[i])
resid = np.where(sig_eff > 0, (mod_eff - obs_eff) / sig_eff, 0.0)

elapsed = time.time() - t0
print(f"\nResult  ({eval_count[0]} evals, {elapsed:.0f}s)")
_theta_str = f"  THETA={THETA:.4f}" if THETA is not None else ""
K_tot = K_res + K_biz + K_sch
phi_biz_out = K_biz / K_tot if K_tot > 0 else 0.0
phi_sch_out = K_sch / K_tot if K_tot > 0 else 0.0
print(f"  K_res={K_res:.4e}  K_biz={K_biz:.4e}  K_sch={K_sch:.4e}  (K={K:.4e})")
print(f"  phi_biz={phi_biz_out:.3f}  phi_sch={phi_sch_out:.3f}")
print(f"  W_BIZ={W_BIZ:.4f}  P={P:.2f}s  ALPHA={ALPHA:.4f}  BETA={BETA:.4f}"
      f"  P_biz={P_biz:.2f}s  ALPHA_biz={ALPHA_biz:.4f}{_theta_str}")
if _has_school:
    print(f"  P_school={P_school:.2f}s  ALPHA_school={ALPHA_school:.4f}  (school magnitude = K_sch)")
print(f"  χ²={chi2:.2f}  χ²/N={chi2_per_n:.4f}  χ²/N_eff={chi2/n_eff:.3f}  (N={n_obs}, N_eff={n_eff})")
if prev_chi2_per_n is not None:
    delta = chi2_per_n - prev_chi2_per_n
    direction = "improvement" if delta < 0 else "regression"
    print(f"  vs previous ({prev_id}):  Δχ²/N={delta:+.4f}  ({direction})")

# External zone values are census-derived (fixed) — no city delta table needed.

# ── Per-slot fraction table ───────────────────────────────────────────────────

_DT_NAMES = {0: "Wkday", 1: "Sat", 2: "Sun"}
print(f"\n  Per-slot hourly fractions (res / biz / school vs NTS priors):")
print(f"  {'Type':<5}  {'Hr':>2}  {'PriorAgg':>9}  {'f_res':>9}  {'f_biz':>9}  {'f_sch':>9}"
      f"  {'Δres%':>6}  {'Δbiz%':>6}  {'Δsch%':>6}  N")
for sk in sorted(slot_fracs_res):
    dt, h                    = sk
    mfa, std_f, mfr, mfb, mfs = slot_prior[sk]
    f_r                      = slot_fracs_res[sk]
    f_b                      = slot_fracs_biz[sk]
    f_s                      = slot_fracs_school[sk]
    n_in_slot                = len(slot_groups.get(sk, []))
    dr  = 100.0 * (f_r - mfr) / mfr if mfr > 0 else 0.0
    db  = 100.0 * (f_b - mfb) / mfb if mfb > 0 else 0.0
    ds  = 100.0 * (f_s - mfs) / mfs if mfs > 0 else 0.0
    print(f"  {_DT_NAMES[dt]:<5}  {h:>2d}  {mfa:>9.6f}  {f_r:>9.6f}  {f_b:>9.6f}  {f_s:>9.6f}"
          f"  {dr:>+5.1f}%  {db:>+5.1f}%  {ds:>+5.1f}%  {n_in_slot}")

# ── Goodness-of-fit table (sorted by |z|) ────────────────────────────────────

fit_rows = []
for i_obs in range(len(observations)):
    sk       = obs_slot_keys[i_obs]
    time_str = format_slot_time(*sk) if sk is not None else ""
    meta     = obs_meta[i_obs]
    fit_rows.append((meta["label"], time_str, obs_eff[i_obs], sig_eff[i_obs],
                     mod_eff[i_obs], resid[i_obs], meta["link"]))

fit_rows.sort(key=lambda r: abs(r[5]), reverse=True)
print_chi2_table(fit_rows, chi2, len(fit_rows), n_eff=n_eff)

# ── Save tuned_params.json ────────────────────────────────────────────────────

tuned = {
    "kernel":    "rational",
    "K":         float(K),         # K_res + K_biz + K_sch (backward compat); full precision
    "K_res":     float(K_res),     # K spans many orders of magnitude — round(x,6) zeroes sub-µ values
    "K_biz":     float(K_biz),
    "K_sch":     float(K_sch),
    "W_BIZ":     round(W_BIZ, 6),
    "P":         round(P, 4),
    "ALPHA":     round(ALPHA, 6),
    "BETA":      round(BETA, 6),
    "P_biz":     round(P_biz, 4),
    "ALPHA_biz": round(ALPHA_biz, 6),
    **( {"THETA": round(THETA, 6)} if THETA is not None else {} ),
    **( {"P_school":     round(P_school,     4),
         "ALPHA_school": round(ALPHA_school, 6)} if _has_school else {} ),
    "slot_fracs_res":    {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_res.items()},
    "slot_fracs_biz":    {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_biz.items()},
    "slot_fracs_school": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_school.items()},
    "chi2":       round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "n_obs":      n_obs,
    "n_slots":    n_slots,
    "n_eff":      n_eff,
    "stage":      "gravity",
    # f is always pinned at the NTS profile now (never tuned); kept for clarity.
    "temporal_profile": "nts_pinned",
}

with open(TUNED_PARAMS, "w") as f:
    json.dump(tuned, f, indent=2)
print(f"\nSaved → {TUNED_PARAMS}")

# ── Gravity model kernel plot ─────────────────────────────────────────────────

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d_sec    = np.logspace(np.log10(30), np.log10(7200), 500)  # 30s – 120 min
    d_min    = d_sec / 60.0
    u_res    = d_sec / P
    u_biz    = d_sec / P_biz
    k_res_c  = K_res * (ALPHA + BETA)     * u_res**BETA / (ALPHA     + BETA * u_res**(ALPHA     + BETA))
    k_biz_c  = K_biz * (ALPHA_biz + BETA) * u_biz**BETA / (ALPHA_biz + BETA * u_biz**(ALPHA_biz + BETA))

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(d_min, k_res_c, linewidth=1.8, label=f"Residential  P={P/60:.1f}min  ALPHA={ALPHA:.3f}  K_res={K_res:.3e}")
    ax.plot(d_min, k_biz_c, linewidth=1.8, linestyle="--",
            label=f"Business     P_biz={P_biz/60:.1f}min  ALPHA_biz={ALPHA_biz:.3f}  K_biz={K_biz:.3e}")
    ax.axvline(P     / 60.0, color="C0", linestyle=":", alpha=0.6)
    ax.axvline(P_biz / 60.0, color="C1", linestyle=":", alpha=0.6)
    if _has_school and P_school is not None:
        u_sch   = d_sec / P_school
        k_sch_c = K_sch * (ALPHA_school + BETA) * u_sch**BETA / (ALPHA_school + BETA * u_sch**(ALPHA_school + BETA))
        ax.plot(d_min, k_sch_c, linewidth=1.8, linestyle=":",
                label=f"School       P_sch={P_school/60:.1f}min  ALPHA_sch={ALPHA_school:.3f}  K_sch={K_sch:.3e}")
        ax.axvline(P_school / 60.0, color="C2", linestyle=":", alpha=0.6)
    ax.set_xlabel("Travel time (minutes)")
    ax.set_ylabel("K_c · kernel(d; P_c, ALPHA_c, BETA)")
    ax.set_title(
        f"Gravity kernels (rational)  BETA={BETA:.3f}\n"
        f"χ²/N={chi2_per_n:.4f}  id={run_id}  git={git_hash}"
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
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
    "kernel":    "rational",
    "K":         float(K),
    "K_res":     float(K_res),
    "K_biz":     float(K_biz),
    "K_sch":     float(K_sch),
    "W_BIZ":     round(W_BIZ, 6),
    "P":         round(P, 4),
    "ALPHA":     round(ALPHA, 6),
    "BETA":      round(BETA, 6),
    "P_biz":     round(P_biz, 4),
    "ALPHA_biz": round(ALPHA_biz, 6),
    **( {"THETA": round(THETA, 6)} if THETA is not None else {} ),
    **( {"P_school":     round(P_school,     4),
         "ALPHA_school": round(ALPHA_school, 6)} if _has_school else {} ),
    "slot_fracs_res":    {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_res.items()},
    "slot_fracs_biz":    {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_biz.items()},
    "slot_fracs_school": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_school.items()},
}

history_entry = {
    "id":        run_id,
    "git_hash":  git_hash,
    "timestamp": datetime.now(timezone.utc).isoformat(),
    "stage":     "gravity",
    "n_evals":   eval_count[0],
    "n_obs":     n_obs,
    "n_slots":   n_slots,
    "n_eff":     n_eff,
    "n_params":  len(log_p0),
    "chi2":      round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "params":    params,
    "tuner_hyperparams": {
        # phi_* priors anchor the scale-share inside solve_scales (degeneracy break).
        # gamma_coupling_scale is no longer used (f is pinned, no per-slot coupling).
        "phi_biz_prior":        phi_prior,
        "phi_biz_std":          phi_std,
        "phi_school_prior":     phi_school_prior,
        "phi_school_std":       phi_school_std,
        "gravity_lambda":       grav_lam_raw,
        "lambda":               lam,
        "fast":                 fast,
        "temporal_profile":     "nts_pinned",
    },
    "initial_gravity": {k: round(v, 6) for k, v in initial_gravity.items()},
    "slot_prior": {
        f"{dt},{h}": [round(mfa, 8), round(std_f, 8), round(mfr, 8), round(mfb, 8), round(mfs, 8)]
        for (dt, h), (mfa, std_f, mfr, mfb, mfs) in slot_prior.items()
    },
    "observations": [
        {
            "kind":     observations[i_obs][0],
            "label":    obs_meta[i_obs]["label"],
            "time":     (format_slot_time(*obs_slot_keys[i_obs])
                         if obs_slot_keys[i_obs] is not None else None),
            "link":     obs_meta[i_obs]["link"],
            "observed": round(float(obs_eff[i_obs]), 1),
            "sigma":    round(float(sig_eff[i_obs]), 1),
            "model":    round(float(mod_eff[i_obs]), 1),
            "z":        round(float(resid[i_obs]), 3),
        }
        for i_obs in range(n_obs)
    ],
}
history_entry["objective"] = "poisson_deviance_walking"
history_entry["temporal_profile"] = "nts_pinned"
if note:
    history_entry["note"] = note

with open(HISTORY_FILE, "a") as f:
    f.write(json.dumps(history_entry) + "\n")
print(f"Appended → {HISTORY_FILE}")
