"""
Powell's-method tuning of gravity model parameters against combined chi-squared
objective across all traffic count observations.

Four-component model: residential (pop↔pop), commute (commuters↔workplace), retail
(pop↔retail), and school (students↔school) flows each carry their own temporal
profile, distance kernel, and global scale.  At each evaluation the four component
scales (K_res, K_commute, K_retail, K_sch) are calibrated directly by a convex
damped-Newton solve (solve_scales); the temporal fractions f are pinned at the NTS
profile and never tuned.

Observations:
  Official sites: hourly count obs from data/official_hourly.json (24 h × 3 day-types
    × 3 sites = 216 obs), with Gaussian error (sigma from between-day variance).
  Walking obs: per-session count obs from data/link_aadt.json, Poisson error (n_eff).
Both types are in count space, unified in _slot_data with per-obs weights and rhs.

Optimizer parameters are stored in unconstrained internal coords (log for the τ's, logit for w).

Tunes a fully-tunable DOUBLE-EXP willingness per component (18 params — {τs, τl, w} × 6 components
incl. 3 INDEPENDENT school levels; 9 with no school demand) for the mode-substitution × willingness
kernel f(c)=driveshare(equiv_miles(c))·[w·exp(−c/τs)+(1−w)·exp(−c/τl)] — the rise + speed are
shared/data-derived (model._modesub_kernel), only {τs, τl, w} are tuned, anchored to the TLD/n_Ire
divide (analysis/fit_kernel.py → tuner_config gravity_ref, seeded by sync_kernel_anchor.py) via a
light L2 pull.  Commute and retail are independent clones of the school component — each a symmetric two-leg,
per-origin-normalised split with NO weight parameter and NO self/cross term (the old
single business component's W_BIZ and biz×biz are gone).
Production-constrained per component: each origin's trip production is fixed by its
producing weight, independent of accessibility.
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
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS, SCHOOL_LEVELS,
                   GENERATION_RATES, MOBILISATION_FILE,
                   constrained_od_flows, scatter_od_to_links, load_self_terms,
                   load_generation_rates, compute_generation_scales,
                   print_chi2_table, assert_paths_cache_fresh,
                   format_slot_time, nice_official, willingness_from_flat)

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

# Temporal fractions f_res/f_commute/f_retail/f_school are pinned at the NTS profile
# (mean_fraction_*) and never tuned — every residual is purely spatial. The inner
# calibration therefore solves only the four component scales (K_res, K_commute,
# K_retail, K_sch) directly via a convex Newton solve (see solve_scales); no φ reparam,
# no f-steps, no alternation, no best-iterate bookkeeping.

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

# External intra-zonal self-term (denominator-only; per-component, from build_intra_times.py).
# Constant across evals — only the kernel F(self_dist) recomputes inside run_assignment.
_self_terms = load_self_terms(node_ids)

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
node_pop_full       = {_pnid(k): v for k, v in wdata["node_population"].items()}
# Commute attractor = car-only jobs (node_commute_attractor); the all-jobs node_workplace layer
# is no longer the commute attractor.
node_commute_attr_full = {_pnid(k): v for k, v in wdata.get("node_commute_attractor", {}).items()}
node_retail_full    = {_pnid(k): v for k, v in wdata.get("node_retail_spaces", {}).items()}
node_commute_prod_full = {_pnid(k): v for k, v in wdata.get("node_commute_producers", {}).items()}
# School attractor + producer are per-level (primary/post-primary/tertiary), each its own layer.
node_school_full = {lvl: {_pnid(k): v for k, v in wdata.get(f"node_school_demand_{lvl}", {}).items()}
                    for lvl in SCHOOL_LEVELS}
node_school_prod_full = {lvl: {_pnid(k): v for k, v in wdata.get(f"node_school_producers_{lvl}", {}).items()}
                         for lvl in SCHOOL_LEVELS}

# Precomputed base weight arrays (from census + OSM demand; external zones fixed).
# Clean separate layers: car-commute jobs (commute attractor, node_commute_attractor) and
# retail spaces (retail attractor) — the all-jobs node_workplace is no longer an attractor.
base_w_pop       = np.array([node_pop_full.get(nid, 0.0)       for nid in node_ids], dtype=np.float64)
base_w_commute_attr = np.array([node_commute_attr_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
base_w_retail    = np.array([node_retail_full.get(nid, 0.0)    for nid in node_ids], dtype=np.float64)
base_w_commute_prod = np.array([node_commute_prod_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)

# Per-level school attractor + producer arrays (each level fully independent).
base_w_school_levels = {lvl: np.array([node_school_full[lvl].get(nid, 0.0) for nid in node_ids],
                                      dtype=np.float64) for lvl in SCHOOL_LEVELS}
base_w_school_prod_levels = {lvl: np.array([node_school_prod_full[lvl].get(nid, 0.0) for nid in node_ids],
                                           dtype=np.float64) for lvl in SCHOOL_LEVELS}
_active_school = [lvl for lvl in SCHOOL_LEVELS if base_w_school_levels[lvl].sum() > 0]
_has_school = len(_active_school) > 0
if not _has_school:
    print("  Warning: no per-level node_school_demand in weights — school components disabled")

# Generation pinning: per-leg producer scales (NTS0409a vehicle-driver trips/day) so
# each K_c should land at ≈ 1.0.  Island per-capita anchors are summed from wdata (the
# external nodes tile the whole island).  Absent rates file ⇒ None ⇒ unpinned (legacy).
_gen_rates = load_generation_rates()
_GEN_SCALE = (compute_generation_scales(wdata, _gen_rates, verbose=True)
              if _gen_rates is not None else None)

# Doubly-constrained (Furness) components: read from tuned_params.json if present.  The listed
# components' legs are ALSO attraction-constrained (Σ_i T_ij ∝ attractor_j; see
# model.constrained_od_flows); absent/empty ⇒ singly (production) constrained everywhere.
# Residential is held singly by design and is not a valid entry.
_DBLC = None
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as _f:
        _dblc_raw = json.load(_f).get("doubly_constrained")
    if _dblc_raw:
        _valid_dblc = {"commute", "retail"} | {f"school_{lvl}" for lvl in SCHOOL_LEVELS}
        _bad = set(_dblc_raw) - _valid_dblc
        if _bad:
            raise SystemExit(f"tuned_params.json doubly_constrained has unknown components "
                             f"{sorted(_bad)} (valid: {sorted(_valid_dblc)})")
        _DBLC = set(_dblc_raw)
        print(f"  Doubly-constrained (Furness) components: {sorted(_DBLC)}")

# Approximate-balancing sweep budget for the doubly-constrained legs (see
# model.constrained_od_flows / project note): the first eval per leg converges to seed the
# cache; every later eval runs this many warm sweeps (production stays exact, attraction
# <1% at k≈10).  From tuner_config.json via tuned_params.json; default 12.  `_FURNESS_STATE`
# persists the balancing factors across evals (Powell objective runs single-process).
_FURNESS_SWEEPS = None
_FURNESS_STATE = {}
if _DBLC and os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as _f:
        _FURNESS_SWEEPS = json.load(_f).get("furness_max_sweeps", 12)
    print(f"  Furness approximate-balancing sweeps/eval (warm): {_FURNESS_SWEEPS}")

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
# ── K-normalisation prior: Λ_K = pinv(Cov_K), anchor 1 (PROPERLY-DERIVED) ─────────────
# Generation pinning + the island m_island rescale (model.load_generation_rates) put every
# producer weight in island-level vehicle-driver trips/day, so K_c ≈ 1 is the value the data
# predict.  The prior around anchor 1 is NOT the old flat diag(1/K_prior_std²); it is a
# Mahalanobis form (K−1)ᵀΛ_K(K−1) whose SHAPE is derived, separating the two epistemically
# distinct directions and using each data source only for what it measures well:
#     Cov_K = P⊥·Cov_boot·P⊥      SPLIT  — NTS cluster-bootstrap sampling cov (generation_rates
#                                          _meta.stat_cov), common mode PROJECTED OUT.
#           + σ_mob²·u·uᵀ         LEVEL  — island mobilisation spatial dispersion
#                                          (mobilisation.json sigma_mob); loose leash the counts
#                                          refine (K departing 1 ⇒ local-vs-island mobilisation).
#           + split_σ²·P⊥         FLOOR  — differential England→Ireland transfer uncertainty on
#                                          the split (the bootstrap CV is ~1.6%, far too confident
#                                          to transfer the purpose MIX at; this floors it) + PSD
#                                          backstop.  (+ optional mob_σ on the common mode.)
# u = 1/√n·1 (common mode), P⊥ = I−uuᵀ (differential subspace).  Everything is in LOG/relative
# space (Var(K_c)≈Var(log ρ̂_c)) but applied as a linear-K cov about anchor 1 — first-order
# identical, keeps the convex linear-K Newton solve.  Components follow solve_scales' comp_keys
# EXACTLY: res, commute, retail, then the ACTIVE school levels — so Λ_K drops straight in.
K_PRIOR_ANCHOR = 1.0
_kcomps = ["res", "commute", "retail"] + [f"school_{lvl}" for lvl in _active_school]
_gen_meta = json.load(open(GENERATION_RATES)).get("_meta", {}) if os.path.exists(GENERATION_RATES) else {}
_stat_cov_k = _gen_meta.get("stat_cov")
_mob = json.load(open(MOBILISATION_FILE)) if os.path.exists(MOBILISATION_FILE) else {}
_k_floor = config.get("K_anchor_floor", {})
_split_sigma = float(_k_floor.get("split_sigma", 0.10))     # differential-transfer floor (relative)
_mob_floor   = float(_k_floor.get("mob_sigma", 0.0))        # extra common-mode floor (σ_mob usually owns it)
_sigma_mob = math.hypot(float(_mob.get("sigma_mob", 0.0)), _mob_floor)
if _sigma_mob < 1e-6:                                        # no mobilisation.json ⇒ level unknown, keep loose
    _sigma_mob = float(_k_floor.get("mob_default", 0.5))
    print(f"  [K-prior] no sigma_mob (run derive_mobilisation.py) — common mode loose at {_sigma_mob}")

def _build_lambda_k(kcomps, stat_cov_meta, sigma_mob, split_sigma):
    n = len(kcomps)
    Cb = np.zeros((n, n))                                    # bootstrap sampling cov, subselected
    if stat_cov_meta:
        order = stat_cov_meta["components"]
        full = np.asarray(stat_cov_meta["cov"], float)
        ix = [order.index(c) for c in kcomps]
        Cb = full[np.ix_(ix, ix)]
    u = np.ones(n) / math.sqrt(n)
    P = np.eye(n) - np.outer(u, u)                           # differential-subspace projector
    Cov = (P @ Cb @ P) + (sigma_mob ** 2) * np.outer(u, u) + (split_sigma ** 2) * P
    L = np.linalg.pinv(Cov)
    return 0.5 * (L + L.T)                                   # symmetrize (numerical)

Lambda_K = _build_lambda_k(_kcomps, _stat_cov_k, _sigma_mob, _split_sigma)
print(f"  K-prior: {_kcomps} | σ_mob={_sigma_mob:.3f} split_floor={_split_sigma:.3f}"
      f" | eff K-prior λ=diag={np.round(np.diag(Lambda_K), 1).tolist()}"
      + ("" if _stat_cov_k else "  [WARNING: no stat_cov in generation_rates.json → split = floor only;"
                                " run derive_generation_rates.py]"))

grav_ref = config.get("gravity_ref", {})
# Gravity params: a fully-tunable DOUBLE-EXP willingness per component for the mode-substitution
# × willingness kernel f(c)=driveshare(equiv_miles(c))·[w·exp(−c/τs)+(1−w)·exp(−c/τl)].  The rise +
# speed are shared/data-derived (model._modesub_kernel); the 3 willingness params {τs, τl, w} per
# component are tuned — 6 components incl. 3 INDEPENDENT school levels ⇒ 18 params (9 with no school
# demand), ALL free (no hard pins).  Anchored to the constrained kernel fit (analysis/iterate_kernel.py
# → tuner_config gravity_ref + gravity_stat_cov, seeded by sync_kernel_anchor.py) and regularized by a
# per-component Gaussian prior (block precision Λ_c, assembled below).
_w_components = ["res", "commute", "retail"]
if _has_school:
    _w_components += [f"school_{lvl}" for lvl in SCHOOL_LEVELS]
_all_wkeys = [f"{c}_{s}" for c in _w_components for s in ("taus", "taul", "w")]
_grav_param_names = list(_all_wkeys)          # ALL free (gravity_fixed retired) — the Powell vector

# Sane fallbacks for any key missing from gravity_ref (the anchor normally supplies all 18).
_WDEFAULTS = {"taus": 600.0, "taul": 3000.0, "w": 0.9}
def _wdefault(key):
    return _WDEFAULTS[key.rsplit("_", 1)[1]]

# Per-key transform between natural units and the unconstrained internal coords Powell optimizes:
# log for the two τ's (positivity), logit for w (∈(0,1)).  τl≥τs is enforced in
# model.willingness_from_flat (a clamp), not by coupling coords.
def _to_internal(key, val):
    if key.endswith("_w"):
        v = min(max(float(val), 1e-6), 1.0 - 1e-6)
        return math.log(v / (1.0 - v))
    return math.log(max(float(val), 1e-4))

def _from_internal(key, x):
    if key.endswith("_w"):
        return 1.0 / (1.0 + math.exp(-max(min(x, 100.0), -100.0)))
    return math.exp(x)

# Anchor in internal coords over all willingness keys.
log_grav_ref = np.array([_to_internal(k, grav_ref.get(k, _wdefault(k))) for k in _grav_param_names])
n_gravity = len(_grav_param_names)

# ── Anchor prior: per-component precision block Λ_c = (Cov_stat + Cov_epi)⁻¹ ───────────────────────
# Replaces the old scalar/per-key gravity_lambda.  Cov_stat (gravity_stat_cov — internal coords
# log τs/log τl/logit w, from the fit's Jacobian or trace-spread) is the DERIVED width; Cov_epi is the
# owned head-tight/tail-loose floor diag(head_σ², tail_σ², tail_σ²) — the two anchor_floor knobs.
# Adding covariances is a PSD floor (width ≥ floor in every direction); the w↔τl fit correlation makes
# each Λ_c a full 3×3, so the penalty is a per-component quadratic form (block-diagonal Λ_full).  Net:
# head floor-dominated (σ_stat tiny ⇒ λ≈30, firm); tail data-dominated (σ_stat large ⇒ looser, floor a
# backstop against the τl→∞ / w→1 runaway).  _grav_param_names is component-major (3 contiguous keys
# per component in _w_components order), so component i occupies slice [3i:3i+3].
_stat_cov = config.get("gravity_stat_cov", {})
_floor    = config.get("anchor_floor", {"head_sigma": 0.182, "tail_sigma": 0.693})
_cov_epi  = np.diag([_floor["head_sigma"] ** 2, _floor["tail_sigma"] ** 2, _floor["tail_sigma"] ** 2])
Lambda_full = np.zeros((n_gravity, n_gravity))
for _i, _c in enumerate(_w_components):
    _cs  = np.asarray(_stat_cov.get(_c, np.zeros((3, 3))), float)          # internal (τs, τl, w)
    _Lc  = np.linalg.pinv(_cs + _cov_epi)
    Lambda_full[3 * _i:3 * _i + 3, 3 * _i:3 * _i + 3] = 0.5 * (_Lc + _Lc.T)  # symmetrize (numerical)

print(f"  gravity: {n_gravity} free willingness params (all free), block anchor prior "
      f"(head_σ={_floor['head_sigma']}, tail_σ={_floor['tail_sigma']}); "
      f"eff. per-key λ=diag={np.round(np.diag(Lambda_full), 1).tolist()}")

# External node weights come from node_weights.json (census data, fixed — not tuned).

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

_raw_fracs = {}  # {(dow, hour): (mean_f, std_f, mf_res, mf_commute, mf_retail, mf_school)}
with open(HOURLY_FRACS_FILE, newline="") as _fh:
    for _row in csv.DictReader(_fh):
        _dow  = int(_row["day_of_week"])
        _hour = int(_row["hour"].split(":")[0])
        _mfr  = float(_row["mean_fraction_res"])      if "mean_fraction_res"      in _row else None
        _mfc  = float(_row["mean_fraction_commute"])  if "mean_fraction_commute"  in _row else None
        _mft  = float(_row["mean_fraction_retail"])   if "mean_fraction_retail"   in _row else None
        _mfs  = {lvl: (float(_row[f"mean_fraction_school_{lvl}"])
                       if f"mean_fraction_school_{lvl}" in _row else None) for lvl in SCHOOL_LEVELS}
        _raw_fracs[(_dow, _hour)] = (float(_row["mean_fraction"]), float(_row["std_fraction"]),
                                     _mfr, _mfc, _mft, _mfs)

# slot_prior[key] = (mean_f_agg, std_f_agg, mf_res, mf_commute, mf_retail, mf_school)
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
        _mfc = (sum(e[3] for e in _entries) / len(_entries)
                if _entries[0][3] is not None else _mf * 0.25)
        _mft = (sum(e[4] for e in _entries) / len(_entries)
                if _entries[0][4] is not None else _mf * 0.10)
        _mfs = {lvl: (sum(e[5][lvl] for e in _entries) / len(_entries)
                      if _entries[0][5][lvl] is not None else 0.0) for lvl in SCHOOL_LEVELS}
        slot_prior[(_dt, _h)] = (_mf, _std, _mfr, _mfc, _mft, _mfs)

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
# model_obs_4c reads modelled flow on only the OBSERVED links (the official-site +
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
_kept_entries  = int(_scatter_keep.sum())
_total_entries = len(link_idx_arr)
print(f"  Observed-link scatter: {_N_obs_links} links carry "
      f"{_kept_entries:,}/{_total_entries:,} incidence entries "
      f"({100*_kept_entries/_total_entries:.0f}%) — scatter restricted to these")

# Free the full-incidence arrays. They are only needed (a) above, to build the
# observed-link subset, and (b) in the _has_stoch CSR path — which is dead for the
# probit cache (_has_stoch is False).  The eval loop touches only the compact *_obs
# arrays, so on a memory-constrained host these full 66.8M-entry arrays (~800 MB) plus
# the per-entry compact/keep masks (~600 MB) would otherwise sit resident for the whole
# optimization and tip the working set into swap (peak RSS ~2.3 GB → ~1.5 GB).
del _compact_link_idx, _scatter_keep, pair_idx, link_idx_arr
if _link_weight is not None:
    del _link_weight
import gc as _gc; _gc.collect()

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

# Per-slot data: builds the constant NTS f vectors (below) and the chi² weight arrays.
# f is pinned at the NTS profile (never tuned), so only the per-component mean
# fractions are needed (the old inv_var/coupling columns of the f-step machinery are
# gone).  Each entry: (slot_key, ia, weights, rhs, Ths,
#                      mf_res, mf_commute, mf_retail, mf_school, mf_agg)
_slot_data = []
for sk, idxs in slot_list:
    ia       = np.array(idxs, dtype=np.int64)
    mfa, std_f, mfr, mfc, mft, mfs = slot_prior[sk]
    _slot_data.append((
        sk, ia,
        obs_w_arr[ia],     # weights: 1/sigma² or 1/n_eff per obs
        obs_rhs_arr[ia],   # rhs: count or n_eff per obs
        obs_Th[ia],        # T/3600 per obs
        mfr, mfc, mft, mfs, mfa,
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

_slot_mfr = np.array([e[5] for e in _slot_data], dtype=np.float64)
_slot_mfc = np.array([e[6] for e in _slot_data], dtype=np.float64)
_slot_mft = np.array([e[7] for e in _slot_data], dtype=np.float64)
_slot_mfs = {lvl: np.array([e[8][lvl] for e in _slot_data], dtype=np.float64) for lvl in SCHOOL_LEVELS}
_slot_mfa = np.array([e[9] for e in _slot_data], dtype=np.float64)

# ── Frozen temporal fractions (NTS profile) ───────────────────────────────────
# f is pinned at the NTS profile (mean_fraction_res/commute/retail/school) and never
# tuned.  Build the constant per-obs fraction vectors ONCE (the per-eval f-steps are
# gone).  Unslotted obs map to the sentinel slot index n_slots → fraction 0; they
# carry zero objective weight anyway (_gauss_w_arr / _walk_slotted are 0 there).
_f_r_by_slot = np.append(_slot_mfr, 0.0)
_f_c_by_slot = np.append(_slot_mfc, 0.0)
_f_t_by_slot = np.append(_slot_mft, 0.0)
_obs_f_r = _f_r_by_slot[_obs_slot_id]
_obs_f_c = _f_c_by_slot[_obs_slot_id]
_obs_f_t = _f_t_by_slot[_obs_slot_id]
# Per-level school per-obs fraction arrays (each level its own temporal shape).
_obs_f_s = {lvl: np.append(_slot_mfs[lvl], 0.0)[_obs_slot_id] for lvl in SCHOOL_LEVELS}

# Constant NTS slot-fraction dicts, returned by solve_scales for persistence to
# tuned_params.json / tuning_history.jsonl (downstream consumers unchanged).
slot_fracs_res_nts     = {sk: float(_slot_mfr[si]) for si, (sk, _) in enumerate(slot_list)}
slot_fracs_commute_nts = {sk: float(_slot_mfc[si]) for si, (sk, _) in enumerate(slot_list)}
slot_fracs_retail_nts  = {sk: float(_slot_mft[si]) for si, (sk, _) in enumerate(slot_list)}
slot_fracs_school_nts  = {lvl: {sk: float(_slot_mfs[lvl][si]) for si, (sk, _) in enumerate(slot_list)}
                          for lvl in SCHOOL_LEVELS}

# ── Assignment and chi-squared helpers ───────────────────────────────────────

def run_assignment(willingness):
    """Production-constrained gravity assignment.

    Returns (flow_res, flow_commute, flow_retail, flow_school).  Each component is
    singly (production) constrained: T^c_ij = K_c·p^c_i·a^c_j·F_c/D^c_i (see
    model.constrained_od_flows).  This returns the PRE-K per-link flows; the
    K_res/K_commute/K_retail/K_sch scaling is applied analytically by solve_scales
    (D_i has no K, so flow stays linear in K — with f pinned at NTS the scale solve
    is convex).

    `willingness` = {component: (w, τs, τl)} — each component's fully-independent
    double-exp willingness for the kernel f(c)=driveshare(equiv_miles(c))·
    [w·exp(−c/τs)+(1−w)·exp(−c/τl)]: flow_res = pop→pop; flow_commute =
    commuters↔workplace; flow_retail = pop↔retail; flow_school = students↔school
    (3 independent levels).  No weight parameter, no self/cross term.
    """
    t_res, t_commute, t_retail, t_sch_by_level = constrained_od_flows(
        od_src, od_dst, od_dist, N_nodes,
        base_w_pop, base_w_commute_attr, base_w_retail,
        willingness,
        with_school=_has_school,
        w_school_levels=base_w_school_levels, w_school_prod_levels=base_w_school_prod_levels,
        self_terms=_self_terms,
        w_commute_prod=base_w_commute_prod,
        gen_scale=_GEN_SCALE,
        doubly_constrained=_DBLC,
        furness_max_sweeps=_FURNESS_SWEEPS,
        furness_state=_FURNESS_STATE)
    # Scatter only onto observed links (compact space) — see the observed-link
    # scatter-restriction block above. flow_* are indexed by COMPACT link id
    # (0.._N_obs_links-1); model_obs_4c reads them via the compact remap.
    def _scatter(t):
        return scatter_od_to_links(t, _pair_idx_obs, _link_idx_obs, _link_weight_obs, _N_obs_links)
    flow_res     = _scatter(t_res)
    flow_commute = _scatter(t_commute)
    flow_retail  = _scatter(t_retail)
    flow_school  = {lvl: _scatter(t_sch_by_level[lvl]) for lvl in SCHOOL_LEVELS}   # zeros if inactive
    return flow_res, flow_commute, flow_retail, flow_school


def model_obs_4c(flow_res, flow_commute, flow_retail, flow_school):
    """Extract per-observation modelled flows. flow_school is {level: array}; returns
    (m_res, m_commute, m_retail, m_school) with m_school a {level: per-obs array} dict."""
    m_r = np.empty(n_obs); m_c = np.empty(n_obs); m_t = np.empty(n_obs)
    m_s = {lvl: np.empty(n_obs) for lvl in SCHOOL_LEVELS}
    # flow_* are indexed by compact observed-link id (see run_assignment).
    for i, idxs in enumerate(_obs_link_idxs_compact[:n_official_hourly]):
        m_r[i] = flow_res[idxs].sum()     if idxs else 0.0
        m_c[i] = flow_commute[idxs].sum() if idxs else 0.0
        m_t[i] = flow_retail[idxs].sum()  if idxs else 0.0
        for lvl in SCHOOL_LEVELS:
            m_s[lvl][i] = flow_school[lvl][idxs].sum() if idxs else 0.0
    m_r[n_official_hourly:] = np.where(_walk_valid, flow_res[_walk_link_safe_compact],     0.0)
    m_c[n_official_hourly:] = np.where(_walk_valid, flow_commute[_walk_link_safe_compact], 0.0)
    m_t[n_official_hourly:] = np.where(_walk_valid, flow_retail[_walk_link_safe_compact],  0.0)
    for lvl in SCHOOL_LEVELS:
        m_s[lvl][n_official_hourly:] = np.where(_walk_valid, flow_school[lvl][_walk_link_safe_compact], 0.0)
    return m_r, m_c, m_t, m_s


def solve_scales(m_res, m_commute, m_retail, m_school):
    """Direct convex solve for the four component scales
    (K_res, K_commute, K_retail, K_sch).

    Temporal fractions f are pinned at the NTS profile (never tuned), so each
    observation prediction is LINEAR in the scales:

        pred_i = K_res·a_i + K_commute·b_i + K_retail·c_i + K_sch·d_i
        a_i = m_res_i·Th_i·f_res[s_i]   (b, c, d analogously) — f_* constant (_obs_f_*).

    The objective is therefore CONVEX over K ≥ 0:
      • Gaussian WLS over the official-hourly obs (convex quadratic),
      • Poisson identity-link deviance 2·Σ(n·log(n/pred)+pred−n) over the slotted
        walking obs (convex in the mean pred>0),
      • a generation-anchored K-prior Σ_c iv_c·(K_c − 1)² over every active
        component — softly pulls each scale toward the generation value 1 (anchor
        + degeneracy break in one); regularises the inner K-solve only and is NOT
        part of the reported χ².

    Solved by damped (Levenberg) Newton with a backtracking line search on the
    full objective — monotone by construction, so there is no K-collapse and no
    best-iterate bookkeeping.

    Returns (K_res, K_commute, K_retail, K_sch, slot_fracs_res, slot_fracs_commute,
    slot_fracs_retail, slot_fracs_school) — the slot_fracs are the constant NTS
    profile dicts (for persistence).
    """
    # Component design columns: pred_i = Σ_comp K_comp·(m_comp_i·Th_i·f_comp[s_i]).  Fixed
    # components res/commute/retail, then one column per ACTIVE school level (each its own
    # temporal f and K-prior) — so this is a 3..6-component solve with no dead/masked columns.
    cols       = [m_res * obs_Th * _obs_f_r, m_commute * obs_Th * _obs_f_c, m_retail * obs_Th * _obs_f_t]
    comp_keys  = ["res", "commute", "retail"]
    for lvl in _active_school:
        cols.append(m_school[lvl] * obs_Th * _obs_f_s[lvl])
        comp_keys.append(f"school_{lvl}")

    C  = np.column_stack(cols)               # (n_obs, n_comp): pred = C @ K
    n_comp = C.shape[1]
    y  = obs_rhs_arr
    wg = _gauss_w_arr                        # Gaussian weights (walking zeroed)
    wm = _walk_slotted                       # slotted-walking (Poisson) mask
    Cw = C[wm]
    nw = _walk_sl_n
    nw_pos = np.maximum(nw, 1e-300)

    # Properly-derived K-prior (anchor = 1): penalty = (K−1)ᵀ Λ_K (K−1), where Λ_K (built at
    # setup over these exact comp_keys) is loose along the common mode (mobilisation, σ_mob),
    # tight on the split (bootstrap sampling cov + transfer floor).  Replaces the old diagonal
    # Σ_c (K_c−1)²/σ_c².  Λ_K is PSD ⇒ the inner objective stays convex.
    LK = Lambda_K
    anchor = K_PRIOR_ANCHOR

    def objval(K):
        pred = C @ K
        r  = pred - y
        Lg = float((wg * r) @ r)
        pw = np.maximum(pred[wm], 1e-30)
        Lp = float(np.sum(2.0 * np.where(
            nw > 0, nw * np.log(nw_pos / pw) + (pw - nw), pw)))
        dK = K - anchor
        Lpr = float(dK @ LK @ dK)
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
        # K-prior L = (K−anchor)ᵀ Λ_K (K−anchor) ⇒ g += 2·Λ_K·(K−anchor), H += 2·Λ_K
        # (full quadratic form; Λ_K PSD keeps H PSD).
        g += 2.0 * (LK @ (K - anchor))
        H += 2.0 * LK
        return g, H

    # ── Init: single global scale from a moment fit, applied uniformly to every
    # component (pred = s0·cc matches the moment fit; s0 ≈ 1 ⇒ each K starts at the
    # generation anchor when the model magnitude is calibrated).  K spans many orders
    # of magnitude as the kernel params change, so a data-driven magnitude is
    # essential to land in the right basin before Newton refines.
    cc    = C.sum(axis=1)
    den_g = float((wg * cc) @ cc)
    num_g = float((wg * cc) @ y)
    if den_g > 0 and num_g > 0:
        s0 = num_g / den_g                   # Gaussian (official) least-squares scale
    else:
        cw_sum = float(cc[wm].sum())
        s0 = (float(nw.sum()) / cw_sum) if cw_sum > 0 else 1.0   # Poisson moment
    s0   = max(s0, 1e-30)
    K = np.maximum(np.full(n_comp, s0), 1e-30)

    f_cur = objval(K)
    for _ in range(60):
        g, H = grad_hess(K)
        lam_lm = 1e-9 * (np.trace(H) / max(n_comp, 1) + 1.0)   # Levenberg damping
        try:
            step = np.linalg.solve(H + lam_lm * np.eye(n_comp), -g)
        except np.linalg.LinAlgError:
            step = -g
        # Backtracking line search on the FULL objective → guaranteed descent.
        t, improved = 1.0, False
        for _ls in range(40):
            Kn = np.maximum(K + t * step, 1e-30)
            fn = objval(Kn)
            if fn < f_cur - 1e-12 * abs(f_cur):
                K, f_cur, improved = Kn, fn, True
                break
            t *= 0.5
        if not improved:
            break
        if np.max(np.abs(t * step)) < 1e-10 * (np.max(np.abs(K)) + 1e-30):
            break

    K_map = {comp_keys[j]: float(K[j]) for j in range(n_comp)}
    K_school = {lvl: K_map.get(f"school_{lvl}", 0.0) for lvl in SCHOOL_LEVELS}
    return (K_map["res"], K_map["commute"], K_map["retail"], K_school,
            slot_fracs_res_nts, slot_fracs_commute_nts,
            slot_fracs_retail_nts, slot_fracs_school_nts)


# ── Objective function ────────────────────────────────────────────────────────

eval_count = [0]
best       = {"chi2": float("inf"), "log_params": None,
              "K_res": 1.0, "K_commute": 1.0, "K_retail": 1.0,
              "K_school": {lvl: 0.0 for lvl in SCHOOL_LEVELS}}
t0         = time.time()


def _unpack_gravity(log_params):
    """Unpack the internal param vector → the willingness dict {component: (w, τs, τl)}
    (all 18 params free; single source of index layout, shared by objective / probe / final eval;
    model.willingness_from_flat clamps τl≥τs)."""
    flat = {k: _from_internal(k, log_params[i]) for i, k in enumerate(_grav_param_names)}
    return willingness_from_flat(flat)


def objective(log_params, log_ref=None):
    log_params = np.clip(log_params, -100, 100)
    willingness = _unpack_gravity(log_params)

    flow_res, flow_commute, flow_retail, flow_school = run_assignment(willingness)
    m_res, m_commute, m_retail, m_school = model_obs_4c(
        flow_res, flow_commute, flow_retail, flow_school)
    K_res, K_commute, K_retail, K_school = solve_scales(
        m_res, m_commute, m_retail, m_school)[:4]
    # f is pinned at the NTS profile, so the per-obs fraction vectors _obs_f_* are
    # module-level constants and the old f-prior/coupling penalty is identically zero; dropped.
    _pred = (K_res * m_res * obs_Th * _obs_f_r
             + K_commute * m_commute * obs_Th * _obs_f_c
             + K_retail * m_retail * obs_Th * _obs_f_t)
    for lvl in SCHOOL_LEVELS:
        _pred = _pred + K_school[lvl] * m_school[lvl] * obs_Th * _obs_f_s[lvl]
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

    # Per-component Gaussian anchor prior: Σ_c Δφ_cᵀ Λ_c Δφ_c = Δφᵀ Λ_full Δφ (block-diagonal).
    _dphi = log_params[:n_gravity] - log_grav_ref
    chi2 += float(_dphi @ Lambda_full @ _dphi)

    eval_count[0] += 1
    if chi2 < best["chi2"]:
        best["chi2"]       = chi2
        best["log_params"] = log_params.copy()
        best["K_res"]      = K_res
        best["K_commute"]  = K_commute
        best["K_retail"]   = K_retail
        best["K_school"]   = K_school

    if eval_count[0] % 20 == 0:
        elapsed = time.time() - t0
        print(f"  {eval_count[0]:4d}  χ²={chi2:.2f}  χ²/N={chi2/n_obs:.3f}"
              f"  best={best['chi2']/n_obs:.3f}  ({elapsed:.0f}s)")

    return chi2


# ── Build initial parameter vector ────────────────────────────────────────────

# Gravity start: from tuned_params.json if available, else the tuner_config refs.
grav_start = dict(grav_ref)
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        tp = json.load(f)
    for k in _grav_param_names:
        if k in tp:
            grav_start[k] = tp[k]
    print(f"Starting gravity params from {TUNED_PARAMS}")

# Clamp to a safe minimum before log-transform; derive the start vector from the param
# name list (no parallel hardcoding — can't drift).
log_p0 = np.array([_to_internal(k, grav_start.get(k, _wdefault(k)))
                   for k in _grav_param_names], dtype=np.float64)

# Guard against param-vector / index drift (loud failure, not a silent mismatch).
assert len(log_p0) == n_gravity == len(_grav_param_names) == len(log_grav_ref), (
    f"param vector length mismatch: log_p0={len(log_p0)} n_gravity={n_gravity} "
    f"names={len(_grav_param_names)} ref={len(log_grav_ref)}")

# Capture starting gravity params for history (before any optimization)
initial_gravity = {k: grav_start[k] for k in _grav_param_names if k in grav_start}

log_ref = None

_grav_note = "  [" + ", ".join(_grav_param_names) + "]"
print(f"Gravity: {len(log_p0)} params{_grav_note}")

# ── Calibration convergence probe (CALIBRATE_PROBE=1) ─────────────────────────
# Diagnostic only: at the start params, run the direct-K solve and report the
# residual global scale λ that would further reduce the data χ² (λ≈1 ⇒ K is at its
# optimum — the convex solver has converged). No optimization, no writes.
if os.environ.get("CALIBRATE_PROBE"):
    print("\n=== CALIBRATE PROBE (no optimization, no writes) ===")
    _will = _unpack_gravity(log_p0)
    _fr,_fc,_ft,_fs = run_assignment(_will)
    _mr,_mc,_mt,_ms = model_obs_4c(_fr,_fc,_ft,_fs)   # _ms is {level: array}
    _fs_lvl = {lvl: np.array([slot_fracs_school_nts[lvl].get(sk,0.0) if sk else 0.0
                              for sk in obs_slot_keys]) for lvl in SCHOOL_LEVELS}
    def _probe_eval(Kr,Kc,Kt,K_school):
        fr=np.array([slot_fracs_res_nts.get(sk,1/24)     if sk else 1/24 for sk in obs_slot_keys])
        fc=np.array([slot_fracs_commute_nts.get(sk,1/24) if sk else 1/24 for sk in obs_slot_keys])
        ft=np.array([slot_fracs_retail_nts.get(sk,1/24)  if sk else 1/24 for sk in obs_slot_keys])
        pred = Kr*_mr*obs_Th*fr + Kc*_mc*obs_Th*fc + Kt*_mt*obs_Th*ft
        for lvl in SCHOOL_LEVELS:
            pred = pred + K_school[lvl]*_ms[lvl]*obs_Th*_fs_lvl[lvl]
        def ch(lam):
            p=np.maximum(lam*pred,1e-30)
            cc=float((_gauss_w_arr*(p-obs_rhs_arr))@(p-obs_rhs_arr))
            pw=p[_walk_slotted]; n=_walk_sl_n; pos=np.maximum(n,1e-300)
            cc+=float(np.sum(2*np.where(n>0,n*np.log(pos/pw)+(pw-n),pw)))
            return cc
        ls=np.logspace(-1,2,2000); cs=[ch(l) for l in ls]; lo=ls[int(np.argmin(cs))]
        return ch(1.0),lo,ch(lo)
    print("  double-exp willingness (w, τs, τl):")
    for _c in _w_components:
        _w, _ts, _tl = _will[_c]
        print(f"    {_c:20s} w={_w:.3f}  τs={_ts:.0f}s  τl={_tl:.0f}s")
    Kr,Kc,Kt,K_school = solve_scales(_mr,_mc,_mt,_ms)[:4]
    c1,lo,clo=_probe_eval(Kr,Kc,Kt,K_school)
    _ksch = " ".join(f"{lvl[:4]}={K_school[lvl]:.3e}" for lvl in SCHOOL_LEVELS)
    print(f"  direct solve  K_res={Kr:.4e} K_commute={Kc:.4e} K_retail={Kt:.4e}  K_school[{_ksch}]"
          f"  data χ²/N={c1/n_obs:7.3f}  opt_λ={lo:6.3f}  χ²/N@optλ={clo/n_obs:7.3f}")
    sys.exit(0)

# ══ DIAGNOSTIC — per-component kernel sweep (SWEEP=res|commute|retail|school) ══
# NOT part of tuning.  Env-gated; runs one assignment + convex K-solve per grid
# cell, does NO optimization and NO writes, then sys.exit(0).  Like CALIBRATE_PROBE.
#   Usage:  SWEEP=res python3 analysis/tune_assignment.py
#
# Varies ONLY the target component's tail scale τl over a sane grid while ALL other params
# (including that component's w and τs) stay FIXED at the start params, and reports that
# component's solved K and the resulting χ²/N per cell.  (τl is the diagnostic-relevant scale —
# the heavy tail the double-exp adds; sweep another sub-param by editing _sweep_key.)
#
# ⚠ READ WITH CARE — the K and χ²/N are CONDITIONAL on everything else being frozen at the start
# values.  They are NOT a joint best fit and must not be quoted as an achievable χ²/N.  The sweep
# answers one narrow question: does THIS component's K collapse, and which way does its preferred
# tail τl lean?  For a real fit, run the full tune (no SWEEP env var).
_sweep_target = os.environ.get("SWEEP")
if _sweep_target:
    if _sweep_target not in _w_components:
        print(f"SWEEP={_sweep_target!r} unknown; use one of {_w_components}"); sys.exit(1)
    _sweep_key = f"{_sweep_target}_taul"          # the sub-param swept (tail scale)
    def _willingness_with(key, val):
        """Full willingness with `key` overridden to `val` (all params free)."""
        flat = {k: _from_internal(k, log_p0[i]) for i, k in enumerate(_grav_param_names)}
        flat[key] = val
        return willingness_from_flat(flat)
    print(f"\n=== SWEEP {_sweep_target} τl (others fixed at start; no writes) ===")

    def _sweep_chi(K4, ms4):
        Kr, Kc, Kt, K_school = K4; mr, mc, mt, ms = ms4   # ms, K_school are {level: ...}
        fr = np.array([slot_fracs_res_nts.get(sk, 1/24)     if sk else 1/24 for sk in obs_slot_keys])
        fc = np.array([slot_fracs_commute_nts.get(sk, 1/24) if sk else 1/24 for sk in obs_slot_keys])
        ft = np.array([slot_fracs_retail_nts.get(sk, 1/24)  if sk else 1/24 for sk in obs_slot_keys])
        p = Kr*mr*obs_Th*fr + Kc*mc*obs_Th*fc + Kt*mt*obs_Th*ft
        for lvl in SCHOOL_LEVELS:
            fs = np.array([slot_fracs_school_nts[lvl].get(sk, 0.0) if sk else 0.0 for sk in obs_slot_keys])
            p = p + K_school[lvl]*ms[lvl]*obs_Th*fs
        p = np.maximum(p, 1e-30)
        cc = float((_gauss_w_arr*(p-obs_rhs_arr)) @ (p-obs_rhs_arr))
        pw = p[_walk_slotted]; n = _walk_sl_n; pos = np.maximum(n, 1e-300)
        cc += float(np.sum(2*np.where(n > 0, n*np.log(pos/pw)+(pw-n), pw)))
        return cc

    _idx3 = {"res": 0, "commute": 1, "retail": 2}
    def _k_total(K4):  return K4[0] + K4[1] + K4[2] + sum(K4[3].values())
    def _k_target(K4): return (K4[_idx3[_sweep_target]] if _sweep_target in _idx3
                               else K4[3][_sweep_target.split("school_")[1]])

    TAU_grid = [300, 600, 900, 1500, 2400, 3600, 5400, 7200, 10800]   # tail τl (s)
    _sw, _sts, _stl = _unpack_gravity(log_p0)[_sweep_target]
    print(f"  start (w={_sw:.3f}, τs={_sts:.0f}s, τl={_stl:.0f}s); sweeping τl")
    print(f"  {'τl_s':>6} {'τl_min':>7} {'K_'+_sweep_target:>13} {'phi':>8} {'chi2/N':>8}")
    for Tv in TAU_grid:
        ms4 = model_obs_4c(*run_assignment(_willingness_with(_sweep_key, Tv)))
        K4 = solve_scales(*ms4)[:4]
        tot = _k_total(K4); kv = _k_target(K4)
        phi = kv / tot if tot > 0 else 0.0
        print(f"  {Tv:6.0f} {Tv/60:7.1f} {kv:13.4e} {phi:8.3f} {_sweep_chi(K4, ms4)/n_obs:8.3f}")
    sys.exit(0)

# ── Run optimization ──────────────────────────────────────────────────────────

_tol = 5e-5 if fast else 1e-5
print(f"\nRunning Powell's method (λ={lam}"
      + (f"  fast: ftol/xtol={_tol:.0e}" if fast else "") + ") …")
print(f"  {'eval':>4s}  χ²/N(curr)  χ²/N(best)  elapsed")

# Evaluate initial point (populates `best` before Powell starts, so a Ctrl-C at any
# point below has a best-seen param vector to write out).
objective(log_p0, log_ref)

try:
    result = scipy.optimize.minimize(
        lambda p: objective(p, log_ref),
        log_p0,
        method="powell",
        options={"maxiter": 5000, "ftol": _tol, "xtol": _tol},
    )
except KeyboardInterrupt:
    # Graceful stop: fall through to the normal write path (final eval + tuned_params +
    # history + plot) using the best params seen so far.  `result` is never read again
    # (log_best comes from `best`), so a placeholder is fine.
    print(f"\n⚠ Interrupted (Ctrl-C) after {eval_count[0]} evals — writing the best-seen "
          f"params (χ²/N={best['chi2']/n_obs:.4f}) and exiting cleanly …")
    result = None

# Use best params seen (Powell may backtrack at convergence, or we were interrupted)
log_best = best["log_params"]

# ── Unpack best params ────────────────────────────────────────────────────────

willingness_best = _unpack_gravity(log_best)   # {component: (w, τs, τl)}

# Final evaluation for clean chi2, the four K's, and slot_fracs (no L2 term)
flow_res, flow_commute, flow_retail, flow_school = run_assignment(willingness_best)
m_res, m_commute, m_retail, m_school = model_obs_4c(
    flow_res, flow_commute, flow_retail, flow_school)
(K_res, K_commute, K_retail, K_school,
 slot_fracs_res, slot_fracs_commute, slot_fracs_retail, slot_fracs_school) = \
    solve_scales(m_res, m_commute, m_retail, m_school)
K_sch = sum(K_school.values())             # total school scale (display / compat)
K = K_res + K_commute + K_retail + K_sch   # total scale (display / backward compat)
# Per-obs combined-school helpers for the display/fit table (each level its own K, m, temporal f).
def _sch_fs(lvl, sk):
    return slot_fracs_school[lvl].get(sk, slot_prior[sk][5][lvl]) if sk is not None else 0.0
def _sch_hourly(i, sk): return sum(K_school[l] * m_school[l][i] * _sch_fs(l, sk) for l in SCHOOL_LEVELS)
def _sch_aadt(i):       return sum(K_school[l] * m_school[l][i] for l in SCHOOL_LEVELS)
def _sch_mf(i, sk):     return sum(m_school[l][i] * _sch_fs(l, sk) for l in SCHOOL_LEVELS)
def _sch_m(i):          return sum(m_school[l][i] for l in SCHOOL_LEVELS)
# f pinned at the NTS profile → constant fraction vectors; the f-prior/coupling
# penalty is identically zero and is dropped (see objective()).
_pred = (K_res * m_res * obs_Th * _obs_f_r
         + K_commute * m_commute * obs_Th * _obs_f_c
         + K_retail * m_retail * obs_Th * _obs_f_t)
for lvl in SCHOOL_LEVELS:
    _pred = _pred + K_school[lvl] * m_school[lvl] * obs_Th * _obs_f_s[lvl]
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
# The total modelled hourly count is Σ_c K_c·m_c·f_c over res/commute/retail/school.
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
            f_r = slot_fracs_res.get(sk,     slot_prior[sk][2])
            f_c = slot_fracs_commute.get(sk, slot_prior[sk][3])
            f_t = slot_fracs_retail.get(sk,  slot_prior[sk][4])
            mod_eff[i] = (K_res * m_res[i] * Th * f_r
                          + K_commute * m_commute[i] * Th * f_c
                          + K_retail * m_retail[i] * Th * f_t
                          + _sch_hourly(i, sk) * Th)
        else:
            mod_eff[i] = (K_res * m_res[i] + K_commute * m_commute[i]
                          + K_retail * m_retail[i] + _sch_aadt(i)) * Th
    else:  # walking: show as vehicles/hour using slot fraction
        n_eff_i = obs_rhs_arr[i]
        if sk is not None and Th > 0 and n_eff_i > 0:
            f_r = slot_fracs_res.get(sk,     slot_prior[sk][2])
            f_c = slot_fracs_commute.get(sk, slot_prior[sk][3])
            f_t = slot_fracs_retail.get(sk,  slot_prior[sk][4])
            m_r_i, m_c_i, m_t_i = m_res[i], m_commute[i], m_retail[i]
            raw_wtd = m_r_i * f_r + m_c_i * f_c + m_t_i * f_t + _sch_mf(i, sk)
            denom   = m_r_i + m_c_i + m_t_i + _sch_m(i)
            f_eff   = raw_wtd / denom if denom > 0 else (f_r + f_c + f_t) / 3
            obs_eff[i] = n_eff_i / (Th * f_eff)
            sig_eff[i] = math.sqrt(n_eff_i) / (Th * f_eff)
            mod_eff[i] = (K_res * m_r_i + K_commute * m_c_i
                          + K_retail * m_t_i + _sch_aadt(i))   # combined AADT
        else:
            obs_eff[i] = n_eff_i / max(Th, 1e-9)
            sig_eff[i] = math.sqrt(n_eff_i) / max(Th, 1e-9)
            mod_eff[i] = (K_res * m_res[i] + K_commute * m_commute[i]
                          + K_retail * m_retail[i] + _sch_aadt(i))
resid = np.where(sig_eff > 0, (mod_eff - obs_eff) / sig_eff, 0.0)

elapsed = time.time() - t0
print(f"\nResult  ({eval_count[0]} evals, {elapsed:.0f}s)")
K_tot = K_res + K_commute + K_retail + K_sch
phi_com_out = K_commute / K_tot if K_tot > 0 else 0.0
phi_ret_out = K_retail  / K_tot if K_tot > 0 else 0.0
phi_sch_out = K_sch     / K_tot if K_tot > 0 else 0.0
print(f"  K_res={K_res:.4e}  K_commute={K_commute:.4e}  K_retail={K_retail:.4e}"
      f"  K_sch={K_sch:.4e}  (K={K:.4e})")
print("  K_school: " + "  ".join(f"{lvl}={K_school[lvl]:.4e}" for lvl in SCHOOL_LEVELS))
print(f"  phi_commute={phi_com_out:.3f}  phi_retail={phi_ret_out:.3f}  phi_sch={phi_sch_out:.3f}")
print("  double-exp willingness (w, τs, τl):")
for _c in _w_components:
    _w, _ts, _tl = willingness_best[_c]
    print(f"    {_c:20s} w={_w:.3f}  τs={_ts:.0f}s  τl={_tl:.0f}s")
print(f"  χ²={chi2:.2f}  χ²/N={chi2_per_n:.4f}  χ²/N_eff={chi2/n_eff:.3f}  (N={n_obs}, N_eff={n_eff})")
if prev_chi2_per_n is not None:
    delta = chi2_per_n - prev_chi2_per_n
    direction = "improvement" if delta < 0 else "regression"
    print(f"  vs previous ({prev_id}):  Δχ²/N={delta:+.4f}  ({direction})")

# External zone values are census-derived (fixed) — no city delta table needed.

# ── Per-slot fraction table ───────────────────────────────────────────────────

_DT_NAMES = {0: "Wkday", 1: "Sat", 2: "Sun"}
print(f"\n  Per-slot hourly fractions (res / commute / retail / school[prim/postp/tert] vs NTS):")
print(f"  {'Type':<5}  {'Hr':>2}  {'PriorAgg':>9}  {'f_res':>9}  {'f_com':>9}  {'f_ret':>9}"
      f"  {'f_scP':>8}  {'f_scPP':>8}  {'f_scT':>8}  N")
for sk in sorted(slot_fracs_res):
    dt, h     = sk
    mfa       = slot_prior[sk][0]
    f_r       = slot_fracs_res[sk]; f_c = slot_fracs_commute[sk]; f_t = slot_fracs_retail[sk]
    fsc       = [slot_fracs_school[lvl].get(sk, 0.0) for lvl in SCHOOL_LEVELS]
    n_in_slot = len(slot_groups.get(sk, []))
    print(f"  {_DT_NAMES[dt]:<5}  {h:>2d}  {mfa:>9.6f}  {f_r:>9.6f}  {f_c:>9.6f}  {f_t:>9.6f}"
          f"  {fsc[0]:>8.5f}  {fsc[1]:>8.5f}  {fsc[2]:>8.5f}  {n_in_slot}")

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

# Flat willingness keys for persistence (natural units): "<comp>_taus/_taul/_w".
_will_flat = {}
for _c in _w_components:
    _w, _ts, _tl = willingness_best[_c]
    _will_flat[f"{_c}_taus"] = round(_ts, 4)
    _will_flat[f"{_c}_taul"] = round(_tl, 4)
    _will_flat[f"{_c}_w"]    = round(_w, 6)

tuned = {
    "kernel":    "modesub_double",
    "K":         float(K),         # K_res+K_commute+K_retail+K_sch (compat); full precision
    "K_res":     float(K_res),     # K spans many orders of magnitude — round(x,6) zeroes sub-µ values
    "K_commute": float(K_commute),
    "K_retail":  float(K_retail),
    "K_sch":     float(K_sch),                       # total (display/compat)
    **{f"K_{lvl}": float(K_school[lvl]) for lvl in SCHOOL_LEVELS},
    **_will_flat,
    "slot_fracs_res":     {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_res.items()},
    "slot_fracs_commute": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_commute.items()},
    "slot_fracs_retail":  {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_retail.items()},
    **{f"slot_fracs_school_{lvl}":
       {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_school[lvl].items()}
       for lvl in SCHOOL_LEVELS},
    "chi2":       round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "n_obs":      n_obs,
    "n_slots":    n_slots,
    "n_eff":      n_eff,
    "stage":      "gravity",
    # f is always pinned at the NTS profile now (never tuned); kept for clarity.
    "temporal_profile": "nts_pinned",
    # Carry the doubly-constrained config through so it survives a tune (build_assignment /
    # diagnose_imbalance read it from here); only written when active.
    **({"doubly_constrained": sorted(_DBLC), "furness_max_sweeps": _FURNESS_SWEEPS}
       if _DBLC else {}),
}

with open(TUNED_PARAMS, "w") as f:
    json.dump(tuned, f, indent=2)
print(f"\nSaved → {TUNED_PARAMS}")

# ── Gravity model kernel plot ─────────────────────────────────────────────────

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from model import _modesub_kernel
    d_sec    = np.logspace(np.log10(30), np.log10(7200), 500)  # 30s – 120 min
    d_min    = d_sec / 60.0
    def _kern(K_c, component):
        return K_c * _modesub_kernel(d_sec, willingness_best[component], component)
    _lbl = {c: f"τs={willingness_best[c][1]/60:.1f}m τl={willingness_best[c][2]/60:.1f}m "
               f"w={willingness_best[c][0]:.2f}" for c in _w_components}

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(d_min, _kern(K_res, "res"), linewidth=1.8,
            label=f"Residential  {_lbl['res']}  K={K_res:.2e}")
    ax.plot(d_min, _kern(K_commute, "commute"), linewidth=1.8, linestyle="--",
            label=f"Commute      {_lbl['commute']}  K={K_commute:.2e}")
    ax.plot(d_min, _kern(K_retail, "retail"), linewidth=1.8, linestyle="-.",
            label=f"Retail       {_lbl['retail']}  K={K_retail:.2e}")
    if _has_school:
        for lvl in SCHOOL_LEVELS:                       # each level its own driveshare + willingness
            if K_school.get(lvl, 0.0) <= 0:
                continue
            _c = f"school_{lvl}"
            ax.plot(d_min, _kern(K_school[lvl], _c),
                    linewidth=1.6, linestyle=":",
                    label=f"School {lvl}  {_lbl[_c]}  K={K_school[lvl]:.2e}")
    ax.set_xlabel("Travel time (minutes)")
    ax.set_ylabel("K_c · driveshare·[w·exp(−c/τs)+(1−w)·exp(−c/τl)]")
    ax.set_title(
        f"Gravity kernels (mode-substitution × double-exp willingness)\n"
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
    "kernel":    "modesub_double",
    "K":         float(K),
    "K_res":     float(K_res),
    "K_commute": float(K_commute),
    "K_retail":  float(K_retail),
    "K_sch":     float(K_sch),                       # total (display/compat)
    **{f"K_{lvl}": float(K_school[lvl]) for lvl in SCHOOL_LEVELS},
    **_will_flat,
    "slot_fracs_res":     {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_res.items()},
    "slot_fracs_commute": {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_commute.items()},
    "slot_fracs_retail":  {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_retail.items()},
    **{f"slot_fracs_school_{lvl}":
       {f"{dt},{h}": round(f, 8) for (dt, h), f in slot_fracs_school[lvl].items()}
       for lvl in SCHOOL_LEVELS},
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
        # Properly-derived K-prior (anchor 1) regularises solve_scales: Λ_K = pinv(
        # P⊥·Cov_boot·P⊥ + σ_mob²·uuᵀ + split_σ²·P⊥) — loose common mode (mobilisation),
        # tight derived split.  gamma_coupling_scale is no longer used (f is pinned).
        "K_prior_anchor":       K_PRIOR_ANCHOR,
        "K_prior_components":   _kcomps,
        "K_prior_sigma_mob":    _sigma_mob,
        "K_prior_split_sigma":  _split_sigma,
        "K_prior_m_island":     float(_mob.get("m_island", 1.0)),
        "K_prior_stat_cov_seed": (_stat_cov_k or {}).get("seed"),
        "K_prior_eff_lambda":   dict(zip(_kcomps, np.round(np.diag(Lambda_K), 4).tolist())),
        # Anchor prior: the two floor knobs + the effective per-key λ (diag of the block precision).
        "anchor_floor":         _floor,
        "gravity_eff_lambda":   dict(zip(_grav_param_names, np.round(np.diag(Lambda_full), 4).tolist())),
        "lambda":               lam,
        "fast":                 fast,
        "temporal_profile":     "nts_pinned",
    },
    "initial_gravity": {k: round(v, 6) for k, v in initial_gravity.items()},
    "slot_prior": {
        f"{dt},{h}": [round(mfa, 8), round(std_f, 8), round(mfr, 8), round(mfc, 8),
                      round(mft, 8), {lvl: round(mfs[lvl], 8) for lvl in SCHOOL_LEVELS}]
        for (dt, h), (mfa, std_f, mfr, mfc, mft, mfs) in slot_prior.items()
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
# Constraint mode, recorded UNCONDITIONALLY so every run's history is unambiguous
# (empty list ⇒ singly-constrained everywhere; non-empty ⇒ those components are also
# attraction-constrained via Furness). Mirrors the tuned_params.json fields, which are
# written only when active — history always records it so runs can be compared.
# _DBLC is None when singly (never set) ⇒ record [] / null.
history_entry["doubly_constrained"] = sorted(_DBLC) if _DBLC else []
history_entry["furness_max_sweeps"] = _FURNESS_SWEEPS if _DBLC else None
if note:
    history_entry["note"] = note

with open(HISTORY_FILE, "a") as f:
    f.write(json.dumps(history_entry) + "\n")
print(f"Appended → {HISTORY_FILE}")
