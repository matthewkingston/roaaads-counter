"""
Gravity model OD matrix assignment.

Requires the precomputed paths cache (simulation/newtownards_paths.npz, built by
build_paths.py). Re-run build_paths.py whenever the road network changes or
through_route_pairs in tuner_config.json changes.

Outputs newtownards_flows.json; run build_demographics.py --map-only afterwards to
refresh the map.

Usage:
  python3 simulation/build_assignment.py
"""

import json, time, os
import numpy as np
import osmnx as ox
from model import (COUNT_SITES, EXCLUDE_LINKS, PATHS_CACHE, WEIGHTS_FILE,
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS, OFFICIAL_HOURLY, SCHOOL_LEVELS,
                   constrained_od_flows, scatter_od_to_links,
                   load_self_terms, aadt_weights,
                   load_generation_rates, compute_generation_scales,
                   site_flow, compute_chi2, print_chi2_table,
                   assert_paths_cache_fresh,
                   willingness_keys, willingness_from_flat)

# ── Config ────────────────────────────────────────────────────────────────────

K      = 1.73   # global flow scale factor

OUT_DIR    = "simulation"
CONS_GRAPH = "simulation/newtownards_reduced.graphml"  # dead-end-reduced routing graph

# ── Require paths cache ───────────────────────────────────────────────────────

if not os.path.exists(PATHS_CACHE):
    print(f"ERROR: paths cache not found: {PATHS_CACHE}")
    print("  Run:  python3 simulation/build_paths.py")
    print("  (build time ~6 min; re-run whenever the road network or")
    print("   through_route_pairs in tuner_config.json changes)")
    raise SystemExit(1)

# ── Load node weights ─────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)

_pnid = lambda k: (int(k) if k.lstrip("-").isdigit() else k)
node_population        = {_pnid(k): v for k, v in weights["node_population"].items()}
node_workplace         = {_pnid(k): v for k, v in weights["node_workplace"].items()}
node_commute_attractor = {_pnid(k): v for k, v in weights.get("node_commute_attractor", {}).items()}
node_retail_spaces     = {_pnid(k): v for k, v in weights.get("node_retail_spaces", {}).items()}
node_commute_producers = {_pnid(k): v for k, v in weights.get("node_commute_producers", {}).items()}
# Per-level school attractor + producer layers.
node_school_demand    = {lvl: {_pnid(k): v for k, v in weights.get(f"node_school_demand_{lvl}", {}).items()}
                         for lvl in SCHOOL_LEVELS}
node_school_producers = {lvl: {_pnid(k): v for k, v in weights.get(f"node_school_producers_{lvl}", {}).items()}
                         for lvl in SCHOOL_LEVELS}

K_res          = None
K_commute      = None
K_retail       = None
K_sch          = None
K_school       = {lvl: 0.0 for lvl in SCHOOL_LEVELS}
willingness    = None   # {component: (w, τs, τl)} — the 6 double-exp kernels
slot_fracs_res     = {}
slot_fracs_commute = {}
slot_fracs_retail  = {}
slot_fracs_school  = {lvl: {} for lvl in SCHOOL_LEVELS}

if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        _tp = json.load(f)
    K             = _tp.get("K",             K)
    if all(k in _tp for k in willingness_keys()):
        willingness = willingness_from_flat(_tp)     # 6 double-exp kernels (natural units)
    if "K_res" in _tp and "K_commute" in _tp and "K_retail" in _tp:
        K_res     = _tp["K_res"]
        K_commute = _tp["K_commute"]
        K_retail  = _tp["K_retail"]
        K_sch     = _tp.get("K_sch", 0.0)
        K_school  = {lvl: _tp.get(f"K_{lvl}", 0.0) for lvl in SCHOOL_LEVELS}
        _parse_sf = lambda key: {tuple(int(x) for x in k.split(",")): v
                                 for k, v in _tp.get(key, {}).items()}
        slot_fracs_res     = _parse_sf("slot_fracs_res")
        slot_fracs_commute = _parse_sf("slot_fracs_commute")
        slot_fracs_retail  = _parse_sf("slot_fracs_retail")
        slot_fracs_school  = {lvl: _parse_sf(f"slot_fracs_school_{lvl}") for lvl in SCHOOL_LEVELS}
    print(f"  [tuned: stage={_tp.get('stage','?')}  χ²/N={_tp.get('chi2_per_n','?')}"
          + (f"  K_res={K_res:.3e}  K_commute={K_commute:.3e}  K_retail={K_retail:.3e}"
             f"  K_sch={K_sch:.3e}" if K_res is not None else "")
          + "]")

# ── Assignment ────────────────────────────────────────────────────────────────

print(f"Loading paths cache ({PATHS_CACHE}) …")
t0    = time.time()
cache = np.load(PATHS_CACHE, allow_pickle=True)
assert_paths_cache_fresh(cache)

node_ids_arr = cache["node_ids"]
od_src       = cache["od_src"]
od_dst       = cache["od_dst"]
od_dist      = cache["od_dist"].astype(np.float64)
pair_idx     = cache["pair_idx"]
link_idx     = cache["link_idx"]
link_u       = cache["link_u"]
link_v       = cache["link_v"]

# Probit stochastic loading (link_weight = fraction of passes routing through each entry).
if "link_weight" not in cache:
    raise SystemExit("paths cache is not in probit format — rebuild: "
                     "python3 simulation/build_paths.py")
link_weight = cache["link_weight"].astype(np.float64)
n_passes    = int(cache["probit_n_passes"])
cv          = float(cache["probit_cv"])
print(f"  Probit loading: {n_passes} passes  CV={cv:.2f}  "
      f"({len(link_idx):,} weighted entries for {len(od_src):,} OD pairs)")

w_pop       = np.array([node_population.get(nid, 0)        for nid in node_ids_arr], dtype=np.float64)
w_workplace = np.array([node_workplace.get(nid, 0)         for nid in node_ids_arr], dtype=np.float64)
# Commute attractor = car-only jobs (node_commute_attractor); w_workplace (all jobs) is kept
# for the map output only.
w_commute_attr = np.array([node_commute_attractor.get(nid, 0) for nid in node_ids_arr], dtype=np.float64)
w_retail    = np.array([node_retail_spaces.get(nid, 0)     for nid in node_ids_arr], dtype=np.float64)
w_commute_prod = np.array([node_commute_producers.get(nid, 0) for nid in node_ids_arr], dtype=np.float64)
if w_commute_prod.sum() == 0:
    w_commute_prod = None   # fall back to population producer (legacy weights)
# Per-level school attractor + producer arrays.
w_school_levels = {lvl: np.array([node_school_demand[lvl].get(nid, 0) for nid in node_ids_arr],
                                 dtype=np.float64) for lvl in SCHOOL_LEVELS}
w_school_prod_levels = {lvl: np.array([node_school_producers[lvl].get(nid, 0) for nid in node_ids_arr],
                                      dtype=np.float64) for lvl in SCHOOL_LEVELS}
for lvl in SCHOOL_LEVELS:
    if w_school_prod_levels[lvl].sum() == 0:
        w_school_prod_levels[lvl] = None    # fall back to population producer
print(f"  {len(node_ids_arr)} nodes  pop {w_pop.sum():,.0f}  workplace {w_workplace.sum():,.0f}"
      f"  commute_attr {w_commute_attr.sum():,.0f}  retail {w_retail.sum():,.0f}")

# Generation pinning: per-leg producer scales (NTS0409a vehicle-driver trips/day) so the
# tuned K_c are interpreted against ≈ 1.0.  Island anchors summed from the weights dict.
# Absent rates file ⇒ None ⇒ unpinned (legacy behaviour).
_gen_rates = load_generation_rates()
_GEN_SCALE = (compute_generation_scales(weights, _gen_rates, verbose=True)
              if _gen_rates is not None else None)

N_links  = len(link_u)
N_nodes  = len(node_ids_arr)
# External intra-zonal self-term (denominator-only; from build_intra_times.py).
self_src, self_dist, self_w = load_self_terms(list(node_ids_arr))
# Require the multi-component params (the six K's + the six double-exp willingness kernels).
if not (K_res is not None and K_commute is not None and K_retail is not None
        and willingness is not None):
    raise SystemExit("tuned_params.json lacks the multi-component double-exp params — "
                     "run reset_gravity_params.py (after sync_kernel_anchor.py) then tune_assignment.py")
# School levels active when their K>0 and demand exists (each level has its own kernel in `willingness`).
_active_school = [lvl for lvl in SCHOOL_LEVELS
                  if K_school.get(lvl, 0.0) > 0 and w_school_levels[lvl].sum() > 0]
_use_school = len(_active_school) > 0

# Production-constrained assignment (singly-constrained per component).
# Per-OD-pair pre-K component flows, then scatter onto links via the probit
# routing incidence and apply K_res/K_commute/K_retail/K_sch.
t_res, t_commute, t_retail, t_sch_by_level = constrained_od_flows(
    od_src, od_dst, od_dist, N_nodes,
    w_pop, w_commute_attr, w_retail,
    willingness,
    with_school=_use_school,
    w_school_levels=w_school_levels, w_school_prod_levels=w_school_prod_levels,
    self_src=self_src, self_dist=self_dist, self_w=self_w,
    w_commute_prod=w_commute_prod,
    gen_scale=_GEN_SCALE)
raw_res     = scatter_od_to_links(t_res,     pair_idx, link_idx, link_weight, N_links)
raw_commute = scatter_od_to_links(t_commute, pair_idx, link_idx, link_weight, N_links)
raw_retail  = scatter_od_to_links(t_retail,  pair_idx, link_idx, link_weight, N_links)
raw_sch     = {lvl: scatter_od_to_links(t_sch_by_level[lvl], pair_idx, link_idx, link_weight, N_links)
               for lvl in SCHOOL_LEVELS}
_nonzero = (raw_res + raw_commute + raw_retail + sum(raw_sch.values())) > 0
_mk = lambda raw, K_c: {(int(link_u[k]), int(link_v[k])): raw[k] * K_c
                        for k in range(N_links) if _nonzero[k]}
link_flow_res     = _mk(raw_res,     K_res)
link_flow_commute = _mk(raw_commute, K_commute)
link_flow_retail  = _mk(raw_retail,  K_retail)
link_flow_school  = {lvl: _mk(raw_sch[lvl], K_school[lvl]) for lvl in SCHOOL_LEVELS}   # {level: dict}
_all_keys = set(link_flow_res) | set(link_flow_commute) | set(link_flow_retail)
for lvl in SCHOOL_LEVELS:
    _all_keys |= set(link_flow_school[lvl])
def _sch_at(lnk): return sum(link_flow_school[lvl].get(lnk, 0.0) for lvl in SCHOOL_LEVELS)
link_flow = {lnk: (link_flow_res.get(lnk, 0.0) + link_flow_commute.get(lnk, 0.0)
                   + link_flow_retail.get(lnk, 0.0) + _sch_at(lnk))
             for lnk in _all_keys}

# ── True AADT (daily) link flows ──────────────────────────────────────────
# link_flow_* above are K_c·m_c — calibrated so K_c·m_c·f_c[slot] matches the
# HOURLY count, so they are NOT daily totals and must NOT be compared to AADT.
# The annual-average daily contribution is K_c·m_c·W_c (W_c from aadt_weights).
# These weighted dicts are what is REPORTED/SERIALISED as AADT; the unweighted
# link_flow_* still feed compute_chi2 (which applies f_c itself — do not double
# weight).  W_res+W_commute+W_retail+W_sch ≈ 1.
W_res, W_commute, W_retail, W_school = aadt_weights(
    slot_fracs_res, slot_fracs_commute, slot_fracs_retail, slot_fracs_school)
W_sch = sum(W_school.values())   # display total
aadt_res     = {lnk: v * W_res     for lnk, v in link_flow_res.items()}
aadt_commute = {lnk: v * W_commute for lnk, v in link_flow_commute.items()}
aadt_retail  = {lnk: v * W_retail  for lnk, v in link_flow_retail.items()}
aadt_school  = {lvl: {lnk: v * W_school[lvl] for lnk, v in link_flow_school[lvl].items()}
                for lvl in SCHOOL_LEVELS}
def _aadt_sch_at(lnk): return sum(aadt_school[lvl].get(lnk, 0.0) for lvl in SCHOOL_LEVELS)
aadt_combined = {lnk: (aadt_res.get(lnk, 0.0) + aadt_commute.get(lnk, 0.0)
                       + aadt_retail.get(lnk, 0.0) + _aadt_sch_at(lnk))
                 for lnk in _all_keys}
print(f"  AADT weights: W_res={W_res:.3f} W_commute={W_commute:.3f} W_retail={W_retail:.3f} "
      + "W_school[" + " ".join(f"{lvl[:4]}={W_school[lvl]:.3f}" for lvl in SCHOOL_LEVELS) + "]")

# Per-external-node trip totals (routed pairs only: through = transiting ext→ext).
# AADT-weighted per component so these are true daily trips.
_n_routed = int(cache.get("n_routed_pairs", len(od_src)))
_is_ext   = np.array([isinstance(nid, str) for nid in node_ids_arr])
_t_total  = (t_res * K_res * W_res + t_commute * K_commute * W_commute
             + t_retail * K_retail * W_retail)
for lvl in SCHOOL_LEVELS:
    _t_total = _t_total + t_sch_by_level[lvl] * K_school[lvl] * W_school[lvl]
_src_r    = od_src[:_n_routed]
_dst_r    = od_dst[:_n_routed]
_t_r      = _t_total[:_n_routed]
_dst_ext  = _is_ext[_dst_r]
ext_node_trips = {}
for _ei in (i for i, nid in enumerate(node_ids_arr) if isinstance(nid, str)):
    _m = _src_r == _ei
    if not _m.any():
        continue
    _t = _t_r[_m]
    _de = _dst_ext[_m]
    ext_node_trips[str(node_ids_arr[_ei])] = {
        "trips_through":  round(float(_t[ _de].sum()), 1),
        "trips_internal": round(float(_t[~_de].sum()), 1),
    }

print(f"  Assignment complete in {time.time()-t0:.2f}s  ({len(link_flow)} loaded links)")

# ── Street name lookup ────────────────────────────────────────────────────────

node_ids    = list(node_ids_arr)
G           = ox.load_graphml(CONS_GRAPH)
node_weight = {nid: float(wp + wk + wr)
               for nid, wp, wk, wr in zip(node_ids_arr, w_pop, w_workplace, w_retail)}

_link_name = {(int(u), int(v)): d["name"]
              for u, v, d in G.edges(data=True) if d.get("name")}


def _link_label(u, v):
    name = _link_name.get((u, v), "")
    return f"{u}→{v}  {name}" if name else f"{u}→{v}"

# ── Report ────────────────────────────────────────────────────────────────────

print(f"\nOfficial count sites  (K = {K:.4e}"
      f"  K_res={K_res:.3e}  K_commute={K_commute:.3e}  K_retail={K_retail:.3e}"
      + (f"  K_sch={K_sch:.3e}" if _use_school else "") + ")")
print(f"  {'Site':<45s}  {'Modelled':>9s}  {'Observed':>9s}  {'Ratio':>6s}")
# Compare TRUE AADT (component-weighted) to the observed AADT — not the unweighted
# K·m link_flow, which is ~1/ΣW ≈ 2.6× larger and is hourly-calibrated, not daily.
for s in COUNT_SITES:
    f = site_flow(aadt_combined, s)
    print(f"  {s['label']:<45s}  {f:>9,.0f}  {s['observed']:>9,}  {f/s['observed']:>6.2f}")

rows, chi2, n_obs, n_eff = compute_chi2(
    link_flow_res,
    label_fn=_link_label,
    link_aadt_file=LINK_AADT,
    exclude_links=EXCLUDE_LINKS,
    link_flow_commute_dict=link_flow_commute,
    link_flow_retail_dict=link_flow_retail,
    link_flow_school_dicts=link_flow_school if _use_school else None,
    slot_fracs_res=slot_fracs_res,
    slot_fracs_commute=slot_fracs_commute,
    slot_fracs_retail=slot_fracs_retail,
    slot_fracs_school_levels=slot_fracs_school if _use_school else None,
)
print_chi2_table(rows, chi2, n_obs, n_eff=n_eff)

# ── Serialise flows ───────────────────────────────────────────────────────────

flows_path = f"{OUT_DIR}/newtownards_flows.json"
# Serialise TRUE AADT (component-weighted) flows — consumed by build_map.py as AADT.
out = {
    "kernel": "modesub_double",
    "K": K,
    "flows": {f"{u},{v}": flow for (u, v), flow in aadt_combined.items()},
    "K_res": K_res, "K_commute": K_commute, "K_retail": K_retail,
    "aadt_weights": {"res": W_res, "commute": W_commute, "retail": W_retail,
                     **{f"school_{lvl}": W_school[lvl] for lvl in SCHOOL_LEVELS}},
    "flows_res":     {f"{u},{v}": flow for (u, v), flow in aadt_res.items()},
    "flows_commute": {f"{u},{v}": flow for (u, v), flow in aadt_commute.items()},
    "flows_retail":  {f"{u},{v}": flow for (u, v), flow in aadt_retail.items()},
    "ext_node_trips": ext_node_trips,
}
if _use_school:
    out["K_sch"] = K_sch
    for lvl in SCHOOL_LEVELS:
        out[f"K_{lvl}"] = K_school[lvl]
        out[f"flows_school_{lvl}"] = {f"{u},{v}": flow for (u, v), flow in aadt_school[lvl].items()}
with open(flows_path, "w") as f:
    json.dump(out, f)
_comp_str = "res/commute/retail/school" if _use_school else "res/commute/retail"
print(f"\nSaved {len(link_flow)} link flows → {flows_path}  (+ {_comp_str} components)")
print(f"Parameters: K={K}  double-exp willingness (w, τs, τl):")
for _c in ("res", "commute", "retail", "school_primary", "school_postprimary", "school_tertiary"):
    _w, _ts, _tl = willingness[_c]
    print(f"    {_c:20s} w={_w:.3f}  τs={_ts:.0f}s  τl={_tl:.0f}s")
