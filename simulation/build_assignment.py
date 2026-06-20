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
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS, OFFICIAL_HOURLY,
                   gravity_assign, site_flow, compute_chi2, print_chi2_table)

# ── Config ────────────────────────────────────────────────────────────────────

K      = 1.73   # global flow scale factor
W_BIZ  = 1.0    # workplace demand weight relative to residential population
P      = 300.0  # peak travel time (seconds); flow peaks at d = P
ALPHA  = 2.0    # tail decay exponent; flow ~ 1/d^ALPHA for large d
BETA   = 1.0    # rise exponent; u^BETA approach to peak from origin

OUT_DIR    = "simulation"
CONS_GRAPH = "simulation/newtownards_consolidated.graphml"

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

node_population      = {int(k): v for k, v in weights["node_population"].items()}
node_business_demand = {int(k): v for k, v in weights["node_business_demand"].items()}
node_school_demand   = {int(k): v for k, v in weights.get("node_school_demand", {}).items()}

THETA          = None
K_res          = None
K_biz          = None
K_sch          = None
P_biz          = None
ALPHA_biz      = None
W_SCHOOL       = None
P_school       = None
ALPHA_school   = None
slot_fracs_res    = {}
slot_fracs_biz    = {}
slot_fracs_school = {}

if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        _tp = json.load(f)
    K         = _tp.get("K",         K)
    W_BIZ     = _tp.get("W_BIZ",     W_BIZ)
    P         = _tp.get("P",         P)
    ALPHA     = _tp.get("ALPHA",     ALPHA)
    BETA      = _tp.get("BETA",      BETA)
    THETA     = _tp.get("THETA",     None)
    P_biz     = _tp.get("P_biz",     None)
    ALPHA_biz = _tp.get("ALPHA_biz", None)
    W_SCHOOL     = _tp.get("W_SCHOOL",     None)
    P_school     = _tp.get("P_school",     None)
    ALPHA_school = _tp.get("ALPHA_school", None)
    if "K_res" in _tp and "K_biz" in _tp:
        K_res = _tp["K_res"]
        K_biz = _tp["K_biz"]
        K_sch = _tp.get("K_sch", 0.0)
        slot_fracs_res = {tuple(int(x) for x in k.split(",")): v
                          for k, v in _tp.get("slot_fracs_res", {}).items()}
        slot_fracs_biz = {tuple(int(x) for x in k.split(",")): v
                          for k, v in _tp.get("slot_fracs_biz", {}).items()}
        slot_fracs_school = {tuple(int(x) for x in k.split(",")): v
                             for k, v in _tp.get("slot_fracs_school", {}).items()}
    print(f"  [tuned: stage={_tp.get('stage','?')}  χ²/N={_tp.get('chi2_per_n','?')}"
          + (f"  THETA={THETA:.4f}"             if THETA is not None     else "")
          + (f"  K_res={K_res:.3e}  K_biz={K_biz:.3e}  K_sch={K_sch:.3e}" if K_res is not None else "")
          + (f"  P_biz={P_biz:.1f}s  ALPHA_biz={ALPHA_biz:.4f}" if P_biz is not None else "")
          + (f"  W_SCHOOL={W_SCHOOL:.4f}  P_school={P_school:.1f}s" if W_SCHOOL is not None else "")
          + "]")

# ── Assignment ────────────────────────────────────────────────────────────────

print(f"Loading paths cache ({PATHS_CACHE}) …")
t0    = time.time()
cache = np.load(PATHS_CACHE)

node_ids_arr = cache["node_ids"]
od_src       = cache["od_src"]
od_dst       = cache["od_dst"]
od_dist      = cache["od_dist"].astype(np.float64)
pair_idx     = cache["pair_idx"]
link_idx     = cache["link_idx"]
link_u       = cache["link_u"]
link_v       = cache["link_v"]

# Probit stochastic loading (new cache format)
if "link_weight" in cache:
    link_weight = cache["link_weight"].astype(np.float64)
    n_passes    = int(cache["probit_n_passes"])
    cv          = float(cache["probit_cv"])
    print(f"  Probit loading: {n_passes} passes  CV={cv:.2f}  "
          f"({len(link_idx):,} weighted entries for {len(od_src):,} OD pairs)")
    od_dist_2 = pair_idx_2 = link_idx_2 = None
    od_dist_3 = pair_idx_3 = link_idx_3 = None
    THETA = None   # not used with probit cache
else:
    # Legacy k=2/k=3 logit cache
    link_weight = None
    _has_stoch  = "pair_idx_2" in cache and THETA is not None
    if _has_stoch:
        od_dist_2  = cache["od_dist_2"].astype(np.float64)
        pair_idx_2 = cache["pair_idx_2"]
        link_idx_2 = cache["link_idx_2"]
        od_dist_3  = cache["od_dist_3"].astype(np.float64)
        pair_idx_3 = cache["pair_idx_3"]
        link_idx_3 = cache["link_idx_3"]
        print(f"  Legacy k=3 paths loaded  THETA={THETA:.4f}")
    else:
        od_dist_2 = pair_idx_2 = link_idx_2 = None
        od_dist_3 = pair_idx_3 = link_idx_3 = None
        if THETA is not None:
            print("  Warning: THETA in params but no stochastic paths in cache — using all-or-nothing")
            THETA = None

w_pop    = np.array([node_population.get(int(nid), 0)      for nid in node_ids_arr], dtype=np.float64)
w_biz    = np.array([node_business_demand.get(int(nid), 0) for nid in node_ids_arr], dtype=np.float64)
w_school = np.array([node_school_demand.get(int(nid), 0)   for nid in node_ids_arr], dtype=np.float64)
print(f"  {len(node_ids_arr)} nodes  total weight {(w_pop + W_BIZ * w_biz).sum():,.0f}  (W_BIZ={W_BIZ})")

N_links  = len(link_u)
_use_3c  = (K_res is not None and K_sch is not None and K_sch > 0
            and W_SCHOOL is not None and w_school.sum() > 0)
_use_2c  = (K_res is not None and not _use_3c)

_kw = dict(BETA=BETA, THETA=THETA,
           P_biz=P_biz, ALPHA_biz=ALPHA_biz,
           od_dist_2=od_dist_2, pair_idx_2=pair_idx_2, link_idx_2=link_idx_2,
           od_dist_3=od_dist_3, pair_idx_3=pair_idx_3, link_idx_3=link_idx_3,
           link_weight=link_weight)

if _use_3c:
    _kw_3c = dict(**_kw, W_SCHOOL=W_SCHOOL, P_school=P_school,
                  ALPHA_school=ALPHA_school, w_school=w_school)
    raw_res, raw_biz, raw_sch = gravity_assign(
        od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
        W_BIZ, P, ALPHA, w_pop, w_biz, return_components=True, **_kw_3c)
    _nonzero = (raw_res + raw_biz + raw_sch) > 0
    link_flow_res    = {(int(link_u[k]), int(link_v[k])): raw_res[k] * K_res
                        for k in range(N_links) if _nonzero[k]}
    link_flow_biz    = {(int(link_u[k]), int(link_v[k])): raw_biz[k] * K_biz
                        for k in range(N_links) if _nonzero[k]}
    link_flow_school = {(int(link_u[k]), int(link_v[k])): raw_sch[k] * K_sch
                        for k in range(N_links) if _nonzero[k]}
    link_flow = {lnk: (link_flow_res.get(lnk, 0.0) + link_flow_biz.get(lnk, 0.0)
                       + link_flow_school.get(lnk, 0.0))
                 for lnk in set(link_flow_res) | set(link_flow_biz) | set(link_flow_school)}
elif _use_2c:
    raw_res, raw_biz = gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
                                      W_BIZ, P, ALPHA, w_pop, w_biz,
                                      return_components=True, **_kw)
    link_flow_res = {(int(link_u[k]), int(link_v[k])): raw_res[k] * K_res
                     for k in range(N_links) if raw_res[k] + raw_biz[k] > 0}
    link_flow_biz = {(int(link_u[k]), int(link_v[k])): raw_biz[k] * K_biz
                     for k in range(N_links) if raw_res[k] + raw_biz[k] > 0}
    link_flow = {lnk: link_flow_res.get(lnk, 0.0) + link_flow_biz.get(lnk, 0.0)
                 for lnk in set(link_flow_res) | set(link_flow_biz)}
    link_flow_school = None
else:
    raw_flow_arr = gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
                                  W_BIZ, P, ALPHA, w_pop, w_biz, **_kw)
    link_flow = {(int(link_u[k]), int(link_v[k])): raw_flow_arr[k] * K
                 for k in range(N_links) if raw_flow_arr[k] > 0}
    link_flow_res = link_flow_biz = link_flow_school = None

print(f"  Assignment complete in {time.time()-t0:.2f}s  ({len(link_flow)} loaded links)")

# ── Street name lookup ────────────────────────────────────────────────────────

node_ids    = [int(nid) for nid in node_ids_arr]
G           = ox.load_graphml(CONS_GRAPH)
node_weight = {int(nid): float(wp + W_BIZ * wb)
               for nid, wp, wb in zip(node_ids_arr, w_pop, w_biz)}

_link_name = {(int(u), int(v)): d["name"]
              for u, v, d in G.edges(data=True) if d.get("name")}


def _link_label(u, v):
    name = _link_name.get((u, v), "")
    return f"{u}→{v}  {name}" if name else f"{u}→{v}"

# ── Report ────────────────────────────────────────────────────────────────────

print(f"\nOfficial count sites  (K = {K:.4e}"
      + (f"  K_res={K_res:.3e}  K_biz={K_biz:.3e}"
         + (f"  K_sch={K_sch:.3e}" if _use_3c else "")
         if (_use_2c or _use_3c) else "") + ")")
print(f"  {'Site':<45s}  {'Modelled':>9s}  {'Observed':>9s}  {'Ratio':>6s}")
for s in COUNT_SITES:
    f = site_flow(link_flow, s)
    print(f"  {s['label']:<45s}  {f:>9,.0f}  {s['observed']:>9,}  {f/s['observed']:>6.2f}")

rows, chi2, n_obs, n_eff = compute_chi2(
    link_flow_res if (_use_2c or _use_3c) else link_flow,
    label_fn=_link_label,
    link_aadt_file=LINK_AADT,
    exclude_links=EXCLUDE_LINKS,
    link_flow_biz_dict=link_flow_biz if (_use_2c or _use_3c) else None,
    link_flow_school_dict=link_flow_school if _use_3c else None,
    slot_fracs_res=slot_fracs_res if (_use_2c or _use_3c) else None,
    slot_fracs_biz=slot_fracs_biz if (_use_2c or _use_3c) else None,
    slot_fracs_school=slot_fracs_school if _use_3c else None,
)
print_chi2_table(rows, chi2, n_obs, n_eff=n_eff)

# ── Serialise flows ───────────────────────────────────────────────────────────

flows_path = f"{OUT_DIR}/newtownards_flows.json"
out = {
    "kernel": "rational", "W_BIZ": W_BIZ, "P": P, "ALPHA": ALPHA, "BETA": BETA, "K": K,
    "flows": {f"{u},{v}": flow for (u, v), flow in link_flow.items()},
}
if _use_3c or _use_2c:
    out["K_res"] = K_res
    out["K_biz"] = K_biz
    out["flows_res"] = {f"{u},{v}": flow for (u, v), flow in link_flow_res.items()}
    out["flows_biz"] = {f"{u},{v}": flow for (u, v), flow in link_flow_biz.items()}
if _use_3c:
    out["K_sch"] = K_sch
    out["flows_school"] = {f"{u},{v}": flow for (u, v), flow in link_flow_school.items()}
with open(flows_path, "w") as f:
    json.dump(out, f)
_comp_str = ("+ res/biz/school" if _use_3c else ("+ res/biz" if _use_2c else ""))
print(f"\nSaved {len(link_flow)} link flows → {flows_path}"
      + (f"  ({_comp_str} components)" if _comp_str else ""))
print(f"Parameters: K={K}  W_BIZ={W_BIZ}  P={P}  ALPHA={ALPHA}  BETA={BETA}"
      + (f"  P_biz={P_biz}  ALPHA_biz={ALPHA_biz}" if P_biz is not None else "")
      + (f"  W_SCHOOL={W_SCHOOL}  P_school={P_school}" if W_SCHOOL is not None else ""))
