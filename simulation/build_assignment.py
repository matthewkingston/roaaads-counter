"""
Gravity model OD matrix + all-or-nothing shortest-path assignment.

Loads node weights from simulation/node_weights.json (produced by
build_demographics.py) and assigns traffic to links using a gravity model.

Fast path (recommended): load precomputed paths from simulation/newtownards_paths.npz
(produced by build_paths.py). Reduces parameter-tuning runs to < 1s.

Slow fallback: if the paths cache is absent, runs NetworkX Dijkstra (~10s).
Re-run build_paths.py whenever the road network or external zone coordinates change.

Outputs newtownards_flows.json; run build_demographics.py afterwards to see
flows on the map.

Usage:
  python3 simulation/build_assignment.py

Tunable parameters: K, W_BIZ, P, ALPHA (see Config section).
Chi²/N uses per-session observations (same data as tuner). No Woodbury correction
is applied, so the value will be slightly higher than the tuner's figure.
"""

import json, math, time, os
import numpy as np
import osmnx as ox
from collections import defaultdict
from model import (COUNT_SITES, EXCLUDE_LINKS, PATHS_CACHE, WEIGHTS_FILE,
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS,
                   gravity_assign, site_flow, compute_chi2, print_chi2_table)

# ── Config ────────────────────────────────────────────────────────────────────

K      = 1.73   # global flow scale factor
W_BIZ  = 1.0    # workplace demand weight relative to residential population
P      = 300.0  # peak travel time (seconds); flow peaks at d = P
ALPHA  = 2.0    # tail decay exponent; flow ~ 1/d^ALPHA for large d
OFFSCREEN_SPEED_MS = 80_000 / 3600   # 80 km/h — assumed speed for off-network boundary legs

OUT_DIR    = "simulation"
CONS_GRAPH = "simulation/newtownards_consolidated.graphml"

# ── Load node weights ─────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)

node_population      = {int(k): v for k, v in weights["node_population"].items()}
node_business_demand = {int(k): v for k, v in weights["node_business_demand"].items()}
node_effective_utm   = {int(k): (v[0], v[1]) for k, v in weights["node_effective_utm"].items()}
boundary_node_ids    = set(weights["boundary_node_ids"])

with open(TUNER_CONFIG) as f:
    _tuner_cfg = json.load(f)
_city_nodes = {name: cfg["nodes"] for name, cfg in _tuner_cfg["cities"].items()}
allowed_through_pairs = set()
for _city_a, _city_b in _tuner_cfg.get("through_route_pairs", []):
    for _na in _city_nodes[_city_a]:
        for _nb in _city_nodes[_city_b]:
            allowed_through_pairs.add((_na, _nb))
            allowed_through_pairs.add((_nb, _na))

if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        _tp = json.load(f)
    K     = _tp.get("K",     K)
    W_BIZ = _tp.get("W_BIZ", W_BIZ)
    P     = _tp.get("P",     P)
    ALPHA = _tp.get("ALPHA", ALPHA)
    for _nid, _val in _tp.get("external_node_pop", {}).items():
        node_population[int(_nid)] = _val
    for _nid, _val in _tp.get("external_node_biz", {}).items():
        node_business_demand[int(_nid)] = _val
    print(f"  [tuned: stage={_tp.get('stage','?')}  χ²/N={_tp.get('chi2_per_n','?')}]")

# ── Assignment ────────────────────────────────────────────────────────────────

link_flow = {}   # (u, v) → raw (pre-K) flow

if os.path.exists(PATHS_CACHE):
    # ── Fast path: vectorised numpy using precomputed paths ───────────────────
    print(f"Loading paths cache ({PATHS_CACHE}) …")
    t0 = time.time()
    cache = np.load(PATHS_CACHE)

    node_ids_arr = cache["node_ids"]          # int32 (N,)
    od_src       = cache["od_src"]            # int32 (P,)  source node index
    od_dst       = cache["od_dst"]            # int32 (P,)  target node index
    od_dist      = cache["od_dist"].astype(np.float64)  # (P,) effective travel time (seconds)
    pair_idx     = cache["pair_idx"]          # int32 (E,)  which OD pair
    link_idx     = cache["link_idx"]          # int32 (E,)  which link
    link_u       = cache["link_u"]            # int32 (L,)
    link_v       = cache["link_v"]            # int32 (L,)

    # Node weight vector in the same order as node_ids_arr
    w_pop = np.array([node_population.get(int(nid), 0)      for nid in node_ids_arr], dtype=np.float64)
    w_biz = np.array([node_business_demand.get(int(nid), 0) for nid in node_ids_arr], dtype=np.float64)
    total_weight = (w_pop + W_BIZ * w_biz).sum()
    print(f"  {len(node_ids_arr)} nodes  total weight {total_weight:,.0f}  (W_BIZ={W_BIZ})")

    N_links      = len(link_u)
    raw_flow_arr = gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
                                  W_BIZ, P, ALPHA, w_pop, w_biz)

    link_flow = {
        (int(link_u[k]), int(link_v[k])): raw_flow_arr[k]
        for k in range(N_links) if raw_flow_arr[k] > 0
    }
    print(f"  Assignment complete in {time.time()-t0:.2f}s  ({len(link_flow)} loaded links)")

    # Reconstruct node_ids list and weights (needed for map / boundary markers)
    node_ids    = [int(nid) for nid in node_ids_arr]
    G           = ox.load_graphml(CONS_GRAPH)
    node_weight = {int(nid): float(wp + W_BIZ * wb)
                   for nid, wp, wb in zip(node_ids_arr, w_pop, w_biz)}

else:
    # ── Slow fallback: per-source NetworkX Dijkstra (~10s) ────────────────────
    print("WARNING: paths cache not found — run build_paths.py for faster iterations")
    import networkx as nx

    print("Loading consolidated graph …")
    G = ox.load_graphml(CONS_GRAPH)
    G = ox.speed.add_edge_speeds(G)
    G = ox.speed.add_edge_travel_times(G)
    node_ids = list(G.nodes())
    print(f"  {len(node_ids)} nodes  {G.number_of_edges()} edges")

    node_weight = {}
    for nid in node_ids:
        pop = node_population.get(nid, 0)
        biz = node_business_demand.get(nid, 0)
        node_weight[nid] = pop + W_BIZ * biz

    boundary_offscreen = {}
    for nid in boundary_node_ids:
        if nid not in G.nodes:
            continue
        nd = G.nodes[nid]
        nx_pos, ny_pos = float(nd["x"]), float(nd["y"])
        cx, cy = node_effective_utm[nid]
        boundary_offscreen[nid] = math.sqrt((nx_pos - cx) ** 2 + (ny_pos - cy) ** 2) / OFFSCREEN_SPEED_MS

    for nid, extra in boundary_offscreen.items():
        print(f"  Node {nid:4d}  offscreen leg: {extra:.0f} s")

    print("Running assignment …")
    t0 = time.time()
    lf = defaultdict(float)

    for idx, source in enumerate(node_ids):
        w_i = node_weight[source]
        if w_i <= 0:
            continue
        src_is_boundary = source in boundary_node_ids
        try:
            lengths, paths = nx.single_source_dijkstra(G, source, weight="travel_time")
        except Exception:
            continue
        for target, path in paths.items():
            if target == source or len(path) < 2:
                continue
            if src_is_boundary and target in boundary_node_ids:
                if (source, target) not in allowed_through_pairs:
                    continue
            w_j = node_weight[target]
            if w_j <= 0:
                continue
            dist = (lengths[target]
                    + boundary_offscreen.get(target, 0)
                    + boundary_offscreen.get(source, 0))
            if dist < 1.0:
                continue
            _u   = dist / P
            t_ij = w_i * w_j * (ALPHA + 1) * _u / (ALPHA + _u ** (ALPHA + 1))
            for u, v in zip(path[:-1], path[1:]):
                lf[(u, v)] += t_ij
        if (idx + 1) % 100 == 0:
            print(f"  {idx + 1}/{len(node_ids)} nodes  ({time.time()-t0:.1f}s)")

    link_flow = dict(lf)
    print(f"  Assignment complete in {time.time()-t0:.1f}s  ({len(link_flow)} loaded links)")

# ── Street name lookup (from already-loaded graph) ────────────────────────────

_link_name = {(int(u), int(v)): d["name"]
              for u, v, d in G.edges(data=True) if d.get("name")}


def _link_label(u, v):
    name = _link_name.get((u, v), "")
    return f"{u}→{v}  {name}" if name else f"{u}→{v}"

# ── Scale by K and report ─────────────────────────────────────────────────────

link_flow = {k: v * K for k, v in link_flow.items()}

print(f"\nOfficial count sites  (K = {K})")
print(f"  {'Site':<45s}  {'Modelled':>9s}  {'Observed':>9s}  {'Ratio':>6s}")
for s in COUNT_SITES:
    f = site_flow(link_flow, s)
    print(f"  {s['label']:<45s}  {f:>9,.0f}  {s['observed']:>9,}  {f/s['observed']:>6.2f}")

rows, chi2, n = compute_chi2(link_flow, label_fn=_link_label,
                             link_aadt_file=LINK_AADT, exclude_links=EXCLUDE_LINKS)
print_chi2_table(rows, chi2, n)

# ── Serialise flows ───────────────────────────────────────────────────────────

flows_path = f"{OUT_DIR}/newtownards_flows.json"
with open(flows_path, "w") as f:
    json.dump({
        "kernel": "rational", "W_BIZ": W_BIZ, "P": P, "ALPHA": ALPHA, "K": K,
        "flows": {f"{u},{v}": flow for (u, v), flow in link_flow.items()},
    }, f)
print(f"\nSaved {len(link_flow)} link flows → {flows_path}")
print(f"Parameters: K={K}  W_BIZ={W_BIZ}  P={P}  ALPHA={ALPHA}")
