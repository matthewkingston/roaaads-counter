"""
Precompute all-pairs shortest-time paths for the Newtownards road network.

Run this once after any network change (build_network.py or edit_network.py).
Produces simulation/newtownards_paths.npz, which build_assignment.py loads
to avoid re-running Dijkstra on every parameter-tuning iteration.

od_dist in the cache stores generalised cost (seconds × road-class factor):
edge travel time multiplied by HIGHWAY_COST_FACTOR, plus the boundary
offscreen leg converted at OFFSCREEN_SPEED_MS. This means major roads appear
shorter than minor roads of equal length, biasing route choice toward them.

Three path alternatives are stored per OD pair (k=1, k=2, k=3):
  k=1: original shortest paths
  k=2: shortest paths after penalising all k=1 edges ×ALT_COST_PENALTY
  k=3: shortest paths after penalising all k=1+k=2 edges ×ALT_COST_PENALTY
All three use original (non-penalised) edge costs for od_dist.
Pairs with no alternative path fall back to the k=1 path (same links, same dist).
These are used for logit stochastic routing (THETA parameter in tune_assignment.py).

Cache must be regenerated when:
  - The road network changes (newtownards_consolidated.graphml is updated)
  - External zone coordinates (lat/lon) change in build_demographics.py EXTERNAL_ZONES
  - HIGHWAY_COST_FACTOR values change

Does NOT need regenerating for changes to W_BIZ, P, ALPHA, THETA, node
populations, damping factors, or count site values.
"""

import json, math, time, os
import numpy as np
import osmnx as ox
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

CONS_GRAPH          = "simulation/newtownards_consolidated.graphml"
WEIGHTS_FILE        = "simulation/node_weights.json"
TUNER_CONFIG        = "simulation/tuner_config.json"
PATHS_CACHE         = "simulation/newtownards_paths.npz"
OFFSCREEN_SPEED_MS  = 80_000 / 3600   # 80 km/h in m/s — assumed speed for off-network boundary legs
ALT_COST_PENALTY    = 10.0             # multiplier applied to k=1 (and k=1+k=2) edges when finding alternatives

# Multipliers applied to edge travel time before Dijkstra.
# < 1 favours that road class; > 1 penalises it.
# Changing these values requires re-running build_paths.py and re-tuning.
HIGHWAY_COST_FACTOR = {
    "trunk":         0.67,
    "trunk_link":    0.67,
    "primary":       0.67,
    "primary_link":  0.67,
    "secondary":     1.0,
    "tertiary":      1.0,
    "tertiary_link": 1.0,
    "residential":   1.2,
    "unclassified":  1.2,
    "living_street": 1.2,
}

# ── Load graph and weights ───────────────────────────────────────────────────────

print("Loading graph …")
G = ox.load_graphml(CONS_GRAPH)
G = ox.routing.add_edge_speeds(G)
G = ox.routing.add_edge_travel_times(G)
node_ids = list(G.nodes())
n = len(node_ids)
node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
print(f"  {n} nodes  {G.number_of_edges()} edges")

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)
node_effective_utm = {int(k): (v[0], v[1]) for k, v in weights["node_effective_utm"].items()}
boundary_node_ids  = set(weights["boundary_node_ids"])

with open(TUNER_CONFIG) as f:
    tuner_cfg = json.load(f)
_city_nodes = {name: cfg["nodes"] for name, cfg in tuner_cfg["cities"].items()}
allowed_through_pairs = set()
for city_a, city_b in tuner_cfg.get("through_route_pairs", []):
    for na in _city_nodes[city_a]:
        for nb in _city_nodes[city_b]:
            allowed_through_pairs.add((na, nb))
            allowed_through_pairs.add((nb, na))
print(f"  Through-route pairs: {len(tuner_cfg.get('through_route_pairs', []))} city pairs → {len(allowed_through_pairs)} node pairs")

# ── Offscreen legs ───────────────────────────────────────────────────────────────

boundary_offscreen = {}
for nid in boundary_node_ids:
    if nid not in G.nodes:
        continue
    nd = G.nodes[nid]
    nx_pos, ny_pos = float(nd["x"]), float(nd["y"])
    cx, cy = node_effective_utm[nid]
    boundary_offscreen[nid] = math.sqrt((nx_pos - cx) ** 2 + (ny_pos - cy) ** 2) / OFFSCREEN_SPEED_MS

# ── Build sparse adjacency matrix ────────────────────────────────────────────────

print("Building sparse adjacency matrix …")
rows, cols, data = [], [], []
link_list    = []
link_to_idx  = {}
row_link_idx = []     # parallel to rows/cols/data: which link index each entry belongs to
edge_orig_cost = {}   # (matrix_i, matrix_j) → original cost (for retracing alternative paths)

for u, v, edata in G.edges(data=True):
    ht = edata.get("highway", "unclassified")
    if isinstance(ht, list):
        ht = ht[0]
    factor = HIGHWAY_COST_FACTOR.get(ht, 1.0)
    cost = float(edata.get("travel_time", 1.0)) * factor
    i, j = node_to_idx[u], node_to_idx[v]
    rows.append(i); cols.append(j); data.append(cost)
    lnk = (u, v)
    if lnk not in link_to_idx:
        link_to_idx[lnk] = len(link_list)
        link_list.append(lnk)
    li = link_to_idx[lnk]
    row_link_idx.append(li)
    # Keep minimum cost for each (i,j) pair (relevant if parallel edges exist)
    if (i, j) not in edge_orig_cost or cost < edge_orig_cost[(i, j)]:
        edge_orig_cost[(i, j)] = cost

adj = csr_matrix((data, (rows, cols)), shape=(n, n))
print(f"  {len(link_list)} directed links")

# ── Helper: build penalised adjacency ────────────────────────────────────────────

def _build_penalized_adj(penalize_link_set):
    """Multiply costs of edges whose link index is in penalize_link_set by ALT_COST_PENALTY."""
    data_pen = [
        cost * ALT_COST_PENALTY if li in penalize_link_set else cost
        for cost, li in zip(data, row_link_idx)
    ]
    return csr_matrix((data_pen, (rows, cols)), shape=(n, n))

# ── Helper: retrace alternative paths for all OD pairs ───────────────────────────

def _retrace_alt(pred, od_src_list, od_dst_list, pair1_starts, link_idx_arr_1, od_dist_list, t_ref):
    """
    Retrace paths using pred (from a penalised Dijkstra) for each OD pair.
    Distances are computed using original (non-penalised) edge costs.
    Falls back to k=1 path for any pair where pred has no route.
    Returns (od_dist_out, pair_idx_out, link_idx_out).
    """
    pair_idx_out = []
    link_idx_out = []
    od_dist_out  = []

    for k, (src_i, dst_i) in enumerate(zip(od_src_list, od_dst_list)):
        src_nid = node_ids[src_i]
        dst_nid = node_ids[dst_i]

        path_links = []
        path_cost  = 0.0
        cur = dst_i
        while cur != src_i:
            prev = pred[src_i, cur]
            if prev < 0:
                path_links = []
                break
            path_cost += edge_orig_cost.get((prev, cur), 0.0)
            lnk = (node_ids[prev], node_ids[cur])
            if lnk in link_to_idx:
                path_links.append(link_to_idx[lnk])
            cur = prev

        if path_links:
            path_links.reverse()
            eff_dist = path_cost + boundary_offscreen.get(src_nid, 0.0) + boundary_offscreen.get(dst_nid, 0.0)
            od_dist_out.append(max(float(eff_dist), 1.0))
            for li in path_links:
                pair_idx_out.append(k)
                link_idx_out.append(li)
        else:
            # No alternative found — duplicate k=1 path
            od_dist_out.append(od_dist_list[k])
            s, e = int(pair1_starts[k]), int(pair1_starts[k + 1])
            for li in link_idx_arr_1[s:e]:
                pair_idx_out.append(k)
                link_idx_out.append(int(li))

        if (k + 1) % 100_000 == 0:
            print(f"    {k+1:,}/{len(od_src_list):,}  ({time.time()-t_ref:.0f}s)")

    return od_dist_out, pair_idx_out, link_idx_out

# ── k=1: all-pairs Dijkstra (scipy) ─────────────────────────────────────────────

print("Running k=1 Dijkstra (scipy) …")
t0 = time.time()
dist_matrix, predecessors = dijkstra(adj, directed=True, return_predecessors=True)
print(f"  Done in {time.time()-t0:.1f}s")

# ── Reconstruct k=1 paths ────────────────────────────────────────────────────────

print("Reconstructing k=1 paths …")
t0 = time.time()

od_src_list   = []
od_dst_list   = []
od_dist_list  = []
pair_idx_list = []
link_idx_list = []

boundary_idx = {node_to_idx[nid] for nid in boundary_node_ids if nid in node_to_idx}

for src_i, src_nid in enumerate(node_ids):
    src_is_boundary = src_i in boundary_idx
    src_offscreen   = boundary_offscreen.get(src_nid, 0.0)

    for dst_i, dst_nid in enumerate(node_ids):
        if dst_i == src_i:
            continue
        if src_is_boundary and dst_i in boundary_idx:
            if (src_nid, dst_nid) not in allowed_through_pairs:
                continue
        d_net = dist_matrix[src_i, dst_i]
        if not math.isfinite(d_net):
            continue

        eff_dist = d_net + src_offscreen + boundary_offscreen.get(dst_nid, 0.0)
        if eff_dist < 1.0:
            continue

        path_links = []
        cur = dst_i
        while cur != src_i:
            prev = predecessors[src_i, cur]
            if prev < 0:
                break
            lnk = (node_ids[prev], node_ids[cur])
            if lnk in link_to_idx:
                path_links.append(link_to_idx[lnk])
            cur = prev
        if not path_links:
            continue
        path_links.reverse()

        pair_k = len(od_src_list)
        od_src_list.append(src_i)
        od_dst_list.append(dst_i)
        od_dist_list.append(eff_dist)
        for li in path_links:
            pair_idx_list.append(pair_k)
            link_idx_list.append(li)

    if (src_i + 1) % 100 == 0:
        print(f"  {src_i + 1}/{n} nodes  ({time.time()-t0:.1f}s)  "
              f"{len(od_src_list):,} pairs  {len(pair_idx_list):,} entries")

print(f"  {len(od_src_list):,} OD pairs  {len(pair_idx_list):,} pair-link entries  "
      f"in {time.time()-t0:.1f}s")

# ── Build per-pair lookup for k=1 fallback ───────────────────────────────────────

pair_idx_arr_1 = np.array(pair_idx_list, dtype=np.int32)
link_idx_arr_1 = np.array(link_idx_list, dtype=np.int32)
pair1_counts   = np.bincount(pair_idx_arr_1, minlength=len(od_src_list))
pair1_starts   = np.concatenate([[0], np.cumsum(pair1_counts)]).astype(np.int64)
used_links_1   = set(link_idx_arr_1.tolist())
print(f"  k=1 uses {len(used_links_1):,} of {len(link_list):,} links")

# ── k=2: penalise k=1 edges, re-run Dijkstra ────────────────────────────────────

print(f"Building k=2 penalised adjacency ({len(used_links_1):,} links ×{ALT_COST_PENALTY}) …")
adj_2 = _build_penalized_adj(used_links_1)

print("Running k=2 Dijkstra …")
t0 = time.time()
_, pred_2 = dijkstra(adj_2, directed=True, return_predecessors=True)
print(f"  Done in {time.time()-t0:.1f}s")

print("Reconstructing k=2 paths …")
t0 = time.time()
od_dist_list_2, pair_idx_list_2, link_idx_list_2 = _retrace_alt(
    pred_2, od_src_list, od_dst_list, pair1_starts, link_idx_arr_1, od_dist_list, t0)
used_links_2 = set(link_idx_list_2)
n_same_2 = sum(1 for k in range(len(od_src_list))
               if od_dist_list_2[k] == od_dist_list[k])
print(f"  {len(pair_idx_list_2):,} entries  {n_same_2:,} pairs fell back to k=1  "
      f"({time.time()-t0:.1f}s)")

# ── k=3: penalise k=1+k=2 edges, re-run Dijkstra ────────────────────────────────

used_links_12 = used_links_1 | used_links_2
print(f"Building k=3 penalised adjacency ({len(used_links_12):,} links ×{ALT_COST_PENALTY}) …")
adj_3 = _build_penalized_adj(used_links_12)

print("Running k=3 Dijkstra …")
t0 = time.time()
_, pred_3 = dijkstra(adj_3, directed=True, return_predecessors=True)
print(f"  Done in {time.time()-t0:.1f}s")

print("Reconstructing k=3 paths …")
t0 = time.time()
od_dist_list_3, pair_idx_list_3, link_idx_list_3 = _retrace_alt(
    pred_3, od_src_list, od_dst_list, pair1_starts, link_idx_arr_1, od_dist_list, t0)
n_same_3 = sum(1 for k in range(len(od_src_list))
               if od_dist_list_3[k] == od_dist_list[k])
print(f"  {len(pair_idx_list_3):,} entries  {n_same_3:,} pairs fell back to k=1  "
      f"({time.time()-t0:.1f}s)")

# ── Save ─────────────────────────────────────────────────────────────────────────

print("Saving cache …")
np.savez_compressed(
    PATHS_CACHE,
    node_ids    = np.array(node_ids,         dtype=np.int32),
    od_src      = np.array(od_src_list,      dtype=np.int32),
    od_dst      = np.array(od_dst_list,      dtype=np.int32),
    od_dist     = np.array(od_dist_list,     dtype=np.float32),
    pair_idx    = np.array(pair_idx_list,    dtype=np.int32),
    link_idx    = np.array(link_idx_list,    dtype=np.int32),
    link_u      = np.array([u for u, v in link_list], dtype=np.int32),
    link_v      = np.array([v for u, v in link_list], dtype=np.int32),
    # Alternative paths (k=2 and k=3) for stochastic logit routing
    od_dist_2   = np.array(od_dist_list_2,   dtype=np.float32),
    pair_idx_2  = np.array(pair_idx_list_2,  dtype=np.int32),
    link_idx_2  = np.array(link_idx_list_2,  dtype=np.int32),
    od_dist_3   = np.array(od_dist_list_3,   dtype=np.float32),
    pair_idx_3  = np.array(pair_idx_list_3,  dtype=np.int32),
    link_idx_3  = np.array(link_idx_list_3,  dtype=np.int32),
)
size_mb = os.path.getsize(PATHS_CACHE) / 1e6
print(f"Saved: {PATHS_CACHE}  ({size_mb:.1f} MB)")
print(f"  k=1: {len(od_src_list):,} pairs  {len(pair_idx_list):,} entries")
print(f"  k=2: {len(od_src_list):,} pairs  {len(pair_idx_list_2):,} entries  ({n_same_2:,} fallbacks)")
print(f"  k=3: {len(od_src_list):,} pairs  {len(pair_idx_list_3):,} entries  ({n_same_3:,} fallbacks)")
