"""
Precompute all-pairs shortest-time paths for the Newtownards road network.

Run this once after any network change (build_network.py or edit_network.py).
Produces simulation/newtownards_paths.npz, which build_assignment.py loads
to avoid re-running Dijkstra on every parameter-tuning iteration.

od_dist in the cache stores effective travel time in seconds: network travel
time (from edge travel_time attributes) plus the boundary offscreen leg
converted at OFFSCREEN_SPEED_MS.

Cache must be regenerated when:
  - The road network changes (newtownards_consolidated.graphml is updated)
  - External zone coordinates (lat/lon) change in build_demographics.py EXTERNAL_ZONES

Does NOT need regenerating for changes to W_BIZ, MU, SIGMA, ALPHA, node
populations, damping factors, or count site values.
"""

import json, math, time, os
import numpy as np
import osmnx as ox
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

CONS_GRAPH          = "simulation/newtownards_consolidated.graphml"
WEIGHTS_FILE        = "simulation/node_weights.json"
PATHS_CACHE         = "simulation/newtownards_paths.npz"
OFFSCREEN_SPEED_MS  = 80_000 / 3600   # 80 km/h in m/s — assumed speed for off-network boundary legs

# ── Load graph and weights ───────────────────────────────────────────────────────

print("Loading graph …")
G = ox.load_graphml(CONS_GRAPH)
G = ox.speed.add_edge_speeds(G)
G = ox.speed.add_edge_travel_times(G)
node_ids = list(G.nodes())
n = len(node_ids)
node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
print(f"  {n} nodes  {G.number_of_edges()} edges")

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)
node_effective_utm = {int(k): (v[0], v[1]) for k, v in weights["node_effective_utm"].items()}
boundary_node_ids  = set(weights["boundary_node_ids"])

# ── Offscreen legs ───────────────────────────────────────────────────────────────

boundary_offscreen = {}
for nid in boundary_node_ids:
    if nid not in G.nodes:
        continue
    nd = G.nodes[nid]
    nx_pos, ny_pos = float(nd["x"]), float(nd["y"])
    cx, cy = node_effective_utm[nid]
    boundary_offscreen[nid] = math.sqrt((nx_pos - cx) ** 2 + (ny_pos - cy) ** 2) / OFFSCREEN_SPEED_MS

# ── Build scipy sparse adjacency matrix ──────────────────────────────────────────

print("Building sparse adjacency matrix …")
rows, cols, data = [], [], []
link_list  = []
link_to_idx = {}

for u, v, edata in G.edges(data=True):
    cost = float(edata.get("travel_time", 1.0))
    i, j = node_to_idx[u], node_to_idx[v]
    rows.append(i); cols.append(j); data.append(cost)
    lnk = (u, v)
    if lnk not in link_to_idx:
        link_to_idx[lnk] = len(link_list)
        link_list.append(lnk)

adj = csr_matrix((data, (rows, cols)), shape=(n, n))
print(f"  {len(link_list)} directed links")

# ── All-pairs Dijkstra (scipy) ───────────────────────────────────────────────────

print("Running all-pairs Dijkstra (scipy) …")
t0 = time.time()
dist_matrix, predecessors = dijkstra(adj, directed=True, return_predecessors=True)
print(f"  Done in {time.time()-t0:.1f}s")

# ── Reconstruct paths and build OD pair + link arrays ───────────────────────────

print("Reconstructing paths and building cache arrays …")
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
            continue
        d_net = dist_matrix[src_i, dst_i]
        if not math.isfinite(d_net):
            continue

        eff_dist = d_net + src_offscreen + boundary_offscreen.get(dst_nid, 0.0)
        if eff_dist < 1.0:
            continue

        # Reconstruct path by backtracking through predecessors
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

# ── Save ─────────────────────────────────────────────────────────────────────────

np.savez_compressed(
    PATHS_CACHE,
    node_ids  = np.array(node_ids,       dtype=np.int32),
    od_src    = np.array(od_src_list,    dtype=np.int32),
    od_dst    = np.array(od_dst_list,    dtype=np.int32),
    od_dist   = np.array(od_dist_list,   dtype=np.float32),
    pair_idx  = np.array(pair_idx_list,  dtype=np.int32),
    link_idx  = np.array(link_idx_list,  dtype=np.int32),
    link_u    = np.array([u for u, v in link_list], dtype=np.int32),
    link_v    = np.array([v for u, v in link_list], dtype=np.int32),
)
size_mb = os.path.getsize(PATHS_CACHE) / 1e6
print(f"Saved: {PATHS_CACHE}  ({size_mb:.1f} MB)")
