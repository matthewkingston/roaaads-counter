"""
Precompute all-pairs shortest-time paths for the Newtownards road network.

Run this once after any network change (build_network.py or edit_network.py).
Produces simulation/newtownards_paths.npz, which build_assignment.py loads
to avoid re-running Dijkstra on every parameter-tuning iteration.

od_dist stores the mean effective path distance (seconds) averaged across
N_PASSES stochastic Dijkstra runs with log-normal edge-cost perturbations
(CV = PROBIT_CV).  link_weight[entry] is the fraction of passes that routed
through that link for the corresponding OD pair.  Weights < 1 represent
genuine route diversity; weight = 1 means that link was always chosen.

Cache must be regenerated when:
  - The road network changes (newtownards_consolidated.graphml is updated)
  - External zone coordinates (lat/lon) change in build_demographics.py EXTERNAL_ZONES
  - HIGHWAY_COST_FACTOR values change
  - N_PASSES or PROBIT_CV change (affects route diversity)

Does NOT need regenerating for changes to W_BIZ, P, ALPHA, node
populations, damping factors, or count site values.
"""

import json, math, time, os
import numpy as np
import osmnx as ox
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

CONS_GRAPH         = "simulation/newtownards_consolidated.graphml"
WEIGHTS_FILE       = "simulation/node_weights.json"
TUNER_CONFIG       = "simulation/tuner_config.json"
PATHS_CACHE        = "simulation/newtownards_paths.npz"
OFFSCREEN_SPEED_MS = 80_000 / 3600   # 80 km/h in m/s
N_PASSES           = 25              # stochastic Dijkstra passes for probit loading
PROBIT_CV          = 0.25            # log-normal noise CV applied to edge costs
RANDOM_SEED        = 42

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

adj    = csr_matrix((data, (rows, cols)), shape=(n, n))
N_links = len(link_list)
print(f"  {N_links} directed links")

# ── Dense lookup tables for fast path tracing ────────────────────────────────────
# edge_to_link[i, j] = link index for directed edge (node_i → node_j), or -1
# edge_cost[i, j]    = original (unperturbed) generalised cost

print("Building edge lookup tables …")
edge_to_link = np.full((n, n), -1, dtype=np.int32)
edge_cost    = np.zeros((n, n), dtype=np.float32)

for (u, v), li in link_to_idx.items():
    i, j = node_to_idx[u], node_to_idx[v]
    edge_to_link[i, j] = li

for row_i, col_i, cost_i in zip(rows, cols, data):
    # Use minimum cost for (i,j) in case of parallel edges
    if edge_cost[row_i, col_i] == 0.0 or cost_i < edge_cost[row_i, col_i]:
        edge_cost[row_i, col_i] = float(cost_i)

# ── k=1: all-pairs Dijkstra ───────────────────────────────────────────────────────

print("Running k=1 Dijkstra …")
t0 = time.time()
dist_matrix, predecessors = dijkstra(adj, directed=True, return_predecessors=True)
print(f"  Done in {time.time()-t0:.1f}s")

# ── Build OD pair list from k=1 paths ────────────────────────────────────────────

print("Building OD pair list …")
t0 = time.time()

od_src_list  = []
od_dst_list  = []
od_k1_links  = []   # k=1 link-index list per pair (for fallback)
od_k1_dist   = []   # k=1 effective distance per pair (for fallback)

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
            li = edge_to_link[prev, cur]
            if li >= 0:
                path_links.append(li)
            cur = prev
        if not path_links:
            continue
        path_links.reverse()

        od_src_list.append(src_i)
        od_dst_list.append(dst_i)
        od_k1_links.append(path_links)
        od_k1_dist.append(eff_dist)

    if (src_i + 1) % 100 == 0:
        print(f"  {src_i + 1}/{n} nodes  ({time.time()-t0:.1f}s)  {len(od_src_list):,} pairs")

n_pairs = len(od_src_list)
print(f"  {n_pairs:,} OD pairs  in {time.time()-t0:.1f}s")

# ── Group OD pairs by source for vectorised stochastic tracing ───────────────────

src_groups = {}   # src_i → {'pks': [...], 'dsts': [...], 'src_off': float}
for k, (si, di) in enumerate(zip(od_src_list, od_dst_list)):
    if si not in src_groups:
        src_groups[si] = {
            'pks':     [],
            'dsts':    [],
            'src_off': boundary_offscreen.get(node_ids[si], 0.0),
        }
    src_groups[si]['pks'].append(k)
    src_groups[si]['dsts'].append(di)

# ── Stochastic probit passes ─────────────────────────────────────────────────────
# For each pass: perturb edge costs with log-normal noise, re-run Dijkstra,
# trace all OD paths using the perturbed predecessor matrix.  Accumulate:
#   hit_accum[pair_k * N_links + link_i] = number of passes using that link
#   dist_accum[pair_k]                   = sum of effective distances across passes
# Pairs whose perturbed path cannot be traced fall back to their k=1 path.

print(f"\nRunning {N_PASSES} stochastic probit passes (CV={PROBIT_CV}, seed={RANDOM_SEED}) …")
t0  = time.time()
rng = np.random.default_rng(seed=RANDOM_SEED)

data_arr = np.array(data, dtype=np.float64)
rows_arr = np.array(rows, dtype=np.int32)
cols_arr = np.array(cols, dtype=np.int32)

hit_accum      = {}                                   # packed int → count
dist_accum     = np.zeros(n_pairs, dtype=np.float64)
n_fallbacks    = 0

for pass_idx in range(N_PASSES):
    t_pass = time.time()

    eps    = rng.normal(0.0, PROBIT_CV, size=len(data_arr))
    data_p = data_arr * np.exp(eps)
    adj_p  = csr_matrix((data_p, (rows_arr, cols_arr)), shape=(n, n))
    _, pred_p = dijkstra(adj_p, directed=True, return_predecessors=True)

    for src_i, grp in src_groups.items():
        pair_ks = grp['pks']
        dsts    = grp['dsts']
        src_off = grp['src_off']
        n_s     = len(pair_ks)

        cur    = np.array(dsts, dtype=np.int32)
        done   = np.zeros(n_s, dtype=bool)
        fallbk = np.zeros(n_s, dtype=bool)
        pc     = np.zeros(n_s, dtype=np.float64)

        for _ in range(300):
            at_src = ~done & (cur == src_i)
            done  |= at_src
            if done.all():
                break

            active = ~done
            prev   = pred_p[src_i, cur]

            dead    = active & (prev < 0)
            fallbk |= dead
            done   |= dead

            active_ok = active & ~dead
            if not active_ok.any():
                break

            aok      = np.where(active_ok)[0]
            prev_aok = prev[aok]
            cur_aok  = cur[aok]
            li_aok   = edge_to_link[prev_aok, cur_aok]
            cost_aok = edge_cost[prev_aok, cur_aok]

            for j in range(len(aok)):
                m  = aok[j]
                li = int(li_aok[j])
                pc[m] += float(cost_aok[j])
                if li >= 0:
                    key = pair_ks[m] * N_links + li
                    hit_accum[key] = hit_accum.get(key, 0) + 1

            cur[active_ok] = prev[active_ok]

        for m in range(n_s):
            pk    = pair_ks[m]
            dst_i = dsts[m]
            if fallbk[m]:
                n_fallbacks += 1
                dist_accum[pk] += od_k1_dist[pk]
                for li in od_k1_links[pk]:
                    key = pk * N_links + li
                    hit_accum[key] = hit_accum.get(key, 0) + 1
            else:
                dst_off = boundary_offscreen.get(node_ids[dst_i], 0.0)
                dist_accum[pk] += pc[m] + src_off + dst_off

    print(f"  Pass {pass_idx + 1}/{N_PASSES} done in {time.time() - t_pass:.1f}s")

print(f"  All passes done in {time.time()-t0:.1f}s  ({n_fallbacks:,} pair-pass fallbacks)")

# ── Build output arrays ───────────────────────────────────────────────────────────

print("Building output arrays …")
od_dist_out = (dist_accum / N_PASSES).astype(np.float32)

pair_idx_out    = []
link_idx_out    = []
link_weight_out = []
for key, count in hit_accum.items():
    pk = key // N_links
    li = key %  N_links
    pair_idx_out.append(pk)
    link_idx_out.append(li)
    link_weight_out.append(count / N_PASSES)

pair_idx_arr    = np.array(pair_idx_out,    dtype=np.int32)
link_idx_arr    = np.array(link_idx_out,    dtype=np.int32)
link_weight_arr = np.array(link_weight_out, dtype=np.float32)

sort_order      = np.argsort(pair_idx_arr, kind='stable')
pair_idx_arr    = pair_idx_arr[sort_order]
link_idx_arr    = link_idx_arr[sort_order]
link_weight_arr = link_weight_arr[sort_order]

n_entries     = len(pair_idx_arr)
n_always_k1   = int((link_weight_arr == 1.0).sum())
mean_lpp      = n_entries / max(n_pairs, 1)
print(f"  {n_pairs:,} OD pairs  {n_entries:,} (pair,link) entries  mean {mean_lpp:.1f} links/pair")
print(f"  weight=1.0 entries: {n_always_k1:,} ({100*n_always_k1/max(n_entries,1):.1f}%)")

# ── Save ─────────────────────────────────────────────────────────────────────────

print("Saving cache …")
np.savez_compressed(
    PATHS_CACHE,
    node_ids        = np.array(node_ids,    dtype=np.int32),
    od_src          = np.array(od_src_list, dtype=np.int32),
    od_dst          = np.array(od_dst_list, dtype=np.int32),
    od_dist         = od_dist_out,
    pair_idx        = pair_idx_arr,
    link_idx        = link_idx_arr,
    link_weight     = link_weight_arr,
    link_u          = np.array([u for u, v in link_list], dtype=np.int32),
    link_v          = np.array([v for u, v in link_list], dtype=np.int32),
    probit_n_passes = np.int32(N_PASSES),
    probit_cv       = np.float32(PROBIT_CV),
)
size_mb = os.path.getsize(PATHS_CACHE) / 1e6
print(f"Saved: {PATHS_CACHE}  ({size_mb:.1f} MB)")
print(f"  {n_pairs:,} OD pairs  {n_entries:,} entries  {N_PASSES} passes  CV={PROBIT_CV}")
