"""
Gravity model OD matrix + all-or-nothing shortest-path assignment.

Loads node weights from simulation/node_weights.json (produced by
build_demographics.py) and assigns traffic to links using a gravity model.

Fast path (recommended): load precomputed paths from simulation/newtownards_paths.npz
(produced by build_paths.py). Reduces parameter-tuning runs to < 1s.

Slow fallback: if the paths cache is absent, runs NetworkX Dijkstra (~10s).
Re-run build_paths.py whenever the road network or external zone coordinates change.

Usage:
  python3 simulation/build_assignment.py           # full run with map
  python3 simulation/build_assignment.py --no-map  # skip map, < 1s with cache

Tunable parameters: W_BIZ, MU, SIGMA, ALPHA, COUNT_SITES (see Config section).
"""

import json, math, time, os, sys
import numpy as np
import osmnx as ox
import folium
from collections import defaultdict

# ── Config ──────────────────────────────────────────────────────────────────────

W_BIZ  = 1.0    # workplace demand weight relative to residential population
MU     = 9.6    # lognormal shift; peak at exp(MU − ALPHA×SIGMA²) = exp(7.6) ≈ 2000m
SIGMA  = 1.0    # lognormal spread in log-distance space
ALPHA  = 2.0    # power-law tail exponent (far-field decay, matches former 1/d² model)

# Traffic count sites used for least-squares calibration of K.
# links: list of directed (u,v) pairs to sum for AADT; None → use cordon_flow(node).
# Bangor Road is a dual carriageway so we sum the two named directed links explicitly.
COUNT_SITES = [
    {"label": "site 507, A21 Bangor Road",    "node": 731, "links": [(731, 730), (730, 731)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road","node":  47, "links": None,                     "observed": 10_792},
    {"label": "site 444, A20 Portaferry Road","node":  92, "links": None,                     "observed":  7_282},
]

OUT_DIR      = "simulation"
WEIGHTS_FILE = "simulation/node_weights.json"
CONS_GRAPH   = "simulation/newtownards_consolidated.graphml"
RAW_GRAPH    = "simulation/newtownards_network.graphml"
PATHS_CACHE  = "simulation/newtownards_paths.npz"

build_map = "--no-map" not in sys.argv

# ── Load node weights ────────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)

node_population      = {int(k): v for k, v in weights["node_population"].items()}
node_business_demand = {int(k): v for k, v in weights["node_business_demand"].items()}
node_effective_utm   = {int(k): (v[0], v[1]) for k, v in weights["node_effective_utm"].items()}
boundary_node_ids    = set(weights["boundary_node_ids"])

# ── Assignment ───────────────────────────────────────────────────────────────────

link_flow = {}   # (u, v) → raw (pre-K) flow

if os.path.exists(PATHS_CACHE):
    # ── Fast path: vectorised numpy using precomputed paths ──────────────────────
    print(f"Loading paths cache ({PATHS_CACHE}) …")
    t0 = time.time()
    cache = np.load(PATHS_CACHE)

    node_ids_arr = cache["node_ids"]          # int32 (N,)
    od_src       = cache["od_src"]            # int32 (P,)  source node index
    od_dst       = cache["od_dst"]            # int32 (P,)  target node index
    od_dist      = cache["od_dist"].astype(np.float64)  # (P,) effective distance
    pair_idx     = cache["pair_idx"]          # int32 (E,)  which OD pair
    link_idx     = cache["link_idx"]          # int32 (E,)  which link
    link_u       = cache["link_u"]            # int32 (L,)
    link_v       = cache["link_v"]            # int32 (L,)

    # Node weight vector in the same order as node_ids_arr
    w_vec = np.array([
        node_population.get(int(nid), 0) + W_BIZ * node_business_demand.get(int(nid), 0)
        for nid in node_ids_arr
    ], dtype=np.float64)

    total_weight = w_vec.sum()
    print(f"  {len(node_ids_arr)} nodes  total weight {total_weight:,.0f}  (W_BIZ={W_BIZ})")

    # Gravity flow for every OD pair — lognormal × power-law decay
    ln_d = np.log(od_dist)
    t_ij = w_vec[od_src] * w_vec[od_dst] * np.exp(-0.5 * ((ln_d - MU) / SIGMA) ** 2) / od_dist ** ALPHA

    # Accumulate onto links
    N_links = len(link_u)
    raw_flow_arr = np.bincount(link_idx, weights=t_ij[pair_idx], minlength=N_links)

    link_flow = {
        (int(link_u[k]), int(link_v[k])): raw_flow_arr[k]
        for k in range(N_links) if raw_flow_arr[k] > 0
    }
    print(f"  Assignment complete in {time.time()-t0:.2f}s  ({len(link_flow)} loaded links)")

    # Reconstruct node_ids list (needed for map / boundary markers)
    node_ids  = [int(nid) for nid in node_ids_arr]
    G = ox.load_graphml(CONS_GRAPH)
    node_weight = {int(nid): float(w) for nid, w in zip(node_ids_arr, w_vec)}

else:
    # ── Slow fallback: per-source NetworkX Dijkstra (~10s) ──────────────────────
    print("WARNING: paths cache not found — run build_paths.py for faster iterations")
    import networkx as nx

    print("Loading consolidated graph …")
    G = ox.load_graphml(CONS_GRAPH)
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
        boundary_offscreen[nid] = math.sqrt((nx_pos - cx) ** 2 + (ny_pos - cy) ** 2)

    for nid, extra in boundary_offscreen.items():
        print(f"  Node {nid:4d}  offscreen leg: {extra/1000:.1f} km")

    print("Running assignment …")
    t0 = time.time()
    lf = defaultdict(float)

    for idx, source in enumerate(node_ids):
        w_i = node_weight[source]
        if w_i <= 0:
            continue
        src_is_boundary = source in boundary_node_ids
        try:
            lengths, paths = nx.single_source_dijkstra(G, source, weight="length")
        except Exception:
            continue
        for target, path in paths.items():
            if target == source or len(path) < 2:
                continue
            if src_is_boundary and target in boundary_node_ids:
                continue
            w_j = node_weight[target]
            if w_j <= 0:
                continue
            dist = (lengths[target]
                    + boundary_offscreen.get(target, 0)
                    + boundary_offscreen.get(source, 0))
            if dist < 1.0:
                continue
            t_ij = w_i * w_j * math.exp(-0.5 * ((math.log(dist) - MU) / SIGMA) ** 2) / dist ** ALPHA
            for u, v in zip(path[:-1], path[1:]):
                lf[(u, v)] += t_ij
        if (idx + 1) % 100 == 0:
            print(f"  {idx + 1}/{len(node_ids)} nodes  ({time.time()-t0:.1f}s)")

    link_flow = dict(lf)
    print(f"  Assignment complete in {time.time()-t0:.1f}s  ({len(link_flow)} loaded links)")

# ── Calibration ──────────────────────────────────────────────────────────────────

def cordon_flow(node):
    return sum(f for (u, v), f in link_flow.items() if u == node or v == node)

def site_raw_flow(site):
    if site["links"]:
        return sum(link_flow.get(lnk, 0) for lnk in site["links"])
    return cordon_flow(site["node"])

raw_flows   = [site_raw_flow(s) for s in COUNT_SITES]
numerator   = sum(f * s["observed"] for f, s in zip(raw_flows, COUNT_SITES))
denominator = sum(f * f             for f      in raw_flows)
K = numerator / denominator if denominator > 0 else 1.0

print(f"\nLeast-squares calibration across {len(COUNT_SITES)} sites  →  K = {K:.6f}")
print(f"  {'Site':<45s}  {'Modelled':>9s}  {'Observed':>9s}  {'Ratio':>6s}")
for f, s in zip(raw_flows, COUNT_SITES):
    modelled = f * K
    print(f"  {s['label']:<45s}  {modelled:>9,.0f}  {s['observed']:>9,}  {modelled/s['observed']:>6.2f}")

# Apply K
link_flow = {k: v * K for k, v in link_flow.items()}

# ── Serialise flows ───────────────────────────────────────────────────────────────

flows_path = f"{OUT_DIR}/newtownards_flows.json"
with open(flows_path, "w") as f:
    json.dump({
        "W_BIZ": W_BIZ, "MU": MU, "SIGMA": SIGMA, "ALPHA": ALPHA, "K": K,
        "flows": {f"{u},{v}": flow for (u, v), flow in link_flow.items()},
    }, f)
print(f"\nSaved {len(link_flow)} link flows → {flows_path}")
print(f"Parameters: W_BIZ={W_BIZ}  MU={MU}  SIGMA={SIGMA}  ALPHA={ALPHA}  K={K:.6f}")

if not build_map:
    sys.exit(0)

# ── Build assignment map ──────────────────────────────────────────────────────────

print("\nBuilding assignment map …")
G_raw = ox.load_graphml(RAW_GRAPH)
CENTRE = (54.5933779, -5.6960935)
m = folium.Map(location=list(CENTRE), zoom_start=14, tiles="CartoDB positron")

all_flows = sorted(f for f in link_flow.values() if f > 0)
if all_flows:
    p10_idx = max(0, int(len(all_flows) * 0.10))
    p90_idx = min(len(all_flows) - 1, int(len(all_flows) * 0.90))
    log_min = math.log10(max(all_flows[p10_idx], 1))
    log_max = math.log10(max(all_flows[p90_idx], 1))
else:
    log_min, log_max = 0.0, 1.0

def flow_color(flow):
    if flow <= 0:
        return "#cccccc"
    t = (math.log10(max(flow, 1)) - log_min) / max(log_max - log_min, 1e-6)
    t = max(0.0, min(1.0, t))
    if t < 0.33:
        r = 0; g = int(180 * (t / 0.33)); b = int(200 * (1 - t / 0.33))
    elif t < 0.66:
        s = (t - 0.33) / 0.33; r = int(220 * s); g = 180; b = 0
    else:
        s = (t - 0.66) / 0.34; r = 220 + int(35 * s); g = int(180 * (1 - s)); b = 0
    return f"#{r:02x}{g:02x}{b:02x}"

def flow_weight(flow):
    if flow <= 0:
        return 1
    t = (math.log10(max(flow, 1)) - log_min) / max(log_max - log_min, 1e-6)
    return 1 + 7 * max(0.0, min(1.0, t))

import pyproj
transformer = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

flow_fg = folium.FeatureGroup(name="Estimated link flows (AADT)", show=True)
for u, v, data in G.edges(data=True):
    flow = link_flow.get((u, v), 0) + link_flow.get((v, u), 0)
    geom = data.get("geometry")
    if geom and hasattr(geom, "coords"):
        coords = [transformer.transform(x, y)[::-1] for x, y in geom.coords]
    else:
        ud, vd = G.nodes[u], G.nodes[v]
        lon_u, lat_u = transformer.transform(float(ud["x"]), float(ud["y"]))
        lon_v, lat_v = transformer.transform(float(vd["x"]), float(vd["y"]))
        coords = [(lat_u, lon_u), (lat_v, lon_v)]
    name = data.get("name", "")
    length = float(data.get("length", 0))
    folium.PolyLine(
        coords, color=flow_color(flow), weight=flow_weight(flow), opacity=0.85,
        tooltip=f"{name or 'link'} ({u}→{v})<br>est. AADT: {flow:,.0f}<br>length: {length:.0f}m",
    ).add_to(flow_fg)
flow_fg.add_to(m)

bn_fg = folium.FeatureGroup(name="Boundary nodes", show=True)
for node_id in boundary_node_ids:
    if node_id not in G.nodes:
        continue
    nd = G.nodes[node_id]
    lon, lat = transformer.transform(float(nd["x"]), float(nd["y"]))
    cf = cordon_flow(node_id) * K
    folium.RegularPolygonMarker(
        location=[lat, lon], number_of_sides=4, radius=7, rotation=45,
        color="#e05c00", fill=True, fill_color="#ff7c20", fill_opacity=0.9, weight=2,
        tooltip=f"<b>Node {node_id}</b> [boundary]<br>est. cordon AADT: {cf:,.0f}",
    ).add_to(bn_fg)
bn_fg.add_to(m)

folium.LayerControl(collapsed=False).add_to(m)
out_path = f"{OUT_DIR}/newtownards_assignment_map.html"
m.save(out_path)
print(f"Saved: {out_path}")
