"""
Gravity model OD matrix + all-or-nothing shortest-path assignment.

Loads node weights from model/node_weights.json (produced by build_demographics.py)
and the consolidated road network, then:
  1. Computes node weight w_i = population_i + W_BIZ × workplace_i
  2. Gravity OD: T[i][j] = w_i × w_j / d_ij^ALPHA  (straight-line centroid distance)
  3. Assigns T[i][j] to shortest directed path for every (i, j) pair
  4. Calibrates global scalar K so the A21 Bangor Road flow matches site 507 AADT
  5. Outputs newtownards_assignment_map.html with links coloured by flow

Tunable parameters: W_BIZ, ALPHA (see Config section below).
"""

import json, math, time
import networkx as nx
import osmnx as ox
import folium
from collections import defaultdict

# ── Config ──────────────────────────────────────────────────────────────────────

W_BIZ  = 1.0    # workplace demand weight relative to residential population
ALPHA  = 2.0    # gravity distance decay exponent

# Traffic count sites used for least-squares calibration of K.
# links: list of directed (u,v) pairs to sum for AADT; None → use cordon_flow(node).
# Bangor Road is a dual carriageway so we sum the two named directed links explicitly.
COUNT_SITES = [
    {"label": "site 507, A21 Bangor Road",    "node": 731, "links": [(731, 730), (730, 731)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road","node":  47, "links": None,                     "observed": 10_792},
    {"label": "site 444, A20 Portaferry Road","node":  92, "links": None,                     "observed":  7_282},
]

OUT_DIR        = "simulation"
WEIGHTS_FILE   = "simulation/node_weights.json"
CONS_GRAPH     = "simulation/newtownards_consolidated.graphml"
RAW_GRAPH      = "simulation/newtownards_network.graphml"

# ── Load node weights ────────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)

node_population      = {int(k): v for k, v in weights["node_population"].items()}
node_business_demand = {int(k): v for k, v in weights["node_business_demand"].items()}
node_effective_utm   = {int(k): (v[0], v[1]) for k, v in weights["node_effective_utm"].items()}
boundary_node_ids    = set(weights["boundary_node_ids"])

# ── Load graph ───────────────────────────────────────────────────────────────────

print("Loading consolidated graph …")
G = ox.load_graphml(CONS_GRAPH)
node_ids = list(G.nodes())
print(f"  {len(node_ids)} nodes  {G.number_of_edges()} edges")

# ── Node weights ─────────────────────────────────────────────────────────────────

node_weight = {}
for n in node_ids:
    pop = node_population.get(n, 0)
    biz = node_business_demand.get(n, 0)
    node_weight[n] = pop + W_BIZ * biz

total_weight = sum(node_weight.values())
print(f"  Total node weight: {total_weight:,.0f}  (W_BIZ={W_BIZ})")

# ── Precompute offscreen legs for boundary nodes ──────────────────────────────────
# For internal→boundary pairs, gravity distance = network distance to the boundary
# node + straight-line from that boundary node to its external zone centroid.
# This adds the "offscreen" leg that Dijkstra can't see, so Belfast feels ~15km
# away rather than just ~2km to the edge of the study area.
# Internal→internal pairs use pure network distance.
# Node 180 (local access, centroid = own position) has offscreen leg = 0.

boundary_offscreen = {}   # node_id → extra straight-line metres to external centroid
for n in boundary_node_ids:
    if n not in G.nodes:
        continue
    nd = G.nodes[n]
    nx_pos, ny_pos = float(nd["x"]), float(nd["y"])
    cx, cy = node_effective_utm[n]
    boundary_offscreen[n] = math.sqrt((nx_pos - cx) ** 2 + (ny_pos - cy) ** 2)

for n, extra in boundary_offscreen.items():
    print(f"  Node {n:4d}  offscreen leg: {extra/1000:.1f} km")

# ── All-or-nothing assignment ────────────────────────────────────────────────────
# Run one Dijkstra per source node (shortest path by edge length).
# For each (source, target) pair with nonzero weights, compute gravity flow
# T[i][j] = w_i × w_j / d_ij^ALPHA and accumulate onto every directed link
# in the shortest path.
#
# Boundary→boundary pairs are skipped: those trips don't traverse the study area.

print("Running assignment …")
t0 = time.time()

link_flow = defaultdict(float)   # (u, v) → accumulated daily flow

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
            continue   # skip external-to-external pairs

        w_j = node_weight[target]
        if w_j <= 0:
            continue

        # Add offscreen leg for both source and target boundary nodes so that
        # T[i,97] and T[97,i] use the same effective distance and are symmetric.
        dist = (lengths[target]
                + boundary_offscreen.get(target, 0)
                + boundary_offscreen.get(source, 0))
        if dist < 1.0:
            continue

        t_ij = w_i * w_j / (dist ** ALPHA)

        for u, v in zip(path[:-1], path[1:]):
            link_flow[(u, v)] += t_ij

    if (idx + 1) % 100 == 0:
        print(f"  {idx + 1}/{len(node_ids)} nodes  ({time.time()-t0:.1f}s)")

print(f"  Assignment complete in {time.time()-t0:.1f}s  ({len(link_flow)} loaded links)")

# ── Calibration ──────────────────────────────────────────────────────────────────
# Sum flows on all directed links incident to the calibration boundary node.
# The total represents all traffic crossing the cordon at that point.
# Divide by 2 because each trip contributes once to the inbound link and
# once to the outbound link (both directions are modelled), so raw sum = 2×AADT.

def cordon_flow(node):
    # Sum both directed links incident to the boundary node.
    # This equals the two-way AADT (what a roadside counter measures).
    # No /2 — inbound + outbound = total daily crossing count.
    total = 0.0
    for (u, v), f in link_flow.items():
        if u == node or v == node:
            total += f
    return total

# Least-squares calibration: K = Σ(f_i · o_i) / Σ(f_i²)
# This minimises Σ(K·f_i − o_i)² across all count sites simultaneously.
def site_raw_flow(site):
    if site["links"]:
        return sum(link_flow.get(lnk, 0) for lnk in site["links"])
    return cordon_flow(site["node"])

raw_flows = [site_raw_flow(s) for s in COUNT_SITES]
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
        "W_BIZ": W_BIZ, "ALPHA": ALPHA, "K": K,
        "flows": {f"{u},{v}": flow for (u, v), flow in link_flow.items()},
    }, f)
print(f"\nSaved {len(link_flow)} link flows → {flows_path}")

# ── Build assignment map ──────────────────────────────────────────────────────────

print("Building assignment map …")
G_raw = ox.load_graphml(RAW_GRAPH)
CENTRE = (54.5933779, -5.6960935)
m = folium.Map(location=list(CENTRE), zoom_start=14, tiles="CartoDB positron")

# Flow colour scale: log scale anchored to 10th/90th percentile of non-zero flows.
# Using global min/max causes the skewed all-or-nothing distribution to compress
# most links into the top of the range (all red). Percentile anchoring spreads
# visual contrast across the busy-to-very-busy range where it matters.
import statistics as _stats
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
    # blue → green → yellow → red
    if t < 0.33:
        r = 0
        g = int(180 * (t / 0.33))
        b = int(200 * (1 - t / 0.33))
    elif t < 0.66:
        s = (t - 0.33) / 0.33
        r = int(220 * s)
        g = 180
        b = 0
    else:
        s = (t - 0.66) / 0.34
        r = 220 + int(35 * s)
        g = int(180 * (1 - s))
        b = 0
    return f"#{r:02x}{g:02x}{b:02x}"

def flow_weight(flow):
    if flow <= 0:
        return 1
    t = (math.log10(max(flow, 1)) - log_min) / max(log_max - log_min, 1e-6)
    t = max(0.0, min(1.0, t))
    return 1 + 7 * t

# Build a lookup from raw graph edges to consolidated node flows.
# Raw graph uses OSM node IDs; consolidated graph uses integer cluster IDs.
# We map raw edges by spatial proximity: find the consolidated link whose
# endpoints are closest to the raw edge endpoints.
# Simpler approach: use the consolidated graph geometry directly.

flow_fg = folium.FeatureGroup(name="Estimated link flows (AADT)", show=True)
for u, v, data in G.edges(data=True):
    flow = link_flow.get((u, v), 0) + link_flow.get((v, u), 0)
    geom = data.get("geometry")
    if geom and hasattr(geom, "coords"):
        import pyproj
        transformer = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)
        coords = [transformer.transform(x, y)[::-1] for x, y in geom.coords]
    else:
        ud = G.nodes[u]
        vd = G.nodes[v]
        lon_u, lat_u = transformer.transform(float(ud["x"]), float(ud["y"]))
        lon_v, lat_v = transformer.transform(float(vd["x"]), float(vd["y"]))
        coords = [(lat_u, lon_u), (lat_v, lon_v)]

    name = data.get("name", "")
    length = float(data.get("length", 0))
    folium.PolyLine(
        coords,
        color=flow_color(flow),
        weight=flow_weight(flow),
        opacity=0.85,
        tooltip=f"{name or 'link'} ({u}→{v})<br>est. AADT: {flow:,.0f}<br>length: {length:.0f}m",
    ).add_to(flow_fg)
flow_fg.add_to(m)

# Boundary node markers
import pyproj as _pyproj
_tr = _pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)
import math as _math2

bn_fg = folium.FeatureGroup(name="Boundary nodes", show=True)
for node_id in boundary_node_ids:
    if node_id not in G.nodes:
        continue
    nd = G.nodes[node_id]
    lon, lat = _tr.transform(float(nd["x"]), float(nd["y"]))
    cf = cordon_flow(node_id) * K
    folium.RegularPolygonMarker(
        location=[lat, lon],
        number_of_sides=4, radius=7, rotation=45,
        color="#e05c00", fill=True, fill_color="#ff7c20", fill_opacity=0.9, weight=2,
        tooltip=f"<b>Node {node_id}</b> [boundary]<br>est. cordon AADT: {cf:,.0f}",
    ).add_to(bn_fg)
bn_fg.add_to(m)

folium.LayerControl(collapsed=False).add_to(m)

out_path = f"{OUT_DIR}/newtownards_assignment_map.html"
m.save(out_path)
print(f"Saved: {out_path}")
print(f"\nParameters: W_BIZ={W_BIZ}  ALPHA={ALPHA}  K={K:.6f}")
