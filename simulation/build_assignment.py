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

Tunable parameters: K, W_BIZ, P, ALPHA, COUNT_SITES (see Config section).
"""

import json, math, time, os
import numpy as np
import osmnx as ox
from collections import defaultdict

# ── Config ──────────────────────────────────────────────────────────────────────


K      = 1.73   # global flow scale factor
W_BIZ  = 1.0    # workplace demand weight relative to residential population
P      = 300.0  # peak travel time (seconds); flow peaks at d = P
ALPHA  = 2.0    # tail decay exponent; flow ~ 1/d^ALPHA for large d
OFFSCREEN_SPEED_MS = 80_000 / 3600   # 80 km/h — assumed speed for off-network boundary legs

# Traffic count sites used in goodness_of_fit().
# links: list of directed (u,v) pairs to sum for AADT; None → use cordon_flow(node).
# Bangor Road is a dual carriageway so we sum the two named directed links explicitly.
COUNT_SITES = [
    {"label": "site 507, A21 Bangor Road",    "node": 731, "links": [(731, 730), (730, 731)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road","node":  47, "links": None,                     "observed": 10_792},
    {"label": "site 444, A20 Portaferry Road","node":  92, "links": None,                     "observed":  7_282},
]

OUT_DIR        = "simulation"
WEIGHTS_FILE   = "simulation/node_weights.json"
TUNER_CONFIG   = "simulation/tuner_config.json"
CONS_GRAPH     = "simulation/newtownards_consolidated.graphml"
PATHS_CACHE    = "simulation/newtownards_paths.npz"
LINK_AADT_FILE = "data/link_aadt.json"

# ── Load node weights ────────────────────────────────────────────────────────────

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

TUNED_PARAMS = "simulation/tuned_params.json"
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
    od_dist      = cache["od_dist"].astype(np.float64)  # (P,) effective travel time (seconds)
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

    # Gravity flow for every OD pair — rational kernel: peak at d=P, tail ~ 1/d^ALPHA
    u    = od_dist / P
    t_ij = w_vec[od_src] * w_vec[od_dst] * (ALPHA + 1) * u / (ALPHA + u ** (ALPHA + 1))

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
            u    = dist / P
            t_ij = w_i * w_j * (ALPHA + 1) * u / (ALPHA + u ** (ALPHA + 1))
            for u, v in zip(path[:-1], path[1:]):
                lf[(u, v)] += t_ij
        if (idx + 1) % 100 == 0:
            print(f"  {idx + 1}/{len(node_ids)} nodes  ({time.time()-t0:.1f}s)")

    link_flow = dict(lf)
    print(f"  Assignment complete in {time.time()-t0:.1f}s  ({len(link_flow)} loaded links)")

# ── Street name lookup (from already-loaded graph) ───────────────────────────────

_link_name = {(int(u), int(v)): d["name"]
              for u, v, d in G.edges(data=True) if d.get("name")}


def _link_label(u, v):
    name = _link_name.get((u, v), "")
    return f"{u}→{v}  {name}" if name else f"{u}→{v}"

# ── Calibration ──────────────────────────────────────────────────────────────────

def cordon_flow(node):
    return sum(f for (u, v), f in link_flow.items() if u == node or v == node)

def site_raw_flow(site):
    if site["links"]:
        return sum(link_flow.get(lnk, 0) for lnk in site["links"])
    return cordon_flow(site["node"])

def goodness_of_fit():
    """Chi-squared goodness of fit against all count data, sorted by |z|."""
    rows = []
    chi2 = 0.0

    for s in COUNT_SITES:
        mod = site_raw_flow(s)
        obs = s["observed"]
        sig = 0.10 * obs
        z   = (mod - obs) / sig
        chi2 += z * z
        rows.append(("official", s["label"], obs, sig, mod, z))

    if os.path.exists(LINK_AADT_FILE):
        with open(LINK_AADT_FILE) as f:
            link_aadt = json.load(f)["links"]
        for key, entry in sorted(link_aadt.items()):
            u, v = map(int, key.split(","))
            mod  = link_flow.get((u, v), 0.0)
            obs  = entry["aadt"]
            sig  = entry["aadt_uncertainty"]
            z    = (mod - obs) / sig
            chi2 += z * z
            rows.append(("walking", _link_label(u, v), obs, sig, mod, z))

    rows.sort(key=lambda r: abs(r[5]), reverse=True)
    n = len(rows)
    chi2_per_n = chi2 / n

    LABEL_W = 52
    print(f"\nGoodness of fit  χ²={chi2:.2f}  n={n}  χ²/N={chi2_per_n:.4f}")
    print(f"  {'':1s}  {'Src':<8}  {'Label':<{LABEL_W}}  {'Obs':>8}  {'σ':>7}  {'Model':>8}  {'z':>6}")
    for kind, lbl, obs, sig, mod, z in rows:
        marker = "*" if abs(z) > 2 else " "
        print(f"  {marker} {kind:<8}  {lbl:<{LABEL_W}}  {obs:>8,.0f}  {sig:>7,.0f}  {mod:>8,.0f}  {z:>+.2f}")

    abs_z      = [abs(r[5]) for r in rows]
    mean_abs_z = sum(abs_z) / len(abs_z)
    n_out2     = sum(1 for a in abs_z if a > 2)
    n_out3     = sum(1 for a in abs_z if a > 3)
    print(f"\n  n={n}  χ²/N={chi2_per_n:.4f}  mean|z|={mean_abs_z:.2f}"
          f"  |z|>2: {n_out2}  |z|>3: {n_out3}")
    return chi2, n


link_flow = {k: v * K for k, v in link_flow.items()}

raw_flows = [site_raw_flow(s) for s in COUNT_SITES]
print(f"\nOfficial count sites  (K = {K})")
print(f"  {'Site':<45s}  {'Modelled':>9s}  {'Observed':>9s}  {'Ratio':>6s}")
for f, s in zip(raw_flows, COUNT_SITES):
    print(f"  {s['label']:<45s}  {f:>9,.0f}  {s['observed']:>9,}  {f/s['observed']:>6.2f}")

goodness_of_fit()

# ── Serialise flows ───────────────────────────────────────────────────────────────

flows_path = f"{OUT_DIR}/newtownards_flows.json"
with open(flows_path, "w") as f:
    json.dump({
        "kernel": "rational", "W_BIZ": W_BIZ, "P": P, "ALPHA": ALPHA, "K": K,
        "flows": {f"{u},{v}": flow for (u, v), flow in link_flow.items()},
    }, f)
print(f"\nSaved {len(link_flow)} link flows → {flows_path}")
print(f"Parameters: K={K}  W_BIZ={W_BIZ}  P={P}  ALPHA={ALPHA}")
