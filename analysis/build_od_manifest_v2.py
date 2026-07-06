"""
Build a second batch of 1000 OD pairs for Google routing calibration.

Excludes all (origin_label, dest_label) pairs already in od_manifest.json so no
query is repeated across the two batches.  Outputs od_manifest_v2.json with
od_ids prefixed "v2_" to avoid filename collisions in data/google_cache/raw/.

Leg-type distribution is skewed toward B2X relative to v1 (B2X had the highest
observed violation rate in v1: ~8% vs ~4% for X2B).

  python3 analysis/build_od_manifest_v2.py            # write manifest_v2 (default 1000)
  python3 analysis/build_od_manifest_v2.py --n 500    # smaller batch

After writing, the user runs:
  GOOGLE_MAPS_API_KEY=... python3 analysis/google_query_routes.py \\
      --manifest data/google_cache/od_manifest_v2.json

Then match (OSRM must be running):
  python3 simulation/build_edge_index.py --match \\
      --manifest data/google_cache/od_manifest_v2.json

Then rebuild skeletons (picks up both v1 and v2 entries in match_cache.jsonl):
  python3 simulation/build_skeleton_index.py
"""

import argparse, json, os, random, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from google_routing_common import haversine_m

import networkx as nx

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CENSUS_FILE    = os.path.join(REPO_ROOT, "data", "census_zones.json")
LINKS_FILE     = os.path.join(REPO_ROOT, "data", "external_links.json")
WEIGHTS_FILE   = os.path.join(REPO_ROOT, "simulation", "node_weights.json")
RAW_GRAPH      = os.path.join(REPO_ROOT, "simulation", "newtownards_network.graphml")
V1_MANIFEST    = os.path.join(REPO_ROOT, "data", "google_cache", "od_manifest.json")
OUT_MANIFEST   = os.path.join(REPO_ROOT, "data", "google_cache", "od_manifest_v2.json")

SEED = 20260626   # different from v1 seed (20260622)

# B2X weighted up relative to v1 (observed ~8% violation rate vs ~4% for X2B).
# INT reduced (lowest violation rate, ~1%).
LEG_FRACTIONS = {"X2B": 0.40, "B2X": 0.25, "X2X": 0.25, "INT": 0.10}

BAND_WEIGHTS = [0.15, 0.20, 0.30, 0.35]
INT_POOL     = 8000


def band_allocate(candidates, key, quota, rng):
    """Same stratified length-band allocation as build_od_manifest.py."""
    if not candidates or quota <= 0:
        return []
    ordered = sorted(candidates, key=key)
    n = len(ordered)
    edges = [0, n // 4, n // 2, 3 * n // 4, n]
    bands  = [ordered[edges[i]:edges[i + 1]] for i in range(4)]
    chosen, carry = [], 0
    for bi, band in enumerate(bands):
        want = round(quota * BAND_WEIGHTS[bi]) + carry
        want = min(want, len(band))
        if want > 0:
            picks = rng.sample(band, want)
            chosen.extend((bi, c) for c in picks)
        carry = round(quota * BAND_WEIGHTS[bi]) + carry - want
    if len(chosen) < quota:
        chosen_set = {id(c) for _, c in chosen}
        leftover   = [c for c in ordered if id(c) not in chosen_set]
        extra      = rng.sample(leftover, min(quota - len(chosen), len(leftover)))
        for c in extra:
            pos = ordered.index(c)
            bi  = min(3, next(k for k in range(4) if pos < edges[k + 1]))
            chosen.append((bi, c))
    return chosen[:quota]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n",   type=int, default=1000, help="target total OD pairs")
    ap.add_argument("--out", default=OUT_MANIFEST)
    args = ap.parse_args()
    rng  = random.Random(SEED)

    # ── Load exclusion set from v1 manifest ──────────────────────────────────
    if not os.path.exists(V1_MANIFEST):
        sys.exit(f"ERROR: {V1_MANIFEST} not found — run build_od_manifest.py first.")
    v1 = json.load(open(V1_MANIFEST))
    existing = {(p["o"]["label"], p["d"]["label"]) for p in v1["od_pairs"]}
    print(f"Loaded v1 manifest: {len(v1['od_pairs'])} pairs → {len(existing)} exclusions")

    # ── Load model data ───────────────────────────────────────────────────────
    print("Loading model data …")
    cz        = json.load(open(CENSUS_FILE))
    centroids = {e["id"]: (e["centroid_lat"], e["centroid_lon"])
                 for e in cz["external_nodes"]}
    el        = json.load(open(LINKS_FILE))
    nw        = json.load(open(WEIGHTS_FILE))
    internal_ids = [str(i) for i in nw["internal_node_ids"]]

    G = nx.read_graphml(RAW_GRAPH)
    def nodexy(osmid):
        nd = G.nodes[str(osmid)]
        return (float(nd["y"]), float(nd["x"]))

    # ── Build candidate pools, excluding v1 pairs ─────────────────────────────
    def fresh(pool):
        return [c for c in pool
                if (c["o_label"], c["d_label"]) not in existing]

    x2b_all = [
        {"o_label": l["from_ext"],      "o": centroids[l["from_ext"]],
         "d_label": str(l["to_boundary"]), "d": nodexy(l["to_boundary"]),
         "len_s":   float(l["duration_s"])}
        for l in el["ext_boundary_links"] if l["from_ext"] in centroids
    ]
    b2x_all = [
        {"o_label": str(l["from_boundary"]), "o": nodexy(l["from_boundary"]),
         "d_label": l["to_ext"],              "d": centroids[l["to_ext"]],
         "len_s":   float(l["duration_s"])}
        for l in el["bnd_external_links"] if l["to_ext"] in centroids
    ]
    x2x_all = [
        {"o_label": src, "o": centroids[src],
         "d_label": dst, "d": centroids[dst],
         "len_s":   haversine_m(*centroids[src], *centroids[dst])}
        for src, dsts in el["allowed_through_pairs"].items()
        if src in centroids
        for dst in dsts if dst in centroids
    ]

    # INT: fresh random pairs from internal nodes (never in v1 by construction
    # since the v1 INT pool was drawn from a different seed with the same universe).
    int_pool = []
    seen = set()
    attempts = 0
    while len(int_pool) < INT_POOL and attempts < INT_POOL * 4:
        attempts += 1
        a, b = rng.sample(internal_ids, 2)
        key  = (a, b) if a < b else (b, a)
        if key in seen or (a, b) in existing:
            continue
        seen.add(key)
        oa, ob = nodexy(a), nodexy(b)
        int_pool.append({"o_label": a, "o": oa, "d_label": b, "d": ob,
                         "len_s": haversine_m(*oa, *ob)})

    pools = {
        "X2B": fresh(x2b_all),
        "B2X": fresh(b2x_all),
        "X2X": fresh(x2x_all),
        "INT": int_pool,
    }
    print("Fresh candidate pools (after v1 exclusions):",
          {k: len(v) for k, v in pools.items()})

    # ── Stratified, length-skewed draw ────────────────────────────────────────
    manifest, counters = [], {}
    for leg, frac in LEG_FRACTIONS.items():
        quota  = round(args.n * frac)
        chosen = band_allocate(pools[leg], key=lambda c: c["len_s"],
                               quota=quota, rng=rng)
        counters[leg] = {"quota": quota, "got": len(chosen), "bands": [0, 0, 0, 0]}
        for bi, c in chosen:
            counters[leg]["bands"][bi] += 1
            i = len(manifest)
            manifest.append({
                "od_id":    f"v2_{leg}_{i:05d}",
                "leg_type": leg, "len_band": bi, "len_s": round(c["len_s"], 1),
                "o": {"label": c["o_label"], "lat": c["o"][0], "lon": c["o"][1]},
                "d": {"label": c["d_label"], "lat": c["d"][0], "lon": c["d"][1]},
            })

    # Sanity: no od_id collision with v1
    v1_ids = {p["od_id"] for p in v1["od_pairs"]}
    v2_ids = {p["od_id"] for p in manifest}
    assert not (v1_ids & v2_ids), "od_id collision with v1 — should not happen"

    out = {
        "meta": {
            "created":       time.strftime("%Y-%m-%dT%H:%M:%S"),
            "seed":          SEED,
            "n_target":      args.n,
            "n_total":       len(manifest),
            "time_basis":    "free_flow",
            "leg_fractions": LEG_FRACTIONS,
            "band_weights":  BAND_WEIGHTS,
            "centre":        cz.get("centre"),
            "excludes":      V1_MANIFEST,
            "n_excluded":    len(existing),
        },
        "od_pairs": manifest,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)

    print(f"\nWrote {args.out}  ({len(manifest)} OD pairs, all fresh vs v1)")
    print("Realised per leg type (quota → got, band counts shortest→longest):")
    for leg, c in counters.items():
        print(f"  {leg}: {c['quota']:4d} → {c['got']:4d}   bands {c['bands']}")
    est = len(manifest) * 5.0 / 1000.0
    print(f"\nEst. live Google cost: ~${est:.2f} ({len(manifest)} queries)")
    print("\nNext steps:")
    print("  1. User runs:")
    print(f"       GOOGLE_MAPS_API_KEY=... python3 analysis/google_query_routes.py \\")
    print(f"           --manifest {args.out}")
    print("  2. Match (OSRM must be running):")
    print(f"       python3 simulation/build_edge_index.py --match \\")
    print(f"           --manifest {args.out}")
    print("  3. Rebuild skeletons (picks up v1 + v2 entries):")
    print("       python3 simulation/build_skeleton_index.py")


if __name__ == "__main__":
    main()
