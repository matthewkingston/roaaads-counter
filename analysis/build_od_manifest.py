"""
Build the fixed OD manifest for the Google routing-time calibration sample.

Writes a deterministic (seeded), model-aligned, length-skewed sample of OD pairs
to data/google_cache/od_manifest.json. Makes NO API calls and does not touch OSRM
— it only reads the model's own link/graph data, so it is safe and repeatable. The
runner (google_query_routes.py) consumes this manifest; the manifest is the single
source of truth for *which* routes get queried, so a crashed run resumes against the
same fixed list.

Model-aligned leg types (the routes the model's impedance actually uses):
  X2B  external centroid → boundary entry node   (external approaches; len = OSRM duration_s)
  B2X  boundary node → external centroid          (reverse)
  X2X  external → external through-routes          (allowlisted transit; len = haversine)
  INT  internal node → internal node               (in-town junction realism; len = haversine)

Within each leg type the sample is stratified into 4 length quartile-bands and
allocated with a deliberate skew toward longer legs (BAND_WEIGHTS), because longer
routes traverse more distinct road-class / junction mixes per paid query. Road-class
× speed-band coverage is treated as emergent here and should be checked after a first
batch (a future refinement could pre-route candidates with local OSRM to stratify on
the realised road-class mix directly).

Usage:
  python3 analysis/build_od_manifest.py            # write manifest (default ~1000)
  python3 analysis/build_od_manifest.py --n 400    # smaller sample
"""

import argparse, json, os, random, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from google_routing_common import haversine_m

import networkx as nx

REPO_ROOT     = "/home/matthew/Documents/CodingFun/roaaads"
CENSUS_FILE   = os.path.join(REPO_ROOT, "data", "census_zones.json")
LINKS_FILE    = os.path.join(REPO_ROOT, "data", "external_links.json")
WEIGHTS_FILE  = os.path.join(REPO_ROOT, "simulation", "node_weights.json")
RAW_GRAPH     = os.path.join(REPO_ROOT, "simulation", "newtownards_network.graphml")
MANIFEST_FILE = os.path.join(REPO_ROOT, "data", "google_cache", "od_manifest.json")

SEED = 20260622

# Leg-type quotas (fractions of total). Tunable; realised counts capped by availability.
LEG_FRACTIONS = {"X2B": 0.45, "X2X": 0.25, "B2X": 0.15, "INT": 0.15}

# Length-band allocation within each leg type — skewed toward longer legs.
# Bands are quartiles of the leg-type's length distribution (band 0 = shortest).
BAND_WEIGHTS = [0.15, 0.20, 0.30, 0.35]

# Size of the random internal-pair candidate pool to draw INT legs from.
INT_POOL = 8000


def band_allocate(candidates, key, quota, rng):
    """Split candidates into 4 length quartile-bands and draw `quota` items with
    BAND_WEIGHTS skew, spilling any shortfall to the next band. Returns chosen list."""
    if not candidates or quota <= 0:
        return []
    ordered = sorted(candidates, key=key)
    n = len(ordered)
    edges = [0, n // 4, n // 2, 3 * n // 4, n]
    bands = [ordered[edges[i]:edges[i + 1]] for i in range(4)]
    chosen, carry = [], 0
    for bi, band in enumerate(bands):
        want = round(quota * BAND_WEIGHTS[bi]) + carry
        want = min(want, len(band))
        if want > 0:
            picks = rng.sample(band, want)
            chosen.extend((bi, c) for c in picks)
        carry = round(quota * BAND_WEIGHTS[bi]) + carry - want  # unfilled spills forward
    # top up if rounding/availability left us short
    if len(chosen) < quota:
        chosenset = {id(c) for _, c in chosen}
        leftover = [c for c in ordered if id(c) not in chosenset]
        extra = rng.sample(leftover, min(quota - len(chosen), len(leftover)))
        # assign each leftover its true band index
        for c in extra:
            pos = ordered.index(c)
            bi = min(3, next(k for k in range(4) if pos < edges[k + 1]))
            chosen.append((bi, c))
    return chosen[:quota]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1000, help="target total OD pairs")
    ap.add_argument("--out", default=MANIFEST_FILE)
    args = ap.parse_args()
    rng = random.Random(SEED)

    print("Loading model data …")
    cz = json.load(open(CENSUS_FILE))
    centroids = {e["id"]: (e["centroid_lat"], e["centroid_lon"]) for e in cz["external_nodes"]}
    el = json.load(open(LINKS_FILE))
    nw = json.load(open(WEIGHTS_FILE))
    internal_ids = [str(i) for i in nw["internal_node_ids"]]

    G = nx.read_graphml(RAW_GRAPH)
    def nodexy(osmid):
        nd = G.nodes[str(osmid)]
        return (float(nd["y"]), float(nd["x"]))   # (lat, lon)

    # ── Build candidate pools per leg type ────────────────────────────────────
    # X2B: external centroid → boundary node; length = model OSRM duration_s
    x2b = [{"o_label": l["from_ext"], "o": centroids[l["from_ext"]],
            "d_label": str(l["to_boundary"]), "d": nodexy(l["to_boundary"]),
            "len_s": float(l["duration_s"])}
           for l in el["ext_boundary_links"] if l["from_ext"] in centroids]

    # B2X: boundary node → external centroid
    b2x = [{"o_label": str(l["from_boundary"]), "o": nodexy(l["from_boundary"]),
            "d_label": l["to_ext"], "d": centroids[l["to_ext"]],
            "len_s": float(l["duration_s"])}
           for l in el["bnd_external_links"] if l["to_ext"] in centroids]

    # X2X: allowlisted external→external through-routes; length = haversine (no model dur)
    x2x = []
    for src, dsts in el["allowed_through_pairs"].items():
        if src not in centroids:
            continue
        for dst in dsts:
            if dst in centroids:
                o, d = centroids[src], centroids[dst]
                x2x.append({"o_label": src, "o": o, "d_label": dst, "d": d,
                            "len_s": haversine_m(*o, *d)})  # metres as length proxy

    # INT: random internal→internal pairs; length = haversine
    int_pool = []
    seen = set()
    attempts = 0
    while len(int_pool) < INT_POOL and attempts < INT_POOL * 4:
        attempts += 1
        a, b = rng.sample(internal_ids, 2)
        key = (a, b) if a < b else (b, a)
        if key in seen:
            continue
        seen.add(key)
        oa, ob = nodexy(a), nodexy(b)
        int_pool.append({"o_label": a, "o": oa, "d_label": b, "d": ob,
                         "len_s": haversine_m(*oa, *ob)})

    pools = {"X2B": x2b, "B2X": b2x, "X2X": x2x, "INT": int_pool}
    print("Candidate pools:", {k: len(v) for k, v in pools.items()})

    # ── Stratified, length-skewed draw per leg type ──────────────────────────
    manifest, counters = [], {}
    for leg, frac in LEG_FRACTIONS.items():
        quota = round(args.n * frac)
        chosen = band_allocate(pools[leg], key=lambda c: c["len_s"], quota=quota, rng=rng)
        counters[leg] = {"quota": quota, "got": len(chosen), "bands": [0, 0, 0, 0]}
        for bi, c in chosen:
            counters[leg]["bands"][bi] += 1
            i = len(manifest)
            manifest.append({
                "od_id": f"{leg}_{i:05d}",
                "leg_type": leg, "len_band": bi, "len_s": round(c["len_s"], 1),
                "o": {"label": c["o_label"], "lat": c["o"][0], "lon": c["o"][1]},
                "d": {"label": c["d_label"], "lat": c["d"][0], "lon": c["d"][1]},
            })

    out = {
        "meta": {
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "seed": SEED, "n_target": args.n, "n_total": len(manifest),
            "time_basis": "free_flow",
            "leg_fractions": LEG_FRACTIONS, "band_weights": BAND_WEIGHTS,
            "centre": cz.get("centre"),
        },
        "od_pairs": manifest,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)

    print(f"\nWrote {args.out}  ({len(manifest)} OD pairs)")
    print("Realised per leg type (quota → got, band counts shortest→longest):")
    for leg, c in counters.items():
        print(f"  {leg}: {c['quota']:4d} → {c['got']:4d}   bands {c['bands']}")
    est = len(manifest) * 5.0 / 1000.0
    print(f"\nEst. live Google cost when run: ~${est:.2f} ({len(manifest)} queries)")
    print("Next: review the manifest, then run google_query_routes.py (needs approval + key).")


if __name__ == "__main__":
    main()
