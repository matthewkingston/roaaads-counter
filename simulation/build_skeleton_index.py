"""
Build the profile-independent *skeleton* cache from the Google route cache —
the one-time step that makes profile benchmarking fast.

For every cached Google route, we map-match its geometry once and record a
representation that does NOT depend on any speed profile:

  * length_by_bucket : metres of matched road in each (class x speed-band) bucket
  * turns            : per-maneuver (angle, junction degree, u-turn flag)
  * n_signals        : traffic-signal nodes on the matched path
  * g_dur            : Google's free-flow duration (the calibration target)

Bucket labels come from an **exact OSM node-id lookup**, not the old probe
profile: build_edge_index.py caches the full tag dict of every way/node along the
matched routes, so each matched segment (node-pair) is resolved to its way's
`highway`/`maxspeed` tags and bucketed via profile_spec.bucket_of. This replaces
the probe's `annotation.speed` readout, which was `distance / round(duration,
0.1s)` and corrupted short-edge buckets (see memory project_google_routing_-
calibration, INT diagnosis). Segment lengths come from the `distance` annotation
(geometry, profile-independent); turns from parsed steps; signals from the
cached node tags.

Run order (build_edge_index first):
  python3 simulation/build_edge_index.py            # matched_nodes.json + osm_ways/osm_nodes.jsonl
  python3 simulation/build_skeleton_index.py        # -> skeletons.jsonl   (real OSRM :5000)
  python3 simulation/build_skeleton_index.py --base-speeds   # -> base_speeds.json

All /match calls hit the real deployed OSRM (:5000) — no probe instance, no
:5001. The map-match is geometry-driven, so skeletons are (near) profile-
independent and only need rebuilding if the cached Google routes change.
"""

import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import profile_spec as ps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "analysis"))
from google_routing_common import (        # noqa: E402
    CONF_MIN, decode_polyline, downsample_by_distance, osrm_match_detail)

REPO_ROOT = "/home/matthew/Documents/CodingFun/roaaads"

CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
RAW_DIR = os.path.join(CACHE_DIR, "raw")
SKELETONS_FILE = os.path.join(CACHE_DIR, "skeletons.jsonl")
BASE_SPEEDS_FILE = os.path.join(CACHE_DIR, "base_speeds.json")
MANIFEST_FILE = os.path.join(CACHE_DIR, "od_manifest.json")
MATCH_CACHE_FILE = os.path.join(CACHE_DIR, "match_cache.jsonl")
OSM_WAYS_FILE = os.path.join(CACHE_DIR, "osm_ways.jsonl")
OSM_NODES_FILE = os.path.join(CACHE_DIR, "osm_nodes.jsonl")

COVERAGE_LO, COVERAGE_HI = 0.8, 1.2       # matched_dist / google_dist sanity band


# ── Edge-index lookup (replaces the probe decode) ────────────────────────────

def _read_match_cache():
    """{(od_id, route_idx): record} from match_cache.jsonl (the single OSRM pass
    cached by build_edge_index.py), tolerant of a partial trailing line."""
    recs = {}
    if not os.path.exists(MATCH_CACHE_FILE):
        return recs
    with open(MATCH_CACHE_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
                recs[(r["od_id"], r["route_idx"])] = r
            except (json.JSONDecodeError, KeyError):
                continue
    return recs


def load_edge_index(match_cache):
    """Load the raw OSM cache from build_edge_index.py into fast lookups:

      seg_tags : (node_u, node_v) -> way tag dict, both directions, restricted to
                 segments that actually appear in a matched route (so the map
                 stays small).
      signals  : set of node ids tagged highway=traffic_signals.

    Returns (seg_tags, signals).
    """
    for f in (OSM_WAYS_FILE, OSM_NODES_FILE):
        if not os.path.exists(f):
            sys.exit(f"ERROR: {f} not found — run build_edge_index.py --extract first.")

    matched_pairs = set()
    for r in match_cache.values():
        nodes = r["nodes"]
        for i in range(len(nodes) - 1):
            a, b = nodes[i], nodes[i + 1]
            if a != b:
                matched_pairs.add((a, b) if a < b else (b, a))

    seg_tags = {}
    with open(OSM_WAYS_FILE) as wf:
        for line in wf:
            w = json.loads(line)
            tags, nds = w["tags"], w["nodes"]
            for i in range(len(nds) - 1):
                a, b = nds[i], nds[i + 1]
                key = (a, b) if a < b else (b, a)
                if key in matched_pairs:
                    seg_tags[(a, b)] = tags
                    seg_tags[(b, a)] = tags

    signals = set()
    with open(OSM_NODES_FILE) as nf:
        for line in nf:
            n = json.loads(line)
            if n["tags"].get("highway") == "traffic_signals":
                signals.add(n["id"])

    print(f"Edge index: {len(seg_tags)//2} matched segments resolved, "
          f"{len(signals)} signal nodes")
    return seg_tags, signals


def _bucket_lengths(nodes, distances, seg_tags):
    """Aggregate matched-segment lengths into (class|band) buckets via the exact
    node-id -> way-tag lookup. nodes has one more entry than distances per leg;
    segment i runs nodes[i]->nodes[i+1] and carries distances[i]."""
    out = {}
    unbucketed = 0.0
    n = min(len(nodes) - 1, len(distances))
    for i in range(n):
        d = distances[i]
        if not d:
            continue
        tags = seg_tags.get((nodes[i], nodes[i + 1]))
        if tags is None:
            unbucketed += d
            continue
        key = ps.bucket_key(*ps.bucket_of(tags))
        out[key] = out.get(key, 0.0) + d
    return out, unbucketed


# ── Phase: map-match the cache -> skeletons.jsonl ────────────────────────────

def _turns(maneuvers):
    turns = []
    for mv in maneuvers:
        ang = mv.get("angle")
        turns.append({
            "angle": ang,
            "degree": mv.get("degree", 0),
            "uturn": mv.get("type") == "uturn" or mv.get("modifier") == "uturn",
        })
    return turns


def build_skeletons():
    """Rebuild skeletons.jsonl from the cached matches + the edge index. No OSRM
    calls — pure recompute, so it is free to re-run after any profile_spec/bucket
    change. Always regenerates from scratch (the slow work is the cached match)."""
    match_cache = _read_match_cache()
    if not match_cache:
        sys.exit(f"ERROR: {MATCH_CACHE_FILE} empty/missing — run build_edge_index.py --match.")
    seg_tags, signals = load_edge_index(match_cache)
    print(f"Building skeletons from {len(match_cache)} cached matched routes …")

    recs = []
    for (od_id, ridx), r in match_cache.items():
        nodes, dist = r["nodes"], r["distances"]
        gdist = r.get("g_dist", 0) or 0
        lbb, unbk = _bucket_lengths(nodes, dist, seg_tags)
        matched_m = sum(lbb.values()) + unbk
        coverage = (matched_m / gdist) if gdist else None
        n_sig = sum(1 for nid in nodes if nid in signals)
        valid = (r.get("conf", 0) >= CONF_MIN and coverage is not None
                 and COVERAGE_LO <= coverage <= COVERAGE_HI)
        recs.append({
            "od_id": od_id, "route_idx": ridx,
            "leg_type": r["leg_type"], "len_band": r["len_band"],
            "g_dur": r["g_dur"], "g_dist": gdist,
            "conf": r.get("conf"),
            "coverage": round(coverage, 4) if coverage is not None else None,
            "valid": valid,
            "length_by_bucket": {kk: round(vv, 1) for kk, vv in lbb.items()},
            "unbucketed_m": round(unbk, 1),
            "n_signals": n_sig,
            "turns": _turns(r["maneuvers"]),
        })

    tmp = SKELETONS_FILE + ".tmp"
    with open(tmp, "w") as f:
        for rec in recs:
            f.write(json.dumps(rec) + "\n")
    os.replace(tmp, SKELETONS_FILE)

    n_valid = sum(1 for r in recs if r["valid"])
    tot_unbk = sum(r["unbucketed_m"] for r in recs)
    tot_len = sum(sum(r["length_by_bucket"].values()) + r["unbucketed_m"] for r in recs)
    print(f"\nWrote {SKELETONS_FILE}")
    print(f"  {len(recs)} routes, {n_valid} valid "
          f"(conf>={CONF_MIN} & coverage in [{COVERAGE_LO},{COVERAGE_HI}])")
    if tot_len:
        print(f"  unbucketed length: {100.0 * tot_unbk / tot_len:.2f}% "
              f"(segments not resolved by the edge index)")


def estimate_base_speeds(speed_url, defactor, sample_n, min_obs):
    """Measure OSRM's realised factor-free base speed per (class x band) bucket.

    Per segment we need its bucket (exact node-id lookup, edge index) and OSRM's
    *realised* speed (annotation.speed on the speed source):
      * factor-free STOCK car.lua -> annotation.speed is the base speed directly
        (--no-defactor), or
      * deployed legacy :5000 -> annotation.speed has the class factor baked in;
        recover the factor-free value as speed * HIGHWAY_COST_FACTOR[class]
        (--defactor, the default).

    Aggregated **length-weighted** (harmonic mean weighted by segment length:
    base = Σ len / Σ(len/speed)) so short-edge speed noise is drowned out rather
    than equally weighted. Writes {bucket_key: km/h} to base_speeds.json.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from routing_config import HIGHWAY_COST_FACTOR

    seg_tags, _signals = load_edge_index(_read_match_cache())

    ods = json.load(open(MANIFEST_FILE))["od_pairs"]
    cached = [o for o in ods if os.path.exists(os.path.join(RAW_DIR, f"{o['od_id']}.json"))]
    stride = max(1, len(cached) // sample_n)
    sample = cached[::stride][:sample_n]
    print(f"Estimating base speeds from {len(sample)} sampled routes  "
          f"(speed={speed_url}, defactor={defactor})")

    # bucket_key -> [sum_len, sum_len_over_speed]
    acc = collections.defaultdict(lambda: [0.0, 0.0])
    n_routes = n_seg = n_skip = 0
    for o in sample:
        raw = json.load(open(os.path.join(RAW_DIR, f"{o['od_id']}.json")))
        for route in raw.get("routes", [])[:1]:
            coords = downsample_by_distance(
                decode_polyline(route["polyline"]["encodedPolyline"]))
            det = osrm_match_detail(speed_url, coords)
            if det is None:
                n_skip += 1
                continue
            nodes, dist, spd = det["nodes"], det["distances"], det["speeds"]
            m = min(len(nodes) - 1, len(dist), len(spd))
            for i in range(m):
                tags = seg_tags.get((nodes[i], nodes[i + 1]))
                d = dist[i]
                v = spd[i]
                if tags is None or not d or v is None or v <= 0:
                    continue
                cls, band = ps.bucket_of(tags)
                v_kmh = v * 3.6
                if defactor:
                    v_kmh *= HIGHWAY_COST_FACTOR.get(cls, 1.0)
                if v_kmh <= 0:
                    continue
                a = acc[ps.bucket_key(cls, band)]
                a[0] += d
                a[1] += d / v_kmh
                n_seg += 1
            n_routes += 1

    base = {k: round(s_len / s_lov, 2)
            for k, (s_len, s_lov) in acc.items() if s_lov > 0 and s_len > 0}
    # require a minimum number of contributing segments (tracked separately)
    with open(BASE_SPEEDS_FILE, "w") as f:
        json.dump(base, f, indent=2, sort_keys=True)
    print(f"\nWrote {BASE_SPEEDS_FILE}")
    print(f"  {n_routes} routes, {n_seg} segment samples, {n_skip} skipped")
    print(f"  {len(base)} buckets covered")
    top = sorted(acc.items(), key=lambda kv: -kv[1][0])[:10]
    print("  bucket                 emp km/h  analytic   km")
    for k, (s_len, s_lov) in top:
        cls, band = k.split("|", 1)
        emp = s_len / s_lov if s_lov else 0
        print(f"  {k:22} {emp:7.1f}  {ps.base_speed_for(cls, band):7.1f}  "
              f"{s_len/1000:6.1f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-speeds", action="store_true",
                    help="measure empirical per-bucket base speeds (single match "
                         "on the speed source)")
    ap.add_argument("--speed-url", default="http://localhost:5000",
                    help="real-speed OSRM for --base-speeds (deployed legacy :5000, "
                         "or a factor-free stock instance with --no-defactor)")
    ap.add_argument("--no-defactor", dest="defactor", action="store_false",
                    help="speed-url is already factor-free (stock); skip the "
                         "HIGHWAY_COST_FACTOR division")
    ap.add_argument("--sample", type=int, default=800,
                    help="routes to sample for --base-speeds")
    ap.add_argument("--min-obs", type=int, default=5,
                    help="(reserved) min segment samples to trust a bucket")
    ap.set_defaults(defactor=True)
    args = ap.parse_args()

    if args.base_speeds:
        estimate_base_speeds(args.speed_url, args.defactor, args.sample, args.min_obs)
    else:
        build_skeletons()


if __name__ == "__main__":
    main()
