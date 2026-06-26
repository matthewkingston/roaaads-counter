"""
Build the raw OSM edge/node index for the Google route set.

This is the rounding-free replacement for the OSRM *probe* profile. The probe
encoded each way's bucket id as its speed and read it back from /match's `speed`
annotation, but that annotation is `distance / round(duration, 0.1s)`, which is
garbage on short urban edges — corrupting per-segment bucket labels and, through
them, the empirical base speeds (see memory project_google_routing_calibration,
the INT diagnosis).

Instead we use the one /match channel that is exact: the matched **OSM node
sequence**. For every cached Google route we map-match its geometry on the real
deployed OSRM (:5000, geometry-driven so profile-independent) and collect the
union of matched node ids and consecutive node-pairs. Then one streaming pass
over the NI .osm.pbf records, for every OSM way/node touching that route set,
**its complete tag dictionary** (plus way node-lists and node coordinates).

Design decision (deliberately raw): we cache EVERYTHING along the route set and
do no bucketing/classification/filtering here. `highway`/`maxspeed` are what the
current profile model reads (via profile_spec.bucket_of), but lanes, surface,
oneway, lit, junction, crossing, … are all kept so a future model version can
use them without re-querying or re-parsing the pbf. Bucketing/filtering is a
pure downstream transform over this cache.

No new pip dependency: the pbf is streamed via the existing osmctools-roaaads
Docker image (osmconvert/osmfilter — same as build_network.py and the old
build_skeleton_index --signals phase) into an all-`highway` .osm, then read with
stdlib xml.etree.iterparse + a matched-node target filter (low RAM — only data
touching the route set is retained). pyosmium would be the streaming-native tool
but ships no wheel for this environment's Python 3.8.

Phases:
  python3 simulation/build_edge_index.py --match      # /match the cache -> matched_nodes.json (free, local OSRM)
  python3 simulation/build_edge_index.py --extract     # osmctools + iterparse -> osm_ways/osm_nodes.jsonl
  python3 simulation/build_edge_index.py                # both, in order
"""

import argparse
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "analysis"))
from google_routing_common import (        # noqa: E402
    decode_polyline, downsample_by_distance, osrm_match_detail)

REPO_ROOT = "/home/matthew/Documents/CodingFun/roaaads"
OSRM_ROOT = "/home/matthew/Documents/CodingFun/osrm"
PBF_NAME = "ireland-and-northern-ireland-latest.osm.pbf"
WORK_DIR = os.path.join(OSRM_ROOT, "edge_index")     # scratch (o5m / .osm extract)

CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
RAW_DIR = os.path.join(CACHE_DIR, "raw")
MANIFEST_FILE = os.path.join(CACHE_DIR, "od_manifest.json")

MATCH_CACHE_FILE = os.path.join(CACHE_DIR, "match_cache.jsonl")
OSM_WAYS_FILE = os.path.join(CACHE_DIR, "osm_ways.jsonl")
OSM_NODES_FILE = os.path.join(CACHE_DIR, "osm_nodes.jsonl")
EDGE_INDEX_META = os.path.join(CACHE_DIR, "edge_index_meta.json")

OSMCTOOLS_IMAGE = "osmctools-roaaads"     # built by build_network.py


# ── Phase 1: /match the cache -> match_cache.jsonl (the single OSRM pass) ─────
# This is the ONLY map-match over the route set. We cache the full match detail
# (node sequence + per-segment distances + step maneuvers + Google duration) per
# route so that (a) the edge-index extract and (b) build_skeleton_index both read
# this cache instead of re-/matching — and re-bucketing after any profile_spec
# change is free. /match is geometry-driven, so the cache is (near) profile-
# independent and only needs rebuilding if the cached Google routes change.

def _read_match_cache():
    """{(od_id, route_idx): record}, tolerant of a partial trailing line."""
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


def collect_matches(osrm_url, limit=None, manifest_file=None):
    """Map-match every cached route on the real OSRM, caching full detail.
    Resumable: skips routes already in match_cache.jsonl; `limit` caps ODs this
    run (for quick validation / batching the slow ~1.7 s/match pass).
    manifest_file overrides MANIFEST_FILE — use for a second batch (v2)."""
    mf = manifest_file or MANIFEST_FILE
    if not os.path.exists(mf):
        sys.exit(f"ERROR: {mf} not found — run build_od_manifest.py first.")
    ods = json.load(open(mf))["od_pairs"]
    cached = [o for o in ods if os.path.exists(os.path.join(RAW_DIR, f"{o['od_id']}.json"))]
    if limit:
        cached = cached[:limit]
    done = set(_read_match_cache())
    print(f"Map-matching {len(cached)} cached ODs on {osrm_url} "
          f"({len(done)} routes already cached) …")

    n_new = n_fail = 0
    out = open(MATCH_CACHE_FILE, "a", buffering=1)
    try:
        for k, o in enumerate(cached):
            raw = json.load(open(os.path.join(RAW_DIR, f"{o['od_id']}.json")))
            for j, route in enumerate(raw.get("routes", [])):
                if (o["od_id"], j) in done:
                    continue
                coords = downsample_by_distance(
                    decode_polyline(route["polyline"]["encodedPolyline"]))
                det = osrm_match_detail(osrm_url, coords)
                if det is None:
                    n_fail += 1
                    continue
                out.write(json.dumps({
                    "od_id": o["od_id"], "route_idx": j,
                    "leg_type": o["leg_type"], "len_band": o["len_band"],
                    "g_dur": float(str(route["duration"]).rstrip("s")),
                    "g_dist": route.get("distanceMeters", 0) or 0,
                    "conf": round(det["conf"], 4),
                    "nodes": det["nodes"],
                    "distances": [round(d, 2) for d in det["distances"]],
                    "maneuvers": det["maneuvers"],
                }) + "\n")
                done.add((o["od_id"], j))
                n_new += 1
            if (k + 1) % 100 == 0:
                print(f"  {k + 1}/{len(cached)} ODs ({n_new} new this run)")
    finally:
        out.close()
    print(f"\nWrote {MATCH_CACHE_FILE}  ({n_new} new this run, {n_fail} match-fail, "
          f"{len(done)} total routes)")


def _matched_nodes_and_pairs():
    """Derive the matched node set + undirected node-pair set from the cache."""
    node_ids, pairs = set(), set()
    for r in _read_match_cache().values():
        nodes = r["nodes"]
        node_ids.update(nodes)
        for i in range(len(nodes) - 1):
            a, b = nodes[i], nodes[i + 1]
            if a != b:
                pairs.add((a, b) if a < b else (b, a))
    return node_ids, pairs


# ── Phase 2: stream the pbf -> raw way/node tag cache ────────────────────────

def _osmctools_extract():
    """pbf -> all-highway .osm (ways + dependent nodes, with tags + coords) via
    the osmctools-roaaads Docker image. Returns the host path of the .osm."""
    os.makedirs(WORK_DIR, exist_ok=True)
    o5m = os.path.join(WORK_DIR, "ni.o5m")
    osm = os.path.join(WORK_DIR, "highways.osm")
    uidgid = f"{os.getuid()}:{os.getgid()}"

    def _docker(cmd):
        subprocess.run(["docker", "run", "--rm", "--user", uidgid,
                        "-v", f"{OSRM_ROOT}:/data", "-v", f"{WORK_DIR}:/out",
                        OSMCTOOLS_IMAGE, "sh", "-c", cmd], check=True)

    if not os.path.exists(o5m):
        print("Converting pbf -> o5m (osmconvert, streaming) …")
        _docker(f"osmconvert /data/{PBF_NAME} -t=/out/_osmconvert_tmp -o=/out/ni.o5m")
    else:
        print(f"Reusing existing {o5m}")
    # Keep every way carrying a `highway` tag (any value); osmfilter follows
    # dependencies by default (the opposite of the --signals phase's
    # --ignore-dependencies), so each kept way's referenced nodes — with their
    # coords + tags — come along too. This is the full road network along/around
    # the routes; the matched-node filter in the iterparse pass narrows it to the
    # actual route set.
    print("Filtering highway=* ways + dependent nodes (osmfilter) …")
    _docker('osmfilter /out/ni.o5m -t=/out/_osmfilter_tmp '
            '--keep="highway=" -o=/out/highways.osm')
    sz = os.path.getsize(osm) / 1e6
    print(f"  Wrote {osm}  ({sz:.0f} MB)")
    return osm


def _iter_clear(osm_path, tag):
    """Stream <tag> elements from an .osm with bounded memory.

    OSM XML lists every <node> before any <way>, all as direct children of the
    <osm> root, so we must drop each parsed top-level element from the root (not
    just clear the target tag) or the root's child list grows to the full file.
    The yielded element is consumed during the `yield` (before the generator
    resumes), so clearing the root on resume is safe."""
    context = ET.iterparse(osm_path, events=("start", "end"))
    _event, root = next(context)
    for event, elem in context:
        if event != "end":
            continue
        if elem.tag == tag:
            yield elem
        if elem.tag in ("node", "way", "relation"):
            root.clear()      # ~1 child at a time -> O(1) amortised, bounded RAM


def extract():
    if not os.path.exists(MATCH_CACHE_FILE):
        sys.exit(f"ERROR: {MATCH_CACHE_FILE} not found — run --match first.")
    matched_nodes, _pairs = _matched_nodes_and_pairs()
    print(f"{len(matched_nodes)} matched nodes (from match cache) to anchor the extract")

    osm = _osmctools_extract()

    # Pass 1: ways. Keep any highway way touching >=1 matched node; write its full
    # tag dict + node list. Accumulate the union of their node refs (needed_nodes)
    # so pass 2 can store those nodes' coords/tags.
    print("Pass 1/2: scanning ways …")
    needed_nodes = set()
    n_ways_seen = n_ways_kept = 0
    with open(OSM_WAYS_FILE, "w") as wf:
        for w in _iter_clear(osm, "way"):
            n_ways_seen += 1
            nds = [int(nd.get("ref")) for nd in w.findall("nd")]
            if not any(n in matched_nodes for n in nds):
                continue
            tags = {t.get("k"): t.get("v") for t in w.findall("tag")}
            wf.write(json.dumps({"id": int(w.get("id")), "tags": tags,
                                 "nodes": nds}) + "\n")
            needed_nodes.update(nds)
            n_ways_kept += 1
            if n_ways_kept % 50000 == 0:
                print(f"    {n_ways_kept} ways kept …")
    print(f"  {n_ways_kept}/{n_ways_seen} ways kept; "
          f"{len(needed_nodes)} referenced nodes")

    # Pass 2: nodes. Store coords + full tags for every node referenced by a kept
    # way (a superset of the matched nodes — gives geometry for turn/angle work
    # and full junction/crossing/signal tags).
    print("Pass 2/2: scanning nodes …")
    n_nodes_seen = n_nodes_kept = 0
    with open(OSM_NODES_FILE, "w") as nf:
        for nd in _iter_clear(osm, "node"):
            n_nodes_seen += 1
            nid = int(nd.get("id"))
            if nid not in needed_nodes:
                continue
            tags = {t.get("k"): t.get("v") for t in nd.findall("tag")}
            nf.write(json.dumps({"id": nid,
                                 "lat": float(nd.get("lat")),
                                 "lon": float(nd.get("lon")),
                                 "tags": tags}) + "\n")
            n_nodes_kept += 1
    print(f"  {n_nodes_kept}/{n_nodes_seen} nodes kept")

    meta = {
        "n_matched_nodes": len(matched_nodes),
        "n_ways_kept": n_ways_kept, "n_ways_seen": n_ways_seen,
        "n_nodes_kept": n_nodes_kept, "n_nodes_seen": n_nodes_seen,
        "n_referenced_nodes": len(needed_nodes),
        "matched_nodes_resolved": sum(1 for n in matched_nodes if n in needed_nodes),
    }
    with open(EDGE_INDEX_META, "w") as f:
        json.dump(meta, f, indent=2)
    resolved = meta["matched_nodes_resolved"]
    print(f"\nWrote {OSM_WAYS_FILE}, {OSM_NODES_FILE}, {EDGE_INDEX_META}")
    print(f"  matched nodes resolved in extract: {resolved}/{len(matched_nodes)} "
          f"({100.0 * resolved / max(1, len(matched_nodes)):.1f}%)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--match", action="store_true",
                    help="phase 1 only: /match the cache -> match_cache.jsonl")
    ap.add_argument("--extract", action="store_true",
                    help="phase 2 only: osmctools + iterparse -> osm_ways/osm_nodes.jsonl")
    ap.add_argument("--osrm-url", default="http://localhost:5000",
                    help="deployed OSRM for /match (geometry-driven; profile-independent)")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap ODs processed in --match this run (resumable; for "
                         "quick validation or batching the slow match pass)")
    ap.add_argument("--manifest", default=None,
                    help="override manifest file for --match (e.g. od_manifest_v2.json); "
                         "appends to the shared match_cache.jsonl; use when processing "
                         "a second batch without touching the v1 manifest")
    args = ap.parse_args()

    if args.match and not args.extract:
        collect_matches(args.osrm_url, args.limit, args.manifest)
    elif args.extract and not args.match:
        extract()
    else:
        collect_matches(args.osrm_url, args.limit, args.manifest)
        extract()


if __name__ == "__main__":
    main()
