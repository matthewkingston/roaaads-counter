"""
Build the profile-independent *skeleton* cache from the Google route cache —
the one-time, expensive step that makes profile benchmarking fast.

For every cached Google route, we map-match its geometry once and record a
representation that does NOT depend on any speed profile:

  * length_by_bucket : metres of matched road in each (class x speed-band) bucket
  * turns            : per-maneuver (angle, junction degree, u-turn flag)
  * n_signals        : traffic-signal nodes on the matched path
  * g_dur            : Google's free-flow duration (the calibration target)

Bucket labels come from OSRM itself via a one-time *probe profile* (osrm_lua.
emit_probe_block) that encodes each way's integer bucket id as its speed; a
single /match then reads the bucket of every matched segment from the `speed`
annotation. Segment lengths come from the `distance` annotation (geometry, so
profile-independent); turns from parsed steps; signals from a traffic-signal
node-id set extracted with osmctools. No pyosmium, no NI graph index.

Three phases (run in order; the heavy OSRM build is a manual Docker step):

  # 1. Generate the probe profile into a SEPARATE osrm data dir, then build it.
  python3 simulation/build_skeleton_index.py --gen-probe
  #    ... run the printed osrm-extract/partition/customize + serve on :5001 ...

  # 2. Extract the traffic-signal node-id set from the pbf (osmctools; ~1-2 min).
  python3 simulation/build_skeleton_index.py --signals

  # 3. Map-match the whole cache through the probe server -> skeletons.jsonl.
  python3 simulation/build_skeleton_index.py --osrm-url http://localhost:5001

The probe instance is intentionally separate from the deployed :5000 OSRM so
this never clobbers production routing data.
"""

import argparse
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import osrm_lua
import profile_spec as ps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "analysis"))
from google_routing_common import (        # noqa: E402
    CONF_MIN, decode_polyline, downsample_by_distance, osrm_match_detail)

REPO_ROOT = "/home/matthew/Documents/CodingFun/roaaads"
OSRM_ROOT = "/home/matthew/Documents/CodingFun/osrm"
PROBE_DIR = os.path.join(OSRM_ROOT, "probe")
PBF_NAME = "ireland-and-northern-ireland-latest.osm.pbf"

CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
RAW_DIR = os.path.join(CACHE_DIR, "raw")
SKELETONS_FILE = os.path.join(CACHE_DIR, "skeletons.jsonl")
SIGNALS_FILE = os.path.join(CACHE_DIR, "signal_nodes.json")
MANIFEST_FILE = os.path.join(CACHE_DIR, "od_manifest.json")

OSMCTOOLS_IMAGE = "osmctools-roaaads"     # built by build_network.py
COVERAGE_LO, COVERAGE_HI = 0.8, 1.2       # matched_dist / google_dist sanity band


# ── Phase 1: generate the probe profile ──────────────────────────────────────

def gen_probe():
    print(f"Pulling stock car.lua from {osrm_lua.OSRM_IMAGE} …")
    base_lua, found_path = osrm_lua.pull_base_lua()
    patched, strategy = osrm_lua.inject(base_lua, osrm_lua.emit_probe_block())
    print(f"  Injection strategy: {strategy}")

    os.makedirs(PROBE_DIR, exist_ok=True)
    out_lua = os.path.join(PROBE_DIR, "car_probe.lua")
    with open(out_lua, "w") as f:
        f.write(patched)
    print(f"  Wrote {out_lua}")
    ok, lib_src = osrm_lua.copy_lib(found_path, PROBE_DIR)
    print("  lib/ copied." if ok else f"  WARNING: copy {lib_src}/ into {PROBE_DIR} manually.")

    base = PBF_NAME.replace(".osm.pbf", "")
    print(f"""
Probe profile ready. Build a SEPARATE probe OSRM instance (does not touch :5000):

  cp {os.path.join(OSRM_ROOT, PBF_NAME)} {PROBE_DIR}/
  cd {PROBE_DIR}
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-extract  -p /data/car_probe.lua /data/{PBF_NAME}
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-partition /data/{base}.osrm
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-customize /data/{base}.osrm
  docker run -t -i -p 5001:5000 -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-routed --algorithm mld /data/{base}.osrm

Then:  python3 simulation/build_skeleton_index.py --signals
       python3 simulation/build_skeleton_index.py --osrm-url http://localhost:5001
""")


# ── Phase 2: traffic-signal node-id set ──────────────────────────────────────

def extract_signals():
    """Write the set of OSM node ids tagged highway=traffic_signals (all NI)."""
    pbf = os.path.join(OSRM_ROOT, PBF_NAME)
    o5m = os.path.join(PROBE_DIR, "ni.o5m")
    sig_osm = os.path.join(PROBE_DIR, "signals.osm")
    os.makedirs(PROBE_DIR, exist_ok=True)

    def _docker(cmd):
        subprocess.run(["docker", "run", "--rm", "-v", f"{OSRM_ROOT}:/data",
                        "-v", f"{PROBE_DIR}:/out", OSMCTOOLS_IMAGE, "sh", "-c", cmd],
                       check=True)

    if not os.path.exists(o5m):
        print("Converting pbf -> o5m (osmconvert, streaming) …")
        _docker(f"osmconvert /data/{PBF_NAME} -o=/out/ni.o5m")
    print("Filtering highway=traffic_signals nodes (osmfilter) …")
    _docker('osmfilter /out/ni.o5m --keep="highway=traffic_signals" '
            '--ignore-dependencies -o=/out/signals.osm')

    ids = set()
    node_re = re.compile(r'<node[^>]*\bid="(\d+)"')
    with open(sig_osm, "r", errors="replace") as f:
        for line in f:
            m = node_re.search(line)
            if m:
                ids.add(int(m.group(1)))
    with open(SIGNALS_FILE, "w") as f:
        json.dump(sorted(ids), f)
    print(f"  Wrote {SIGNALS_FILE}  ({len(ids)} signal nodes)")


# ── Phase 3: map-match the cache -> skeletons.jsonl ──────────────────────────

def _bucket_lengths(speeds_mps, distances):
    """Aggregate matched-segment lengths into (class|band) buckets."""
    out, unbucketed = {}, 0.0
    n = min(len(speeds_mps), len(distances))
    for i in range(n):
        d = distances[i]
        if not d:
            continue
        bk = ps.bucket_from_probe_speed(speeds_mps[i] * 3.6)
        if bk is None:
            unbucketed += d
            continue
        key = ps.bucket_key(*bk)
        out[key] = out.get(key, 0.0) + d
    return out, unbucketed


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


def build_skeletons(osrm_url):
    if not os.path.exists(MANIFEST_FILE):
        sys.exit(f"ERROR: {MANIFEST_FILE} not found — run build_od_manifest.py first.")
    ods = json.load(open(MANIFEST_FILE))["od_pairs"]

    signals = set()
    if os.path.exists(SIGNALS_FILE):
        signals = set(json.load(open(SIGNALS_FILE)))
        print(f"Loaded {len(signals)} traffic-signal nodes")
    else:
        print(f"WARNING: {SIGNALS_FILE} missing — n_signals will be 0. Run --signals.")

    cached = [o for o in ods if os.path.exists(os.path.join(RAW_DIR, f"{o['od_id']}.json"))]
    print(f"Map-matching {len(cached)} cached ODs through probe at {osrm_url} …")

    n_routes = n_valid = n_matchfail = 0
    tmp = SKELETONS_FILE + ".tmp"
    with open(tmp, "w") as out:
        for k, o in enumerate(cached):
            raw = json.load(open(os.path.join(RAW_DIR, f"{o['od_id']}.json")))
            for j, route in enumerate(raw.get("routes", [])):
                gdur = float(str(route["duration"]).rstrip("s"))
                gdist = route.get("distanceMeters", 0) or 0
                coords = downsample_by_distance(
                    decode_polyline(route["polyline"]["encodedPolyline"]))
                det = osrm_match_detail(osrm_url, coords)
                n_routes += 1
                if det is None:
                    n_matchfail += 1
                    continue
                lbb, unbk = _bucket_lengths(det["speeds"], det["distances"])
                matched_m = sum(lbb.values()) + unbk
                coverage = (matched_m / gdist) if gdist else None
                n_sig = sum(1 for nid in det["nodes"] if nid in signals)
                valid = (det["conf"] >= CONF_MIN and coverage is not None
                         and COVERAGE_LO <= coverage <= COVERAGE_HI)
                n_valid += valid
                out.write(json.dumps({
                    "od_id": o["od_id"], "route_idx": j,
                    "leg_type": o["leg_type"], "len_band": o["len_band"],
                    "g_dur": gdur, "g_dist": gdist,
                    "conf": round(det["conf"], 4),
                    "coverage": round(coverage, 4) if coverage is not None else None,
                    "valid": valid,
                    "length_by_bucket": {kk: round(vv, 1) for kk, vv in lbb.items()},
                    "unbucketed_m": round(unbk, 1),
                    "n_signals": n_sig,
                    "turns": _turns(det["maneuvers"]),
                }) + "\n")
            if (k + 1) % 100 == 0:
                print(f"  {k + 1}/{len(cached)} ODs")
    os.replace(tmp, SKELETONS_FILE)
    print(f"\nWrote {SKELETONS_FILE}")
    print(f"  {n_routes} routes, {n_valid} valid (conf>={CONF_MIN} & coverage in "
          f"[{COVERAGE_LO},{COVERAGE_HI}]), {n_matchfail} match-fail")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen-probe", action="store_true",
                    help="generate car_probe.lua + print the probe-build commands")
    ap.add_argument("--signals", action="store_true",
                    help="extract the traffic-signal node-id set from the pbf")
    ap.add_argument("--osrm-url", default="http://localhost:5001",
                    help="probe OSRM server (default :5001, separate from deployed :5000)")
    args = ap.parse_args()

    if args.gen_probe:
        gen_probe()
    elif args.signals:
        extract_signals()
    else:
        build_skeletons(args.osrm_url)


if __name__ == "__main__":
    main()
