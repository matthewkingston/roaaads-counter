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
import collections
import json
import os
import re
import statistics
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
BASE_SPEEDS_FILE = os.path.join(CACHE_DIR, "base_speeds.json")
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

    uidgid = f"{os.getuid()}:{os.getgid()}"

    def _docker(cmd):
        # Run as the host user (like build_network.py) so osmconvert/osmfilter
        # outputs are owned by us and readable back; without --user the container
        # writes them root-owned mode 600.
        subprocess.run(["docker", "run", "--rm", "--user", uidgid,
                        "-v", f"{OSRM_ROOT}:/data", "-v", f"{PROBE_DIR}:/out",
                        OSMCTOOLS_IMAGE, "sh", "-c", cmd], check=True)

    if not os.path.exists(o5m):
        print("Converting pbf -> o5m (osmconvert, streaming) …")
        _docker(f"osmconvert /data/{PBF_NAME} -t=/out/_osmconvert_tmp -o=/out/ni.o5m")
    print("Filtering highway=traffic_signals nodes (osmfilter) …")
    # -t= points osmfilter's temp files at the writable /out mount; the container
    # cwd is not writable as a non-root user (matches build_network.py).
    _docker('osmfilter /out/ni.o5m -t=/out/_osmfilter_tmp '
            '--keep="highway=traffic_signals" --ignore-dependencies -o=/out/signals.osm')

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


def _read_jsonl(path):
    """Load a skeletons jsonl into {(od_id, route_idx): record}, tolerating a
    partial/corrupt trailing line left by an interrupted run."""
    recs = {}
    if not os.path.exists(path):
        return recs
    with open(path) as f:
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


def _consolidate_existing():
    """Merge any prior skeletons.jsonl + a leftover .tmp into one clean file and
    return the set of (od_id, route_idx) already done. Salvages progress from an
    interrupted earlier run — including the old non-resumable .tmp, which the
    previous code wrote to line-by-line before its final atomic rename."""
    tmp_old = SKELETONS_FILE + ".tmp"
    merged = _read_jsonl(SKELETONS_FILE)
    salvaged = 0
    if os.path.exists(tmp_old):
        for key, rec in _read_jsonl(tmp_old).items():
            if key not in merged:
                merged[key] = rec
                salvaged += 1
    if merged:
        reb = SKELETONS_FILE + ".rebuild"
        with open(reb, "w") as f:
            for rec in merged.values():
                f.write(json.dumps(rec) + "\n")
        os.replace(reb, SKELETONS_FILE)
    if os.path.exists(tmp_old):
        os.remove(tmp_old)          # consumed; a new resumable run never recreates it
    if merged:
        msg = f"Resuming: {len(merged)} routes already done"
        if salvaged:
            msg += f" ({salvaged} salvaged from an interrupted .tmp)"
        print(msg)
    return set(merged)


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

    # Resumable: skip routes already cached; append each completed route durably
    # (line-buffered) so a crash/kill loses at most the in-flight route.
    done = _consolidate_existing()
    n_new = n_matchfail = 0
    out = open(SKELETONS_FILE, "a", buffering=1)
    try:
        for k, o in enumerate(cached):
            raw = json.load(open(os.path.join(RAW_DIR, f"{o['od_id']}.json")))
            for j, route in enumerate(raw.get("routes", [])):
                if (o["od_id"], j) in done:
                    continue
                gdur = float(str(route["duration"]).rstrip("s"))
                gdist = route.get("distanceMeters", 0) or 0
                coords = downsample_by_distance(
                    decode_polyline(route["polyline"]["encodedPolyline"]))
                det = osrm_match_detail(osrm_url, coords)
                if det is None:
                    n_matchfail += 1          # not recorded -> retried next run
                    continue
                lbb, unbk = _bucket_lengths(det["speeds"], det["distances"])
                matched_m = sum(lbb.values()) + unbk
                coverage = (matched_m / gdist) if gdist else None
                n_sig = sum(1 for nid in det["nodes"] if nid in signals)
                valid = (det["conf"] >= CONF_MIN and coverage is not None
                         and COVERAGE_LO <= coverage <= COVERAGE_HI)
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
                done.add((o["od_id"], j))
                n_new += 1
            if (k + 1) % 100 == 0:
                print(f"  {k + 1}/{len(cached)} ODs ({n_new} new this run)")
    finally:
        out.close()

    final = _read_jsonl(SKELETONS_FILE)
    n_valid = sum(1 for r in final.values() if r.get("valid"))
    print(f"\nWrote {SKELETONS_FILE}")
    print(f"  {len(final)} routes total ({n_new} new this run, {n_matchfail} match-fail), "
          f"{n_valid} valid (conf>={CONF_MIN} & coverage in [{COVERAGE_LO},{COVERAGE_HI}])")


def estimate_base_speeds(probe_url, speed_url, defactor, sample_n, min_obs):
    """Measure OSRM's realised factor-free base speed per (class x band) bucket.

    Per segment we need two facts: its bucket (from the probe's bucket-id speed
    annotation) and OSRM's *realised* speed. The probe carries the bucket, so the
    realised speed must come from a second instance:
      * a factor-free STOCK car.lua -> annotation.speed is the base speed directly
        (--no-defactor), or
      * the deployed legacy :5000 -> annotation.speed has the class factor baked
        in; recover the factor-free value as speed * HIGHWAY_COST_FACTOR[class]
        (--defactor, the default).
    Base speed is near-deterministic per bucket, so a sample of routes suffices.
    Segments are paired by OSM edge (node pair), robust to small match
    differences. Writes {bucket_key: km/h} (median per bucket) to base_speeds.json.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from routing_config import HIGHWAY_COST_FACTOR

    ods = json.load(open(MANIFEST_FILE))["od_pairs"]
    cached = [o for o in ods if os.path.exists(os.path.join(RAW_DIR, f"{o['od_id']}.json"))]
    stride = max(1, len(cached) // sample_n)
    sample = cached[::stride][:sample_n]
    print(f"Estimating base speeds from {len(sample)} sampled routes")
    print(f"  probe={probe_url}  speed={speed_url}  defactor={defactor}")

    acc = collections.defaultdict(list)        # bucket_key -> [km/h]
    n_routes = n_pairs = n_skip = 0
    for o in sample:
        raw = json.load(open(os.path.join(RAW_DIR, f"{o['od_id']}.json")))
        for route in raw.get("routes", [])[:1]:   # best route is plenty for sampling
            coords = downsample_by_distance(
                decode_polyline(route["polyline"]["encodedPolyline"]))
            pdet = osrm_match_detail(probe_url, coords)
            sdet = osrm_match_detail(speed_url, coords)
            if pdet is None or sdet is None:
                n_skip += 1
                continue
            # realised speed keyed by OSM edge (node pair)
            smap = {}
            sn, ss = sdet["nodes"], sdet["speeds"]
            for i in range(min(len(sn) - 1, len(ss))):
                smap[(sn[i], sn[i + 1])] = ss[i]
            pn, pspd = pdet["nodes"], pdet["speeds"]
            for i in range(min(len(pn) - 1, len(pspd))):
                bk = ps.bucket_from_probe_speed(pspd[i] * 3.6)
                if bk is None:
                    continue
                v = smap.get((pn[i], pn[i + 1]))
                if v is None:
                    continue
                v_kmh = v * 3.6
                if defactor:
                    v_kmh *= HIGHWAY_COST_FACTOR.get(bk[0], 1.0)
                if v_kmh <= 0:
                    continue
                acc[ps.bucket_key(*bk)].append(v_kmh)
                n_pairs += 1
            n_routes += 1

    base = {k: round(statistics.median(v), 2)
            for k, v in acc.items() if len(v) >= min_obs}
    with open(BASE_SPEEDS_FILE, "w") as f:
        json.dump(base, f, indent=2, sort_keys=True)
    print(f"\nWrote {BASE_SPEEDS_FILE}")
    print(f"  {n_routes} routes paired, {n_pairs} segment samples, {n_skip} skipped")
    print(f"  {len(base)} buckets with >= {min_obs} obs "
          f"({len(acc) - len(base)} buckets too sparse -> analytical fallback)")
    # quick compare to the analytical estimate for a few high-coverage buckets
    top = sorted(acc.items(), key=lambda kv: -len(kv[1]))[:8]
    print("  bucket                 emp km/h  analytic  n")
    for k, v in top:
        cls, band = k.split("|", 1)
        print(f"  {k:22} {statistics.median(v):7.1f}  "
              f"{ps.base_speed_for(cls, band):7.1f}  {len(v)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen-probe", action="store_true",
                    help="generate car_probe.lua + print the probe-build commands")
    ap.add_argument("--signals", action="store_true",
                    help="extract the traffic-signal node-id set from the pbf")
    ap.add_argument("--base-speeds", action="store_true",
                    help="measure empirical per-bucket base speeds (dual match)")
    ap.add_argument("--osrm-url", default="http://localhost:5001",
                    help="probe OSRM server (default :5001, separate from deployed :5000)")
    ap.add_argument("--speed-url", default="http://localhost:5000",
                    help="real-speed OSRM for --base-speeds (deployed legacy :5000, "
                         "or a factor-free stock instance with --no-defactor)")
    ap.add_argument("--no-defactor", dest="defactor", action="store_false",
                    help="speed-url is already factor-free (stock); skip the "
                         "HIGHWAY_COST_FACTOR division")
    ap.add_argument("--sample", type=int, default=400,
                    help="routes to sample for --base-speeds")
    ap.add_argument("--min-obs", type=int, default=5,
                    help="min segment samples to trust a bucket's empirical speed")
    ap.set_defaults(defactor=True)
    args = ap.parse_args()

    if args.gen_probe:
        gen_probe()
    elif args.signals:
        extract_signals()
    elif args.base_speeds:
        estimate_base_speeds(args.osrm_url, args.speed_url, args.defactor,
                             args.sample, args.min_obs)
    else:
        build_skeletons(args.osrm_url)


if __name__ == "__main__":
    main()
