"""
Resumable Google Routes runner for the time-calibration sample.

Consumes the fixed manifest (build_od_manifest.py) and is crash-safe:

  Phase A (spendy, resumable): for each manifest OD with no cached raw response,
    query Google and write the raw JSON to data/google_cache/raw/<od_id>.json
    *immediately*. A crash loses at most the in-flight query; re-running skips
    every OD whose raw file already exists. --limit caps queries per run so spend
    can be released in controlled batches.

  Phase B (free, idempotent): rebuild data/google_cache/results.jsonl by running
    OSRM /match over every cached raw response (best route + all alternatives) and
    OSRM /route for each OD. No API calls; safe to re-run any time (e.g. after
    changing the matching logic). Run standalone with --reprocess-only.

NEVER makes a live Google call without GOOGLE_MAPS_API_KEY set AND will not start
Phase A unless there is at least one uncached OD; see memory
feedback_no_google_api_without_approval — get explicit per-run approval before use.

Usage:
  python3 analysis/google_query_routes.py --dry-run            # counts + cost, no spend
  python3 analysis/google_query_routes.py --limit 100          # query <=100 uncached ODs
  python3 analysis/google_query_routes.py                      # query all remaining
  python3 analysis/google_query_routes.py --reprocess-only     # rebuild results.jsonl only
"""

import argparse, json, os, sys, time, urllib.error

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from google_routing_common import (
    CONF_MIN, decode_polyline, downsample_by_distance,
    osrm_route, osrm_match_detail, google_routes, parse_google_duration)

REPO_ROOT     = "/home/matthew/Documents/CodingFun/roaaads"
CACHE_DIR         = os.path.join(REPO_ROOT, "data", "google_cache")
RAW_DIR           = os.path.join(CACHE_DIR, "raw")
MANIFEST_FILE     = os.path.join(CACHE_DIR, "od_manifest.json")
RESULTS_FILE      = os.path.join(CACHE_DIR, "results.jsonl")
MATCH_CACHE_FILE  = os.path.join(CACHE_DIR, "match_cache.jsonl")


def raw_path(od_id):
    return os.path.join(RAW_DIR, f"{od_id}.json")


def phase_a_query(ods, args):
    """Query uncached ODs from Google, writing each raw response immediately."""
    remaining = [o for o in ods if not os.path.exists(raw_path(o["od_id"]))]
    done = len(ods) - len(remaining)
    print(f"Phase A: {done}/{len(ods)} already cached; {len(remaining)} remaining")
    if args.limit:
        remaining = remaining[:args.limit]
        print(f"  --limit {args.limit}: querying {len(remaining)} this run")

    if args.dry_run:
        est = len(remaining) * 5.0 / 1000.0
        print(f"  DRY RUN — would make {len(remaining)} live Google calls (~${est:.2f}). "
              f"No calls made.")
        return 0
    if not remaining:
        print("  Nothing to query.")
        return 0

    api_key = os.environ.get("GOOGLE_MAPS_API_KEY")
    if not api_key:
        print("  ERROR: GOOGLE_MAPS_API_KEY not set — refusing to start Phase A.")
        sys.exit(1)

    dep_iso = None  # free-flow; manifest meta records time_basis
    n = 0
    t0 = time.time()
    for o in remaining:
        try:
            gdata = google_routes(api_key, o["o"]["lat"], o["o"]["lon"],
                                  o["d"]["lat"], o["d"]["lon"], traffic=False,
                                  departure_iso=dep_iso)
        except urllib.error.HTTPError as e:
            print(f"  {o['od_id']}: Google HTTP {e.code}: "
                  f"{e.read()[:160].decode(errors='replace')}")
            continue
        except Exception as e:
            print(f"  {o['od_id']}: Google error: {e}")
            continue
        # Write raw immediately (crash-safe): temp then atomic rename.
        tmp = raw_path(o["od_id"]) + ".tmp"
        with open(tmp, "w") as f:
            json.dump(gdata, f)
        os.replace(tmp, raw_path(o["od_id"]))
        n += 1
        if n % 25 == 0:
            rate = n / (time.time() - t0)
            print(f"  {n}/{len(remaining)} queried ({rate:.1f}/s)")
        time.sleep(0.05)
    print(f"  Phase A done: {n} new queries, {n * 5.0 / 1000.0:.2f} USD this run")
    return n


def _read_match_cache():
    """Return {(od_id, route_idx): record} from match_cache.jsonl."""
    cache = {}
    if not os.path.exists(MATCH_CACHE_FILE):
        return cache
    with open(MATCH_CACHE_FILE) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
                cache[(r["od_id"], r["route_idx"])] = r
            except (json.JSONDecodeError, KeyError):
                continue
    return cache


def phase_b_reprocess(ods, osrm_url):
    """Rebuild results.jsonl from cached raw responses, writing full match detail
    to match_cache.jsonl in the same pass (one /match per route, not two).

    Routes already in match_cache.jsonl are skipped for the detail write but
    their cached detail is used to rebuild results.jsonl entries.
    """
    cached = [o for o in ods if os.path.exists(raw_path(o["od_id"]))]
    print(f"\nPhase B: reprocessing {len(cached)} cached ODs via OSRM /match …")

    match_cache = _read_match_cache()
    n_already = sum(1 for o in cached
                    for j in range(3) if (o["od_id"], j) in match_cache)
    print(f"  {len(match_cache)} routes already in match_cache.jsonl — no second /match for those")

    n_routes = n_valid = n_lowconf = n_matchfail = n_skip = 0
    tmp = RESULTS_FILE + ".tmp"
    mc_out = open(MATCH_CACHE_FILE, "a", buffering=1)
    try:
        with open(tmp, "w") as out:
            for k, o in enumerate(cached):
                gdata = json.load(open(raw_path(o["od_id"])))
                groutes = gdata.get("routes", [])
                if not groutes:
                    continue
                r = osrm_route(osrm_url, o["o"]["lat"], o["o"]["lon"],
                               o["d"]["lat"], o["d"]["lon"])
                if r is None:
                    continue
                osrm_route_dur, osrm_route_dist, route_nodes = r

                per_route = []
                for j, route in enumerate(groutes):
                    gdur  = parse_google_duration(route["duration"])
                    gdist = route.get("distanceMeters", 0)

                    cached_det = match_cache.get((o["od_id"], j))
                    if cached_det is not None:
                        # Already matched — read from cache, no second /match call.
                        match_dur   = cached_det.get("match_dur")   # None for old v1 entries
                        match_nodes = cached_det.get("nodes", [])
                        conf        = cached_det.get("conf", 0.0)
                        n_skip += 1
                    else:
                        gsamp = downsample_by_distance(
                            decode_polyline(route["polyline"]["encodedPolyline"]))
                        det = osrm_match_detail(osrm_url, gsamp)
                        if det is not None:
                            match_dur   = det["duration"]
                            match_nodes = det["nodes"]
                            conf        = det["conf"]
                            mc_out.write(json.dumps({
                                "od_id":      o["od_id"],
                                "route_idx":  j,
                                "leg_type":   o["leg_type"],
                                "len_band":   o["len_band"],
                                "g_dur":      gdur,
                                "g_dist":     gdist,
                                "match_dur":  round(match_dur, 3),
                                "conf":       round(conf, 4),
                                "nodes":      match_nodes,
                                "distances":  [round(d, 2) for d in det["distances"]],
                                "maneuvers":  det["maneuvers"],
                            }) + "\n")
                        else:
                            match_dur = match_nodes = conf = None

                    if match_dur is None:
                        match_nodes = []
                        conf = 0.0
                        n_matchfail += 1

                    valid = (match_dur is not None) and (conf >= CONF_MIN)
                    n_routes += 1
                    n_valid  += valid
                    if (match_dur is not None) and not valid:
                        n_lowconf += 1

                    overlap = None
                    if j == 0 and match_nodes:
                        sm, sr = set(match_nodes), set(route_nodes)
                        overlap = len(sm & sr) / len(sm)
                    per_route.append({
                        "idx": j, "g_dur": gdur, "g_dist": gdist,
                        "match_dur": match_dur, "conf": conf,
                        "te_matched": (match_dur / gdur) if match_dur else None,
                        "valid": valid,
                        "route_overlap": overlap,
                    })

                best = per_route[0]
                rec = {
                    "od_id": o["od_id"], "leg_type": o["leg_type"],
                    "len_band": o["len_band"], "len_s": o["len_s"],
                    "o": o["o"], "d": o["d"],
                    "google_best_dur": best["g_dur"], "n_alts": len(groutes),
                    "osrm_route_dur": osrm_route_dur, "osrm_route_dist": osrm_route_dist,
                    "te_endpoint": osrm_route_dur / best["g_dur"],
                    "route_overlap": best["route_overlap"],
                    "routes": per_route,
                }
                out.write(json.dumps(rec) + "\n")
                if (k + 1) % 100 == 0:
                    print(f"  {k+1}/{len(cached)} reprocessed")
    finally:
        mc_out.close()
    os.replace(tmp, RESULTS_FILE)

    tem = []
    for line in open(RESULTS_FILE):
        rec = json.loads(line)
        tem += [r["te_matched"] for r in rec["routes"] if r["valid"]]
    print(f"\nWrote {RESULTS_FILE}  (also appended new routes to {MATCH_CACHE_FILE})")
    print(f"  {n_routes} routes, {n_valid} valid (conf>={CONF_MIN}), "
          f"{n_lowconf} low-conf, {n_matchfail} match-fail, {n_skip} already cached")
    if tem:
        tem.sort()
        print(f"  te_matched: median {tem[len(tem)//2]:.2f}, mean {sum(tem)/len(tem):.2f}, "
              f"min {min(tem):.2f}, max {max(tem):.2f}, n={len(tem)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default=MANIFEST_FILE)
    ap.add_argument("--osrm-url", default="http://localhost:5000")
    ap.add_argument("--limit", type=int, default=0, help="max Google queries this run (0=all)")
    ap.add_argument("--dry-run", action="store_true", help="report counts/cost, no spend")
    ap.add_argument("--reprocess-only", action="store_true",
                    help="rebuild results.jsonl from cache via OSRM; no Google calls")
    args = ap.parse_args()

    if not os.path.exists(args.manifest):
        print(f"ERROR: manifest {args.manifest} not found. Run build_od_manifest.py first.")
        sys.exit(1)
    man = json.load(open(args.manifest))
    ods = man["od_pairs"]
    os.makedirs(RAW_DIR, exist_ok=True)
    print(f"Manifest: {len(ods)} OD pairs (created {man['meta'].get('created')}, "
          f"basis {man['meta'].get('time_basis')})")

    if args.reprocess_only:
        phase_b_reprocess(ods, args.osrm_url)
        return

    phase_a_query(ods, args)
    if args.dry_run:
        return
    # Only reprocess if OSRM is reachable; Phase B needs it.
    if osrm_route(args.osrm_url, 54.5933, -5.6960, 54.6033, -5.6960) is None:
        print(f"\nWARNING: OSRM not reachable at {args.osrm_url} — skipping Phase B. "
              f"Run --reprocess-only later.")
        return
    phase_b_reprocess(ods, args.osrm_url)


if __name__ == "__main__":
    main()
