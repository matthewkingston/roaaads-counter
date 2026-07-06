"""
Google Maps routing feasibility experiment (see task_google_routing.txt).

Goal: prove, cheaply, that we can use Google as a "source of truth" for routing,
and quantify how wrong the current OSRM profile is — separating *time* error from
*route-choice* error. This is the de-risking pilot before committing to a full
Google-calibration project (see memory project_google_routing_calibration).

For a small, data-driven sample of OD pairs (Ballyrainey/Belfast corridor + a
compass spread + nearest + farthest external zones, all → Newtownards centre) it:

  1. Queries Google Routes API v2 (computeRoutes, with alternatives) → duration,
     distance, encoded polyline for each returned route.
  2. Decodes Google's polyline and runs OSRM /match along *Google's geometry* →
     OSRM's own time on the roads Google chose (isolates time error from route
     choice). This /match step is the technical linchpin for stage-1 calibration.
  3. Runs OSRM /route for the same OD → OSRM's freely-chosen route + duration.
  4. Reports per-OD:
       time_err_matched  = osrm_match_dur / google_dur   (pure time error)
       time_err_endpoint = osrm_route_dur / google_dur    (includes route choice)
       route_overlap     = frac of Google's matched OSM nodes OSRM also uses
       n_alts            = number of Google alternative routes

Raw Google + OSRM responses are cached per OD (gitignored) so re-runs don't
re-spend on the API. Use --dry-run to see the OD list + cost estimate without
calling Google (still pings OSRM /route to confirm reachability).

Free-flow times: defaults to routingPreference=TRAFFIC_UNAWARE (cheapest tier,
free-flow) because the model is daily-average AADT, not a peak snapshot. Use
--traffic to switch to TRAFFIC_AWARE (live/typical traffic; higher cost tier) —
e.g. to reproduce a Google app number that included traffic.

Usage:
  export GOOGLE_MAPS_API_KEY=...        # pay-as-you-go key, Routes API enabled
  python3 analysis/google_feasibility.py --dry-run      # no spend; show plan+cost
  python3 analysis/google_feasibility.py                # live run (~pennies)
"""

import argparse, json, math, os, sys, time, urllib.request, urllib.error
import http.client

# ── Repo paths ────────────────────────────────────────────────────────────────
# This script may run from an isolated worktree, but census_zones.json and the
# cache live in the real checkout so results persist across worktree cleanup.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CENSUS_FILE = os.path.join(REPO_ROOT, "data", "census_zones.json")
CACHE_DIR   = os.path.join(REPO_ROOT, "data", "google_cache")

# Newtownards town centre — mirrors simulation/zones_config.py CENTRE (lat, lon).
CENTRE = (54.5933779, -5.6960935)

# Explicit reproduction of the user's spot-check: South Belfast → Newtownards.
SOUTH_BELFAST = (54.5620, -5.9500)

GOOGLE_URL = "https://routes.googleapis.com/directions/v2:computeRoutes"

# Drop /match results below this confidence — they produce garbage durations
# (one feasibility OD came back at conf 7.6e-5 with a bogus matched time).
CONF_MIN = 0.5


# ── Geometry helpers ────────────────────────────────────────────────────────────

def haversine_m(lat1, lon1, lat2, lon2):
    R = 6371000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


def bearing_deg(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


def decode_polyline(polyline_str):
    """Decode a Google encoded polyline to a list of (lat, lon)."""
    index, lat, lng = 0, 0, 0
    coords = []
    while index < len(polyline_str):
        for is_lat in (True, False):
            shift, result = 0, 0
            while True:
                b = ord(polyline_str[index]) - 63
                index += 1
                result |= (b & 0x1f) << shift
                shift += 5
                if b < 0x20:
                    break
            delta = ~(result >> 1) if (result & 1) else (result >> 1)
            if is_lat:
                lat += delta
            else:
                lng += delta
        coords.append((lat / 1e5, lng / 1e5))
    return coords


def downsample_by_distance(coords, step_m=80.0, max_pts=95):
    """Thin a dense polyline: keep first/last + a point roughly every step_m.
    Caps total at max_pts (OSRM /match coordinate limit) by widening the step."""
    if len(coords) <= max_pts:
        kept = [coords[0]]
        acc = 0.0
        for i in range(1, len(coords)):
            acc += haversine_m(*coords[i - 1], *coords[i])
            if acc >= step_m:
                kept.append(coords[i]); acc = 0.0
        if kept[-1] != coords[-1]:
            kept.append(coords[-1])
        if len(kept) <= max_pts:
            return kept
    # too many even after stepping — fall back to uniform stride
    stride = max(1, math.ceil(len(coords) / (max_pts - 1)))
    kept = coords[::stride]
    if kept[-1] != coords[-1]:
        kept.append(coords[-1])
    return kept


# ── OSRM ────────────────────────────────────────────────────────────────────────

def osrm_get(base_url, path):
    url = base_url.rstrip("/") + path
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.loads(r.read())


def osrm_route(base_url, lat1, lon1, lat2, lon2):
    """OSRM's freely-chosen route. Returns (duration_s, distance_m, node_list) or None."""
    path = (f"/route/v1/driving/{lon1},{lat1};{lon2},{lat2}"
            f"?annotations=nodes,duration&overview=false")
    try:
        data = osrm_get(base_url, path)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
        return None
    if data.get("code") != "Ok":
        return None
    rt = data["routes"][0]
    nodes = rt["legs"][0]["annotation"]["nodes"]
    return rt["duration"], rt["distance"], nodes


def osrm_match(base_url, coords, radius_m=40):
    """Map-match a (lat,lon) trace onto OSM via OSRM /match, returning OSRM's own
    duration + node sequence along that geometry. Returns
    (duration_s, node_list, confidence, n_matchings) or None on failure."""
    pts = ";".join(f"{lon},{lat}" for lat, lon in coords)
    radii = ";".join(str(radius_m) for _ in coords)
    path = (f"/match/v1/driving/{pts}"
            f"?annotations=nodes,duration&overview=false&geometries=geojson"
            f"&gaps=ignore&tidy=true&radiuses={radii}")
    try:
        data = osrm_get(base_url, path)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
        return None
    if data.get("code") != "Ok":
        return None
    total_dur = 0.0
    nodes = []
    confs = []
    for m in data.get("matchings", []):
        confs.append(m.get("confidence", 0.0))
        for leg in m.get("legs", []):
            ann = leg.get("annotation", {})
            total_dur += sum(ann.get("duration", []))
            nodes.extend(ann.get("nodes", []))
    if not data.get("matchings"):
        return None
    conf = sum(confs) / len(confs) if confs else 0.0
    return total_dur, nodes, conf, len(data["matchings"])


# ── Google Routes API ─────────────────────────────────────────────────────────

def google_routes(api_key, lat1, lon1, lat2, lon2, traffic=False):
    """Query Google Routes API v2 computeRoutes (with alternatives).
    Returns parsed JSON dict, or raises on HTTP error."""
    body = {
        "origin":      {"location": {"latLng": {"latitude": lat1, "longitude": lon1}}},
        "destination": {"location": {"latLng": {"latitude": lat2, "longitude": lon2}}},
        "travelMode":  "DRIVE",
        "computeAlternativeRoutes": True,
        "routingPreference": "TRAFFIC_AWARE" if traffic else "TRAFFIC_UNAWARE",
        "units": "METRIC",
    }
    if traffic:
        # TRAFFIC_AWARE needs a departure time; use a representative weekday daytime.
        body["departureTime"] = time.strftime(
            "%Y-%m-%dT%H:%M:%SZ", time.gmtime(time.time() + 3600))
    req = urllib.request.Request(
        GOOGLE_URL, data=json.dumps(body).encode(),
        headers={
            "Content-Type": "application/json",
            "X-Goog-Api-Key": api_key,
            "X-Goog-FieldMask": ("routes.duration,routes.distanceMeters,"
                                 "routes.polyline.encodedPolyline"),
        }, method="POST")
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())


def parse_google_duration(dur_str):
    """Google returns duration like '1234s'."""
    return float(str(dur_str).rstrip("s"))


# ── OD sample selection (data-driven, reproducible) ───────────────────────────

def build_od_sample(n_west=5, n_spread=6, n_near=4, n_far=2):
    with open(CENSUS_FILE) as f:
        zones = json.load(f)
    ext = zones["external_nodes"]
    clat, clon = CENTRE
    for e in ext:
        e["_dist"] = haversine_m(clat, clon, e["centroid_lat"], e["centroid_lon"])
        e["_brg"]  = bearing_deg(clat, clon, e["centroid_lat"], e["centroid_lon"])

    chosen = {}  # id -> (label, lat, lon)

    def add(e, tag):
        if e["id"] not in chosen:
            chosen[e["id"]] = (tag, e["centroid_lat"], e["centroid_lon"], e["_dist"])

    # Westernmost (Belfast / Dundonald / Ballyrainey corridor)
    for e in sorted(ext, key=lambda z: z["centroid_lon"])[:n_west]:
        add(e, "west")
    # Compass spread: pick the node nearest each of n_spread evenly-spaced bearings
    for k in range(n_spread):
        target = 360.0 * k / n_spread
        best = min(ext, key=lambda z: min(abs(z["_brg"] - target),
                                          360 - abs(z["_brg"] - target)))
        add(best, "spread")
    # Nearest (short journeys) and farthest (long journeys)
    for e in sorted(ext, key=lambda z: z["_dist"])[:n_near]:
        add(e, "near")
    for e in sorted(ext, key=lambda z: -z["_dist"])[:n_far]:
        add(e, "far")

    ods = []
    # Explicit Ballyrainey reproduction first
    ods.append(("ballyrainey_spotcheck", SOUTH_BELFAST[0], SOUTH_BELFAST[1],
                CENTRE[0], CENTRE[1]))
    for nid, (tag, lat, lon, dist) in chosen.items():
        ods.append((f"{tag}:{nid}", lat, lon, CENTRE[0], CENTRE[1]))
    return ods


# ── Cache ─────────────────────────────────────────────────────────────────────

def cache_key(lat1, lon1, lat2, lon2, traffic):
    t = "T" if traffic else "F"
    return f"{lat1:.5f}_{lon1:.5f}__{lat2:.5f}_{lon2:.5f}__{t}.json"


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="show OD list + cost estimate, ping OSRM, do NOT call Google")
    ap.add_argument("--traffic", action="store_true",
                    help="use TRAFFIC_AWARE (live/typical traffic) instead of free-flow")
    ap.add_argument("--osrm-url", default="http://localhost:5000")
    ap.add_argument("--refresh", action="store_true",
                    help="ignore cached Google responses and re-query (re-spends)")
    args = ap.parse_args()

    ods = build_od_sample()
    os.makedirs(CACHE_DIR, exist_ok=True)

    # Cost estimate: Routes API "Compute Routes" basic ~ $5 / 1000 requests.
    est_cost = len(ods) * 5.0 / 1000.0
    print(f"Sample: {len(ods)} OD pairs  (→ Newtownards centre)")
    print(f"Mode:   {'TRAFFIC_AWARE (traffic)' if args.traffic else 'TRAFFIC_UNAWARE (free-flow)'}")
    print(f"Est. Google cost: ~${est_cost:.2f}  (1 request per OD; alternatives free)")
    print(f"Cache dir: {CACHE_DIR}")
    print()

    # Confirm OSRM reachable.
    test = osrm_route(args.osrm_url, *CENTRE, CENTRE[0] + 0.01, CENTRE[1])
    if test is None:
        print(f"ERROR: OSRM not reachable at {args.osrm_url}")
        sys.exit(1)
    print(f"OSRM reachable at {args.osrm_url}\n")

    if args.dry_run:
        print("DRY RUN — planned OD pairs:")
        for label, lat1, lon1, lat2, lon2 in ods:
            d = haversine_m(lat1, lon1, lat2, lon2) / 1000
            print(f"  {label:28s} ({lat1:.4f},{lon1:.4f}) → centre   ~{d:.1f} km straight")
        print("\nNo Google calls made. Re-run without --dry-run to execute.")
        return

    # Key is only needed for a live (cache-miss) call. Fully-cached re-runs (e.g.
    # re-analysing alternatives) make zero API calls and need no key.
    api_key = os.environ.get("GOOGLE_MAPS_API_KEY")

    rows = []          # one record per Google route (best + alternatives)
    n_api_calls = 0
    for label, lat1, lon1, lat2, lon2 in ods:
        ck = os.path.join(CACHE_DIR, cache_key(lat1, lon1, lat2, lon2, args.traffic))
        if os.path.exists(ck) and not args.refresh:
            with open(ck) as f:
                gdata = json.load(f)
            cached = True
        else:
            if not api_key:
                print(f"  {label}: not cached and GOOGLE_MAPS_API_KEY not set — "
                      f"skipping (no live call made).")
                continue
            try:
                gdata = google_routes(api_key, lat1, lon1, lat2, lon2, args.traffic)
            except urllib.error.HTTPError as e:
                print(f"  {label}: Google HTTP {e.code}: {e.read()[:200].decode(errors='replace')}")
                continue
            except Exception as e:
                print(f"  {label}: Google error: {e}")
                continue
            with open(ck, "w") as f:
                json.dump(gdata, f)
            n_api_calls += 1
            cached = False
            time.sleep(0.1)

        groutes = gdata.get("routes", [])
        if not groutes:
            print(f"  {label}: no Google route ({gdata})")
            continue

        # OSRM's own chosen route (once per OD) — basis for route-overlap / TEends.
        r = osrm_route(args.osrm_url, lat1, lon1, lat2, lon2)
        if r is None:
            print(f"  {label}: OSRM /route failed")
            continue
        osrm_route_dur, osrm_route_dist, route_nodes = r

        # Match EVERY Google route (best + alternatives) — each is a free
        # (geometry, duration) timing-calibration pair at no extra API cost.
        for j, route in enumerate(groutes):
            gdur = parse_google_duration(route["duration"])
            gdist = route.get("distanceMeters", 0)
            gcoords = decode_polyline(route["polyline"]["encodedPolyline"])
            gsamp = downsample_by_distance(gcoords)

            m = osrm_match(args.osrm_url, gsamp)
            if m is None:
                match_dur, match_nodes, conf, n_m = None, [], 0.0, 0
            else:
                match_dur, match_nodes, conf, n_m = m

            te_matched = (match_dur / gdur) if match_dur else None
            valid = (match_dur is not None) and (conf >= CONF_MIN)

            # route-overlap & TEends only meaningful for the best route vs OSRM's own route
            overlap = te_endpoint = None
            if j == 0:
                te_endpoint = osrm_route_dur / gdur
                if match_nodes:
                    sm, sr = set(match_nodes), set(route_nodes)
                    overlap = len(sm & sr) / len(sm)

            rows.append({
                "label": label, "route_idx": j, "is_best": j == 0,
                "cached": cached, "n_alts": len(groutes),
                "google_dur": gdur, "google_dist": gdist,
                "osrm_route_dur": osrm_route_dur if j == 0 else None,
                "osrm_match_dur": match_dur, "match_conf": conf, "n_matchings": n_m,
                "te_matched": te_matched, "te_matched_valid": valid,
                "te_endpoint": te_endpoint, "route_overlap": overlap,
            })

    # ── Report ──────────────────────────────────────────────────────────────
    n_routes = len(rows)
    n_valid = sum(1 for r in rows if r["te_matched_valid"])
    n_ods = len({r["label"] for r in rows})
    print(f"\n{n_api_calls} live Google calls; {n_routes} Google routes across "
          f"{n_ods} OD pairs ({n_valid} passed conf>={CONF_MIN})\n")
    hdr = (f"{'label':22s} {'rt':>2s} {'gGoog':>6s} {'oMatch':>6s} {'oRoute':>6s} "
           f"{'TEmatch':>7s} {'TEends':>7s} {'overlap':>7s} {'conf':>5s}")
    print(hdr)
    print("-" * len(hdr))
    def fmt(v, n=2):
        return "   -  " if v is None else f"{v:.{n}f}"
    for r in rows:
        flag = "" if r["te_matched_valid"] else " <lowconf>"
        print(f"{r['label']:22s} {r['route_idx']:2d} "
              f"{r['google_dur']/60:6.1f} "
              f"{(r['osrm_match_dur'] or 0)/60:6.1f} "
              f"{(r['osrm_route_dur']/60) if r['osrm_route_dur'] else 0:6.1f} "
              f"{fmt(r['te_matched']):>7s} {fmt(r['te_endpoint']):>7s} "
              f"{fmt(r['route_overlap']):>7s} {fmt(r['match_conf']):>5s}{flag}")
    print("\nrt = Google route index (0=best, 1+=alternative). gGoog/oMatch/oRoute = minutes.")
    print("  TEmatch = oMatch/gGoog (pure time error per route; <1 = OSRM too fast).")
    print("  TEends/overlap shown for best route only (vs OSRM's own route).")
    print(f"  <lowconf> = /match confidence < {CONF_MIN}; excluded from aggregate.")

    # Aggregate over all confidence-valid routes (best + alternatives)
    tem = sorted(r["te_matched"] for r in rows if r["te_matched_valid"])
    if tem:
        med = tem[len(tem) // 2]
        mean = sum(tem) / len(tem)
        print(f"\nMatched-geometry time error over {len(tem)} valid routes (best+alts): "
              f"median {med:.2f}, mean {mean:.2f}, min {min(tem):.2f}, max {max(tem):.2f}")
    n_match_fail = sum(1 for r in rows if r["osrm_match_dur"] is None)
    n_lowconf = n_routes - n_valid - n_match_fail
    print(f"/match: {n_match_fail} hard failures, {n_lowconf} low-confidence "
          f"(<{CONF_MIN}) out of {n_routes} routes")

    summ = os.path.join(CACHE_DIR, "feasibility_summary.json")
    with open(summ, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"\nWrote {summ}")


if __name__ == "__main__":
    main()
