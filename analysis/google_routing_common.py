"""
Shared helpers for the Google-routing calibration tooling:
  - build_od_manifest.py   (writes the fixed OD sample; no API spend)
  - google_query_routes.py (resumable runner that queries Google + OSRM /match)

Pure stdlib (urllib/http) so there is no third-party dependency beyond what the
runner's caller already has. See memory project_google_routing_calibration and
the feasibility pilot analysis/google_feasibility.py.
"""

import json, math, urllib.request, urllib.error

GOOGLE_URL = "https://routes.googleapis.com/directions/v2:computeRoutes"

# Drop /match results below this confidence — low confidence yields garbage
# durations (see feasibility: one OD came back at conf 7.6e-5 / te 1.30).
CONF_MIN = 0.5


# ── Geometry ────────────────────────────────────────────────────────────────

def haversine_m(lat1, lon1, lat2, lon2):
    R = 6371000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


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
    """Thin a dense polyline: keep first/last + a point ~every step_m, capped at
    max_pts (OSRM /match coordinate limit) by widening the stride if needed."""
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
    stride = max(1, math.ceil(len(coords) / (max_pts - 1)))
    kept = coords[::stride]
    if kept[-1] != coords[-1]:
        kept.append(coords[-1])
    return kept


# ── OSRM ────────────────────────────────────────────────────────────────────

def osrm_get(base_url, path):
    url = base_url.rstrip("/") + path
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.loads(r.read())


def osrm_route(base_url, lat1, lon1, lat2, lon2):
    """OSRM's freely-chosen route → (duration_s, distance_m, node_list) or None."""
    path = (f"/route/v1/driving/{lon1},{lat1};{lon2},{lat2}"
            f"?annotations=nodes,duration&overview=false")
    try:
        data = osrm_get(base_url, path)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
        return None
    if data.get("code") != "Ok":
        return None
    rt = data["routes"][0]
    return rt["duration"], rt["distance"], rt["legs"][0]["annotation"]["nodes"]


def osrm_match(base_url, coords, radius_m=40):
    """Map-match a (lat,lon) trace onto OSM via OSRM /match → OSRM's own duration +
    node sequence along that geometry. Returns (duration_s, node_list, confidence,
    n_matchings) or None on failure."""
    pts = ";".join(f"{lon},{lat}" for lat, lon in coords)
    radii = ";".join(str(radius_m) for _ in coords)
    path = (f"/match/v1/driving/{pts}"
            f"?annotations=nodes,duration&overview=false&geometries=geojson"
            f"&gaps=ignore&tidy=true&radiuses={radii}")
    try:
        data = osrm_get(base_url, path)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
        return None
    if data.get("code") != "Ok" or not data.get("matchings"):
        return None
    total_dur, nodes, confs = 0.0, [], []
    for m in data["matchings"]:
        confs.append(m.get("confidence", 0.0))
        for leg in m.get("legs", []):
            ann = leg.get("annotation", {})
            total_dur += sum(ann.get("duration", []))
            nodes.extend(ann.get("nodes", []))
    conf = sum(confs) / len(confs) if confs else 0.0
    return total_dur, nodes, conf, len(data["matchings"])


# ── Google Routes API ─────────────────────────────────────────────────────────

def google_routes(api_key, lat1, lon1, lat2, lon2, traffic=False, departure_iso=None):
    """Query Google Routes API v2 computeRoutes (with alternatives). Returns parsed
    JSON dict, or raises urllib.error.HTTPError / URLError."""
    body = {
        "origin":      {"location": {"latLng": {"latitude": lat1, "longitude": lon1}}},
        "destination": {"location": {"latLng": {"latitude": lat2, "longitude": lon2}}},
        "travelMode":  "DRIVE",
        "computeAlternativeRoutes": True,
        "routingPreference": "TRAFFIC_AWARE" if traffic else "TRAFFIC_UNAWARE",
        "units": "METRIC",
    }
    if traffic and departure_iso:
        body["departureTime"] = departure_iso
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
