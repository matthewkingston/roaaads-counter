"""
Generate external zone links for the Newtownards model via local OSRM routing.

Queries a local OSRM instance to determine:
  1. External→boundary links (X→B): valid directed edges from each external node to
     each boundary node (discarded if another boundary node is hit first en route).
  2. Boundary→external links (B→X): valid directed edges from each boundary node to
     each external node (discarded if the route re-enters the core before exiting).
  3. Boundary→boundary exterior shortcuts (B1→B2): when the shortest B1→B2 route goes
     outside the core, a direct link is added with the OSRM segment duration.
  4. External→external through-route allowlist: (X1, X2) pairs whose OSRM route passes
     through any boundary node (i.e., the route genuinely transits the core area).
  5. External→external non-through times: (X1, X2) pairs whose OSRM route does NOT
     transit the core — recorded with their direct OSRM time as denominator-only
     virtual edges (production-constrained gravity needs these so a distant zone's
     per-capita trip budget isn't entirely dumped into the core; they carry no flow).

Requires:
  - data/census_zones.json          (from build_census_zones.py)
  - simulation/node_weights.json    (from build_demographics.py — needs boundary_node_ids + internal_node_ids)
  - simulation/newtownards_network.graphml (from build_network.py — raw graph, OSM node IDs + WGS84 coords)
  - Local OSRM instance at OSRM_HOST:OSRM_PORT (car profile, NI extract)

Output: data/external_links.json

Start OSRM with:
  docker run -t -i -p 5000:5000 -v $(pwd):/data osrm/osrm-backend \\
    osrm-routed --algorithm mld /data/northern-ireland.osrm
"""

import json, math, time, sys, http.client
import osmnx as ox

CENSUS_ZONES_FILE = "data/census_zones.json"
WEIGHTS_FILE      = "simulation/node_weights.json"
RAW_GRAPH         = "simulation/newtownards_network.graphml"
OUTPUT_FILE       = "data/external_links.json"
OSRM_HOST         = "localhost"
OSRM_PORT         = 5000

# ── Load inputs ────────────────────────────────────────────────────────────────

print("Loading census zones …")
with open(CENSUS_ZONES_FILE) as f:
    census_zones = json.load(f)
external_nodes = census_zones["external_nodes"]
print(f"  {len(external_nodes)} external nodes")

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)
boundary_node_ids = set(weights["boundary_node_ids"])   # OSM node IDs
internal_node_ids = set(weights["internal_node_ids"])   # OSM node IDs
print(f"  {len(boundary_node_ids)} boundary nodes, {len(internal_node_ids)} internal nodes")

# Raw graph has OSM node IDs and WGS84 coordinates (x=lon, y=lat) — no transform needed.
# The consolidated graph renumbers nodes and cannot be used here.
print("Loading raw graph for boundary node positions …")
G_raw = ox.load_graphml(RAW_GRAPH)
boundary_nodes = []
for nid in boundary_node_ids:
    nd = G_raw.nodes[nid]
    boundary_nodes.append({"id": int(nid), "lat": float(nd["y"]), "lon": float(nd["x"])})
boundary_nodes.sort(key=lambda x: x["id"])
print(f"  {len(boundary_nodes)} boundary nodes with WGS84 positions")

# ── OSRM helper ────────────────────────────────────────────────────────────────

_request_count = 0
_t_start = time.time()
_conn = None


def _get_conn():
    global _conn
    if _conn is None:
        _conn = http.client.HTTPConnection(OSRM_HOST, OSRM_PORT, timeout=15)
    return _conn


def osrm_route(lat1, lon1, lat2, lon2, retries=3, want_geometry=False):
    """Query OSRM route from (lat1,lon1) to (lat2,lon2) using a persistent connection.

    Returns (node_sequence, total_duration_s, annotation_durations, snaps, coords)
    or None if no route found.
      node_sequence       — list of OSM node IDs along the route
      total_duration_s    — total trip duration in seconds
      annotation_durations— per-step durations (len = len(node_sequence)-1)
      snaps               — list of {distance: metres from input to snap point}
      coords              — per-node [lon,lat] geometry (1:1 with node_sequence) if
                            want_geometry else None; only the legs that run a crow-flies
                            revisit check (B→X, B→B) request it, to keep X→X/X→B lean.
    """
    global _request_count, _conn
    _ov = "full&geometries=geojson" if want_geometry else "false"
    path = (f"/route/v1/driving/"
            f"{lon1},{lat1};{lon2},{lat2}"
            f"?annotations=nodes,duration&overview={_ov}")
    for attempt in range(retries):
        try:
            conn = _get_conn()
            conn.request("GET", path)
            r = conn.getresponse()
            body = r.read()
            _request_count += 1
            if _request_count % 500 == 0:
                elapsed = time.time() - _t_start
                print(f"  {_request_count} queries in {elapsed:.0f}s "
                      f"({_request_count/elapsed:.0f} q/s)")
            data = json.loads(body)
            if data.get("code") != "Ok":
                return None
            leg = data["routes"][0]["legs"][0]
            nodes = leg["annotation"]["nodes"]
            durs  = leg["annotation"]["duration"]
            total = data["routes"][0]["duration"]
            snaps = [{"distance": w.get("distance", 0)} for w in data["waypoints"]]
            coords = data["routes"][0]["geometry"]["coordinates"] if want_geometry else None
            return nodes, total, durs, snaps, coords
        except (http.client.HTTPException, ConnectionError, json.JSONDecodeError, KeyError):
            _conn = None  # force reconnect on next call
            if attempt < retries - 1:
                time.sleep(0.1)
            else:
                return None


def _first_boundary_in_sequence(node_seq, boundary_ids):
    """Return the index of the first boundary node in node_seq (excluding index 0)."""
    for i, nid in enumerate(node_seq):
        if i == 0:
            continue
        if nid in boundary_ids:
            return i
    return None


def _haversine_m(lat1, lon1, lat2, lon2):
    """Great-circle distance in metres."""
    r = 6371000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = math.radians(lat2 - lat1)
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


# A route may legitimately return to its origin boundary node when that node sits on a
# roundabout/complex junction (OSRM's fine per-arc nodes collapse to a single model node,
# so entering and leaving the junction re-passes the same OSM node) — benign, priced out
# of the leg.  But returning to it after wandering far is a snapping/routing fault → loud.
REVISIT_MAX_LOOP_M = 500.0


def _classify_bb_shortcut(node_seq, ann_durs, total_dur, coords, b1_id, b2_id,
                          boundary_ids, internal_ids):
    """Decide whether an OSRM route B1→B2 is a valid *exterior* boundary→boundary
    shortcut: a route that leaves B1 and reaches B2 as the first boundary node
    encountered WITHOUT passing through the core interior.

    Returns the penalty-inclusive B1→B2 leg duration to KEEP the shortcut, or None to
    DISCARD it.  Raises ValueError on a snapping anomaly that needs human review.

      1. Locate B1.  OSRM snaps the passed coordinate to the nearest *edge* and reports
         both endpoints of that first edge, so the origin boundary node lands at index 0
         OR 1.  Absent, or index > 1 → start-snap anomaly → raise.
      2. B2 must appear somewhere in the route.  If not, the destination coordinate
         snapped elsewhere — most often onto B2's dual-carriageway twin ~14 m away → raise.
      3. Scan forward from just past B1 to the first boundary node, tracking crow-flies
         distance from B1 and whether the route dipped into the core interior:
           - B1 again        → benign junction re-pass IF the route stayed within
                               REVISIT_MAX_LOOP_M of B1 (a roundabout collapses to one
                               model node, so OSRM re-passes it); that loop is priced out
                               of the leg and does not count as a core transit.  A return
                               after a farther excursion is a snapping/routing fault → raise.
           - a different B3   → discard (that trip is covered by B1→B3 + B3→B2 directly);
           - B2               → stop.
      4. If an interior node was seen before B2 → discard (the route went *through* the
         core, not around it).  Else keep.

    Duration = `total_dur` minus the per-edge travel AFTER B2 (strips a destination
    overshoot) minus the per-edge travel up to B1's *last* occurrence (strips the start
    junction loop).  Using `total_dur` (not a raw annotation sum) keeps OSRM turn/signal
    penalties in the leg.

    Known limitation (accepted): the interior test only sees model nodes, so a boundary→
    boundary road crossing the interior on a *junction-free* corridor (no model junction)
    is not detected as through-core.  Such a road has no junctions (nothing measured on
    it); the only cost is a redundant near-duplicate of the real internal route, mildly
    perturbing the probit split there.  Closing it fully needs a route-polyline vs
    core-polygon test.
    """
    # 1. locate the origin boundary node — must be at index 0 or 1
    try:
        b1_pos = node_seq.index(b1_id)
    except ValueError:
        raise ValueError(
            f"[bb-shortcut] B1 {b1_id} absent from OSRM route to B2 {b2_id} "
            f"(start-snap anomaly — review)")
    if b1_pos > 1:
        raise ValueError(
            f"[bb-shortcut] B1 {b1_id} at position {b1_pos} (>1) in route to B2 {b2_id} "
            f"(start-snap anomaly — review)")

    # 2. the destination boundary node must appear at all
    if b2_id not in node_seq:
        raise ValueError(
            f"[bb-shortcut] B2 {b2_id} absent from OSRM route from B1 {b1_id} "
            f"(dest-snap anomaly, likely a carriageway-twin snap — review)")

    # 3-4. scan to the first boundary node; track excursion distance + interior incursion
    # coords fetched lazily by the caller only when B1 re-appears; None on the common path.
    b1_lat, b1_lon = (coords[b1_pos][1], coords[b1_pos][0]) if coords is not None else (None, None)
    last_b1 = b1_pos
    max_excursion = 0.0
    entered_core = False
    for i in range(b1_pos + 1, len(node_seq)):
        nid = node_seq[i]
        if nid in boundary_ids:            # boundary tested first (boundary ⊆ internal)
            if nid == b1_id:
                if max_excursion > REVISIT_MAX_LOOP_M:
                    raise ValueError(
                        f"[bb-shortcut] route B1 {b1_id} → B2 {b2_id} returns to B1 at "
                        f"position {i} after a {max_excursion:.0f} m excursion "
                        f"(> {REVISIT_MAX_LOOP_M:.0f} m — review)")
                last_b1 = i                # benign junction loop — price it out, not a transit
                max_excursion = 0.0
                entered_core = False
                continue
            if nid == b2_id:
                if entered_core:
                    return None
                return total_dur - sum(ann_durs[i:]) - sum(ann_durs[:last_b1])
            return None                    # a different boundary first → B1→B3→B2, drop
        if coords is not None:
            max_excursion = max(max_excursion,
                                _haversine_m(b1_lat, b1_lon, coords[i][1], coords[i][0]))
        if nid in internal_ids:
            entered_core = True

    # defensive: B2 present but never reached as the first boundary after B1
    raise ValueError(
        f"[bb-shortcut] route B1 {b1_id} → B2 {b2_id}: scan ended without reaching B2 "
        f"as the first boundary (review)")


def _classify_xb_link(node_seq, ann_durs, total_dur, b_id, boundary_ids, internal_ids):
    """Decide whether an OSRM route X→B yields a valid external→boundary (inbound) link:
    B is the natural entry into the core from external node X (the first boundary node
    reached).  The criterion is pure boundary-node membership, so off-model shape points
    are transparent — only snapping needs guarding.

    Returns the leg duration (X→B, priced *up to B's first appearance*) to KEEP the link,
    or None to DISCARD (B is present but a *different* boundary is X's natural entry —
    that boundary's own X→B′ query keeps it).  Raises ValueError on an anomaly for review.

      - node_seq[0] is X's snapped origin.  If it is an internal/boundary model node the
        external centroid snapped onto/into the core (boundary ⊆ internal) — a
        data-consistency fault → raise.
      - B must appear in the route; if not, the destination coordinate snapped elsewhere
        (e.g. a carriageway twin) or the route never reached B → raise.
      - Duration is `total_dur` minus the per-edge travel time AFTER B, i.e. the time to
        reach B.  Priced this way (not as a raw sum of per-edge annotation durations, which
        OMITS OSRM turn/signal penalties) so it preserves the penalty-inclusive leg time in
        the normal case (nothing after B) while still stripping a destination snap onto the
        wrong side of a one-way, where OSRM overshoots B and loops back.
    """
    if node_seq[0] in internal_ids:                       # boundary ⊆ internal
        raise ValueError(
            f"[X→B] external origin snapped onto core node {node_seq[0]} en route to "
            f"boundary {b_id} (external centroid inside/on the core — review)")
    if b_id not in node_seq:
        raise ValueError(
            f"[X→B] boundary {b_id} absent from OSRM route from external origin "
            f"(dest-snap anomaly, e.g. carriageway twin — review)")
    first_b_idx = _first_boundary_in_sequence(node_seq, boundary_ids)  # B present ⇒ non-None
    if node_seq[first_b_idx] != b_id:
        return None                                       # a different boundary is the entry
    return total_dur - sum(ann_durs[first_b_idx:])        # penalty-inclusive time to reach B


def _classify_bx_link(node_seq, ann_durs, total_dur, coords, b_id, boundary_ids):
    """Decide whether an OSRM route B→X yields a valid boundary→external (outbound) link:
    B is the last boundary departed on the way to external node X (no *other* boundary
    node appears).  Boundary-membership criterion again, so only snapping needs guarding.

    Returns the penalty-inclusive B→X leg duration to KEEP, or None to DISCARD (another
    boundary B′ appears → the trip is covered by B→B′ + B′→X).  Raises ValueError on a
    snapping anomaly for review.

      - B (the origin) must sit at index 0 or 1 (OSRM snaps the origin coordinate to the
        nearest edge, reporting both endpoints).  Absent, or index > 1 → raise.
      - B may re-appear when it sits on a roundabout/complex junction (which collapses to
        one model node, so OSRM re-passes the same OSM node): benign IF the route stayed
        within REVISIT_MAX_LOOP_M of B, and that start loop is priced out of the leg; a
        return after a farther excursion is a snapping/routing fault → raise.
    Duration = `total_dur` minus the per-edge travel up to B's *last* occurrence (strips
    the start junction loop), keeping OSRM turn/signal penalties.  X, the destination, is
    an external centroid — guarded by the Step-1 X→B check on the same centroid.
    """
    try:
        b_pos = node_seq.index(b_id)
    except ValueError:
        raise ValueError(
            f"[B→X] origin boundary {b_id} absent from its own OSRM route "
            f"(start-snap anomaly — review)")
    if b_pos > 1:
        raise ValueError(
            f"[B→X] origin boundary {b_id} at position {b_pos} (>1) "
            f"(start-snap anomaly — review)")
    # coords is only needed for the crow-flies revisit check; the caller fetches geometry
    # (a second OSRM query) lazily, ONLY when B re-appears, so coords is None on the common
    # no-revisit path (max_excursion stays 0 and is never consulted).
    b_lat, b_lon = (coords[b_pos][1], coords[b_pos][0]) if coords is not None else (None, None)
    last_b = b_pos
    max_excursion = 0.0
    for i in range(b_pos + 1, len(node_seq)):
        nid = node_seq[i]
        if nid == b_id:
            if max_excursion > REVISIT_MAX_LOOP_M:
                raise ValueError(
                    f"[B→X] route from boundary {b_id} returns to it at position {i} "
                    f"after a {max_excursion:.0f} m excursion "
                    f"(> {REVISIT_MAX_LOOP_M:.0f} m — review)")
            last_b = i                         # benign junction loop — price it out
            max_excursion = 0.0
            continue
        if nid in boundary_ids:                # a different boundary B′ appears
            return None                        # B→B′ + B′→X covers it → discard
        if coords is not None:
            max_excursion = max(max_excursion,
                                _haversine_m(b_lat, b_lon, coords[i][1], coords[i][0]))
    return total_dur - sum(ann_durs[:last_b])  # keep; strip the start junction loop


def _ext_ext_transits_core(node_seq, boundary_ids, internal_ids):
    """Classify an external→external OSRM route: True if it transits the core (passes a
    boundary node → a flow-carrying through-pair), False if it stays outside (a
    denominator-only pair).  The boundary-membership criterion is robust to off-model
    shape points, and a pressure test found no through-route that dodges every boundary
    node, so only endpoint snapping is guarded here.

    Raises ValueError if either external endpoint (node_seq[0] origin, node_seq[-1]
    destination) snapped onto a core node (boundary ⊆ internal) — its centroid is
    inside/on the core, a data-consistency fault to review.  Mirrors the X→B origin
    guard; guarding both ends here additionally catches a boundary-straddling snap edge
    whose exterior endpoint let the same centroid pass its Step-1 origin check.
    """
    if node_seq[0] in internal_ids:
        raise ValueError(
            f"[X→X] origin external centroid snapped onto core node {node_seq[0]} "
            f"(centroid inside/on the core — review)")
    if node_seq[-1] in internal_ids:
        raise ValueError(
            f"[X→X] destination external centroid snapped onto core node {node_seq[-1]} "
            f"(centroid inside/on the core — review)")
    return any(nid in boundary_ids for nid in node_seq)


# ── Check OSRM is reachable ────────────────────────────────────────────────────

print(f"\nChecking OSRM at {OSRM_HOST}:{OSRM_PORT} …")
try:
    test = osrm_route(54.5933, -5.6960, 54.5933, -5.6960)
except Exception:
    test = None
if test is None:
    try:
        c = http.client.HTTPConnection(OSRM_HOST, OSRM_PORT, timeout=5)
        c.request("GET", "/route/v1/driving/-5.696,54.593;-5.696,54.593")
        c.getresponse().read()
        c.close()
    except Exception as e:
        print(f"ERROR: Cannot reach OSRM at {OSRM_HOST}:{OSRM_PORT}: {e}")
        print("Start OSRM with:")
        print("  docker run -t -i -p 5000:5000 -v $(pwd):/data osrm/osrm-backend \\")
        print("    osrm-routed --algorithm mld /data/northern-ireland.osrm")
        sys.exit(1)
print("  OSRM reachable")

# ── Step 1: External→boundary links (X→B) ────────────────────────────────────

print(f"\nStep 1: External→boundary links ({len(external_nodes)} × {len(boundary_nodes)} = "
      f"{len(external_nodes)*len(boundary_nodes):,} queries) …")

ext_boundary_links = []   # {"from_ext": id, "to_boundary": osm_id, "duration_s": float}
snap_warnings = []

for xi, ext in enumerate(external_nodes):
    xlat, xlon = ext["centroid_lat"], ext["centroid_lon"]
    xid = ext["id"]

    # Check snap distance on first query to this centroid
    first_snap = None

    for bi, bnode in enumerate(boundary_nodes):
        blat, blon = bnode["lat"], bnode["lon"]
        bid = bnode["id"]

        result = osrm_route(xlat, xlon, blat, blon)
        if result is None:
            continue
        node_seq, total_dur, ann_durs, snaps, coords = result

        # Log snap distance for first boundary query per external node
        if first_snap is None:
            first_snap = snaps[0]["distance"]
            if first_snap > 2000:
                snap_warnings.append((ext["code"], ext["level"], first_snap))

        if len(node_seq) < 2:
            continue

        # Keep only if B is the natural entry (first boundary reached); price up to B.
        # Fails loud on snap faults (external origin snapped into the core, or B absent).
        dur = _classify_xb_link(node_seq, ann_durs, total_dur, bid,
                                 boundary_node_ids, internal_node_ids)
        if dur is None:
            continue  # B present but a different boundary is X's natural entry
        ext_boundary_links.append({
            "from_ext":    xid,
            "to_boundary": bid,
            "duration_s":  round(dur, 2),
        })

    if (xi + 1) % 10 == 0:
        print(f"  {xi+1}/{len(external_nodes)} external nodes  "
              f"({len(ext_boundary_links)} links so far)")

print(f"  {len(ext_boundary_links)} X→B links")

if snap_warnings:
    print(f"  Warning: {len(snap_warnings)} external centroids snapped >2km from road:")
    for code, level, dist in snap_warnings[:5]:
        print(f"    {level} {code}: {dist:.0f}m snap distance")

# ── Step 2: Boundary→external links (B→X) ────────────────────────────────────

print(f"\nStep 2: Boundary→external links ({len(boundary_nodes)} × {len(external_nodes)} = "
      f"{len(boundary_nodes)*len(external_nodes):,} queries) …")

bnd_external_links = []   # {"from_boundary": osm_id, "to_ext": id, "duration_s": float}

for bi, bnode in enumerate(boundary_nodes):
    blat, blon = bnode["lat"], bnode["lon"]
    bid = bnode["id"]

    for xi, ext in enumerate(external_nodes):
        xlat, xlon = ext["centroid_lat"], ext["centroid_lon"]
        xid = ext["id"]

        result = osrm_route(blat, blon, xlat, xlon)
        if result is None:
            continue
        node_seq, total_dur, ann_durs, snaps, coords = result

        if len(node_seq) < 2:
            continue

        # B re-appears (e.g. it sits on a roundabout) → re-query WITH geometry so the
        # crow-flies revisit check has per-node coords.  Rare, so the common no-revisit
        # path stays geometry-free (geometry adds ~1.5× per query).
        if node_seq.count(bid) > 1:
            result = osrm_route(blat, blon, xlat, xlon, want_geometry=True)
            if result is None:
                continue
            node_seq, total_dur, ann_durs, snaps, coords = result

        # Keep only if B is the last boundary departed (no other boundary appears);
        # B→B'+B'→X covers the rest.  Fails loud on origin-snap faults / far revisits.
        dur = _classify_bx_link(node_seq, ann_durs, total_dur, coords, bid, boundary_node_ids)
        if dur is None:
            continue
        bnd_external_links.append({
            "from_boundary": bid,
            "to_ext":        xid,
            "duration_s":    round(dur, 2),
        })

    if (bi + 1) % 5 == 0:
        print(f"  {bi+1}/{len(boundary_nodes)} boundary nodes  "
              f"({len(bnd_external_links)} links so far)")

print(f"  {len(bnd_external_links)} B→X links")

# ── Step 3: Boundary→boundary exterior shortcuts ──────────────────────────────

print(f"\nStep 3: Boundary→boundary shortcuts "
      f"({len(boundary_nodes)*(len(boundary_nodes)-1):,} queries) …")

boundary_boundary_links = []   # {"from": osm_id, "to": osm_id, "duration_s": float}

for b1 in boundary_nodes:
    for b2 in boundary_nodes:
        if b1["id"] == b2["id"]:
            continue

        result = osrm_route(b1["lat"], b1["lon"], b2["lat"], b2["lon"])
        if result is None:
            continue
        node_seq, total_dur, ann_durs, snaps, coords = result

        if len(node_seq) < 2:
            continue

        # B1 re-appears (e.g. it sits on a roundabout) → re-query WITH geometry for the
        # crow-flies revisit check.  Rare, so the common path stays geometry-free.
        if node_seq.count(b1["id"]) > 1:
            result = osrm_route(b1["lat"], b1["lon"], b2["lat"], b2["lon"], want_geometry=True)
            if result is None:
                continue
            node_seq, total_dur, ann_durs, snaps, coords = result

        # Keep only genuine EXTERIOR shortcuts: B2 is the first boundary reached after
        # leaving B1, and the route never dips into the core interior.  Handles the
        # start/dest coordinate-snap cases explicitly and fails loud on anomalies.
        seg_dur = _classify_bb_shortcut(
            node_seq, ann_durs, total_dur, coords, b1["id"], b2["id"],
            boundary_node_ids, internal_node_ids)
        if seg_dur is None:
            continue

        boundary_boundary_links.append({
            "from":       b1["id"],
            "to":         b2["id"],
            "duration_s": round(seg_dur, 2),
        })

print(f"  {len(boundary_boundary_links)} boundary→boundary exterior shortcuts")

# ── Step 4: External→external through-route check + non-through times ─────────
# Each ordered external pair is classified by whether its OSRM route transits the
# core (passes a boundary node):
#   through     → added to allowed_through_pairs; the path distance is computed by
#                 Dijkstra on the augmented graph in build_paths.py (X→B→…→B'→X'),
#                 so it carries flow AND a denominator term. No time stored here.
#   non-through → recorded in external_external_times as a denominator-only virtual
#                 edge with the direct OSRM time. These pairs carry NO flow across
#                 observed/core links (production-constrained denominator only), but
#                 are needed so a distant external zone's per-capita trip budget is
#                 not entirely dumped into the core. See the production-constrained
#                 gravity plan (inter-external links: "necessary, not optional").

n_ext = len(external_nodes)
print(f"\nStep 4: External→external through-route check ({n_ext*(n_ext-1):,} queries) …")

allowed_through_pairs = {}      # {src_id: [dst_id, ...]}        (through-routed, flow + denom)
external_external_times = {}    # {src_id: {dst_id: duration_s}} (non-through, denom-only)

for xi, ext1 in enumerate(external_nodes):
    dsts  = []
    times = {}
    for xj, ext2 in enumerate(external_nodes):
        if xi == xj:
            continue

        result = osrm_route(ext1["centroid_lat"], ext1["centroid_lon"],
                            ext2["centroid_lat"], ext2["centroid_lon"])
        if result is None:
            continue   # no road route at all (e.g. across water) → no denom term
        node_seq, total_dur, ann_durs, snaps, coords = result

        if _ext_ext_transits_core(node_seq, boundary_node_ids, internal_node_ids):
            dsts.append(ext2["id"])
        else:
            # Fastest real route does not enter the core → denominator-only virtual edge.
            times[ext2["id"]] = round(total_dur, 2)

    if dsts:
        allowed_through_pairs[ext1["id"]] = dsts
    if times:
        external_external_times[ext1["id"]] = times

    if (xi + 1) % 10 == 0:
        n_pairs   = sum(len(v) for v in allowed_through_pairs.values())
        n_nonthru = sum(len(v) for v in external_external_times.values())
        print(f"  {xi+1}/{n_ext} external nodes  "
              f"({n_pairs} through, {n_nonthru} non-through pairs so far)")

n_through  = sum(len(v) for v in allowed_through_pairs.values())
n_nonthru  = sum(len(v) for v in external_external_times.values())
print(f"  {n_through} allowed through-route pairs ({len(allowed_through_pairs)} sources)")
print(f"  {n_nonthru} non-through ext→ext virtual edges ({len(external_external_times)} sources)")

# ── Write output ────────────────────────────────────────────────────────────────

total_queries = _request_count
elapsed = time.time() - _t_start

output = {
    "ext_boundary_links":     ext_boundary_links,
    "bnd_external_links":     bnd_external_links,
    "boundary_boundary_links": boundary_boundary_links,
    "allowed_through_pairs":  {str(k): v for k, v in allowed_through_pairs.items()},
    "external_external_times": {str(k): {str(d): t for d, t in v.items()}
                                for k, v in external_external_times.items()},
    "boundary_node_ids":      sorted(boundary_node_ids),
}

with open(OUTPUT_FILE, "w") as f:
    json.dump(output, f)

print(f"\nSaved {OUTPUT_FILE}")
print(f"  {len(ext_boundary_links)} X→B links")
print(f"  {len(bnd_external_links)} B→X links")
print(f"  {len(boundary_boundary_links)} boundary→boundary shortcuts")
print(f"  {sum(len(v) for v in allowed_through_pairs.values())} allowed through-route pairs")
print(f"  {sum(len(v) for v in external_external_times.values())} non-through ext→ext virtual edges")
print(f"  {total_queries:,} OSRM queries in {elapsed:.0f}s ({total_queries/elapsed:.0f} q/s)")
print(f"\nNext: python3 simulation/build_paths.py")
