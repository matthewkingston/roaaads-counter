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


def osrm_route(lat1, lon1, lat2, lon2, retries=3):
    """Query OSRM route from (lat1,lon1) to (lat2,lon2) using a persistent connection.

    Returns (node_sequence, total_duration_s, annotation_durations, snaps)
    or None if no route found.
      node_sequence       — list of OSM node IDs along the route
      total_duration_s    — total trip duration in seconds
      annotation_durations— per-step durations (len = len(node_sequence)-1)
      snaps               — list of {distance: metres from input to snap point}
    """
    global _request_count, _conn
    path = (f"/route/v1/driving/"
            f"{lon1},{lat1};{lon2},{lat2}"
            f"?annotations=nodes,duration&overview=false")
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
            return nodes, total, durs, snaps
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


def _classify_bb_shortcut(node_seq, ann_durs, total_dur, b1_id, b2_id, boundary_ids, internal_ids):
    """Decide whether an OSRM route B1→B2 is a valid *exterior* boundary→boundary
    shortcut: a route that leaves B1 and reaches B2 as the first boundary node
    encountered WITHOUT passing through the core interior.

    Returns the penalty-inclusive B1→B2 leg duration (`total_dur` minus the per-edge
    travel time AFTER B2 — a raw sum of annotation durations would drop OSRM turn/signal
    penalties, and also strips a destination snap that overshoots B2 and loops back) if
    the shortcut should be KEPT, or None to DISCARD it.  Raises ValueError on a snapping
    anomaly that needs human review — deliberately loud, per design.

      1. Locate B1.  OSRM snaps the passed coordinate to the nearest *edge* and reports
         both endpoints of that first edge, so the origin boundary node lands at index 0
         OR 1.  Absent, or index > 1 → start-snap anomaly → raise.
      2. B2 must appear somewhere in the route.  If not, the destination coordinate
         snapped elsewhere — most often onto B2's dual-carriageway twin ~14 m away — which
         needs review rather than a silent drop → raise.  (We'd likely want to KEEP the
         shortcut in that case, but for now surface it.)
      3. Scan forward from just past B1 to the first boundary node, tracking whether the
         route dipped into the core interior:
           - B1 again        → raise (a shortest route must not revisit its origin);
           - a different B3   → discard (that trip is covered by B1→B3 + B3→B2 directly);
           - B2               → stop.
      4. If an interior node was seen before B2 → discard (the route went *through* the
         core, not around it).  Else keep.

    Known limitation (accepted): the interior test only sees nodes present in the model's
    node set.  A boundary→boundary road crossing the core interior on a *junction-free*
    corridor touches only off-model shape points and no model junction, so it is not
    detected as through-core and would be recorded as an exterior shortcut.  Such a road
    has no junctions (nothing is measured on it); the only cost is a redundant
    near-duplicate of the real internal route, which mildly perturbs the probit
    route-choice split there.  Closing this fully would need a geometric
    (route-polyline vs core-polygon) test.
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

    # 3-4. scan to the first boundary node, tracking interior incursion
    entered_core = False
    for i in range(b1_pos + 1, len(node_seq)):
        nid = node_seq[i]
        if nid in boundary_ids:            # boundary tested first (boundary ⊆ internal)
            if nid == b1_id:
                raise ValueError(
                    f"[bb-shortcut] route B1 {b1_id} → B2 {b2_id} revisits B1 at "
                    f"position {i} (review)")
            if nid == b2_id:
                return None if entered_core else total_dur - sum(ann_durs[i:])
            return None                    # a different boundary first → B1→B3→B2, drop
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


def _classify_bx_link(node_seq, b_id, boundary_ids):
    """Decide whether an OSRM route B→X yields a valid boundary→external (outbound) link:
    B is the last boundary departed on the way to external node X (no *other* boundary
    node appears).  Boundary-membership criterion again, so only snapping needs guarding.

    Returns True to KEEP, or None to DISCARD (another boundary B′ appears → the trip is
    covered by B→B′ + B′→X).  Raises ValueError on a snapping anomaly for review.

      - B (the origin) must sit at index 0 or 1 (OSRM snaps the origin coordinate to the
        nearest edge, reporting both endpoints).  Absent, or index > 1 → raise.
      - B must not reappear after its start (a shortest route revisiting its origin) → raise.
    (X, the destination, is an external centroid — not a node — so its snap is guarded by
    the Step-1 X→B check on the same centroid, which runs first.)
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
    for i in range(b_pos + 1, len(node_seq)):
        nid = node_seq[i]
        if nid == b_id:
            raise ValueError(
                f"[B→X] route from boundary {b_id} revisits it at position {i} (review)")
        if nid in boundary_ids:                # a different boundary B′ appears
            return None                        # B→B′ + B′→X covers it → discard
    return True                                # B is the sole/last boundary → keep


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
        node_seq, total_dur, ann_durs, snaps = result

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
        node_seq, total_dur, ann_durs, snaps = result

        if len(node_seq) < 2:
            continue

        # Keep only if B is the last boundary departed (no other boundary appears);
        # B→B'+B'→X covers the rest.  Fails loud on origin-snap faults / B revisits.
        if _classify_bx_link(node_seq, bid, boundary_node_ids) is None:
            continue
        bnd_external_links.append({
            "from_boundary": bid,
            "to_ext":        xid,
            "duration_s":    round(total_dur, 2),
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
        node_seq, total_dur, ann_durs, snaps = result

        if len(node_seq) < 2:
            continue

        # Keep only genuine EXTERIOR shortcuts: B2 is the first boundary reached after
        # leaving B1, and the route never dips into the core interior.  Handles the
        # start/dest coordinate-snap cases explicitly and fails loud on anomalies.
        seg_dur = _classify_bb_shortcut(
            node_seq, ann_durs, total_dur, b1["id"], b2["id"],
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
        node_seq, total_dur, ann_durs, snaps = result

        if any(nid in boundary_node_ids for nid in node_seq):
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
