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


def _duration_up_to_index(annotation_durs, idx):
    """Sum annotation durations for steps 0..idx-1 (travel from node[0] to node[idx])."""
    return sum(annotation_durs[:idx])


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

        # Find first boundary node in sequence (after the start)
        first_b_idx = _first_boundary_in_sequence(node_seq, boundary_node_ids)
        if first_b_idx is None:
            continue  # route doesn't reach core

        first_b_nid = node_seq[first_b_idx]
        if first_b_nid == bid:
            # B is the natural entry point — keep X→B link
            ext_boundary_links.append({
                "from_ext":    xid,
                "to_boundary": bid,
                "duration_s":  round(total_dur, 2),
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

        # Symmetric with X→B: discard if any other boundary node appears in the route.
        # That journey is covered by B→B' + B'→X, so B→X would be redundant.
        # (A route that re-enters the core must also pass another boundary node,
        # so this single check subsumes the old internal-node check.)
        if any(nid in boundary_node_ids and nid != bid for nid in node_seq[1:]):
            continue

        # B is the last boundary node departed on the way to X — keep B→X
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

        # Only keep if the route exits core first (first node after B1 is external)
        if node_seq[1] in internal_node_ids:
            continue

        # Symmetric with X→B: keep only if B2 is the first boundary node hit.
        # If an intermediate B_mid is hit first, discard — B1→B_mid will be
        # found and kept when querying the (B1, B_mid) pair directly.
        first_b_idx = _first_boundary_in_sequence(node_seq, boundary_node_ids)
        if first_b_idx is None:
            continue

        if node_seq[first_b_idx] != b2["id"]:
            continue

        seg_dur = _duration_up_to_index(ann_durs, first_b_idx)
        boundary_boundary_links.append({
            "from":       b1["id"],
            "to":         b2["id"],
            "duration_s": round(seg_dur, 2),
        })

print(f"  {len(boundary_boundary_links)} boundary→boundary exterior shortcuts")

# ── Step 4: External→external through-route allowlist ─────────────────────────

n_ext = len(external_nodes)
print(f"\nStep 4: External→external through-route check ({n_ext*(n_ext-1):,} queries) …")

allowed_through_pairs = {}   # {src_id: [dst_id, ...]}

for xi, ext1 in enumerate(external_nodes):
    dsts = []
    for xj, ext2 in enumerate(external_nodes):
        if xi == xj:
            continue

        result = osrm_route(ext1["centroid_lat"], ext1["centroid_lon"],
                            ext2["centroid_lat"], ext2["centroid_lon"])
        if result is None:
            continue
        node_seq, total_dur, ann_durs, snaps = result

        if any(nid in boundary_node_ids for nid in node_seq):
            dsts.append(ext2["id"])

    if dsts:
        allowed_through_pairs[ext1["id"]] = dsts

    if (xi + 1) % 10 == 0:
        n_pairs = sum(len(v) for v in allowed_through_pairs.values())
        print(f"  {xi+1}/{n_ext} external nodes  ({n_pairs} through pairs so far)")

n_through = sum(len(v) for v in allowed_through_pairs.values())
print(f"  {n_through} allowed through-route pairs ({len(allowed_through_pairs)} sources)")

# ── Write output ────────────────────────────────────────────────────────────────

total_queries = _request_count
elapsed = time.time() - _t_start

output = {
    "ext_boundary_links":     ext_boundary_links,
    "bnd_external_links":     bnd_external_links,
    "boundary_boundary_links": boundary_boundary_links,
    "allowed_through_pairs":  {str(k): v for k, v in allowed_through_pairs.items()},
    "boundary_node_ids":      sorted(boundary_node_ids),
}

with open(OUTPUT_FILE, "w") as f:
    json.dump(output, f)

print(f"\nSaved {OUTPUT_FILE}")
print(f"  {len(ext_boundary_links)} X→B links")
print(f"  {len(bnd_external_links)} B→X links")
print(f"  {len(boundary_boundary_links)} boundary→boundary shortcuts")
print(f"  {len(allowed_through_pairs)} allowed through-route pairs")
print(f"  {total_queries:,} OSRM queries in {elapsed:.0f}s ({total_queries/elapsed:.0f} q/s)")
print(f"\nNext: python3 simulation/build_paths.py")
