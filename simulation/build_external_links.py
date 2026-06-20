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
  - simulation/node_weights.json    (from build_demographics.py — needs boundary_node_ids)
  - simulation/newtownards_consolidated.graphml (from build_network.py)
  - Local OSRM instance at OSRM_URL (car profile, NI extract)

Output: data/external_links.json

Start OSRM with:
  docker run -t -i -p 5000:5000 -v $(pwd):/data osrm/osrm-backend \\
    osrm-routed --algorithm mld /data/northern-ireland.osrm
"""

import json, math, time, sys, urllib.request, urllib.parse, urllib.error
import osmnx as ox

CENSUS_ZONES_FILE = "data/census_zones.json"
WEIGHTS_FILE      = "simulation/node_weights.json"
CONS_GRAPH        = "simulation/newtownards_consolidated.graphml"
OUTPUT_FILE       = "data/external_links.json"
OSRM_URL          = "http://localhost:5000"

# ── Load inputs ────────────────────────────────────────────────────────────────

print("Loading census zones …")
with open(CENSUS_ZONES_FILE) as f:
    census_zones = json.load(f)
external_nodes = census_zones["external_nodes"]
print(f"  {len(external_nodes)} external nodes")

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    weights = json.load(f)
boundary_node_ids = set(weights["boundary_node_ids"])
print(f"  {len(boundary_node_ids)} boundary nodes")

print("Loading consolidated graph …")
G = ox.load_graphml(CONS_GRAPH)
# Get lat/lon positions of boundary nodes
boundary_nodes = []
for nid in boundary_node_ids:
    nd = G.nodes[nid]
    # osmnx graphml stores x=longitude, y=latitude in WGS84 — no transform needed
    lon = float(nd["x"])
    lat = float(nd["y"])
    boundary_nodes.append({"id": int(nid), "lat": lat, "lon": lon})
boundary_nodes.sort(key=lambda x: x["id"])
print(f"  {len(boundary_nodes)} boundary nodes with positions")

# Set of all internal (including boundary) node IDs — used for routing checks
all_node_ids = set(int(n) for n in G.nodes())

# ── OSRM helper ────────────────────────────────────────────────────────────────

_request_count = 0
_t_start = time.time()


def osrm_route(lat1, lon1, lat2, lon2, retries=3):
    """Query OSRM route from (lat1,lon1) to (lat2,lon2).

    Returns (node_sequence, total_duration_s, annotation_durations, snaps)
    or None if no route found.
      node_sequence       — list of OSM node IDs along the route
      total_duration_s    — total trip duration in seconds
      annotation_durations— per-step durations (len = len(node_sequence)-1)
      snaps               — list of {distance: metres from input to snap point}
    """
    global _request_count
    url = (f"{OSRM_URL}/route/v1/driving/"
           f"{lon1},{lat1};{lon2},{lat2}"
           f"?annotations=nodes,duration&overview=false")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as r:
                data = json.loads(r.read())
            _request_count += 1
            if _request_count % 500 == 0:
                elapsed = time.time() - _t_start
                print(f"  {_request_count} queries in {elapsed:.0f}s "
                      f"({_request_count/elapsed:.0f} q/s)")
            if data.get("code") != "Ok":
                return None
            leg = data["routes"][0]["legs"][0]
            nodes = leg["annotation"]["nodes"]
            durs  = leg["annotation"]["duration"]
            total = data["routes"][0]["duration"]
            snaps = [{"distance": w.get("distance", 0)} for w in data["waypoints"]]
            return nodes, total, durs, snaps
        except (urllib.error.URLError, json.JSONDecodeError, KeyError) as e:
            if attempt < retries - 1:
                time.sleep(0.5)
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

print(f"\nChecking OSRM at {OSRM_URL} …")
try:
    test = osrm_route(54.5933, -5.6960, 54.5933, -5.6960)
except Exception:
    test = None
if test is None:
    # Simple connectivity test
    try:
        with urllib.request.urlopen(f"{OSRM_URL}/route/v1/driving/-5.696,54.593;-5.696,54.593",
                                    timeout=5) as r:
            pass
    except Exception as e:
        print(f"ERROR: Cannot reach OSRM at {OSRM_URL}: {e}")
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

        # Check: is the first node after B (index 1) outside the core (internal set)?
        # We check if it's NOT in the full internal node set.
        first_next = node_seq[1] if len(node_seq) > 1 else None
        if first_next is None:
            continue
        if first_next in all_node_ids:
            # Route immediately re-enters core — B is not the natural exit for X.
            continue

        # B opens directly onto the external network toward X — keep B→X
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

        # Check if first node after B1 is outside the internal network
        first_next = node_seq[1]
        if first_next in all_node_ids:
            # Route goes through the core — no external shortcut needed
            continue

        # Route exits core first — find the first boundary node hit (may not be B2)
        first_b_idx = _first_boundary_in_sequence(node_seq, boundary_node_ids)
        if first_b_idx is None:
            continue

        b_mid_nid = node_seq[first_b_idx]
        seg_dur = _duration_up_to_index(ann_durs, first_b_idx)

        boundary_boundary_links.append({
            "from":       b1["id"],
            "to":         b_mid_nid,
            "duration_s": round(seg_dur, 2),
        })

print(f"  {len(boundary_boundary_links)} boundary→boundary exterior shortcuts")

# ── Step 4: External→external through-route allowlist ─────────────────────────

n_ext = len(external_nodes)
print(f"\nStep 4: External→external through-route check ({n_ext*(n_ext-1):,} queries) …")

allowed_through_pairs = []   # [ext_id1, ext_id2]

for xi, ext1 in enumerate(external_nodes):
    for xj, ext2 in enumerate(external_nodes):
        if xi == xj:
            continue

        result = osrm_route(ext1["centroid_lat"], ext1["centroid_lon"],
                            ext2["centroid_lat"], ext2["centroid_lon"])
        if result is None:
            continue
        node_seq, total_dur, ann_durs, snaps = result

        # Does route pass through any boundary node?
        hits_core = any(nid in boundary_node_ids for nid in node_seq)
        if hits_core:
            allowed_through_pairs.append([ext1["id"], ext2["id"]])

    if (xi + 1) % 10 == 0:
        print(f"  {xi+1}/{n_ext} external nodes  ({len(allowed_through_pairs)} through pairs so far)")

print(f"  {len(allowed_through_pairs)} allowed through-route pairs")

# ── Write output ────────────────────────────────────────────────────────────────

total_queries = _request_count
elapsed = time.time() - _t_start

output = {
    "ext_boundary_links":     ext_boundary_links,
    "bnd_external_links":     bnd_external_links,
    "boundary_boundary_links": boundary_boundary_links,
    "allowed_through_pairs":  allowed_through_pairs,
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
