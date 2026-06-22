"""
Build the Newtownards road network from the local NI .osm.pbf (the same snapshot
OSRM is built from) and produce Level-1 segments: one directed link per road
between adjacent intersection nodes, with complex junctions (roundabouts,
staggered crossings) consolidated into single nodes.

A small drivable-highway extract is streamed out of the pbf with osmium-tool
(Docker) before osmnx reads it, so the full ~400 MB island file is never parsed
in process. Requires Docker (image auto-built from simulation/osmium.Dockerfile).

Outputs (in the same directory):
  newtownards_nodes.geojson        — raw intersection nodes
  newtownards_links.geojson        — raw directed links
  newtownards_network.graphml      — raw graph (geographic CRS)
  newtownards_consolidated.graphml — junction-consolidated graph (projected CRS)
"""

import json, os, sys, subprocess
import osmnx as ox
import geopandas as gpd
import pyproj
from shapely.geometry import Polygon as _Poly, Point as _Pt
from shapely.ops import transform as _shp_transform

# Newtownards town centre — CENTRE lives in zones_config.py (single source).
# CORE_RADIUS is also defined there and is written to census_zones.json by
# build_census_zones.py; the network extent is read from that polygon below.
from zones_config import CENTRE
# PBF source + bbox margin live in demographics_config (single source for the
# OSM snapshot path, shared conceptually with the OSRM build).
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from demographics_config import PBF_PATH, BOUNDARY_BBOX_MARGIN_M
CONSOLIDATION_TOLERANCE_M = 15   # merge nodes within this distance
OUT_DIR   = "simulation"

CENSUS_ZONES_FILE = "data/census_zones.json"

# osmium-tool Docker image used to stream a small extract out of the full NI pbf
# (the 400 MB whole-Ireland file OOMs an in-process parse on modest RAM). Built
# from simulation/osmium.Dockerfile; auto-built here on first run if absent.
OSMIUM_IMAGE = "osmium-roaaads"
OSMIUM_DOCKERFILE = os.path.join(os.path.dirname(__file__), "osmium.Dockerfile")
# Drivable highway values kept by the extract — the positive form of osmnx's
# "drive" network filter (which excludes footway/cycleway/path/service/track/…).
DRIVE_HIGHWAYS = ("motorway,motorway_link,trunk,trunk_link,primary,primary_link,"
                  "secondary,secondary_link,tertiary,tertiary_link,unclassified,"
                  "residential,living_street,road")
_BBOX_PBF   = os.path.join(OUT_DIR, "_pbf_bbox_extract.osm.pbf")   # intermediate
_DRIVE_OSM  = os.path.join(OUT_DIR, "_pbf_drive_extract.osm")      # intermediate (XML)

# ── 1. Extract the road network from the local NI pbf ────────────────────────────
# Source the graph from the same .osm.pbf OSRM is built from (one OSM snapshot →
# boundary/internal node IDs match OSRM route node IDs). The extract is bounded by
# the core polygon buffered by BOUNDARY_BBOX_MARGIN_M, which supersedes the old
# Overpass `dist` circle margin (the buffer only needs to reach boundary nodes'
# external neighbours; 5 km is generous). The consolidated routing graph is still
# clipped to the core polygon in step 3, so only the raw graph's extent grows.
#
# Two-step osmium pipeline (streaming, ~30 MB RAM — the whole pbf is never loaded
# in process), then osmnx's native XML reader builds the graph exactly as the old
# graph_from_point("drive") path did:
#   osmium extract --bbox    → bbox.osm.pbf   (geographic clip of the snapshot)
#   osmium tags-filter w/highway=<drive set> → drive.osm  (drivable highways + nodes)
#   ox.graph_from_xml(drive.osm)             → simplified MultiDiGraph

_to_utm = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
_to_wgs = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

if os.path.exists(CENSUS_ZONES_FILE):
    _cz = json.load(open(CENSUS_ZONES_FILE))
    _core_utm = _shp_transform(_to_utm.transform, _Poly(_cz["core_polygon"]))
    _bbox_utm = _core_utm.buffer(BOUNDARY_BBOX_MARGIN_M)
    print(f"Extracting drive network from pbf: core polygon + {BOUNDARY_BBOX_MARGIN_M}m "
          f"(bbox ≈ {_bbox_utm.area / 1e6:.0f} km²) …")
else:
    # Fallback when census_zones.json is absent: a circle around CENTRE.
    _cx, _cy = _to_utm.transform(CENTRE[1], CENTRE[0])   # always_xy → (lon, lat)
    _bbox_utm = _Pt(_cx, _cy).buffer(5000 + BOUNDARY_BBOX_MARGIN_M)
    print(f"No {CENSUS_ZONES_FILE} — falling back to a "
          f"{5000 + BOUNDARY_BBOX_MARGIN_M}m circle around CENTRE "
          f"(run build_census_zones.py first).")

_minx, _miny, _maxx, _maxy = _shp_transform(_to_wgs.transform, _bbox_utm).bounds


def _docker_image_exists(image):
    return subprocess.run(["docker", "image", "inspect", image],
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL).returncode == 0


if not _docker_image_exists(OSMIUM_IMAGE):
    print(f"  Building {OSMIUM_IMAGE} image (one-off) from {OSMIUM_DOCKERFILE} …")
    subprocess.run(["docker", "build", "-f", OSMIUM_DOCKERFILE,
                    "-t", OSMIUM_IMAGE, os.path.dirname(OSMIUM_DOCKERFILE)],
                   check=True)

_pbf_dir  = os.path.abspath(os.path.dirname(PBF_PATH))
_pbf_name = os.path.basename(PBF_PATH)
_out_abs  = os.path.abspath(OUT_DIR)
_uidgid   = f"{os.getuid()}:{os.getgid()}"
_docker_pre = ["docker", "run", "--rm", "--user", _uidgid,
               "-v", f"{_pbf_dir}:/pbf:ro", "-v", f"{_out_abs}:/out",
               OSMIUM_IMAGE]

print(f"  osmium extract --bbox {_minx:.5f},{_miny:.5f},{_maxx:.5f},{_maxy:.5f} …")
subprocess.run(_docker_pre + [
    "extract", "--bbox", f"{_minx},{_miny},{_maxx},{_maxy}",
    "--strategy", "smart", "--overwrite",
    "-o", f"/out/{os.path.basename(_BBOX_PBF)}", f"/pbf/{_pbf_name}"], check=True)

print(f"  osmium tags-filter w/highway=<drive set> …")
subprocess.run(_docker_pre + [
    "tags-filter", "--overwrite",
    "-o", f"/out/{os.path.basename(_DRIVE_OSM)}",
    f"/out/{os.path.basename(_BBOX_PBF)}", f"w/highway={DRIVE_HIGHWAYS}"], check=True)

G = ox.graph_from_xml(_DRIVE_OSM, simplify=True, retain_all=False)
# graph_from_xml does not populate the `street_count` node attribute that the
# Overpass download path adds. consolidate_intersections relies on it to identify
# true intersections; without it the graph badly under-merges (≈1416 vs ≈1002
# core nodes). Add it here so consolidation matches the previous pipeline.
import networkx as nx
nx.set_node_attributes(G, values=ox.stats.count_streets_per_node(G),
                       name="street_count")
print(f"  {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ── 2. Save raw graph ──────────────────────────────────────────────────────────

NODE_COLS = ["geometry", "street_count"]
LINK_COLS = ["geometry", "name", "highway", "oneway", "length",
             "maxspeed", "lanes", "bridge", "tunnel", "osmid"]

nodes_gdf, links_gdf = ox.graph_to_gdfs(G, nodes=True, edges=True)
nodes_out = nodes_gdf[[c for c in NODE_COLS if c in nodes_gdf.columns]].copy()
links_out = links_gdf[[c for c in LINK_COLS if c in links_gdf.columns]].copy()

for col in links_out.columns:
    if col != "geometry":
        links_out[col] = links_out[col].apply(lambda v: v[0] if isinstance(v, list) else v)

nodes_out.to_file(f"{OUT_DIR}/newtownards_nodes.geojson", driver="GeoJSON")
links_out.to_file(f"{OUT_DIR}/newtownards_links.geojson", driver="GeoJSON")
ox.save_graphml(G, f"{OUT_DIR}/newtownards_network.graphml")

# ── 3. Clip to core polygon, then consolidate ─────────────────────────────────
# The full download covers RADIUS_M so boundary nodes have their external
# neighbours in the graph (needed for build_demographics.py OSRM boundary
# detection via newtownards_network.graphml).  The consolidated routing graph
# used by build_paths.py should contain only core-polygon nodes — buffer nodes
# outside the core are not road nodes in the model; external zones connect to
# boundary nodes via OSRM-derived links instead.

if os.path.exists(CENSUS_ZONES_FILE):
    from shapely.geometry import Point as _Pt, Polygon as _Poly
    _core_poly_wgs = _Poly(_cz["core_polygon"])
    _core_nodes = [n for n, d in G.nodes(data=True)
                   if _core_poly_wgs.contains(_Pt(d["x"], d["y"]))]
    G_core = G.subgraph(_core_nodes).copy()
    print(f"\nClipped to core polygon: {G.number_of_nodes()} → {G_core.number_of_nodes()} nodes")
else:
    G_core = G
    print("\nNo census_zones.json — skipping core-polygon clip (run build_census_zones.py first)")

print(f"Consolidating junctions (tolerance={CONSOLIDATION_TOLERANCE_M}m) …")
G_proj = ox.project_graph(G_core)  # auto-selects UTM zone (EPSG:32630 for NI)
G_cons = ox.consolidate_intersections(
    G_proj,
    tolerance=CONSOLIDATION_TOLERANCE_M,
    rebuild_graph=True,
    dead_ends=False,
)
print(f"  {G_proj.number_of_nodes()} -> {G_cons.number_of_nodes()} nodes  "
      f"({G_proj.number_of_nodes() - G_cons.number_of_nodes()} merged)")
print(f"  {G_proj.number_of_edges()} -> {G_cons.number_of_edges()} edges")

# Relabel consolidated nodes to use stable OSM-based IDs so that re-running
# build_network.py doesn't break downstream node references.
# Rule: consolidated node → min(osmid_original); single node → int(osmid_original).
# ox.consolidate_intersections stores osmid_original as a string like
# '[286949408, 13098508098]' (not a real list) on the in-memory graph; parse
# with ast.literal_eval before the GraphML write converts it properly.
import ast, networkx as nx
_relabel = {}
for _n, _d in G_cons.nodes(data=True):
    _orig = _d.get("osmid_original")
    if isinstance(_orig, list):
        _relabel[_n] = min(int(x) for x in _orig)
    elif isinstance(_orig, str):
        _parsed = ast.literal_eval(_orig)
        if isinstance(_parsed, list):
            _relabel[_n] = min(int(x) for x in _parsed)
        else:
            _relabel[_n] = int(_parsed)
    elif _orig is not None:
        _relabel[_n] = int(_orig)
    else:
        _relabel[_n] = _n
G_cons = nx.relabel_nodes(G_cons, _relabel)
print(f"  Relabeled to OSM IDs: {G_cons.number_of_nodes()} nodes, "
      f"{G_cons.number_of_edges()} edges")

ox.save_graphml(G_cons, f"{OUT_DIR}/newtownards_consolidated.graphml")

# ── 4. Summary ────────────────────────────────────────────────────────────────

total_km = links_out["length"].sum() / 1000
print(f"\nNetwork summary:")
print(f"  Total length:  {total_km:.1f} km")
print(f"  Mean link:     {links_out['length'].mean():.0f} m")
print(f"  Median link:   {links_out['length'].median():.0f} m")
print(f"  Centre:        {CENTRE[0]:.6f}, {CENTRE[1]:.6f}  "
      f"bbox margin: {BOUNDARY_BBOX_MARGIN_M}m (from pbf)")

print(f"\nLink breakdown by highway type:")
for htype, count in links_out["highway"].value_counts().items():
    print(f"  {htype:<30} {count:>5}")
