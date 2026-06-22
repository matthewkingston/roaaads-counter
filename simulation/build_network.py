"""
Download the Newtownards road network from OSM and produce Level-1 segments:
one directed link per road between adjacent intersection nodes, with complex
junctions (roundabouts, staggered crossings) consolidated into single nodes.

Outputs (in the same directory):
  newtownards_nodes.geojson        — raw intersection nodes
  newtownards_links.geojson        — raw directed links
  newtownards_network.graphml      — raw graph (geographic CRS)
  newtownards_consolidated.graphml — junction-consolidated graph (projected CRS)
"""

import json, os, sys
import osmnx as ox
import geopandas as gpd
import pyproj
from shapely.geometry import Polygon as _Poly, Point as _Pt
from shapely.ops import transform as _shp_transform
from pyrosm import OSM

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

# ── 1. Read road network from the local NI pbf ──────────────────────────────────
# Source the graph from the same .osm.pbf OSRM is built from (one OSM snapshot →
# boundary/internal node IDs match OSRM route node IDs). The read is bounded by the
# core polygon buffered by BOUNDARY_BBOX_MARGIN_M, which supersedes the old Overpass
# `dist` circle margin (the buffer only needs to reach boundary nodes' external
# neighbours; 5 km is generous). The consolidated routing graph is still clipped to
# the core polygon in step 3, so only the raw graph's extent grows.

_to_utm = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
_to_wgs = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

if os.path.exists(CENSUS_ZONES_FILE):
    _cz = json.load(open(CENSUS_ZONES_FILE))
    _core_utm = _shp_transform(_to_utm.transform, _Poly(_cz["core_polygon"]))
    _bbox_utm = _core_utm.buffer(BOUNDARY_BBOX_MARGIN_M)
    print(f"Reading drive network from pbf: core polygon + {BOUNDARY_BBOX_MARGIN_M}m "
          f"(bbox ≈ {_bbox_utm.area / 1e6:.0f} km²) …")
else:
    # Fallback when census_zones.json is absent: a circle around CENTRE.
    _cx, _cy = _to_utm.transform(CENTRE[1], CENTRE[0])   # always_xy → (lon, lat)
    _bbox_utm = _Pt(_cx, _cy).buffer(5000 + BOUNDARY_BBOX_MARGIN_M)
    print(f"No {CENSUS_ZONES_FILE} — falling back to a "
          f"{5000 + BOUNDARY_BBOX_MARGIN_M}m circle around CENTRE "
          f"(run build_census_zones.py first).")

_bbox_wgs = _shp_transform(_to_wgs.transform, _bbox_utm)

osm = OSM(PBF_PATH, bounding_box=_bbox_wgs)
nodes, edges = osm.get_network(network_type="driving", nodes=True,
                               extra_attributes=["bridge", "tunnel", "lanes"])
G = osm.to_graph(nodes, edges, graph_type="networkx", retain_all=False)
if "crs" not in G.graph:
    G.graph["crs"] = "epsg:4326"
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
