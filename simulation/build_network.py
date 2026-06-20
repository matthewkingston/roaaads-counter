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

import osmnx as ox
import geopandas as gpd

# Newtownards town centre
CENTRE      = (54.5933779, -5.6960935)
CORE_RADIUS = 3000   # metres — defines core area (SDZs intersecting this → DZ-level)
RADIUS_M    = 5000   # metres — OSMnx download bbox; larger than CORE_RADIUS so all
                     # boundary nodes and their external neighbours are in the graph
CONSOLIDATION_TOLERANCE_M = 15   # merge nodes within this distance
OUT_DIR   = "simulation"

# ── 1. Download ────────────────────────────────────────────────────────────────

print(f"Downloading drive network within {RADIUS_M}m of Newtownards town centre …")
G = ox.graph_from_point(CENTRE, dist=RADIUS_M, network_type="drive", retain_all=False)
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

# ── 3. Consolidate junctions ───────────────────────────────────────────────────
# Project to UTM first so the tolerance is in real metres, not degrees

print(f"\nConsolidating junctions (tolerance={CONSOLIDATION_TOLERANCE_M}m) …")
G_proj = ox.project_graph(G)  # auto-selects UTM zone (EPSG:32630 for NI)
G_cons = ox.consolidate_intersections(
    G_proj,
    tolerance=CONSOLIDATION_TOLERANCE_M,
    rebuild_graph=True,
    dead_ends=False,
)
print(f"  {G_proj.number_of_nodes()} -> {G_cons.number_of_nodes()} nodes  "
      f"({G_proj.number_of_nodes() - G_cons.number_of_nodes()} merged)")
print(f"  {G_proj.number_of_edges()} -> {G_cons.number_of_edges()} edges")

ox.save_graphml(G_cons, f"{OUT_DIR}/newtownards_consolidated.graphml")

# ── 4. Summary ────────────────────────────────────────────────────────────────

total_km = links_out["length"].sum() / 1000
print(f"\nNetwork summary:")
print(f"  Total length:  {total_km:.1f} km")
print(f"  Mean link:     {links_out['length'].mean():.0f} m")
print(f"  Median link:   {links_out['length'].median():.0f} m")
print(f"  Centre:        {CENTRE[0]:.6f}, {CENTRE[1]:.6f}  radius: {RADIUS_M}m")

print(f"\nLink breakdown by highway type:")
for htype, count in links_out["highway"].value_counts().items():
    print(f"  {htype:<30} {count:>5}")
