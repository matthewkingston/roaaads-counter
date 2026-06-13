"""
Render the Newtownards road network as an interactive Leaflet map.

Uses the consolidated graph (roundabouts / complex junctions merged).
A toggle layer shows the raw (pre-consolidation) nodes for comparison.

Output: model/newtownards_map.html  (self-contained, open in any browser)
"""

import osmnx as ox
import folium

RAW_GRAPH_PATH  = "simulation/newtownards_network.graphml"
CONS_GRAPH_PATH = "simulation/newtownards_consolidated.graphml"
OUT_PATH        = "simulation/newtownards_map.html"

HIGHWAY_STYLE = {
    "motorway":      {"color": "#e8684a", "weight": 5},
    "trunk":         {"color": "#f5a623", "weight": 4},
    "trunk_link":    {"color": "#f5a623", "weight": 2},
    "primary":       {"color": "#f5d623", "weight": 3},
    "primary_link":  {"color": "#f5d623", "weight": 2},
    "secondary":     {"color": "#a8d08d", "weight": 2},
    "tertiary":      {"color": "#7bafd4", "weight": 2},
    "tertiary_link": {"color": "#7bafd4", "weight": 1},
    "residential":   {"color": "#cccccc", "weight": 1},
    "unclassified":  {"color": "#bbbbbb", "weight": 1},
    "living_street": {"color": "#dddddd", "weight": 1},
}
DEFAULT_STYLE = {"color": "#aaaaaa", "weight": 1}


def add_edges(G, m, use_xy=False):
    """Draw all edges grouped by highway type as FeatureGroups."""
    from collections import defaultdict
    by_type = defaultdict(list)
    for u, v, data in G.edges(data=True):
        htype = data.get("highway", "unclassified")
        if isinstance(htype, list):
            htype = htype[0]
        by_type[htype].append((u, v, data))

    type_order = [
        "living_street", "unclassified", "residential",
        "tertiary_link", "tertiary", "secondary",
        "primary_link", "primary",
        "trunk_link", "trunk", "motorway",
    ]
    all_types = type_order + [t for t in by_type if t not in type_order]

    for htype in all_types:
        edges = by_type.get(htype)
        if not edges:
            continue
        style = HIGHWAY_STYLE.get(htype, DEFAULT_STYLE)
        fg = folium.FeatureGroup(name=f"{htype} ({len(edges)})", show=True)

        for u, v, data in edges:
            u_data = G.nodes[u]
            v_data = G.nodes[v]

            geom = data.get("geometry")
            if geom and hasattr(geom, "coords"):
                if use_xy:
                    # projected CRS — coords are (easting, northing), need lat/lon
                    import pyproj
                    transformer = pyproj.Transformer.from_crs(
                        G.graph.get("crs", "EPSG:32630"), "EPSG:4326", always_xy=True
                    )
                    coords = [
                        transformer.transform(x, y)[::-1]  # (lat, lon)
                        for x, y in geom.coords
                    ]
                else:
                    coords = [(lat, lon) for lon, lat in geom.coords]
            else:
                if use_xy:
                    coords = [(u_data["y"], u_data["x"]), (v_data["y"], v_data["x"])]
                else:
                    coords = [(u_data["y"], u_data["x"]), (v_data["y"], v_data["x"])]

            name   = data.get("name", "")
            speed  = data.get("maxspeed", "?")
            lanes  = data.get("lanes", "?")
            length = float(data.get("length", 0))
            oneway = data.get("oneway", False)
            junc   = data.get("junction", "")

            tooltip = (
                f"<b>{name or htype}</b><br>"
                f"Type: {htype}{' · ' + junc if junc else ''}<br>"
                f"Length: {length:.0f} m<br>"
                f"Speed: {speed} · Lanes: {lanes}<br>"
                f"One-way: {oneway}"
            )

            folium.PolyLine(
                coords,
                color=style["color"],
                weight=style["weight"],
                opacity=0.85,
                tooltip=tooltip,
            ).add_to(fg)

        fg.add_to(m)


def add_nodes(G, m, label, color, show, use_xy=False):
    """Draw intersection nodes as a toggleable dot layer."""
    import pyproj
    transformer = None
    if use_xy:
        transformer = pyproj.Transformer.from_crs(
            G.graph.get("crs", "EPSG:32630"), "EPSG:4326", always_xy=True
        )

    intersections = [(n, d) for n, d in G.nodes(data=True) if G.degree(n) > 2]
    fg = folium.FeatureGroup(name=f"{label} ({len(intersections)})", show=show)

    for node_id, data in intersections:
        if use_xy:
            lon, lat = transformer.transform(data["x"], data["y"])
        else:
            lat, lon = data["y"], data["x"]

        folium.CircleMarker(
            location=[lat, lon],
            radius=4,
            color=color,
            fill=True,
            fill_color=color,
            fill_opacity=0.7,
            weight=1,
            tooltip=f"Node {node_id} · degree {G.degree(node_id)}",
        ).add_to(fg)

    fg.add_to(m)


# ── Build map ──────────────────────────────────────────────────────────────────

print("Loading graphs …")
G_raw  = ox.load_graphml(RAW_GRAPH_PATH)
G_cons = ox.load_graphml(CONS_GRAPH_PATH)

centre_lat = sum(d["y"] for _, d in G_raw.nodes(data=True)) / G_raw.number_of_nodes()
centre_lon = sum(d["x"] for _, d in G_raw.nodes(data=True)) / G_raw.number_of_nodes()

print(f"Building map …")
m = folium.Map(
    location=[centre_lat, centre_lon],
    zoom_start=14,
    tiles="CartoDB positron",
)

# Road links (from raw graph — geometry is in geographic CRS)
add_edges(G_raw, m, use_xy=False)

# Consolidated nodes (white) — default ON
add_nodes(G_cons, m, label="Consolidated intersections", color="#1a73e8", show=True,  use_xy=True)

# Raw nodes (grey) — default OFF, for comparison
add_nodes(G_raw,  m, label="Raw nodes (pre-consolidation)", color="#999999", show=False, use_xy=False)

folium.LayerControl(collapsed=False).add_to(m)

m.save(OUT_PATH)
print(f"Saved: {OUT_PATH}")
print(f"  Raw:          {G_raw.number_of_nodes()} nodes, {G_raw.number_of_edges()} edges")
print(f"  Consolidated: {G_cons.number_of_nodes()} nodes, {G_cons.number_of_edges()} edges")
