"""
Fetch NISRA Census 2021 Data Zone population data, filter to the Newtownards
study area, and generate an interactive map showing:
  - Road network (colour-coded by type)
  - DZ boundaries as a population-density choropleth
  - Road nodes with road-length-weighted population estimates

Outputs:
  model/newtownards_demographics.geojson  — DZ polygons + population for study area
  model/newtownards_map.html              — updated map (overwrites previous)
"""

import json, sys, urllib.request
import geopandas as gpd
import pandas as pd
import osmnx as ox
import folium
import numpy as np
from shapely.geometry import Point
from scipy.spatial import cKDTree
import pyproj

# ── Config ─────────────────────────────────────────────────────────────────────

CENTRE    = (54.5933779, -5.6960935)
RADIUS_M  = 3000
OUT_DIR   = "simulation"

POPULATION_API = (
    "https://ws-data.nisra.gov.uk/public/api.restful/"
    "PxStat.Data.Cube_API.ReadDataset/MYE01T011/CSV/1.0/en/"
)
DZ_BOUNDARY_FILE    = "simulation/dz2021/DZ2021.geojson"
GRAPH_PATH          = "simulation/newtownards_consolidated.graphml"
WORKPLACE_DATA_FILE = "data/census-2021-apwp001.xlsx"
POPULATION_CACHE    = "data/cache_nisra_population.csv"
POI_CACHE           = "data/cache_osm_pois.geojson"

HIGHWAY_STYLE = {
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

# ── External zone weights for boundary nodes ───────────────────────────────────
# Geographic positions (name, lat, lon) are hardcoded here.
# Pop, workplace, and damping are read from tuner_config.json — single source of
# truth shared with the tuner's L2 regularisation anchors.
# Node 180 (local access) has no external centroid and no tuner_config entry;
# its pop/wp are hardcoded as it is too small to tune.

_EXT_GEO = {
    #  node: (name,                    lat,      lon)
     47: ("Donaghadee",           54.6408, -5.5328),
     65: ("Comber",               54.5503, -5.7419),
     92: ("Lower Ards Peninsula", 54.4892, -5.5283),
     97: ("Belfast",              54.5973, -5.9301),
     98: ("Bangor",               54.6536, -5.6697),
     99: ("Holywood",             54.6322, -5.8325),
    119: ("Belfast",              54.5973, -5.9301),
    180: (None,                   None,    None   ),
    617: ("Comber",               54.5503, -5.7419),
    618: ("Comber",               54.5503, -5.7419),
    620: ("Comber",               54.5503, -5.7419),
    731: ("Bangor",               54.6536, -5.6697),
    748: ("Millisle",             54.6015, -5.5031),
    749: ("Millisle",             54.6015, -5.5031),
}

with open("simulation/tuner_config.json") as _f:
    _tuner_cfg_ext = json.load(_f)

_node_cfg = {}  # node_id → (ref_pop, ref_wp, damping)
for _city_cfg in _tuner_cfg_ext["cities"].values():
    for _nid in _city_cfg["nodes"]:
        _node_cfg[_nid] = (
            _city_cfg["ref_pop"],
            _city_cfg["ref_wp"],
            _city_cfg["dampings"][str(_nid)],
        )

EXTERNAL_ZONES = {}
for _nid, (name, lat, lon) in _EXT_GEO.items():
    if _nid == 180:
        EXTERNAL_ZONES[_nid] = (name, lat, lon, 50, 0, 1.0)
    else:
        ref_pop, ref_wp, damp = _node_cfg[_nid]
        EXTERNAL_ZONES[_nid] = (name, lat, lon, ref_pop, ref_wp, damp)

# ── Fast path: --zones-only ────────────────────────────────────────────────────
# Patches only boundary node entries in node_weights.json. Skips all internal
# population processing and map building. Use after editing ref_pop/ref_wp/dampings
# in tuner_config.json (not when lat/lon has changed).

if "--zones-only" in sys.argv:
    import os as _os
    weights_path = f"{OUT_DIR}/node_weights.json"
    with open(weights_path) as f:
        w = json.load(f)
    print("Updating external zone weights …")
    for node_id, (name, lat, lon, pop, workplace, damping) in EXTERNAL_ZONES.items():
        w["node_population"][str(node_id)]      = pop * damping
        w["node_business_demand"][str(node_id)] = workplace * damping
        zone = name or "local access"
        print(f"  Node {node_id:4d}  {zone:<22}  pop={pop * damping:>8.0f}  workplace={workplace * damping:>8.0f}  damping={damping}")
    with open(weights_path, "w") as f:
        json.dump(w, f)
    print(f"Saved {len(EXTERNAL_ZONES)} boundary nodes → {weights_path}")
    print("Next: python3 simulation/build_assignment.py")
    sys.exit(0)

# ── Fast path: --map-only ──────────────────────────────────────────────────────
# Rebuilds only the map HTML, reusing node_weights.json and
# newtownards_demographics.geojson written by a prior full run.
# Use after retuning (build_assignment.py already updates the flows).

if "--map-only" in sys.argv:
    import os as _os
    print("--map-only: loading saved outputs …")
    transformer_to_utm = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
    transformer_to_wgs = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)
    centre_utm_x, centre_utm_y = transformer_to_utm.transform(CENTRE[1], CENTRE[0])
    print("Loading graph …")
    G_cons = ox.load_graphml(GRAPH_PATH)
    node_ids = list(G_cons.nodes())
    node_coords_utm = [(G_cons.nodes[n]["x"], G_cons.nodes[n]["y"]) for n in node_ids]
    print("Loading node weights …")
    with open(f"{OUT_DIR}/node_weights.json") as _f:
        _w = json.load(_f)
    node_population      = {int(k): v for k, v in _w["node_population"].items()}
    node_business_demand = {int(k): v for k, v in _w["node_business_demand"].items()}
    print("Loading DZ boundaries …")
    dz_final = gpd.read_file(f"{OUT_DIR}/newtownards_demographics.geojson")
    print(f"  {len(dz_final)} Data Zones · {len(node_ids)} nodes")

else:
    # ── 1. Download population data ────────────────────────────────────────────

    import os as _os
    if _os.path.exists(POPULATION_CACHE):
        print(f"Loading NISRA population from cache ({POPULATION_CACHE}) …")
        pop_df = pd.read_csv(POPULATION_CACHE)
    else:
        print("Fetching NISRA mid-2021 DZ population …")
        req = urllib.request.Request(POPULATION_API, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            pop_csv = r.read().decode("utf-8-sig")
        from io import StringIO
        pop_df = pd.read_csv(StringIO(pop_csv))
        pop_df.to_csv(POPULATION_CACHE, index=False)
        print(f"  Cached → {POPULATION_CACHE}")
    # Filter to 2021 only; exclude the aggregate NI-total row (code N92000002)
    pop_df = pop_df[
        (pop_df["TLIST(A1)"] == 2021) &
        (pop_df["DZ2021"].str.startswith("N20"))
    ][["DZ2021", "Data Zones", "VALUE"]].rename(
        columns={"DZ2021": "DZ2021_cd", "Data Zones": "DZ2021_nm", "VALUE": "population"}
    )
    print(f"  {len(pop_df)} Data Zones, NI total pop: {pop_df['population'].sum():,}")

    # ── 2. Load boundaries and join population ─────────────────────────────────

    print("Loading DZ boundaries …")
    dz = gpd.read_file(DZ_BOUNDARY_FILE)  # EPSG:4326
    dz = dz.merge(pop_df[["DZ2021_cd", "population"]], on="DZ2021_cd", how="left")

    # ── 3. Clip to study circle & estimate population ──────────────────────────

    print("Clipping Data Zones to study circle …")

    transformer_to_utm = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
    transformer_to_wgs = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

    # Work in UTM throughout for accurate area calculations
    dz_utm = dz.to_crs("EPSG:32630")
    centre_utm_x, centre_utm_y = transformer_to_utm.transform(CENTRE[1], CENTRE[0])
    study_circle = Point(centre_utm_x, centre_utm_y).buffer(RADIUS_M)

    dz_intersect = dz_utm[dz_utm.geometry.intersects(study_circle)].copy()
    dz_intersect["area_original_m2"] = dz_intersect.geometry.area

    # Hard clip to circle
    dz_intersect = dz_intersect.copy()
    dz_intersect["geometry"] = dz_intersect.geometry.intersection(study_circle)
    dz_intersect["area_clipped_m2"] = dz_intersect.geometry.area
    dz_intersect["area_pct"] = dz_intersect["area_clipped_m2"] / dz_intersect["area_original_m2"]

    # Population scaled by area fraction (constant density assumption)
    dz_intersect["pop_estimated"] = (
        dz_intersect["population"] * dz_intersect["area_pct"]
    ).round().astype("Int64")

    print(f"  {len(dz_intersect)} DZs intersect circle before node filter")

    # ── 4. Filter to DZs with road nodes inside ───────────────────────────────

    print("Filtering to DZs containing road nodes …")
    G_cons = ox.load_graphml(GRAPH_PATH)

    # Consolidated graph nodes already have x/y in EPSG:32630 (UTM)
    node_ids = list(G_cons.nodes())
    node_coords_utm = [(G_cons.nodes[n]["x"], G_cons.nodes[n]["y"]) for n in node_ids]

    # Spatial join: which clipped DZs contain at least one junction node?
    node_gdf = gpd.GeoDataFrame(
        {"node_id": node_ids},
        geometry=[Point(x, y) for x, y in node_coords_utm],
        crs="EPSG:32630",
    )
    joined = gpd.sjoin(node_gdf, dz_intersect[["DZ2021_cd", "geometry"]], how="inner", predicate="within")
    dzs_with_nodes = set(joined["DZ2021_cd"])

    dz_final = dz_intersect[dz_intersect["DZ2021_cd"].isin(dzs_with_nodes)].copy()
    n_dropped = len(dz_intersect) - len(dz_final)
    print(f"  Dropped {n_dropped} DZs with no road nodes inside clipped area")
    print(f"  Kept {len(dz_final)} DZs, estimated pop: {dz_final['pop_estimated'].sum():,}")

    # Save clipped DZs (convert back to WGS84 for GeoJSON)
    dz_final.to_crs("EPSG:4326").to_file(f"{OUT_DIR}/newtownards_demographics.geojson", driver="GeoJSON")

    # ── Road-length-weighted population per node ───────────────────────────────
    # For each node within a DZ, sum outgoing edge lengths as its "road weight".
    # Each node receives: DZ_pop_estimated × (node_road_length / total_DZ_road_length)

    print("Computing road-length-weighted node populations …")

    # Sum of outgoing edge lengths for every node in the consolidated graph
    node_road_length = {
        n: sum(float(d.get("length", 0)) for _, _, d in G_cons.out_edges(n, data=True))
        for n in node_ids
    }

    # joined maps each node to its DZ; build dz_code → [node_id, ...] lookup
    dz_to_nodes = joined.groupby("DZ2021_cd")["node_id"].apply(list).to_dict()

    node_population = {}   # node_id → float population share
    pop_lookup = dz_final.set_index("DZ2021_cd")["pop_estimated"]

    for dz_code, nodes_in_dz in dz_to_nodes.items():
        dz_pop = float(pop_lookup.get(dz_code, 0) or 0)
        total_len = sum(node_road_length[n] for n in nodes_in_dz)
        for n in nodes_in_dz:
            if total_len > 0:
                node_population[n] = dz_pop * node_road_length[n] / total_len
            else:
                node_population[n] = dz_pop / len(nodes_in_dz)

    assigned = sum(1 for v in node_population.values() if v > 0)
    print(f"  {assigned} nodes assigned population (of {len(node_ids)} total)")

    # ── OSM POI download + workplace population allocation ─────────────────────
    # OSM POIs (amenity / shop / office) are snapped to their nearest consolidated
    # node and used as spatial weights within each DZ.
    # NISRA Census 2021 workplace population (APWP001) provides the DZ control total;
    # each node's business demand = DZ_workplace_pop × (node_POI_count / DZ_POI_count).
    # Fallback for DZs with no POIs: equal distribution among nodes.

    print("Loading NISRA Census 2021 workplace population (APWP001) …")
    wp_df = pd.read_excel(WORKPLACE_DATA_FILE, sheet_name="DZ", header=5)
    wp_df = wp_df[["Geography Code", "Workplace population"]].rename(
        columns={"Geography Code": "DZ2021_cd", "Workplace population": "workplace_pop"}
    )
    wp_df = wp_df[wp_df["DZ2021_cd"].astype(str).str.startswith("N20")].copy()
    wp_df["workplace_pop"] = pd.to_numeric(wp_df["workplace_pop"], errors="coerce").fillna(0)
    print(f"  {len(wp_df)} DZs · NI total workplace pop: {int(wp_df['workplace_pop'].sum()):,}")

    # Join to dz_final so we only keep study-area DZs
    dz_final = dz_final.merge(wp_df[["DZ2021_cd", "workplace_pop"]], on="DZ2021_cd", how="left")
    dz_final["workplace_pop"] = dz_final["workplace_pop"].fillna(0)
    wp_lookup = dz_final.set_index("DZ2021_cd")["workplace_pop"]
    print(f"  Study area workplace pop: {int(dz_final['workplace_pop'].sum()):,}")


    EXCLUDE_AMENITY = {
        "parking", "parking_space", "parking_entrance", "fuel",
        "atm", "vending_machine", "post_box", "waste_basket",
        "bench", "bicycle_parking", "recycling", "toilets",
        "shelter", "telephone",
    }

    if _os.path.exists(POI_CACHE):
        print(f"Loading OSM POIs from cache ({POI_CACHE}) …")
        pois_raw = gpd.read_file(POI_CACHE)
    else:
        print("Downloading OSM POIs (amenity / shop / office) …")
        pois_raw = ox.features_from_point(CENTRE, tags={"amenity": True, "shop": True, "office": True}, dist=RADIUS_M)
        _save_cols = [c for c in ["amenity", "shop", "office"] if c in pois_raw.columns]
        pois_raw[_save_cols + ["geometry"]].to_crs("EPSG:4326").to_file(POI_CACHE, driver="GeoJSON")
        print(f"  Cached → {POI_CACHE}")

    # Filter out low-trip-generating amenity types
    if "amenity" in pois_raw.columns:
        mask = pois_raw["amenity"].isna() | ~pois_raw["amenity"].isin(EXCLUDE_AMENITY)
        pois_raw = pois_raw[mask]

    # Normalise all geometries to points (polygon/linestring features → centroid)
    pois_utm = pois_raw.to_crs("EPSG:32630").copy()
    pois_utm["geometry"] = pois_utm.geometry.centroid
    print(f"  {len(pois_utm)} POIs after filtering")

    # Snap each POI to its nearest consolidated node
    kdtree = cKDTree(node_coords_utm)
    poi_coords = [(geom.x, geom.y) for geom in pois_utm.geometry]
    _, nearest_indices = kdtree.query(poi_coords)
    pois_utm["nearest_node"] = [node_ids[int(i)] for i in nearest_indices]

    # POI count per node
    node_poi_count = pois_utm.groupby("nearest_node").size().to_dict()

    # Allocate workplace population to nodes within each DZ by POI weight
    node_business_demand = {}
    for dz_code, nodes_in_dz in dz_to_nodes.items():
        dz_wp = float(wp_lookup.get(dz_code, 0) or 0)
        poi_weights = {n: node_poi_count.get(n, 0) for n in nodes_in_dz}
        total_pois = sum(poi_weights.values())
        for n in nodes_in_dz:
            if total_pois > 0:
                node_business_demand[n] = dz_wp * poi_weights[n] / total_pois
            else:
                node_business_demand[n] = dz_wp / len(nodes_in_dz)

    active_biz_nodes = sum(1 for v in node_business_demand.values() if v > 0)
    total_biz = sum(node_business_demand.values())
    print(f"  {active_biz_nodes} nodes with business demand · total {total_biz:.0f} workplace pop attributed")

    # node_id → effective UTM centroid used for gravity-model distances:
    #   internal nodes  → their own network coordinates
    #   boundary nodes  → external destination centroid (or own coords for node 180)
    node_effective_utm = {n: (float(x), float(y)) for n, (x, y) in zip(node_ids, node_coords_utm)}

    print("Assigning external zone weights to boundary nodes …")
    for node_id, (name, lat, lon, pop, workplace, damping) in EXTERNAL_ZONES.items():
        node_population[node_id]      = pop * damping
        node_business_demand[node_id] = workplace * damping
        if lat is not None:
            cx, cy = transformer_to_utm.transform(lon, lat)
            node_effective_utm[node_id] = (cx, cy)
        zone = name or "local access"
        print(f"  Node {node_id:4d}  {zone:<22}  pop={pop * damping:>8.0f}  workplace={workplace * damping:>8.0f}  damping={damping}")

    # ── Serialise node weights for assignment script ───────────────────────────
    weights_path = f"{OUT_DIR}/node_weights.json"
    with open(weights_path, "w") as f:
        json.dump({
            "node_population":      {str(k): v for k, v in node_population.items()},
            "node_business_demand": {str(k): v for k, v in node_business_demand.items()},
            "node_effective_utm":   {str(k): list(v) for k, v in node_effective_utm.items()},
            "boundary_node_ids":    list(EXTERNAL_ZONES.keys()),
        }, f)
    print(f"Saved node weights → {weights_path}")

# ── 5. Build map ───────────────────────────────────────────────────────────────

print("Building map …")
G_raw = ox.load_graphml(f"{OUT_DIR}/newtownards_network.graphml")

m = folium.Map(location=list(CENTRE), zoom_start=14, tiles="CartoDB positron")

# 5a. Road network edges (from raw geographic graph)
from collections import defaultdict
by_type = defaultdict(list)
for u, v, data in G_raw.edges(data=True):
    htype = data.get("highway", "unclassified")
    if isinstance(htype, list): htype = htype[0]
    by_type[htype].append((u, v, data))

type_order = [
    "living_street", "unclassified", "residential",
    "tertiary_link", "tertiary", "secondary",
    "primary_link", "primary", "trunk_link", "trunk", "motorway",
]
all_types = type_order + [t for t in by_type if t not in type_order]

ROAD_TYPE_LABELS = {
    "trunk":         "Roads · trunk",
    "trunk_link":    "Roads · trunk (links)",
    "primary":       "Roads · primary",
    "primary_link":  "Roads · primary (links)",
    "secondary":     "Roads · secondary",
    "tertiary":      "Roads · tertiary",
    "tertiary_link": "Roads · tertiary (links)",
    "residential":   "Roads · residential",
    "unclassified":  "Roads · unclassified",
    "living_street": "Roads · living street",
}
for htype in all_types:
    edges = by_type.get(htype)
    if not edges:
        continue
    style = HIGHWAY_STYLE.get(htype, {"color": "#aaaaaa", "weight": 1})
    label = ROAD_TYPE_LABELS.get(htype, f"Roads · {htype}")
    fg = folium.FeatureGroup(name=label, show=False)
    for u, v, data in edges:
        geom = data.get("geometry")
        if geom and hasattr(geom, "coords"):
            coords = [(lat, lon) for lon, lat in geom.coords]
        else:
            ud, vd = G_raw.nodes[u], G_raw.nodes[v]
            coords = [(ud["y"], ud["x"]), (vd["y"], vd["x"])]
        name = data.get("name", "")
        if isinstance(name, list):
            name = name[0]
        length = float(data.get("length", 0))
        folium.PolyLine(
            coords, color=style["color"], weight=style["weight"],
            opacity=0.7,
            tooltip=f"{name or '(unnamed)'} [{htype}] · {length:.0f}m",
        ).add_to(fg)
    fg.add_to(m)

# 5b. Classify nodes as boundary (in/out flow) vs interior.
#
# Set established by auto-detection followed by manual inspection.
# To update: add/remove node IDs from BOUNDARY_NODE_IDS.

import math

BOUNDARY_NODE_IDS = {
    47, 65, 92, 97, 98, 99,
    119, 180, 617, 618, 620,
    731, 748, 749,
}

boundary_nodes_map = {}   # node_id → (wgs_lat, wgs_lon, dist, degree)
interior_nodes_map = {}

for node_id, (nx_utm, ny_utm) in zip(node_ids, node_coords_utm):
    dist = math.sqrt((nx_utm - centre_utm_x)**2 + (ny_utm - centre_utm_y)**2)
    deg  = G_cons.degree(node_id)
    nlon, nlat = transformer_to_wgs.transform(nx_utm, ny_utm)
    if node_id in BOUNDARY_NODE_IDS:
        boundary_nodes_map[node_id] = (nlat, nlon, dist, deg)
    else:
        interior_nodes_map[node_id] = (nlat, nlon, dist, deg)

print(f"  Node classification: {len(boundary_nodes_map)} boundary, {len(interior_nodes_map)} interior")

# Boundary nodes — ON by default, orange diamonds
boundary_fg = folium.FeatureGroup(name=f"In/Out flow nodes — auto-detected ({len(boundary_nodes_map)})", show=True)
for node_id, (nlat, nlon, dist, deg) in boundary_nodes_map.items():
    node_pop = node_population.get(node_id, 0)
    node_biz = node_business_demand.get(node_id, 0)
    ez = EXTERNAL_ZONES.get(node_id)
    zone_name = ez[0] or "local access" if ez else "—"
    damping = ez[5] if ez else 1.0
    damp_str = f" ×{damping}" if damping < 1.0 else ""
    folium.RegularPolygonMarker(
        location=[nlat, nlon],
        number_of_sides=4, radius=6, rotation=45,
        color="#e05c00", fill=True, fill_color="#ff7c20", fill_opacity=0.9, weight=1.5,
        tooltip=(
            f"<b>Node {node_id}</b> [boundary → {zone_name}{damp_str}]<br>"
            f"degree={deg} · {dist:.0f}m from centre<br>"
            f"ext. pop: {node_pop:,.0f}<br>"
            f"ext. workplace: {node_biz:,.0f}"
        ),
    ).add_to(boundary_fg)
boundary_fg.add_to(m)

# Interior nodes — OFF by default, small blue dots
interior_fg = folium.FeatureGroup(name=f"Interior nodes ({len(interior_nodes_map)})", show=False)
for node_id, (nlat, nlon, dist, deg) in interior_nodes_map.items():
    node_pop = node_population.get(node_id, 0)
    node_biz = node_business_demand.get(node_id, 0)
    folium.CircleMarker(
        location=[nlat, nlon], radius=3,
        color="#1a73e8", fill=True, fill_color="#1a73e8", fill_opacity=0.8, weight=1,
        tooltip=(
            f"Node {node_id} · degree={deg} · {dist:.0f}m from centre<br>"
            f"est. pop: {node_pop:.1f}<br>"
            f"workplace pop: {node_biz:.1f}"
        ),
    ).add_to(interior_fg)
interior_fg.add_to(m)

# 5d. DZ choropleth (estimated population in clipped area)
# Convert clipped DZs back to WGS84 for Leaflet rendering
dz_plot = dz_final.to_crs("EPSG:4326")
pop_vals = dz_plot["pop_estimated"].dropna().astype(float)
pop_min, pop_max = pop_vals.min(), pop_vals.max()

def pop_color(pop):
    if pd.isna(pop): return "#eeeeee"
    t = (float(pop) - pop_min) / max(pop_max - pop_min, 1)
    r = int(255 * (1 - t * 0.8))
    g = int(255 * (1 - t * 0.8))
    b = int(255 * (1 - t * 0.3))
    return f"#{r:02x}{g:02x}{b:02x}"

dz_fg = folium.FeatureGroup(name=f"Data Zones — pop estimate ({len(dz_plot)})", show=True)
for _, row in dz_plot.iterrows():
    geojson_str = json.dumps(row.geometry.__geo_interface__)
    pop_est = int(row["pop_estimated"]) if pd.notna(row["pop_estimated"]) else None
    pop_full = int(row["population"]) if pd.notna(row["population"]) else None
    area_pct = row["area_pct"]
    clipped = area_pct < 0.99
    folium.GeoJson(
        geojson_str,
        style_function=lambda f, pop=row["pop_estimated"]: {
            "fillColor":   pop_color(pop),
            "color":       "#444444",
            "weight":      1.5,
            "fillOpacity": 0.45,
        },
        tooltip=folium.Tooltip(
            f"<b>{row['DZ2021_nm']}</b><br>"
            f"Code: {row['DZ2021_cd']}<br>"
            f"Est. pop (clipped): {pop_est if pop_est else 'N/A'}<br>"
            + (f"Census pop (full DZ): {pop_full} · {area_pct*100:.0f}% of DZ in study area<br>" if clipped else
               f"Census pop: {pop_full}<br>")
            + f"Area in study: {row['area_clipped_m2']/10000:.1f} ha"
        ),
    ).add_to(dz_fg)
dz_fg.add_to(m)

# 5e. Business demand nodes — proportional purple circles, OFF by default
max_biz = max(node_business_demand.values(), default=1)
biz_fg = folium.FeatureGroup(name="Workplace population nodes", show=False)
all_nodes_map = {**boundary_nodes_map, **interior_nodes_map}
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    biz = node_business_demand.get(node_id, 0)
    if biz <= 0:
        continue
    radius = max(3, 14 * (biz / max_biz) ** 0.5)
    node_pop = node_population.get(node_id, 0)
    folium.CircleMarker(
        location=[nlat, nlon], radius=radius,
        color="#7b2d8b", fill=True, fill_color="#b05ec0", fill_opacity=0.65, weight=1,
        tooltip=(
            f"<b>Node {node_id}</b><br>"
            f"workplace pop: {biz:.1f}<br>"
            f"est. pop: {node_pop:.1f}"
        ),
    ).add_to(biz_fg)
biz_fg.add_to(m)

# ── Optional flow layer (loaded from newtownards_flows.json if it exists) ────────
import os as _os, math as _math
_flows_path = f"{OUT_DIR}/newtownards_flows.json"
if _os.path.exists(_flows_path):
    import pyproj as _pyproj
    _tr_flow = _pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

    with open(_flows_path) as _f:
        _flows_data = json.load(_f)
    _link_flow = {tuple(int(x) for x in k.split(",")): v
                  for k, v in _flows_data["flows"].items()}

    # Percentile-anchored log colour scale (same logic as build_assignment.py)
    _sorted = sorted(v for v in _link_flow.values() if v > 0)
    if _sorted:
        _p10 = _sorted[max(0, int(len(_sorted) * 0.10))]
        _p90 = _sorted[min(len(_sorted) - 1, int(len(_sorted) * 0.90))]
        _lmin = _math.log10(max(_p10, 1))
        _lmax = _math.log10(max(_p90, 1))
    else:
        _lmin, _lmax = 0.0, 1.0

    def _flow_color(flow):
        if flow <= 0:
            return "#cccccc"
        t = (_math.log10(max(flow, 1)) - _lmin) / max(_lmax - _lmin, 1e-6)
        t = max(0.0, min(1.0, t))
        if t < 0.33:
            r, g, b = 0, int(180 * (t / 0.33)), int(200 * (1 - t / 0.33))
        elif t < 0.66:
            s = (t - 0.33) / 0.33
            r, g, b = int(220 * s), 180, 0
        else:
            s = (t - 0.66) / 0.34
            r, g, b = 220 + int(35 * s), int(180 * (1 - s)), 0
        return f"#{r:02x}{g:02x}{b:02x}"

    def _flow_weight(flow):
        if flow <= 0:
            return 1
        t = (_math.log10(max(flow, 1)) - _lmin) / max(_lmax - _lmin, 1e-6)
        return 1 + 7 * max(0.0, min(1.0, t))

    flow_fg = folium.FeatureGroup(name="Road flows — est. AADT", show=True)
    for u, v, data in G_cons.edges(data=True):
        flow = _link_flow.get((u, v), 0) + _link_flow.get((v, u), 0)
        geom = data.get("geometry")
        if geom and hasattr(geom, "coords"):
            _coords = [_tr_flow.transform(x, y)[::-1] for x, y in geom.coords]
        else:
            ud, vd = G_cons.nodes[u], G_cons.nodes[v]
            lon_u, lat_u = _tr_flow.transform(float(ud["x"]), float(ud["y"]))
            lon_v, lat_v = _tr_flow.transform(float(vd["x"]), float(vd["y"]))
            _coords = [(lat_u, lon_u), (lat_v, lon_v)]
        name = data.get("name", "")
        length = float(data.get("length", 0))
        folium.PolyLine(
            _coords,
            color=_flow_color(flow),
            weight=_flow_weight(flow),
            opacity=0.85,
            tooltip=f"{name or 'link'}<br>est. AADT: {flow:,.0f}<br>length: {length:.0f}m",
        ).add_to(flow_fg)
    flow_fg.add_to(m)
    print(f"Added flow layer from {_flows_path} ({len(_link_flow)} links)")
else:
    print(f"No flow data found at {_flows_path} — skipping flow layer")

folium.LayerControl(collapsed=False).add_to(m)

out_path = f"{OUT_DIR}/newtownards_map.html"
m.save(out_path)
total_node_pop = sum(node_population.values())
print(f"\nSaved: {out_path}")
print(f"  {len(dz_plot)} Data Zones · est. pop range {pop_min:.0f}–{pop_max:.0f}")
print(f"  {len(node_population)} nodes with population assigned · total {total_node_pop:.0f}")
