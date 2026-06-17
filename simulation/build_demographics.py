"""
Fetch NISRA Census 2021 Data Zone population data, filter to the Newtownards
study area, and generate an interactive map showing:
  - Road network (colour-coded by type)
  - DZ boundaries as a population-density choropleth
  - Road nodes with OSM-building-snapped population estimates (road-length fallback per DZ)

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
from shapely.geometry import Point, LineString
from shapely.strtree import STRtree
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
BUILDING_CACHE      = "data/cache_osm_buildings.geojson"
PARKING_CACHE       = "data/cache_osm_parking.geojson"

EXCLUDE_AMENITY = {
    "parking", "parking_space", "parking_entrance",
    "vending_machine", "post_box", "waste_basket",
    "bench", "bicycle_parking", "recycling",
    "shelter", "telephone", "grit_bin",
}

# Per-tag trip-generation weights relative to baseline (café/small shop = 1.0).
# Parking layer already handles large retail anchors; weights here add signal
# from institutional employers and high-turnover stops without parking polygons.
POI_WEIGHTS = {
    # amenity tag → weight
    "hospital":        5.0,
    "school":          3.0,
    "college":         3.0,
    "university":      3.0,
    "cinema":          3.0,
    "theatre":         3.0,
    "fuel":            2.0,
    "fast_food":       1.5,
    "place_of_worship": 0.5,
    "atm":             0.5,
    "toilets":         0.25,
    # shop tag → weight
    "supermarket":     1.5,
    # any office tag → 2.0 (applied inline; not listed here as the value covers all subtypes)
}

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
  10000: ("Dundonald",           54.5790, -5.8450),
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
    node_parking_equiv   = {int(k): v for k, v in _w.get("node_parking_equiv", {}).items()}
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

    print("Computing road-length-weighted node populations (fallback) …")

    # Sum of outgoing edge lengths for every node in the consolidated graph (fallback)
    node_road_length = {
        n: sum(float(d.get("length", 0)) for _, _, d in G_cons.out_edges(n, data=True))
        for n in node_ids
    }

    # joined maps each node to its DZ; build dz_code → [node_id, ...] lookup
    dz_to_nodes = joined.groupby("DZ2021_cd")["node_id"].apply(list).to_dict()

    # ── Build edge STRtree for rigorous "snap building to nearest road" ──────────
    # STRtree.nearest() finds the true closest point on edge geometry, not midpoint.

    _edge_keys = []
    _edge_geom_list = []
    for _eu, _ev, _edata in G_cons.edges(data=True):
        _egeom = _edata.get("geometry")
        if _egeom is None:
            _ux, _uy = float(G_cons.nodes[_eu]["x"]), float(G_cons.nodes[_eu]["y"])
            _vx, _vy = float(G_cons.nodes[_ev]["x"]), float(G_cons.nodes[_ev]["y"])
            _egeom = LineString([(_ux, _uy), (_vx, _vy)])
        _edge_keys.append((_eu, _ev))
        _edge_geom_list.append(_egeom)
    _edge_strtree = STRtree(_edge_geom_list)
    print(f"  Built edge STRtree ({len(_edge_keys)} edges)")

    # ── Download and cache OSM buildings ─────────────────────────────────────────

    if _os.path.exists(BUILDING_CACHE):
        print(f"Loading OSM buildings from cache ({BUILDING_CACHE}) …")
        _bld_raw = gpd.read_file(BUILDING_CACHE)
    else:
        print("Downloading OSM buildings …")
        _bld1 = ox.features_from_point(
            CENTRE,
            tags={"building": ["house", "detached", "semidetached_house", "terrace",
                               "apartments", "bungalow", "residential", "flat",
                               "dormitory", "cottage", "yes"]},
            dist=RADIUS_M,
        )
        _bld2 = ox.features_from_point(
            CENTRE, tags={"addr:housenumber": True}, dist=RADIUS_M
        )
        _bld_raw = pd.concat([_bld1, _bld2])
        _bld_raw = _bld_raw[~_bld_raw.index.duplicated(keep="first")]
        _bld_raw = _bld_raw[_bld_raw.geometry.notna()].copy()
        _bld_raw = _bld_raw.to_crs("EPSG:32630")
        _bld_raw["geometry"] = _bld_raw.geometry.centroid
        _bld_raw = _bld_raw[_bld_raw.geometry.geom_type == "Point"]
        _bld_raw = _bld_raw[["geometry"]].to_crs("EPSG:4326")
        _bld_raw.to_file(BUILDING_CACHE, driver="GeoJSON")
        print(f"  Cached → {BUILDING_CACHE}")

    # Reset to clean integer index so sjoin indices are valid positional offsets
    buildings_utm = _bld_raw.to_crs("EPSG:32630").reset_index(drop=True)
    print(f"  {len(buildings_utm)} buildings in study area")

    # Spatial join: which buildings fall inside each DZ?
    _bld_dz = gpd.sjoin(
        buildings_utm,
        dz_final[["DZ2021_cd", "geometry"]],
        how="inner", predicate="within",
    )
    dz_bld_idx = _bld_dz.groupby("DZ2021_cd").apply(lambda df: list(df.index)).to_dict()

    # ── Building-snapped or road-length-weighted population per node ─────────────

    FALLBACK_THRESHOLD = 3    # min buildings per DZ to use OSM method
    POP_PER_BLD_WARN   = 10.0 # pop/building above this suggests missing data

    node_population = {}
    pop_lookup = dz_final.set_index("DZ2021_cd")["pop_estimated"]
    _quality_rows = []

    for dz_code, nodes_in_dz in dz_to_nodes.items():
        dz_pop = float(pop_lookup.get(dz_code, 0) or 0)
        bld_indices = dz_bld_idx.get(dz_code, [])
        n_bld = len(bld_indices)
        pop_per_bld = dz_pop / n_bld if n_bld > 0 else float("inf")

        too_few   = n_bld < FALLBACK_THRESHOLD
        too_dense = pop_per_bld >= POP_PER_BLD_WARN
        use_osm   = not too_few and not too_dense

        weights = {}
        if use_osm:
            weights = {n: 0.0 for n in nodes_in_dz}
            per_bld = 1.0 / n_bld
            nodes_in_dz_set = set(nodes_in_dz)
            dz_node_coords = [(G_cons.nodes[n]["x"], G_cons.nodes[n]["y"]) for n in nodes_in_dz]
            _dz_kd = cKDTree(dz_node_coords)
            for _, row in buildings_utm.iloc[bld_indices].iterrows():
                pt = row.geometry
                nearest_edge_idx = _edge_strtree.nearest(pt)
                eu, ev = _edge_keys[nearest_edge_idx]
                line = _edge_geom_list[nearest_edge_idx]
                t = line.project(pt, normalized=True)  # 0=u end, 1=v end
                in_dz = nodes_in_dz_set.intersection({eu, ev})
                if len(in_dz) == 2:
                    weights[eu] += (1 - t) * per_bld
                    weights[ev] += t * per_bld
                elif len(in_dz) == 1:
                    weights[next(iter(in_dz))] += per_bld
                else:
                    # Neither edge endpoint is in this DZ — assign to nearest DZ node
                    _, _nn = _dz_kd.query((pt.x, pt.y))
                    weights[nodes_in_dz[_nn]] += per_bld
            method = f"OSM buildings (n={n_bld}, pop/bld={pop_per_bld:.1f})"

        if not use_osm:
            total_len = sum(node_road_length[n] for n in nodes_in_dz)
            for n in nodes_in_dz:
                weights[n] = (node_road_length[n] / total_len) if total_len > 0 else 1.0 / len(nodes_in_dz)
            if too_few:
                method = f"road-length fallback (only {n_bld} buildings)"
            else:
                method = f"road-length fallback (pop/bld={pop_per_bld:.1f} >= {POP_PER_BLD_WARN})"

        for n, w in weights.items():
            node_population[n] = node_population.get(n, 0.0) + dz_pop * w

        _quality_rows.append((dz_code, n_bld, pop_per_bld, method))

    # Quality report
    _n_osm = sum(1 for *_, m in _quality_rows if "fallback" not in m)
    _n_fb  = len(_quality_rows) - _n_osm
    print(f"\nDZ population distribution ({_n_osm} OSM buildings / {_n_fb} road-length fallback):")
    print(f"  {'DZ code':<15} | {'n_bld':>5} | {'pop/bld':>7} | method")
    for _dz_code, _n_bld, _ppb, _method in sorted(_quality_rows):
        _ppb_str = f"{_ppb:.1f}" if _ppb != float("inf") else "    inf"
        print(f"  {_dz_code:<15} | {_n_bld:>5} | {_ppb_str:>7} | {_method}")
    print(f"  Totals: {len(_quality_rows)} DZs  |  OSM: {_n_osm}  |  Fallback: {_n_fb}")

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
    dz_final["workplace_pop"] = dz_final["workplace_pop"].fillna(0) * dz_final["area_pct"]
    wp_lookup = dz_final.set_index("DZ2021_cd")["workplace_pop"]
    print(f"  Study area workplace pop: {int(dz_final['workplace_pop'].sum()):,}")


    if _os.path.exists(POI_CACHE):
        print(f"Loading OSM POIs from cache ({POI_CACHE}) …")
        pois_raw = gpd.read_file(POI_CACHE)
    else:
        print("Downloading OSM POIs (amenity / shop / office) …")
        pois_raw = ox.features_from_point(CENTRE, tags={"amenity": True, "shop": True, "office": True}, dist=RADIUS_M)
        _save_cols = [c for c in ["amenity", "shop", "office", "name"] if c in pois_raw.columns]
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

    # Snap each POI centroid to nearest road edge; split weight (1-t)/t between endpoints
    node_poi_weight = {}
    for _, _poi_row in pois_utm.iterrows():
        _pt = _poi_row.geometry
        _eidx = _edge_strtree.nearest(_pt)
        _eu, _ev = _edge_keys[_eidx]
        _t = _edge_geom_list[_eidx].project(_pt, normalized=True)
        # Determine trip-generation weight: check amenity and shop tags first,
        # then default to 2.0 for any office, else 1.0.
        _amenity = _poi_row.get("amenity") if "amenity" in pois_utm.columns else None
        _shop    = _poi_row.get("shop")    if "shop"    in pois_utm.columns else None
        _office  = _poi_row.get("office")  if "office"  in pois_utm.columns else None
        if pd.notna(_amenity) and _amenity in POI_WEIGHTS:
            _w = POI_WEIGHTS[_amenity]
        elif pd.notna(_shop) and _shop in POI_WEIGHTS:
            _w = POI_WEIGHTS[_shop]
        elif pd.notna(_office):
            _w = 2.0
        else:
            _w = 1.0
        node_poi_weight[_eu] = node_poi_weight.get(_eu, 0.0) + _w * (1.0 - _t)
        node_poi_weight[_ev] = node_poi_weight.get(_ev, 0.0) + _w * _t

    # Allocate workplace population to nodes within each DZ by edge-snapped POI weight
    node_business_demand = {}
    for dz_code, nodes_in_dz in dz_to_nodes.items():
        dz_wp = float(wp_lookup.get(dz_code, 0) or 0)
        poi_weights_dz = {n: node_poi_weight.get(n, 0.0) for n in nodes_in_dz}
        total_weight = sum(poi_weights_dz.values())
        for n in nodes_in_dz:
            if total_weight > 0:
                node_business_demand[n] = dz_wp * poi_weights_dz[n] / total_weight
            else:
                node_business_demand[n] = dz_wp / len(nodes_in_dz)

    active_biz_nodes = sum(1 for v in node_business_demand.values() if v > 0)
    total_biz = sum(node_business_demand.values())
    print(f"  {active_biz_nodes} nodes with business demand · total {total_biz:.0f} workplace pop attributed")

    # ── Car park area → additional business demand ────────────────────────────────
    # OSM parking polygons proxy vehicle trip attraction (visitor/customer demand).
    # Public: area/25 ≈ equivalent persons.  Private (access=private): area/50 (half weight).

    if _os.path.exists(PARKING_CACHE):
        print(f"Loading OSM parking from cache ({PARKING_CACHE}) …")
        _park_raw = gpd.read_file(PARKING_CACHE)
    else:
        print("Downloading OSM parking polygons …")
        _park_raw = ox.features_from_point(
            CENTRE,
            tags={"amenity": "parking", "landuse": "parking"},
            dist=RADIUS_M,
        )
        _park_raw = _park_raw[_park_raw.geometry.notna()].copy()
        _save_park_cols = [c for c in ["amenity", "landuse", "access", "name"] if c in _park_raw.columns]
        _park_raw[_save_park_cols + ["geometry"]].to_crs("EPSG:4326").to_file(PARKING_CACHE, driver="GeoJSON")
        print(f"  Cached → {PARKING_CACHE}")

    # Polygon features only — excludes point-tagged parking_entrance, parking_space nodes
    _park_utm = _park_raw.to_crs("EPSG:32630").copy()
    _park_utm = _park_utm[_park_utm.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    _park_utm["area_m2"] = _park_utm.geometry.area
    _park_utm["centroid_geom"] = _park_utm.geometry.centroid
    _park_utm["is_private"] = (
        _park_utm["access"].isin(["private"]) if "access" in _park_utm.columns
        else pd.Series(False, index=_park_utm.index)
    )

    # Snap each parking centroid to nearest road edge; split equiv between endpoints by t
    PARKING_SCALE_PUBLIC  = 25.0
    PARKING_SCALE_PRIVATE = 50.0
    node_parking_equiv = {}
    for _, _row in _park_utm.iterrows():
        _pt = _row["centroid_geom"]
        _eidx = _edge_strtree.nearest(_pt)
        _eu, _ev = _edge_keys[_eidx]
        _t = _edge_geom_list[_eidx].project(_pt, normalized=True)
        _scale = PARKING_SCALE_PRIVATE if _row["is_private"] else PARKING_SCALE_PUBLIC
        _equiv = _row["area_m2"] / _scale
        _eu_share, _ev_share = (1.0 - _t) * _equiv, _t * _equiv
        node_business_demand[_eu] = node_business_demand.get(_eu, 0.0) + _eu_share
        node_business_demand[_ev] = node_business_demand.get(_ev, 0.0) + _ev_share
        node_parking_equiv[_eu]   = node_parking_equiv.get(_eu, 0.0) + _eu_share
        node_parking_equiv[_ev]   = node_parking_equiv.get(_ev, 0.0) + _ev_share

    _n_priv  = int(_park_utm["is_private"].sum())
    _n_pub   = len(_park_utm) - _n_priv
    _eq_pub  = _park_utm[~_park_utm["is_private"]]["area_m2"].sum() / PARKING_SCALE_PUBLIC
    _eq_priv = _park_utm[_park_utm["is_private"]]["area_m2"].sum()  / PARKING_SCALE_PRIVATE
    print(f"  {len(_park_utm)} parking polygons "
          f"({_n_pub} public +{_eq_pub:.0f} eq-persons, "
          f"{_n_priv} private +{_eq_priv:.0f} eq-persons)")

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
            "node_parking_equiv":   {str(k): v for k, v in node_parking_equiv.items()},
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
    731, 748, 749, 10000,
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

# Interior nodes — ON by default, small blue dots
interior_fg = folium.FeatureGroup(name=f"Interior nodes ({len(interior_nodes_map)})", show=True)
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

dz_fg = folium.FeatureGroup(name=f"Data Zones — pop estimate ({len(dz_plot)})", show=False)
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
# Scale based on internal nodes only so boundary giants (Belfast) don't collapse the range.
max_biz_internal = max((node_business_demand.get(n, 0) for n in interior_nodes_map), default=1)
biz_fg = folium.FeatureGroup(name="Business demand nodes", show=False)
all_nodes_map = {**boundary_nodes_map, **interior_nodes_map}
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    biz = node_business_demand.get(node_id, 0)
    if biz <= 0:
        continue
    radius = max(3, 14 * (biz / max_biz_internal) ** 0.5)
    node_pop  = node_population.get(node_id, 0)
    park_eq   = node_parking_equiv.get(node_id, 0)
    wp_demand = biz - park_eq
    folium.CircleMarker(
        location=[nlat, nlon], radius=radius,
        color="#7b2d8b", fill=True, fill_color="#b05ec0", fill_opacity=0.65, weight=1,
        tooltip=(
            f"<b>Node {node_id}</b><br>"
            f"workplace population: {wp_demand:.1f}<br>"
            f"parking eq.: {park_eq:.1f}<br>"
            f"business demand: {biz:.1f}<br>"
            f"est. pop: {node_pop:.1f}"
        ),
    ).add_to(biz_fg)
biz_fg.add_to(m)

# 5f. Parking polygons — realism check layer (red = private, blue = public/untagged)
_park_wgs = None
if "_park_utm" in dir():  # full run: variable already in scope
    _park_wgs = _park_utm.to_crs("EPSG:4326")
elif _os.path.exists(PARKING_CACHE):  # --map-only: reload and reprocess from cache
    _p = gpd.read_file(PARKING_CACHE).to_crs("EPSG:32630")
    _p = _p[_p.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    _p["area_m2"] = _p.geometry.area
    _p["is_private"] = _p["access"].isin(["private"]) if "access" in _p.columns else False
    _park_wgs = _p.to_crs("EPSG:4326")

if _park_wgs is not None:
    park_fg = folium.FeatureGroup(name=f"Car parks — OSM ({len(_park_wgs)})", show=False)
    for _, _row in _park_wgs.iterrows():
        _access_col = "access" in _park_wgs.columns
        _access_val = _row["access"] if _access_col else None
        _access_label = str(_access_val) if (_access_col and pd.notna(_access_val)) else "public/untagged"
        _area = _row["area_m2"]
        _name_val = _row["name"] if "name" in _park_wgs.columns and pd.notna(_row.get("name")) else ""
        _is_priv = bool(_row.get("is_private", False))
        _color = "#cc4444" if _is_priv else "#4488cc"
        _tip = (
            f"<b>{_name_val or 'Car park'}</b><br>"
            f"access: {_access_label}<br>"
            f"area: {_area:,.0f} m²"
        )
        folium.GeoJson(
            _row.geometry.__geo_interface__,
            style_function=lambda f, c=_color: {
                "fillColor": c, "color": c, "weight": 1.5, "fillOpacity": 0.4,
            },
            tooltip=folium.Tooltip(_tip),
        ).add_to(park_fg)
    park_fg.add_to(m)
    print(f"Added parking layer ({len(_park_wgs)} polygons)")

# 5g. POI layer — amenity/shop/office features used for workplace allocation
_poi_wgs = None
if "pois_utm" in dir():  # full run: already filtered + centroid geometry in UTM
    _poi_wgs = pois_utm.to_crs("EPSG:4326")
elif _os.path.exists(POI_CACHE):  # --map-only: reload and re-filter from cache
    _pp = gpd.read_file(POI_CACHE)
    if "amenity" in _pp.columns:
        _pp = _pp[_pp["amenity"].isna() | ~_pp["amenity"].isin(EXCLUDE_AMENITY)]
    _pp = _pp[_pp.geometry.notna()].copy()
    _pp["geometry"] = _pp.geometry.centroid
    _poi_wgs = _pp

if _poi_wgs is not None:
    _POI_COLOURS = {"amenity": "#e67e22", "shop": "#27ae60", "office": "#2980b9"}
    poi_fg = folium.FeatureGroup(name=f"POIs — workplace allocation ({len(_poi_wgs)})", show=False)
    for _, _row in _poi_wgs.iterrows():
        _kind, _val, _color = None, None, "#888888"
        for _col in ("amenity", "shop", "office"):
            if _col in _poi_wgs.columns and pd.notna(_row.get(_col)):
                _kind, _val, _color = _col, _row[_col], _POI_COLOURS[_col]
                break
        _name = ""
        if "name" in _poi_wgs.columns and pd.notna(_row.get("name")):
            _name = str(_row["name"])
        _type_str = f"{_kind}: {_val}" if _kind else "unknown"
        _tip = f"<b>{_name}</b><br>{_type_str}" if _name else f"<b>{_type_str}</b>"
        folium.CircleMarker(
            location=[_row.geometry.y, _row.geometry.x],
            radius=5,
            color=_color, fill=True, fill_color=_color, fill_opacity=0.85, weight=1,
            tooltip=folium.Tooltip(_tip),
        ).add_to(poi_fg)
    poi_fg.add_to(m)
    print(f"Added POI layer ({len(_poi_wgs)} POIs: orange=amenity, green=shop, blue=office)")

# ── Optional flow layers (loaded from newtownards_flows.json if it exists) ────
import os as _os, math as _math
_flows_path = f"{OUT_DIR}/newtownards_flows.json"
if _os.path.exists(_flows_path):
    import pyproj as _pyproj
    _tr_flow = _pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

    with open(_flows_path) as _f:
        _flows_data = json.load(_f)

    def _parse_flows(key):
        return {tuple(int(x) for x in k.split(",")): v
                for k, v in _flows_data.get(key, {}).items()}

    _link_flow     = _parse_flows("flows")
    _link_flow_res = _parse_flows("flows_res")
    _link_flow_biz = _parse_flows("flows_biz")
    _has_components = bool(_link_flow_res)

    # ── colour helpers ────────────────────────────────────────────────────────
    def _log_scale(flow_dict):
        """Return (lmin, lmax) from P10/P90 of a flow dict."""
        vals = sorted(v for v in flow_dict.values() if v > 0)
        if not vals:
            return 0.0, 1.0
        p10 = vals[max(0, int(len(vals) * 0.10))]
        p90 = vals[min(len(vals) - 1, int(len(vals) * 0.90))]
        return _math.log10(max(p10, 1)), _math.log10(max(p90, 1))

    def _t(flow, lmin, lmax):
        if flow <= 0:
            return 0.0
        return max(0.0, min(1.0, (_math.log10(max(flow, 1)) - lmin) / max(lmax - lmin, 1e-6)))

    def _weight(t):
        return 1 + 7 * t

    def _color_combined(t):
        """Blue → yellow → red."""
        if t < 0.33:
            r, g, b = 0, int(180 * (t / 0.33)), int(200 * (1 - t / 0.33))
        elif t < 0.66:
            s = (t - 0.33) / 0.33
            r, g, b = int(220 * s), 180, 0
        else:
            s = (t - 0.66) / 0.34
            r, g, b = 220 + int(35 * s), int(180 * (1 - s)), 0
        return f"#{r:02x}{g:02x}{b:02x}"

    def _color_res(t):
        """Light teal → dark forest green."""
        r = int(30  + 80  * (1 - t))
        g = int(140 + 80  * (1 - t * 0.5))
        b = int(80  + 80  * (1 - t))
        return f"#{r:02x}{g:02x}{b:02x}"

    def _color_biz(t):
        """Amber → dark orange-red."""
        r = min(255, int(200 + 55 * t))
        g = int(160 * (1 - t * 0.85))
        b = int(20  * (1 - t))
        return f"#{r:02x}{g:02x}{b:02x}"

    # ── helper to build one FeatureGroup ─────────────────────────────────────
    def _add_flow_fg(flow_dict, name, color_fn, tooltip_fn, show):
        lmin, lmax = _log_scale(flow_dict)
        fg = folium.FeatureGroup(name=name, show=show)
        for u, v, data in G_cons.edges(data=True):
            flow = flow_dict.get((u, v), 0) + flow_dict.get((v, u), 0)
            geom = data.get("geometry")
            if geom and hasattr(geom, "coords"):
                coords = [_tr_flow.transform(x, y)[::-1] for x, y in geom.coords]
            else:
                ud, vd = G_cons.nodes[u], G_cons.nodes[v]
                lon_u, lat_u = _tr_flow.transform(float(ud["x"]), float(ud["y"]))
                lon_v, lat_v = _tr_flow.transform(float(vd["x"]), float(vd["y"]))
                coords = [(lat_u, lon_u), (lat_v, lon_v)]
            ti = _t(flow, lmin, lmax)
            folium.PolyLine(
                coords,
                color=color_fn(ti) if flow > 0 else "#cccccc",
                weight=_weight(ti),
                opacity=0.85,
                tooltip=tooltip_fn(data.get("name", ""), flow,
                                   float(data.get("length", 0))),
            ).add_to(fg)
        fg.add_to(m)

    # ── Combined layer (always shown) ─────────────────────────────────────────
    def _tt_combined(name, flow, length):
        tip = f"{name or 'link'}<br>est. AADT: {flow:,.0f}"
        if _has_components:
            r = _link_flow_res.get((_u, _v), 0) + _link_flow_res.get((_v, _u), 0)
            b = _link_flow_biz.get((_u, _v), 0) + _link_flow_biz.get((_v, _u), 0)
            tot = r + b
            if tot > 0:
                tip += f"<br>  residential: {r:,.0f} ({100*r/tot:.0f}%)"
                tip += f"<br>  business: {b:,.0f} ({100*b/tot:.0f}%)"
        tip += f"<br>length: {length:.0f}m"
        return tip

    # The combined tooltip references (u,v) from the outer loop — build it inline
    lmin_c, lmax_c = _log_scale(_link_flow)
    fg_combined = folium.FeatureGroup(name="Road flows — est. AADT", show=True)
    for _u, _v, _data in G_cons.edges(data=True):
        _flow = _link_flow.get((_u, _v), 0) + _link_flow.get((_v, _u), 0)
        _geom = _data.get("geometry")
        if _geom and hasattr(_geom, "coords"):
            _coords = [_tr_flow.transform(x, y)[::-1] for x, y in _geom.coords]
        else:
            _ud, _vd = G_cons.nodes[_u], G_cons.nodes[_v]
            _lon_u, _lat_u = _tr_flow.transform(float(_ud["x"]), float(_ud["y"]))
            _lon_v, _lat_v = _tr_flow.transform(float(_vd["x"]), float(_vd["y"]))
            _coords = [(_lat_u, _lon_u), (_lat_v, _lon_v)]
        _ti = _t(_flow, lmin_c, lmax_c)
        _tip = f"{_data.get('name', '') or 'link'}<br>est. AADT: {_flow:,.0f}"
        if _has_components:
            _r = _link_flow_res.get((_u, _v), 0) + _link_flow_res.get((_v, _u), 0)
            _b = _link_flow_biz.get((_u, _v), 0) + _link_flow_biz.get((_v, _u), 0)
            _tot = _r + _b
            if _tot > 0:
                _tip += f"<br>&nbsp;&nbsp;residential: {_r:,.0f} ({100*_r/_tot:.0f}%)"
                _tip += f"<br>&nbsp;&nbsp;business: {_b:,.0f} ({100*_b/_tot:.0f}%)"
        _tip += f"<br>length: {float(_data.get('length', 0)):.0f}m"
        folium.PolyLine(
            _coords,
            color=_color_combined(_ti) if _flow > 0 else "#cccccc",
            weight=_weight(_ti), opacity=0.85, tooltip=_tip,
        ).add_to(fg_combined)
    fg_combined.add_to(m)

    # ── Component layers (off by default, only added when data available) ─────
    if _has_components:
        _add_flow_fg(
            _link_flow_res, "Road flows — residential (pop→pop)", _color_res,
            lambda name, flow, length: (
                f"{name or 'link'}<br>residential AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        _add_flow_fg(
            _link_flow_biz, "Road flows — business-adjacent (pb+bb)", _color_biz,
            lambda name, flow, length: (
                f"{name or 'link'}<br>business AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        print(f"Added flow layers from {_flows_path} (combined + res + biz, {len(_link_flow)} links)")
    else:
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
