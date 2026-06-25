"""
Fetch NISRA Census 2021 Data Zone population data, filter to the Newtownards
study area, and compute per-node demand weights:
  - residential population (OSM-building-snapped; road-length fallback per DZ)
  - business demand (workplace population × POI proxy + car-park area)
  - school demand (OSM school POIs)
  - boundary-node detection + external-zone weights from census data

Outputs:
  simulation/newtownards_demographics.geojson  — DZ polygons + population
  simulation/node_weights.json                 — per-node demand + boundary IDs

The interactive map is built separately by build_map.py (run it afterwards).

Usage:
  python3 simulation/build_demographics.py              # full node-weight run
  python3 simulation/build_demographics.py --zones-only # patch external zones only
"""

import json, os, sys, urllib.request, argparse
import geopandas as gpd
import pandas as pd
import osmnx as ox
import numpy as np
from shapely.geometry import Point, LineString
from shapely.strtree import STRtree
from scipy.spatial import cKDTree
import pyproj

from demographics_config import (
    CENTRE, OUT_DIR, NETWORK_MARGIN_M, POPULATION_API, DZ_BOUNDARY_FILE,
    GRAPH_PATH, WORKPLACE_DATA_FILE, POPULATION_CACHE, POI_CACHE,
    BUILDING_CACHE, PARKING_CACHE, CENSUS_ZONES_FILE, TUNER_CONFIG_FILE,
    EXCLUDE_AMENITY, POI_WEIGHTS, SCHOOL_ENROLL_FALLBACK,
)

_SCHOOL_TAGS = set(SCHOOL_ENROLL_FALLBACK)

# ── CLI ──────────────────────────────────────────────────────────────────────
_ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument(
    "--zones-only", action="store_true",
    help="Patch only external-zone weights in node_weights.json (fast). "
         "Use after editing ext_biz_scale in tuner_config.json; not when the "
         "core polygon or internal demographics change.")
args = _ap.parse_args()

# ── Config (constants live in demographics_config.py) ──────────────────────────
# RADIUS_M is derived from census data below, once _census_zones is loaded.
RADIUS_M = None  # set from census_zones.json (see below)

# ── External zone weights from census_zones.json ────────────────────────────────
# Population and workplace population for external nodes come from census data
# computed by build_census_zones.py.  No hand-crafted values, no dampings.

_census_zones = None
if os.path.exists(CENSUS_ZONES_FILE):
    with open(CENSUS_ZONES_FILE) as _f:
        _census_zones = json.load(_f)

# Size the OSM download to cover the whole core polygon (matches build_network.py).
if _census_zones is not None:
    RADIUS_M = round(_census_zones["max_core_vertex_dist_m"]) + NETWORK_MARGIN_M

_ext_biz_scale    = 1.0
_ext_school_per_pop = None
if os.path.exists(TUNER_CONFIG_FILE):
    with open(TUNER_CONFIG_FILE) as _f:
        _tc = json.load(_f)
    _ext_biz_scale      = _tc.get("ext_biz_scale", 1.0)
    _ext_school_per_pop = _tc.get("ext_school_per_pop")

if _ext_school_per_pop is None:
    print("ERROR: 'ext_school_per_pop' missing from tuner_config.json.")
    print("  Run: python3 simulation/compute_ext_school_scale.py")
    print("  Then add the printed value to tuner_config.json.")
    sys.exit(1)

# ── Fast path: --zones-only ────────────────────────────────────────────────────
# Patches only external-zone entries in node_weights.json. Skips all internal
# population processing. Use after editing ext_biz_scale in tuner_config.json
# (not when the core polygon or internal demographics have changed).

if args.zones_only:
    if _census_zones is None:
        print(f"ERROR: {CENSUS_ZONES_FILE} not found. Run build_census_zones.py first.")
        sys.exit(1)
    weights_path = f"{OUT_DIR}/node_weights.json"
    with open(weights_path) as f:
        w = json.load(f)
    print(f"Updating external zone weights from census data (ext_biz_scale={_ext_biz_scale:.4f}) …")
    ext_nodes = _census_zones["external_nodes"]
    for ext in ext_nodes:
        nid = ext["id"]
        w["node_population"][str(nid)]      = float(ext["population"])
        w["node_business_demand"][str(nid)] = float(ext["workplace_pop"]) * _ext_biz_scale
        if "node_school_demand" in w:
            w["node_school_demand"][str(nid)] = float(ext["population"]) * _ext_school_per_pop
        print(f"  {nid}  {ext['level']}  "
              f"pop={ext['population']:>8,}  wp={ext['workplace_pop']:>8,}  "
              f"biz={w['node_business_demand'][str(nid)]:>10,.1f}  "
              f"school={w['node_school_demand'][str(nid)]:>8,.1f}")
    with open(weights_path, "w") as f:
        json.dump(w, f)
    print(f"Saved {len(ext_nodes)} external nodes → {weights_path}")
    print("Next: python3 simulation/build_assignment.py")

else:
    # ── Full run: compute per-node demand weights → node_weights.json ──────────
    # (The interactive map is built separately by build_map.py.)
    if _census_zones is None:
        print(f"ERROR: {CENSUS_ZONES_FILE} not found. Run build_census_zones.py first.")
        sys.exit(1)

    # ── 1. Download population data ────────────────────────────────────────────

    if os.path.exists(POPULATION_CACHE):
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

    # ── 3. Select core-polygon DZs & estimate population ───────────────────────
    # Model the whole core polygon (data/census_zones.json), not a 3 km circle.
    # The core polygon is the union of whole core DZs, so each core DZ lies
    # entirely inside it; select by centroid-within and keep full DZ geometry and
    # full population (no area clipping).  area_pct is held at 1.0 so downstream
    # consumers (workplace scaling, choropleth tooltip) keep working unchanged.

    print("Selecting Data Zones within the core polygon …")

    transformer_to_utm = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
    transformer_to_wgs = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

    # Work in UTM throughout for accurate area calculations
    dz_utm = dz.to_crs("EPSG:32630")
    centre_utm_x, centre_utm_y = transformer_to_utm.transform(CENTRE[1], CENTRE[0])

    from shapely.geometry import Polygon as _Polygon
    from shapely.ops import transform as _shp_transform
    core_poly_wgs = _Polygon(_census_zones["core_polygon"])
    core_poly_utm = _shp_transform(transformer_to_utm.transform, core_poly_wgs)

    dz_intersect = dz_utm[dz_utm.geometry.centroid.within(core_poly_utm)].copy()
    dz_intersect["area_original_m2"] = dz_intersect.geometry.area
    # No clipping: full DZ geometry and population are inside the core polygon.
    dz_intersect["area_clipped_m2"] = dz_intersect["area_original_m2"]
    dz_intersect["area_pct"] = 1.0
    dz_intersect["pop_estimated"] = dz_intersect["population"].round().astype("Int64")

    print(f"  {len(dz_intersect)} core DZs before node filter")

    # ── 4. Filter to DZs with road nodes inside ───────────────────────────────

    print("Filtering to DZs containing road nodes …")
    G_cons = ox.load_graphml(GRAPH_PATH)

    # Consolidated graph nodes already have x/y in EPSG:32630 (UTM)
    node_ids = list(G_cons.nodes())
    node_coords_utm = [(G_cons.nodes[n]["x"], G_cons.nodes[n]["y"]) for n in node_ids]

    # Spatial join: which core DZs contain at least one junction node?
    node_gdf = gpd.GeoDataFrame(
        {"node_id": node_ids},
        geometry=[Point(x, y) for x, y in node_coords_utm],
        crs="EPSG:32630",
    )
    joined = gpd.sjoin(node_gdf, dz_intersect[["DZ2021_cd", "geometry"]], how="inner", predicate="within")
    dzs_with_nodes = set(joined["DZ2021_cd"])

    dz_final = dz_intersect[dz_intersect["DZ2021_cd"].isin(dzs_with_nodes)].copy()
    n_dropped = len(dz_intersect) - len(dz_final)
    print(f"  Dropped {n_dropped} core DZs with no road nodes inside")
    print(f"  Kept {len(dz_final)} DZs, estimated pop: {dz_final['pop_estimated'].sum():,}")

    # Save core DZs (convert back to WGS84 for GeoJSON)
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
    print(f"  Built edge STRtree ({len(_edge_keys)} consolidated edges)")

    # ── Ghost edges for absorbed dead-end streets ─────────────────────────────────
    # OSMnx simplify_graph treats bidirectional dead-end termini as degree-2 nodes
    # (in=1, out=1 in the directed graph) and removes them, causing the dead-end edge
    # to vanish from the consolidated graph. Buildings on those streets then snap to
    # the nearest surviving edge — often but not always the correct junction road.
    # Ghost edges restore the raw dead-end geometry to the STRtree so buildings snap
    # to the right location; all their demand is then routed to the surviving junction
    # consolidated node (since all traffic must enter the network there).
    import ast as _ast
    _G_raw = ox.load_graphml(f"{OUT_DIR}/newtownards_network.graphml")

    # Build raw OSM node id → consolidated node id from osmid_original attribute
    _raw_to_cons = {}
    for _cid, _cdata in G_cons.nodes(data=True):
        _oids = _cdata.get("osmid_original", [])
        if isinstance(_oids, int):
            _oids = [_oids]
        elif isinstance(_oids, str):
            _s = _oids.strip()
            _oids = _ast.literal_eval(_s) if _s.startswith("[") else [int(_s)]
        for _oid in _oids:
            _raw_to_cons[int(_oid)] = _cid

    _cons_node_set = set(G_cons.nodes())
    _ghost_junction = {}   # STRtree index → consolidated junction node id

    for _rn in list(_G_raw.nodes()):
        # Dead-end terminus: directed degree 2 (in=1 out=1 on a bidirectional stub)
        if _G_raw.in_degree(_rn) != 1 or _G_raw.out_degree(_rn) != 1:
            continue
        # Must not have survived into the consolidated graph
        if _rn in _cons_node_set or int(_rn) in _raw_to_cons:
            continue
        # Find the single junction neighbour
        _junc_raw = next(iter(_G_raw.successors(_rn)), None)
        if _junc_raw is None:
            continue
        _junc_cons = _raw_to_cons.get(int(_junc_raw))
        if _junc_cons is None:
            continue
        # Get raw edge geometry (WGS84 lon/lat) and project to UTM
        _raw_edata = next(iter(_G_raw[_rn][_junc_raw].values()), {})
        _raw_geom = _raw_edata.get("geometry")
        if _raw_geom is not None:
            _pts = [transformer_to_utm.transform(x, y) for x, y in _raw_geom.coords]
        else:
            _rnd = _G_raw.nodes[_rn]
            _jnd = _G_raw.nodes[_junc_raw]
            _pts = [
                transformer_to_utm.transform(float(_rnd["x"]), float(_rnd["y"])),
                transformer_to_utm.transform(float(_jnd["x"]), float(_jnd["y"])),
            ]
        if len(_pts) < 2:
            continue
        _ghost_geom = LineString(_pts)
        if _ghost_geom.length < 0.1:
            continue
        _ghost_idx = len(_edge_keys)
        _edge_keys.append(None)
        _edge_geom_list.append(_ghost_geom)
        _ghost_junction[_ghost_idx] = _junc_cons

    # Rebuild STRtree to include ghost edges
    _edge_strtree = STRtree(_edge_geom_list)
    print(f"  Added {len(_ghost_junction)} ghost dead-end edges (absorbed termini)")

    # ── Download and cache OSM buildings ─────────────────────────────────────────

    if os.path.exists(BUILDING_CACHE):
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
                if nearest_edge_idx in _ghost_junction:
                    # Absorbed dead-end: all demand → surviving junction consolidated node
                    _junc = _ghost_junction[nearest_edge_idx]
                    if _junc in nodes_in_dz_set:
                        weights[_junc] += per_bld
                    else:
                        _, _nn = _dz_kd.query((pt.x, pt.y))
                        weights[nodes_in_dz[_nn]] += per_bld
                else:
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


    if os.path.exists(POI_CACHE):
        print(f"Loading OSM POIs from cache ({POI_CACHE}) …")
        pois_raw = gpd.read_file(POI_CACHE)
    else:
        print("Downloading OSM POIs (amenity / shop / office) …")
        pois_raw = ox.features_from_point(CENTRE, tags={"amenity": True, "shop": True, "office": True}, dist=RADIUS_M)
        _save_cols = [c for c in ["amenity", "shop", "office", "name", "capacity"] if c in pois_raw.columns]
        pois_raw[_save_cols + ["geometry"]].to_crs("EPSG:4326").to_file(POI_CACHE, driver="GeoJSON")
        print(f"  Cached → {POI_CACHE}")

    # Filter out low-trip-generating amenity types
    if "amenity" in pois_raw.columns:
        mask = pois_raw["amenity"].isna() | ~pois_raw["amenity"].isin(EXCLUDE_AMENITY)
        pois_raw = pois_raw[mask]

    # Normalise all geometries to points (polygon/linestring features → centroid)
    pois_utm = pois_raw.to_crs("EPSG:32630").copy()
    pois_utm["geometry"] = pois_utm.geometry.centroid
    # Keep only POIs inside the core polygon (the download circle overshoots it).
    pois_utm = pois_utm[pois_utm.geometry.within(core_poly_utm)].copy()
    print(f"  {len(pois_utm)} POIs after filtering (within core polygon)")

    # Snap each POI centroid to nearest road edge; split weight (1-t)/t between endpoints
    node_poi_weight = {}
    for _, _poi_row in pois_utm.iterrows():
        _pt = _poi_row.geometry
        _eidx = _edge_strtree.nearest(_pt)
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
        if _eidx in _ghost_junction:
            _junc = _ghost_junction[_eidx]
            node_poi_weight[_junc] = node_poi_weight.get(_junc, 0.0) + _w
        else:
            _eu, _ev = _edge_keys[_eidx]
            _t = _edge_geom_list[_eidx].project(_pt, normalized=True)
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

    if os.path.exists(PARKING_CACHE):
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
    # Keep only car parks whose centroid is inside the core polygon.
    _park_utm = _park_utm[_park_utm["centroid_geom"].within(core_poly_utm)].copy()
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
        _scale = PARKING_SCALE_PRIVATE if _row["is_private"] else PARKING_SCALE_PUBLIC
        _equiv = _row["area_m2"] / _scale
        if _eidx in _ghost_junction:
            _junc = _ghost_junction[_eidx]
            node_business_demand[_junc] = node_business_demand.get(_junc, 0.0) + _equiv
            node_parking_equiv[_junc]   = node_parking_equiv.get(_junc, 0.0) + _equiv
        else:
            _eu, _ev = _edge_keys[_eidx]
            _t = _edge_geom_list[_eidx].project(_pt, normalized=True)
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

    # ── School demand layer ───────────────────────────────────────────────────────
    # Snap school POIs to road edges; accumulate enrollment (pupils) per node.
    # Enrollment = OSM capacity tag where present, else type-based fallback from
    # SCHOOL_ENROLL_FALLBACK.  Units of node_school_demand are pupils, so W_SCHOOL
    # is interpretable as a trip-production ratio relative to residential population.

    node_school_demand = {}
    _n_school = 0
    _n_capacity_used = 0
    for _, _poi_row in pois_utm.iterrows():
        _amenity = _poi_row.get("amenity") if "amenity" in pois_utm.columns else None
        if not (pd.notna(_amenity) and _amenity in _SCHOOL_TAGS):
            continue
        # Use OSM capacity if present and numeric, else fall back to type estimate
        _cap_raw = _poi_row.get("capacity") if "capacity" in pois_utm.columns else None
        if pd.notna(_cap_raw):
            try:
                _enrollment = float(_cap_raw)
                _n_capacity_used += 1
            except (ValueError, TypeError):
                _enrollment = SCHOOL_ENROLL_FALLBACK[_amenity]
        else:
            _enrollment = SCHOOL_ENROLL_FALLBACK[_amenity]
        _eidx = _edge_strtree.nearest(_poi_row.geometry)
        if _eidx in _ghost_junction:
            _junc = _ghost_junction[_eidx]
            node_school_demand[_junc] = node_school_demand.get(_junc, 0.0) + _enrollment
        else:
            _eu, _ev = _edge_keys[_eidx]
            _t = _edge_geom_list[_eidx].project(_poi_row.geometry, normalized=True)
            node_school_demand[_eu] = node_school_demand.get(_eu, 0.0) + _enrollment * (1.0 - _t)
            node_school_demand[_ev] = node_school_demand.get(_ev, 0.0) + _enrollment * _t
        _n_school += 1

    # External zone nodes have no school demand (schools are internal features).
    # Their school_demand entries will be set to 0 in the external node weight block below.

    _tot_sch = sum(node_school_demand.values())
    print(f"  {_n_school} school POIs → {len(node_school_demand)} nodes"
          f"  total enrollment={_tot_sch:.0f} pupils"
          f"  ({_n_capacity_used} from OSM capacity, {_n_school - _n_capacity_used} from fallback)")

    # ── Auto-detect boundary nodes from core polygon ──────────────────────────
    # Uses the raw (pre-consolidation) graph for OSM node IDs (OSRM-compatible).
    # build_network.py downloads the raw graph to core_max_vertex + ~1 km margin,
    # so boundary nodes' immediate external neighbours are present for the edge check.
    if _census_zones is not None:
        from shapely.geometry import Polygon as _Polygon
        _G_raw_bd = ox.load_graphml(f"{OUT_DIR}/newtownards_network.graphml")
        _core_poly_wgs = _Polygon(_census_zones["core_polygon"])
        _internal_ids = {n for n, d in _G_raw_bd.nodes(data=True)
                         if _core_poly_wgs.contains(Point(float(d["x"]), float(d["y"])))}
        _boundary_ids = {u for u, v in _G_raw_bd.edges()
                         if u in _internal_ids and v not in _internal_ids}

        # Map OSM boundary IDs → consolidated IDs for map display and build_paths.py
        # osmnx stores original OSM IDs in 'osmid_original': string for single nodes,
        # list of ints for merged nodes.
        _osmid_to_cons = {}
        for _cid, _cdata in G_cons.nodes(data=True):
            _osmids = _cdata.get("osmid_original")
            if _osmids is None:
                continue
            if isinstance(_osmids, list):
                for _oid in _osmids:
                    _osmid_to_cons[int(_oid)] = _cid
            else:
                _osmid_to_cons[int(_osmids)] = _cid
        _boundary_ids_cons = {_osmid_to_cons[o] for o in _boundary_ids if o in _osmid_to_cons}
        print(f"Auto-detected {len(_internal_ids)} internal nodes, "
              f"{len(_boundary_ids)} boundary nodes (OSM IDs)")
        print(f"  {len(_boundary_ids_cons)} boundary nodes map to consolidated graph")
    else:
        print(f"WARNING: {CENSUS_ZONES_FILE} not found — boundary_node_ids will be empty.")
        print(f"  Run build_census_zones.py first.")
        _boundary_ids = set()
        _boundary_ids_cons = set()
        _internal_ids = set()

    # ── Add external node weights from census data ────────────────────────────
    if _census_zones is not None:
        ext_nodes = _census_zones["external_nodes"]
        print(f"Adding {len(ext_nodes)} external node weights from census data "
              f"(ext_biz_scale={_ext_biz_scale:.4f}) …")
        for ext in ext_nodes:
            nid = ext["id"]
            node_population[nid]      = float(ext["population"])
            node_business_demand[nid] = float(ext["workplace_pop"]) * _ext_biz_scale
            node_school_demand[nid]   = float(ext["population"]) * _ext_school_per_pop
    else:
        ext_nodes = []

    # ── Serialise node weights for assignment script ───────────────────────────
    weights_path = f"{OUT_DIR}/node_weights.json"
    with open(weights_path, "w") as f:
        json.dump({
            "node_population":      {str(k): v for k, v in node_population.items()},
            "node_business_demand": {str(k): v for k, v in node_business_demand.items()},
            "node_school_demand":   {str(k): v for k, v in node_school_demand.items()},
            "node_parking_equiv":   {str(k): v for k, v in node_parking_equiv.items()},
            "boundary_node_ids":      sorted(_boundary_ids),
            "boundary_node_ids_cons": sorted(_boundary_ids_cons),
            "internal_node_ids":      sorted(_internal_ids),
        }, f)
    print(f"Saved node weights → {weights_path}"
          f"  ({len(node_ids)} internal + {len(ext_nodes)} external nodes)")
    print("Next: python3 simulation/build_map.py  (build the interactive map)")
