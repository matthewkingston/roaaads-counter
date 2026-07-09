"""
Fetch NISRA Census 2021 Data Zone population data, filter to the Newtownards
study area, and compute per-node demand weights:
  - residential population (OSM-building-snapped; road-length fallback per DZ)
  - business demand (workplace population, POI-distributed + retail parking spaces)
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
    BUILDING_CACHE, PARKING_ISLAND_CACHE, SCHOOL_ISLAND_CACHE, CENSUS_ZONES_FILE,
    EXCLUDE_AMENITY, POI_WEIGHTS, PROJECTED_CRS,
)
from parking_demand import parking_spaces
from school_attractor import add_level_enrolments, LEVEL_ENROL_COLS
from census_supply import load_supply
import census_attractor
import census_school_producers

# Phase-2 per-level school node layers (primary/post-primary/tertiary), each paired with the single
# total layer. Attractor: node_school_demand_<lvl>; producer: node_school_producers_<lvl>.
SCHOOL_LEVELS = ("primary", "postprimary", "tertiary")
SCHOOL_LEVEL_LAYERS    = ["node_school_demand_"    + l for l in SCHOOL_LEVELS]  # attractor (aligned w/ LEVEL_ENROL_COLS)
SCHOOL_PRODUCER_LAYERS = ["node_school_producers_" + l for l in SCHOOL_LEVELS]  # producer

# ── CLI ──────────────────────────────────────────────────────────────────────
_ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument(
    "--zones-only", action="store_true",
    help="Patch only external-zone weights in node_weights.json (fast). "
         "Use after re-running build_census_zones.py (external retail/census); "
         "not when the core polygon or internal demographics change.")
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

# External school demand is now measured per zone (retail_spaces-style) from the island
# school cache and written into census_zones.json by build_census_zones.py — the old
# population × ext_school_per_pop approximation is removed (2026-06-28).

# ── Fast path: --zones-only ────────────────────────────────────────────────────
# Patches only external-zone entries in node_weights.json. Skips all internal
# population processing. Use after re-running build_census_zones.py
# (not when the core polygon or internal demographics have changed).

if args.zones_only:
    if _census_zones is None:
        print(f"ERROR: {CENSUS_ZONES_FILE} not found. Run build_census_zones.py first.")
        sys.exit(1)
    weights_path = f"{OUT_DIR}/node_weights.json"
    with open(weights_path) as f:
        w = json.load(f)
    print("Updating external zone weights from census data …")
    ext_nodes = _census_zones["external_nodes"]
    w.setdefault("node_retail_spaces", {})
    w.setdefault("node_workplace", {})
    for ext in ext_nodes:
        nid = ext["id"]
        _retail = float(ext.get("retail_spaces", 0.0))
        w["node_population"][str(nid)]    = float(ext["population"])
        w["node_retail_spaces"][str(nid)] = _retail
        w["node_workplace"][str(nid)]     = float(ext["workplace_pop"])
        w.setdefault("node_commute_attractor", {})[str(nid)] = float(ext.get("commute_attractor", 0.0))
        for _lyr in SCHOOL_LEVEL_LAYERS:
            w.setdefault(_lyr, {})[str(nid)] = float(ext.get("school_demand_" + _lyr.split("_")[-1], 0.0))
        w.setdefault("node_commute_producers", {})[str(nid)] = float(ext.get("commute_producers", 0.0))
        for _lyr in SCHOOL_PRODUCER_LAYERS:
            w.setdefault(_lyr, {})[str(nid)] = float(ext.get("school_producers_" + _lyr.split("_")[-1], 0.0))
        _sch_tot = sum(float(ext.get("school_demand_" + l, 0.0)) for l in SCHOOL_LEVELS)
        print(f"  {nid}  {ext['level']}  "
              f"pop={ext['population']:>8,}  wp={ext['workplace_pop']:>8,}  "
              f"retail_sp={_retail:>8,.0f}  school={_sch_tot:>8,.1f}")
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

    transformer_to_utm = pyproj.Transformer.from_crs("EPSG:4326", PROJECTED_CRS, always_xy=True)
    transformer_to_wgs = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)

    # Work in UTM throughout for accurate area calculations
    dz_utm = dz.to_crs(PROJECTED_CRS)
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

    # Consolidated graph nodes already have x/y in PROJECTED_CRS (ITM)
    node_ids = list(G_cons.nodes())
    node_coords_utm = [(G_cons.nodes[n]["x"], G_cons.nodes[n]["y"]) for n in node_ids]

    # Spatial join: which core DZs contain at least one junction node?
    node_gdf = gpd.GeoDataFrame(
        {"node_id": node_ids},
        geometry=[Point(x, y) for x, y in node_coords_utm],
        crs=PROJECTED_CRS,
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
        _bld_raw = _bld_raw.to_crs(PROJECTED_CRS)
        _bld_raw["geometry"] = _bld_raw.geometry.centroid
        _bld_raw = _bld_raw[_bld_raw.geometry.geom_type == "Point"]
        _bld_raw = _bld_raw[["geometry"]].to_crs("EPSG:4326")
        _bld_raw.to_file(BUILDING_CACHE, driver="GeoJSON")
        print(f"  Cached → {BUILDING_CACHE}")

    # Reset to clean integer index so sjoin indices are valid positional offsets
    buildings_utm = _bld_raw.to_crs(PROJECTED_CRS).reset_index(drop=True)
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
    pois_utm = pois_raw.to_crs(PROJECTED_CRS).copy()
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

    # Allocate workplace population to nodes within each DZ by edge-snapped POI weight.
    # node_workplace is the clean place-of-work jobs layer (the commute component's
    # attractor) — retail parking spaces are kept separately in node_retail_spaces and
    # are NOT folded in here.
    node_workplace = {}
    for dz_code, nodes_in_dz in dz_to_nodes.items():
        dz_wp = float(wp_lookup.get(dz_code, 0) or 0)
        poi_weights_dz = {n: node_poi_weight.get(n, 0.0) for n in nodes_in_dz}
        total_weight = sum(poi_weights_dz.values())
        for n in nodes_in_dz:
            if total_weight > 0:
                node_workplace[n] = dz_wp * poi_weights_dz[n] / total_weight
            else:
                node_workplace[n] = dz_wp / len(nodes_in_dz)

    active_wp_nodes = sum(1 for v in node_workplace.values() if v > 0)
    total_wp = sum(node_workplace.values())
    print(f"  {active_wp_nodes} nodes with workplace jobs · total {total_wp:.0f} workplace pop attributed")

    # Car-commute attractor (jobs reached by car) — the commute component's attractor. Sourced
    # island-wide per census area by census_attractor.load_attractor() (NI: apwp001 × SDZ
    # car-share; RoI: WZ driver count), then distributed across each area's nodes by the SAME
    # edge-snapped POI weights as node_workplace. Portable: no apwp001 read, no NI assumption.
    _commute_attr = census_attractor.load_attractor()
    node_commute_attractor = {}
    for dz_code, nodes_in_dz in dz_to_nodes.items():
        dz_car = float(_commute_attr.get(dz_code, 0.0))
        poi_weights_dz = {n: node_poi_weight.get(n, 0.0) for n in nodes_in_dz}
        total_weight = sum(poi_weights_dz.values())
        for n in nodes_in_dz:
            if total_weight > 0:
                node_commute_attractor[n] = dz_car * poi_weights_dz[n] / total_weight
            else:
                node_commute_attractor[n] = dz_car / len(nodes_in_dz)
    print(f"  internal commute attractor (car jobs) total "
          f"{sum(node_commute_attractor.values()):.0f}")

    # ── Trip producers (commute, school) — census per-DZ counts distributed to nodes ──
    # Residents → distribute each DZ's producer totals across its nodes by population share.
    # Stored as separate layers; the school components use node_school_producers_<level> as their
    # producing weights, node_commute_producers feeds the commute component.
    _supply = load_supply()                                   # commute producers (census_supply)
    _school_prod = census_school_producers.load_school_producers()   # per-level school producers
    node_commute_producers = {}
    node_school_prod_levels = {lyr: {} for lyr in SCHOOL_PRODUCER_LAYERS}
    for dz_code, nodes_in_dz in dz_to_nodes.items():
        _com = float(_supply.get(dz_code, {}).get("commute", 0.0))
        _sp = _school_prod.get(dz_code, {})
        _lvl_vals = [float(_sp.get(l, 0.0)) for l in SCHOOL_LEVELS]  # aligned w/ SCHOOL_PRODUCER_LAYERS
        pops = {n: node_population.get(n, 0.0) for n in nodes_in_dz}
        tot = sum(pops.values())
        for n in nodes_in_dz:
            frac = (pops[n] / tot) if tot > 0 else (1.0 / len(nodes_in_dz))
            node_commute_producers[n] = node_commute_producers.get(n, 0.0) + _com * frac
            for _lyr, _v in zip(SCHOOL_PRODUCER_LAYERS, _lvl_vals):
                node_school_prod_levels[_lyr][n] = node_school_prod_levels[_lyr].get(n, 0.0) + _v * frac
    _sch_prod_tot = sum(sum(d.values()) for d in node_school_prod_levels.values())
    print(f"  internal producers: commute {sum(node_commute_producers.values()):.0f}, "
          f"school {_sch_prod_tot:.0f} "
          f"(" + " / ".join(f"{lyr.split('_')[-1]} {sum(node_school_prod_levels[lyr].values()):.0f}"
                            for lyr in SCHOOL_PRODUCER_LAYERS) + ")")

    # ── Car-park retail demand (estimated parking spaces) ─────────────────────────
    # OSM parking polygons proxy retail/customer vehicle-trip attraction. Each lot is
    # turned into an estimate of retail parking *spaces* by parking_demand.parking_spaces
    # (access filter + capacity plausibility gate + area fallback; see CLAUDE.md /
    # plan_parking_retail) — replacing the old area/25-/50 "equivalent persons" hack.
    # Spaces are the retail attractor's currency (the tuner's K absorbs spaces→trips).
    # Sourced from the island-wide cache (build_parking.py) so internal and external
    # zones use one estimator with identical tag handling.

    if not os.path.exists(PARKING_ISLAND_CACHE):
        print(f"ERROR: {PARKING_ISLAND_CACHE} not found. "
              f"Run: python3 simulation/build_parking.py")
        sys.exit(1)
    print(f"Loading island parking from cache ({PARKING_ISLAND_CACHE}) …")
    _park_raw = gpd.read_file(PARKING_ISLAND_CACHE)

    # Polygon features only; clip to the core polygon (internal retail attractor).
    _park_utm = _park_raw.to_crs(PROJECTED_CRS).copy()
    _park_utm = _park_utm[_park_utm.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    _park_utm["area_m2"] = _park_utm.geometry.area
    _park_utm["centroid_geom"] = _park_utm.geometry.centroid
    _park_utm = _park_utm[_park_utm["centroid_geom"].within(core_poly_utm)].copy()

    # Estimate each lot's retail spaces, snap to the nearest road edge, split by t.
    _park_tag_cols = [c for c in _park_utm.columns
                      if c not in ("geometry", "centroid_geom", "area_m2")]
    node_retail_spaces = {}
    _total_spaces = 0.0
    _n_kept = 0
    for _, _row in _park_utm.iterrows():
        _spaces = parking_spaces({c: _row[c] for c in _park_tag_cols}, _row["area_m2"])
        if _spaces <= 0:
            continue
        _n_kept += 1
        _total_spaces += _spaces
        _pt = _row["centroid_geom"]
        _eidx = _edge_strtree.nearest(_pt)
        if _eidx in _ghost_junction:
            _junc = _ghost_junction[_eidx]
            node_retail_spaces[_junc]   = node_retail_spaces.get(_junc, 0.0) + _spaces
        else:
            _eu, _ev = _edge_keys[_eidx]
            _t = _edge_geom_list[_eidx].project(_pt, normalized=True)
            _eu_share, _ev_share = (1.0 - _t) * _spaces, _t * _spaces
            node_retail_spaces[_eu]   = node_retail_spaces.get(_eu, 0.0) + _eu_share
            node_retail_spaces[_ev]   = node_retail_spaces.get(_ev, 0.0) + _ev_share

    print(f"  {len(_park_utm)} parking polygons in core → "
          f"{_n_kept} retail lots, {_total_spaces:.0f} estimated spaces")

    # ── School demand layer (estimated enrolment) ─────────────────────────────────
    # Per-POI enrolment — primary/secondary (jurisdiction-aware sourced averages),
    # curated/sourced universities split across their POIs, kindergarten, SEN — is
    # precomputed island-wide by build_schools.py (school_demand.assign_enrolments).
    # Here we clip to the core polygon and snap each school's enrolment to a road edge.

    if not os.path.exists(SCHOOL_ISLAND_CACHE):
        print(f"ERROR: {SCHOOL_ISLAND_CACHE} not found. "
              f"Run: python3 simulation/build_schools.py")
        sys.exit(1)
    print(f"Loading island schools from cache ({SCHOOL_ISLAND_CACHE}) …")
    _sch_utm = gpd.read_file(SCHOOL_ISLAND_CACHE).to_crs(PROJECTED_CRS).copy()
    _sch_utm = add_level_enrolments(_sch_utm)          # + enrol_primary/postprimary/tertiary (special split)
    _sch_utm["centroid_geom"] = _sch_utm.geometry.centroid
    _sch_utm = _sch_utm[_sch_utm["centroid_geom"].within(core_poly_utm)].copy()

    # Three per-level school-demand layers, all snapped to the same road node(s) per POI.
    node_school_levels = {lyr: {} for lyr in SCHOOL_LEVEL_LAYERS}
    _n_school = 0
    def _add(_d, _k, _v):
        _d[_k] = _d.get(_k, 0.0) + _v
    for _, _row in _sch_utm.iterrows():
        _enrollment = float(_row["enrolment"])
        if _enrollment <= 0:
            continue
        _pt = _row["centroid_geom"]
        _eidx = _edge_strtree.nearest(_pt)
        _lvl_vals = [float(_row[c]) for c in LEVEL_ENROL_COLS]   # aligned with SCHOOL_LEVEL_LAYERS
        if _eidx in _ghost_junction:
            _junc = _ghost_junction[_eidx]
            for _lyr, _v in zip(SCHOOL_LEVEL_LAYERS, _lvl_vals):
                _add(node_school_levels[_lyr], _junc, _v)
        else:
            _eu, _ev = _edge_keys[_eidx]
            _t = _edge_geom_list[_eidx].project(_pt, normalized=True)
            for _lyr, _v in zip(SCHOOL_LEVEL_LAYERS, _lvl_vals):
                _add(node_school_levels[_lyr], _eu, _v * (1.0 - _t))
                _add(node_school_levels[_lyr], _ev, _v * _t)
        _n_school += 1

    # External zone nodes get their school demand from census_zones.json (set in the
    # external node weight block below); internal core schools are snapped above.
    _tot_sch = sum(sum(d.values()) for d in node_school_levels.values())
    print(f"  {_n_school} school POIs in core → total enrolment={_tot_sch:.0f} pupils "
          f"(" + " / ".join(f"{lyr.split('_')[-1]} {sum(node_school_levels[lyr].values()):.0f}"
                            for lyr in SCHOOL_LEVEL_LAYERS) + ")")

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
        # A node is on the boundary if it is the internal endpoint of an edge that
        # crosses the core polygon — in EITHER direction. Direction matters on one-way
        # dual carriageways split by the core boundary: the outbound lane's node has an
        # internal→external edge, but the inbound lane's entry node only has an
        # external→internal edge, so an outbound-only test silently drops it (leaving
        # that direction of traffic with no clean entry point).
        _boundary_ids = set()
        for _u, _v in _G_raw_bd.edges():
            _ui, _vi = _u in _internal_ids, _v in _internal_ids
            if _ui != _vi:                       # exactly one endpoint inside the core
                _boundary_ids.add(_u if _ui else _v)

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
        print(f"Adding {len(ext_nodes)} external node weights from census data …")
        _missing_retail = 0
        for ext in ext_nodes:
            nid = ext["id"]
            _retail = float(ext.get("retail_spaces", 0.0))
            if "retail_spaces" not in ext:
                _missing_retail += 1
            node_population[nid]      = float(ext["population"])
            node_retail_spaces[nid]   = _retail
            node_workplace[nid]       = float(ext["workplace_pop"])
            node_commute_attractor[nid] = float(ext.get("commute_attractor", 0.0))
            for _lyr in SCHOOL_LEVEL_LAYERS:            # school_demand_<lvl> → node_school_demand_<lvl>
                node_school_levels[_lyr][nid] = float(ext.get("school_demand_" + _lyr.split("_")[-1], 0.0))
            node_commute_producers[nid] = float(ext.get("commute_producers", 0.0))
            for _lyr in SCHOOL_PRODUCER_LAYERS:        # school_producers_<lvl> → node_school_producers_<lvl>
                node_school_prod_levels[_lyr][nid] = float(ext.get("school_producers_" + _lyr.split("_")[-1], 0.0))
        if _missing_retail:
            print(f"  WARNING: {_missing_retail}/{len(ext_nodes)} external nodes lack "
                  f"'retail_spaces' — re-run build_parking.py then build_census_zones.py. "
                  f"Treated as 0 (commute-only) for now.")
    else:
        ext_nodes = []

    # ── Serialise node weights for assignment script ───────────────────────────
    weights_path = f"{OUT_DIR}/node_weights.json"
    with open(weights_path, "w") as f:
        json.dump({
            "node_population":      {str(k): v for k, v in node_population.items()},
            "node_workplace":       {str(k): v for k, v in node_workplace.items()},
            "node_commute_attractor": {str(k): v for k, v in node_commute_attractor.items()},
            **{lyr: {str(k): v for k, v in node_school_levels[lyr].items()} for lyr in SCHOOL_LEVEL_LAYERS},
            "node_retail_spaces":   {str(k): v for k, v in node_retail_spaces.items()},
            "node_commute_producers": {str(k): v for k, v in node_commute_producers.items()},
            **{lyr: {str(k): v for k, v in node_school_prod_levels[lyr].items()} for lyr in SCHOOL_PRODUCER_LAYERS},
            "boundary_node_ids":      sorted(_boundary_ids),
            "boundary_node_ids_cons": sorted(_boundary_ids_cons),
            "internal_node_ids":      sorted(_internal_ids),
        }, f)
    print(f"Saved node weights → {weights_path}"
          f"  ({len(node_ids)} internal + {len(ext_nodes)} external nodes)")
    print("Next: python3 simulation/build_map.py  (build the interactive map)")
