"""
Build the hierarchical census-area external zone network.

Works for any CENTRE on the island of Ireland: NI areas use DZ/SDZ/DEA from NISRA;
RoI areas use SA/ED/LEA from CSO.  Both are loaded and merged before classification,
so changing zones_config.CENTRE to a RoI town produces a valid output automatically.

Hierarchy
---------
  Core area   : small areas (DZ or SA) that intersect CORE_RADIUS
  Intermediate: SDZs/EDs within SDZ_ZONE_RADIUS whose parent outer-zone is "broken"
                and that have no core small areas
  Orphan      : non-core small areas from partially-core intermediate zones
  Outer       : DEAs/LEAs outside SDZ_ZONE_RADIUS (one centroid node each)

Data files required (NI)
    simulation/dz2021/DZ2021.geojson
    simulation/sdz2021/SDZ2021.geojson
    simulation/dea2021/DEA2021.geojson
    data/census-2021-apwp001.xlsx

Data files required (RoI, in data/ireland_data/)
    Small_Area_National_Statistical_Boundaries_2022_*.geojson
    Complete_set_of_Census_2022_SAPs/SAPS_2022_Small_Area_UR_171024.csv
    Workplace_Zones_ITM/Workplace_Zones_ITM.shp   (WZ boundaries in ITM)
    cache_sa_workplace.csv   (pre-computed WZ→SA apportionment)
        → run simulation/build_wz_apportionment.py once to generate this

Output: data/census_zones.json  (same schema as before; external_nodes now include
RoI LEA/ED/SA entries alongside NI DEA/SDZ/DZ entries)
"""

import json, math, os
import geopandas as gpd
import pandas as pd
import numpy as np
import pyproj
from shapely.geometry import Point, mapping
from shapely.ops import unary_union
import shapely.ops as sops

from zones_config import CENTRE, CORE_RADIUS, SDZ_ZONE_RADIUS
from demographics_config import PROJECTED_CRS, PARKING_ISLAND_CACHE, SCHOOL_ISLAND_CACHE
from parking_demand import parking_spaces
import census_supply
from ingest_ni_census import load_ni_census
from ingest_roi_census import load_roi_census

OUTPUT_FILE = "data/census_zones.json"

# ── Coordinate transformers ───────────────────────────────────────────────────

to_utm = pyproj.Transformer.from_crs("EPSG:4326", PROJECTED_CRS, always_xy=True)
to_wgs = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)

centre_utm_x, centre_utm_y = to_utm.transform(CENTRE[1], CENTRE[0])
core_circle = Point(centre_utm_x, centre_utm_y).buffer(CORE_RADIUS)
sdz_circle  = Point(centre_utm_x, centre_utm_y).buffer(SDZ_ZONE_RADIUS)

# ── Load and combine NI + RoI census data ────────────────────────────────────

dz_gdf,  sdz_gdf,  dea_gdf  = load_ni_census()
sa_gdf,  ed_gdf,   lea_gdf  = load_roi_census()

# Combined GeoDataFrames — classification logic operates on these uniformly.
# "dz" = all small areas (NI DZ + RoI SA)
# "sdz" = all intermediate zones (NI SDZ + RoI ED)
# "dea" = all outer zones (NI DEA + RoI LEA)
dz  = pd.concat([dz_gdf,  sa_gdf],  ignore_index=True)
sdz = pd.concat([sdz_gdf, ed_gdf],  ignore_index=True)
dea = pd.concat([dea_gdf, lea_gdf], ignore_index=True)

# Per-small-area trip producers (commute, school) — NI DZ + RoI SA, harmonised (RoI ÷2).
_supply = census_supply.load_supply()
dz["commute_producers"] = dz["area_code"].map(lambda c: _supply.get(c, {}).get("commute", 0.0))
dz["school_producers"]  = dz["area_code"].map(lambda c: _supply.get(c, {}).get("school", 0.0))
_n_missing = int((~dz["area_code"].isin(_supply)).sum())
print(f"Producers: matched {len(dz) - _n_missing}/{len(dz)} small areas to census_supply"
      + (f" — {_n_missing} unmatched (treated as 0)" if _n_missing else ""))

# ── Hierarchy lookups ────────────────────────────────────────────────────────
# small-area code → parent intermediate-zone code
dz_to_sdz  = dz.set_index("area_code")["parent_code"].to_dict()
# intermediate-zone code → parent outer-zone code
sdz_to_dea = sdz.set_index("area_code")["parent_code"].to_dict()

# ── Core classification ───────────────────────────────────────────────────────

dz["is_core"] = dz.geometry.intersects(core_circle)
core_sa_codes = set(dz.loc[dz["is_core"], "area_code"])

sdz["in_sdz_zone"] = sdz.geometry.intersects(sdz_circle)
dea["in_sdz_zone"] = dea.geometry.intersects(sdz_circle)

n_core     = dz["is_core"].sum()
n_sdz_zone = sdz["in_sdz_zone"].sum()
n_dea_broken = dea["in_sdz_zone"].sum()
n_dea_single = (~dea["in_sdz_zone"]).sum()

print(f"\nClassification (CENTRE {CENTRE[0]:.4f}°N, {CENTRE[1]:.4f}°E):")
print(f"  {n_core} small areas intersect CORE_RADIUS ({CORE_RADIUS} m) → core")
print(f"  {n_sdz_zone} intermediate zones intersect SDZ_ZONE_RADIUS ({SDZ_ZONE_RADIUS} m)")
print(f"  {n_dea_broken} outer zones intersect SDZ_ZONE_RADIUS → broken into intermediate nodes")
print(f"  {n_dea_single} outer zones remain as single external nodes")

# Intermediate zones that have at least one core small area
_int_has_core = {sdz_cd: True for sa_cd, sdz_cd in dz_to_sdz.items()
                 if sa_cd in core_sa_codes}
sdz["has_core_sa"] = sdz["area_code"].map(_int_has_core).fillna(False)

# Build core polygon (union of all core small-area geometries)
core_sa_gdf  = dz[dz["is_core"]]
core_polygon = unary_union(core_sa_gdf.geometry.values)

n_core_intermediates = len({dz_to_sdz.get(cd) for cd in core_sa_codes} - {None})
print(f"\nCore area: {len(core_sa_gdf)} small areas (from {n_core_intermediates} partially/fully-core intermediate zones)")
print(f"  Core polygon area: {core_polygon.area / 1e6:.2f} km²")

centre_pt_utm   = Point(centre_utm_x, centre_utm_y)
max_vertex_dist = max(centre_pt_utm.distance(Point(x, y))
                      for x, y in core_polygon.exterior.coords)
print(f"  Max core polygon vertex distance: {max_vertex_dist:.0f} m")

# Outer zones that are broken into intermediate nodes
broken_outer_codes = set(dea.loc[dea["in_sdz_zone"], "area_code"])

# ── Weighted centroid helper ─────────────────────────────────────────────────

def weighted_centroid_wgs(geometries, weights):
    """Population-weighted centroid in WGS84 → (lat, lon, utm_x, utm_y)."""
    weights = np.array(weights, dtype=float)
    if weights.sum() == 0:
        weights = np.ones(len(weights))
    weights /= weights.sum()
    centroids = [g.centroid for g in geometries]
    cx = sum(c.x * w for c, w in zip(centroids, weights))
    cy = sum(c.y * w for c, w in zip(centroids, weights))
    lon, lat = to_wgs.transform(cx, cy)
    return lat, lon, cx, cy

# ── Build external node list ─────────────────────────────────────────────────

external_nodes = []
zone_geom = {}   # external node id → its census polygon (ITM), for parking aggregation
print("\nBuilding external node list …")

# 1. Intermediate external nodes (SDZ/ED): in broken outer zones with no core SAs
sdz["parent_outer"] = sdz["area_code"].map(sdz_to_dea)
sdz_in_broken_outer = sdz["parent_outer"].isin(broken_outer_codes)
sdz_external_mask   = sdz_in_broken_outer & ~sdz["has_core_sa"]
sdz_external        = sdz[sdz_external_mask]

for _, row in sdz_external.iterrows():
    int_code      = row["area_code"]
    child_sa_codes = [sa for sa, p in dz_to_sdz.items() if p == int_code]
    child_sa       = dz[dz["area_code"].isin(child_sa_codes)]
    pop = child_sa["population"].sum()
    wp  = child_sa["workplace_pop"].sum()
    lat, lon, utm_x, utm_y = weighted_centroid_wgs(
        child_sa.geometry.values, child_sa["population"].values)
    zone_geom[int_code] = child_sa.geometry.unary_union
    external_nodes.append({
        "id":             int_code,
        "level":          row["level"],
        "centroid_lat":   round(lat,   6),
        "centroid_lon":   round(lon,   6),
        "centroid_utm_x": round(utm_x, 1),
        "centroid_utm_y": round(utm_y, 1),
        "population":     int(round(pop)),
        "workplace_pop":  int(round(wp)),
        "commute_producers": int(round(child_sa["commute_producers"].sum())),
        "school_producers":  int(round(child_sa["school_producers"].sum())),
    })
print(f"  {len(external_nodes)} intermediate external nodes (SDZ/ED)")

# 2. Orphan small-area nodes: non-core SAs from partially-core intermediate zones
#    (intermediate zone has some core SAs but these SAs themselves are outside the core)
partial_core_int_codes = set(sdz.loc[sdz_in_broken_outer & sdz["has_core_sa"], "area_code"])
orphan_sa_mask = (
    dz["area_code"].map(dz_to_sdz).isin(partial_core_int_codes)
    & ~dz["is_core"]
)
orphan_sa = dz[orphan_sa_mask]
n_int_start = len(external_nodes)

for _, row in orphan_sa.iterrows():
    cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
    lon, lat = to_wgs.transform(cx, cy)
    zone_geom[row["area_code"]] = row.geometry
    external_nodes.append({
        "id":             row["area_code"],
        "level":          row["level"],
        "centroid_lat":   round(lat,   6),
        "centroid_lon":   round(lon,   6),
        "centroid_utm_x": round(cx,    1),
        "centroid_utm_y": round(cy,    1),
        "population":     int(round(row["population"])),
        "workplace_pop":  int(round(row["workplace_pop"])),
        "commute_producers": int(round(row["commute_producers"])),
        "school_producers":  int(round(row["school_producers"])),
    })
print(f"  {len(external_nodes) - n_int_start} orphan small-area external nodes (DZ/SA)")

# 3. Outer external nodes (DEA/LEA): outer zones NOT intersecting SDZ_ZONE_RADIUS
dea_external = dea[~dea["in_sdz_zone"]]
n_outer_start = len(external_nodes)

for _, row in dea_external.iterrows():
    outer_code      = row["area_code"]
    child_int_codes = {s for s, d in sdz_to_dea.items() if d == outer_code}
    child_sa_codes  = {sa for sa, p in dz_to_sdz.items() if p in child_int_codes}
    child_sa        = dz[dz["area_code"].isin(child_sa_codes)]
    pop = child_sa["population"].sum()
    wp  = child_sa["workplace_pop"].sum()
    cprod = child_sa["commute_producers"].sum()
    sprod = child_sa["school_producers"].sum()
    if len(child_sa) == 0:
        cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
        lon, lat = to_wgs.transform(cx, cy)
        utm_x, utm_y = cx, cy
        pop, wp, cprod, sprod = 0, 0, 0, 0
    else:
        lat, lon, utm_x, utm_y = weighted_centroid_wgs(
            child_sa.geometry.values, child_sa["population"].values)
    zone_geom[outer_code] = (child_sa.geometry.unary_union if len(child_sa) else row.geometry)
    external_nodes.append({
        "id":             outer_code,
        "level":          row["level"],
        "centroid_lat":   round(lat,   6),
        "centroid_lon":   round(lon,   6),
        "centroid_utm_x": round(utm_x, 1),
        "centroid_utm_y": round(utm_y, 1),
        "population":     int(round(pop)),
        "workplace_pop":  int(round(wp)),
        "commute_producers": int(round(cprod)),
        "school_producers":  int(round(sprod)),
    })
print(f"  {len(external_nodes) - n_outer_start} outer external nodes (DEA/LEA)")
print(f"  {len(external_nodes)} external nodes total")
print(f"  Total external pop: {sum(n['population'] for n in external_nodes):,}")
print(f"  Total external wp:  {sum(n['workplace_pop'] for n in external_nodes):,}")
print(f"  Total external commute producers: {sum(n['commute_producers'] for n in external_nodes):,}")
print(f"  Total external school producers:  {sum(n['school_producers'] for n in external_nodes):,}")

# ── External retail demand: estimated parking spaces per zone ─────────────────
# Sum parking_demand.parking_spaces over every parking polygon whose centroid falls
# inside each external zone's census polygon — the SAME estimator build_demographics.py
# applies to internal core nodes, so internal and external retail are on one scale
# (this is what lets ext_biz_scale be removed). Island-wide parking from build_parking.py.
print("\nExternal retail demand (parking spaces per zone) …")
if not os.path.exists(PARKING_ISLAND_CACHE):
    raise SystemExit(f"ERROR: {PARKING_ISLAND_CACHE} not found. "
                     f"Run: python3 simulation/build_parking.py")
_park = gpd.read_file(PARKING_ISLAND_CACHE).to_crs(PROJECTED_CRS)
_park = _park[_park.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
_park_tag_cols = [c for c in _park.columns if c != "geometry"]
_park["spaces"] = [
    parking_spaces({c: r[c] for c in _park_tag_cols}, g.area)
    for r, g in zip(_park.to_dict("records"), _park.geometry)
]
_park = _park[_park["spaces"] > 0].copy()
_park["geometry"] = _park.geometry.centroid          # assign lots by centroid

# Spatial join parking centroids → zone polygons, sum spaces per zone.
_zones_gdf = gpd.GeoDataFrame(
    {"id": list(zone_geom.keys())},
    geometry=list(zone_geom.values()), crs=PROJECTED_CRS)
_joined = gpd.sjoin(_park[["spaces", "geometry"]], _zones_gdf,
                    how="inner", predicate="within")
_spaces_by_zone = _joined.groupby("id")["spaces"].sum().to_dict()

# Workplace-derived fallback for zones with no mapped parking: scale the zone's
# workplace_pop by the island-wide median spaces-per-workplace ratio (over matched
# zones). Loud-logged, never a silent 0 (some small zones genuinely lack OSM parking).
_ratios = [_spaces_by_zone[n["id"]] / n["workplace_pop"]
           for n in external_nodes
           if _spaces_by_zone.get(n["id"], 0) > 0 and n["workplace_pop"] > 0]
_fallback_ratio = float(np.median(_ratios)) if _ratios else 0.0
_fallback_ids = []
for n in external_nodes:
    sp = _spaces_by_zone.get(n["id"], 0.0)
    if sp <= 0:
        sp = n["workplace_pop"] * _fallback_ratio
        if sp > 0:
            _fallback_ids.append(n["id"])
    n["retail_spaces"] = round(float(sp), 1)

print(f"  {len(_park)} retail lots, {sum(_spaces_by_zone.values()):,.0f} spaces "
      f"matched into {len(_spaces_by_zone)} zones")
print(f"  median spaces/workplace ratio = {_fallback_ratio:.3f}")
if _fallback_ids:
    print(f"  {len(_fallback_ids)} zones had NO mapped parking → workplace-derived "
          f"fallback: {', '.join(map(str, _fallback_ids))}")
print(f"  Total external retail spaces: {sum(n['retail_spaces'] for n in external_nodes):,.0f}")

# ── External school demand: enrolment per zone ────────────────────────────────
# Sum the per-POI enrolment from the island school cache (build_schools.py — already
# primary/secondary jurisdiction-aware, curated universities, etc.) over each external
# zone's census polygon. Same estimator as internal core nodes → ext_school_per_pop
# removed. Zones with no mapped school get 0 (a village with no school genuinely
# attracts no school trips; its school-age population's trips distribute to other zones).
print("\nExternal school demand (enrolment per zone) …")
if not os.path.exists(SCHOOL_ISLAND_CACHE):
    raise SystemExit(f"ERROR: {SCHOOL_ISLAND_CACHE} not found. "
                     f"Run: python3 simulation/build_schools.py")
_sch = gpd.read_file(SCHOOL_ISLAND_CACHE).to_crs(PROJECTED_CRS)
_sch["geometry"] = _sch.geometry.centroid
_sjoin = gpd.sjoin(_sch[["enrolment", "geometry"]], _zones_gdf, how="inner", predicate="within")
_school_by_zone = _sjoin.groupby("id")["enrolment"].sum().to_dict()
for n in external_nodes:
    n["school_demand"] = round(float(_school_by_zone.get(n["id"], 0.0)), 1)
print(f"  {len(_sch)} school POIs, {sum(_school_by_zone.values()):,.0f} pupils "
      f"matched into {len(_school_by_zone)} zones")
print(f"  {sum(1 for n in external_nodes if n['school_demand'] == 0)} zones with no mapped school (school_demand=0)")
print(f"  Total external school demand: {sum(n['school_demand'] for n in external_nodes):,.0f} pupils")

# ── Serialise core polygon to WGS84 ─────────────────────────────────────────

core_polygon_wgs = sops.transform(lambda x, y: to_wgs.transform(x, y), core_polygon)
core_coords = list(mapping(core_polygon_wgs)["coordinates"][0])

# ── Write output ─────────────────────────────────────────────────────────────

output = {
    "core_radius":             CORE_RADIUS,
    "sdz_zone_radius":         SDZ_ZONE_RADIUS,
    "centre_lat":              CENTRE[0],
    "centre_lon":              CENTRE[1],
    "max_core_vertex_dist_m":  round(max_vertex_dist, 1),
    "core_polygon":            core_coords,
    "n_core_dzs":              int(len(core_sa_gdf)),
    "n_core_sdzs":             n_core_intermediates,
    "external_nodes":          external_nodes,
}

with open(OUTPUT_FILE, "w") as f:
    json.dump(output, f, indent=2)

print(f"\nSaved {OUTPUT_FILE}")
print(f"  core polygon: {len(core_coords)} vertices")
print(f"  external nodes: {len(external_nodes)}")
print(f"\nNext: python3 simulation/build_network.py  (if not already done)")
print(f"      python3 simulation/build_demographics.py")
print(f"      python3 simulation/build_external_links.py")
