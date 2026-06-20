"""
Build the hierarchical census-area external zone network for Newtownards.

Classifies all of NI into a three-level hierarchy:
  - Core area  : union of DZs that directly intersect CORE_RADIUS
  - DZ nodes   : non-core DZs from partially-core SDZs (individual DZ-level external nodes)
  - SDZ nodes  : SDZs within SDZ_ZONE_RADIUS with no core DZs
  - DEA nodes  : DEAs outside SDZ_ZONE_RADIUS (represented by a single centroid)

Population-weighted centroids are computed for each external (SDZ/DEA) node
from DZ-level census population data.

Data files required:
  simulation/dz2021/DZ2021.geojson   — DZ polygon boundaries
    https://www.nisra.gov.uk/support/geography/data-zones-census-2021
  simulation/sdz2021/SDZ2021.geojson  — SDZ polygon boundaries
    https://www.nisra.gov.uk/support/geography/super-data-zones-census-2021
  simulation/dea2021/DEA2021.geojson  — DEA polygon boundaries (OSNI Open Data)
    https://admin.opendatani.gov.uk/dataset/osni-open-data-largescale-boundaries-district-electoral-areas-2012
  data/census-2021-apwp001.xlsx       — DZ-level workplace population

Population is fetched from the NISRA API (cached to data/cache_nisra_population.csv).

Output: data/census_zones.json
"""

import json, math, os, sys, urllib.request
import geopandas as gpd
import pandas as pd
import numpy as np
import pyproj
from shapely.geometry import Point, mapping
from shapely.ops import unary_union

# ── Config ─────────────────────────────────────────────────────────────────────

CENTRE          = (54.5933779, -5.6960935)   # (lat, lon) of Newtownards town centre
CORE_RADIUS     = 3000     # metres — SDZs intersecting this circle → DZs become core
SDZ_ZONE_RADIUS = 10000    # metres — DEAs intersecting this circle → broken into SDZs

DZ_BOUNDARY_FILE  = "simulation/dz2021/DZ2021.geojson"
SDZ_BOUNDARY_FILE = "simulation/sdz2021/SDZ2021.geojson"
DEA_BOUNDARY_FILE = "simulation/dea2021/DEA2021.geojson"
WORKPLACE_FILE    = "data/census-2021-apwp001.xlsx"
POPULATION_CACHE  = "data/cache_nisra_population.csv"
OUTPUT_FILE       = "data/census_zones.json"

POPULATION_API = (
    "https://ws-data.nisra.gov.uk/public/api.restful/"
    "PxStat.Data.Cube_API.ReadDataset/MYE01T011/CSV/1.0/en/"
)

# ── Check required files ───────────────────────────────────────────────────────

missing = []
for f in [DZ_BOUNDARY_FILE, SDZ_BOUNDARY_FILE, DEA_BOUNDARY_FILE, WORKPLACE_FILE]:
    if not os.path.exists(f):
        missing.append(f)
if missing:
    print("ERROR: Missing required data files:")
    for f in missing:
        print(f"  {f}")
    print("\nDownload from:")
    print("  DZ2021.geojson  — https://www.nisra.gov.uk/support/geography/data-zones-census-2021")
    print("  SDZ2021.geojson — https://www.nisra.gov.uk/support/geography/super-data-zones-census-2021")
    print("  DEA2021.geojson — https://admin.opendatani.gov.uk/dataset/osni-open-data-largescale-boundaries-district-electoral-areas-2012")
    print("  APWP001.xlsx    — https://www.nisra.gov.uk/statistics/census/2021-census/workplace-population")
    sys.exit(1)

# ── Coordinate systems ─────────────────────────────────────────────────────────

to_utm = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
to_wgs = pyproj.Transformer.from_crs("EPSG:32630", "EPSG:4326", always_xy=True)

centre_utm_x, centre_utm_y = to_utm.transform(CENTRE[1], CENTRE[0])
core_circle = Point(centre_utm_x, centre_utm_y).buffer(CORE_RADIUS)
sdz_circle  = Point(centre_utm_x, centre_utm_y).buffer(SDZ_ZONE_RADIUS)

# ── Load census population data ────────────────────────────────────────────────

print("Loading DZ population …")
if os.path.exists(POPULATION_CACHE):
    pop_df = pd.read_csv(POPULATION_CACHE)
else:
    print("  Fetching from NISRA API …")
    req = urllib.request.Request(POPULATION_API, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        pop_csv = r.read().decode("utf-8-sig")
    from io import StringIO
    pop_df = pd.read_csv(StringIO(pop_csv))
    pop_df.to_csv(POPULATION_CACHE, index=False)

pop_df = pop_df[
    (pop_df["TLIST(A1)"] == 2021) &
    (pop_df["DZ2021"].str.startswith("N20"))
][["DZ2021", "VALUE"]].rename(columns={"DZ2021": "DZ2021_cd", "VALUE": "population"})
pop_lookup = pop_df.set_index("DZ2021_cd")["population"].to_dict()
print(f"  {len(pop_lookup)} DZs, NI total pop: {sum(pop_lookup.values()):,}")

print("Loading DZ workplace population …")
wp_df = pd.read_excel(WORKPLACE_FILE, sheet_name="DZ", header=5)
wp_df = wp_df[["Geography Code", "Workplace population"]].rename(
    columns={"Geography Code": "DZ2021_cd", "Workplace population": "workplace_pop"}
)
wp_df = wp_df[wp_df["DZ2021_cd"].astype(str).str.startswith("N20")].copy()
wp_df["workplace_pop"] = pd.to_numeric(wp_df["workplace_pop"], errors="coerce").fillna(0)
wp_lookup = wp_df.set_index("DZ2021_cd")["workplace_pop"].to_dict()
print(f"  {len(wp_lookup)} DZs with workplace data, NI total: {int(sum(wp_lookup.values())):,}")

# ── Load boundary files ────────────────────────────────────────────────────────

print("Loading DZ boundaries …")
dz = gpd.read_file(DZ_BOUNDARY_FILE).to_crs("EPSG:32630")
dz["population"]   = dz["DZ2021_cd"].map(pop_lookup).fillna(0)
dz["workplace_pop"] = dz["DZ2021_cd"].map(wp_lookup).fillna(0)
dz["centroid_utm"]  = dz.geometry.centroid
print(f"  {len(dz)} DZs loaded")

print("Loading SDZ boundaries …")
sdz = gpd.read_file(SDZ_BOUNDARY_FILE).to_crs("EPSG:32630")
print(f"  {len(sdz)} SDZs loaded")

print("Loading DEA boundaries …")
dea = gpd.read_file(DEA_BOUNDARY_FILE).to_crs("EPSG:32630")
print(f"  {len(dea)} DEAs loaded")

# ── Detect SDZ code column in DZ file (for hierarchy lookup) ──────────────────
# NISRA files typically contain a parent-area code column.

def _find_col(df, candidates):
    """Return the first column name from candidates that exists in df.columns."""
    for c in candidates:
        if c in df.columns:
            return c
    return None

dz_sdz_col  = _find_col(dz,  ["SDZ2021_cd", "SDZ_cd", "SDZ2021", "SDZCODE"])
dea_code_col = _find_col(dea, ["DEA2021_cd", "DEA_cd",  "DEA2021", "DEACODE"])
sdz_code_col = _find_col(sdz, ["SDZ2021_cd", "SDZ_cd",  "SDZ2021", "SDZCODE"])

if dz_sdz_col:
    print(f"  DZ→SDZ hierarchy via column '{dz_sdz_col}'")
else:
    print("  DZ→SDZ hierarchy via spatial join (no parent code column found)")

if dea_code_col is None:
    # Try to infer a unique code column
    candidate = [c for c in dea.columns if "code" in c.lower() or "cd" in c.lower()]
    dea_code_col = candidate[0] if candidate else dea.columns[0]
    print(f"  DEA code column inferred: '{dea_code_col}'")

if sdz_code_col is None:
    candidate = [c for c in sdz.columns if "code" in c.lower() or "cd" in c.lower()]
    sdz_code_col = candidate[0] if candidate else sdz.columns[0]
    print(f"  SDZ code column inferred: '{sdz_code_col}'")

# ── Core DZ classification (DZ-level, not SDZ-level) ─────────────────────────
# Using DZ polygon intersection avoids pulling in whole SDZs whose boundary
# merely clips the core circle.

dz["is_core"]    = dz.geometry.intersects(core_circle)
core_dz_codes    = set(dz.loc[dz["is_core"], "DZ2021_cd"])

sdz["in_sdz_zone"] = sdz.geometry.intersects(sdz_circle)
n_sdz_zone = sdz["in_sdz_zone"].sum()
print(f"\nDZ/SDZ classification:")
print(f"  {dz['is_core'].sum()} DZs directly intersect CORE_RADIUS ({CORE_RADIUS}m) → core area")
print(f"  {n_sdz_zone} SDZs intersect SDZ_ZONE_RADIUS ({SDZ_ZONE_RADIUS}m)")

# ── Determine which DEAs intersect the SDZ zone ──────────────────────────────

dea["in_sdz_zone"] = dea.geometry.intersects(sdz_circle)
n_dea_broken = dea["in_sdz_zone"].sum()
n_dea_single = (~dea["in_sdz_zone"]).sum()
print(f"\nDEA classification:")
print(f"  {n_dea_broken} DEAs intersect SDZ_ZONE_RADIUS → broken into SDZs")
print(f"  {n_dea_single} DEAs remain as single external nodes")

# ── Map DZs to their parent SDZ ───────────────────────────────────────────────

if dz_sdz_col:
    dz_to_sdz = dz.set_index("DZ2021_cd")[dz_sdz_col].to_dict()
else:
    # Spatial join: assign each DZ to the SDZ whose centroid it falls in
    dz_centroids = dz[["DZ2021_cd", "centroid_utm"]].copy()
    dz_centroids = dz_centroids.rename(columns={"centroid_utm": "geometry"}).set_geometry("geometry")
    joined = gpd.sjoin(dz_centroids, sdz[[sdz_code_col, "geometry"]], how="left", predicate="within")
    dz_to_sdz = joined.set_index("DZ2021_cd")[sdz_code_col].to_dict()

# Map DZs to their parent DEA via SDZ
sdz_to_dea = {}
if "DEA2021_cd" in sdz.columns or dea_code_col in sdz.columns:
    _dea_col_in_sdz = _find_col(sdz, ["DEA2021_cd", "DEA_cd", dea_code_col])
    if _dea_col_in_sdz:
        sdz_to_dea = sdz.set_index(sdz_code_col)[_dea_col_in_sdz].to_dict()

if not sdz_to_dea:
    # Spatial join: assign each SDZ to the DEA it falls in
    sdz_centroids = sdz[[sdz_code_col, "geometry"]].copy()
    sdz_centroids["geometry"] = sdz_centroids.geometry.centroid
    joined2 = gpd.sjoin(sdz_centroids, dea[[dea_code_col, "geometry"]], how="left", predicate="within")
    sdz_to_dea = joined2.set_index(sdz_code_col)[dea_code_col].to_dict()

# Which SDZs have at least one core DZ? (used for external node classification)
_sdz_has_core = {}
for dz_cd, sdz_cd in dz_to_sdz.items():
    if dz_cd in core_dz_codes:
        _sdz_has_core[sdz_cd] = True
sdz["has_core_dz"] = sdz[sdz_code_col].map(_sdz_has_core).fillna(False)

# ── Build core area polygon ────────────────────────────────────────────────────

core_dz_gdf  = dz[dz["is_core"]]
core_polygon  = unary_union(core_dz_gdf.geometry.values)

n_core_sdzs = len({dz_to_sdz.get(cd) for cd in core_dz_codes} - {None})
print(f"\nCore area: {len(core_dz_gdf)} DZs (from {n_core_sdzs} partially/fully-core SDZs)")
core_area_km2 = core_polygon.area / 1e6
print(f"  Core polygon area: {core_area_km2:.2f} km²")

# Max distance from centre to any core polygon vertex — determines minimum RADIUS_M
centre_pt_utm   = Point(centre_utm_x, centre_utm_y)
max_vertex_dist = max(centre_pt_utm.distance(Point(x, y))
                      for x, y in core_polygon.exterior.coords)
print(f"  Max core polygon vertex distance: {max_vertex_dist:.0f}m")

# ── Determine SDZ zone DEA codes ──────────────────────────────────────────────

# DEAs that intersect SDZ_ZONE_RADIUS are broken into SDZs
broken_dea_codes = set(dea.loc[dea["in_sdz_zone"], dea_code_col].tolist())

# ── Compute population-weighted centroids ─────────────────────────────────────

def weighted_centroid_wgs(geometries, weights):
    """Return (lat, lon) population-weighted centroid."""
    weights = np.array(weights, dtype=float)
    if weights.sum() == 0:
        weights = np.ones(len(weights))
    weights /= weights.sum()
    centroids = [g.centroid for g in geometries]
    cx = sum(c.x * w for c, w in zip(centroids, weights))
    cy = sum(c.y * w for c, w in zip(centroids, weights))
    lon, lat = to_wgs.transform(cx, cy)
    return lat, lon, cx, cy

# ── Build external node list ──────────────────────────────────────────────────

external_nodes = []

print("\nBuilding external node list …")

# SDZ external nodes: SDZs in broken DEAs with NO core DZs (fully external)
sdz["parent_dea"] = sdz[sdz_code_col].map(sdz_to_dea)
sdz_in_broken_dea = sdz["parent_dea"].isin(broken_dea_codes)

sdz_external_mask = sdz_in_broken_dea & ~sdz["has_core_dz"]
sdz_external = sdz[sdz_external_mask].copy()

for _, row in sdz_external.iterrows():
    sdz_code = row[sdz_code_col]
    child_dz_codes = [dz_cd for dz_cd, parent in dz_to_sdz.items() if parent == sdz_code]
    child_dz = dz[dz["DZ2021_cd"].isin(child_dz_codes)]
    pop   = child_dz["population"].sum()
    wp    = child_dz["workplace_pop"].sum()
    lat, lon, utm_x, utm_y = weighted_centroid_wgs(child_dz.geometry.values,
                                                     child_dz["population"].values)
    external_nodes.append({
        "id":             sdz_code,
        "level":          "SDZ",
        "centroid_lat":   round(lat, 6),
        "centroid_lon":   round(lon, 6),
        "centroid_utm_x": round(utm_x, 1),
        "centroid_utm_y": round(utm_y, 1),
        "population":     int(round(pop)),
        "workplace_pop":  int(round(wp)),
    })
print(f"  {len(external_nodes)} SDZ external nodes")

# DZ external nodes: non-core DZs from partially-core SDZs
# (their SDZ has some core DZs but these DZs themselves don't intersect the core circle)
partial_core_sdz_codes = set(sdz.loc[sdz_in_broken_dea & sdz["has_core_dz"], sdz_code_col])
orphan_dz_mask = (
    dz["DZ2021_cd"].map(dz_to_sdz).isin(partial_core_sdz_codes)
    & ~dz["is_core"]
)
orphan_dz = dz[orphan_dz_mask].copy()
n_dz_ext_start = len(external_nodes)

for _, row in orphan_dz.iterrows():
    cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
    lon, lat = to_wgs.transform(cx, cy)
    external_nodes.append({
        "id":             row["DZ2021_cd"],
        "level":          "DZ",
        "centroid_lat":   round(lat, 6),
        "centroid_lon":   round(lon, 6),
        "centroid_utm_x": round(cx, 1),
        "centroid_utm_y": round(cy, 1),
        "population":     int(round(row["population"])),
        "workplace_pop":  int(round(row["workplace_pop"])),
    })
print(f"  {len(external_nodes) - n_dz_ext_start} DZ external nodes (orphan DZs from partially-core SDZs)")

# DEA external nodes: DEAs NOT intersecting SDZ_ZONE_RADIUS
dea_external = dea[~dea["in_sdz_zone"]].copy()
n_dea_ext_start = len(external_nodes)

for _, row in dea_external.iterrows():
    dea_code = row[dea_code_col]
    # Child SDZs → child DZs
    child_sdz_codes = [s for s, d in sdz_to_dea.items() if d == dea_code]
    child_dz_codes  = [dz_cd for dz_cd, s in dz_to_sdz.items() if s in child_sdz_codes]
    child_dz = dz[dz["DZ2021_cd"].isin(child_dz_codes)]
    pop   = child_dz["population"].sum()
    wp    = child_dz["workplace_pop"].sum()
    if len(child_dz) == 0:
        # Fallback: use DEA polygon centroid
        cx, cy = row.geometry.centroid.x, row.geometry.centroid.y
        lon, lat = to_wgs.transform(cx, cy)
        utm_x, utm_y = cx, cy
        pop, wp = 0, 0
    else:
        lat, lon, utm_x, utm_y = weighted_centroid_wgs(child_dz.geometry.values,
                                                         child_dz["population"].values)
    external_nodes.append({
        "id":             dea_code,
        "level":          "DEA",
        "centroid_lat":   round(lat, 6),
        "centroid_lon":   round(lon, 6),
        "centroid_utm_x": round(utm_x, 1),
        "centroid_utm_y": round(utm_y, 1),
        "population":     int(round(pop)),
        "workplace_pop":  int(round(wp)),
    })
print(f"  {len(external_nodes) - n_dea_ext_start} DEA external nodes")
print(f"  {len(external_nodes)} external nodes total")
print(f"  Total external pop: {sum(n['population'] for n in external_nodes):,}")
print(f"  Total external wp:  {sum(n['workplace_pop'] for n in external_nodes):,}")

# ── Serialise core polygon to WGS84 ──────────────────────────────────────────

# Convert core polygon from UTM back to WGS84 for JSON storage
import shapely.ops as sops
core_polygon_wgs = sops.transform(lambda x, y: to_wgs.transform(x, y), core_polygon)
core_coords = list(mapping(core_polygon_wgs)["coordinates"][0])  # exterior ring

# ── Write output ──────────────────────────────────────────────────────────────

output = {
    "core_radius":             CORE_RADIUS,
    "sdz_zone_radius":         SDZ_ZONE_RADIUS,
    "centre_lat":              CENTRE[0],
    "centre_lon":              CENTRE[1],
    "max_core_vertex_dist_m":  round(max_vertex_dist, 1),
    "core_polygon":            core_coords,
    "n_core_dzs":              int(len(core_dz_gdf)),
    "n_core_sdzs":             n_core_sdzs,
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
