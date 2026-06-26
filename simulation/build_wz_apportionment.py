"""
Pre-compute WZ→SA workplace apportionment using POI-weighted geometric intersection.

The 2016 CSO Workplace Zones (WZs) do not align with 2022 SA boundaries.  This
script intersects WZ polygons with 2022 SA polygons geometrically — bypassing any
2016→2022 SA change-code tracking — and splits each WZ's headcount across the SAs
it overlaps, weighted by the count of commercial POIs in each intersection piece.
Area-proportional fallback is used for pieces that contain no POIs.

Output
------
data/ireland_data/cache_sa_workplace.csv   (sa_code, workplace_pop)
    One row per 2022 SA; workplace_pop is the sum of apportioned WZ headcounts.
    Loaded by ingest_roi_census.load_roi_census() at every pipeline run.

Run once.  Re-run only when:
  - WZ boundaries change  (Workplace_Zones_ITM.shp)
  - SA boundaries change  (Small_Area_National_Statistical_Boundaries_2022_*.geojson)
  - WZ SAPS headcounts change  (T1_T column already in the shapefile)
  - OSM POI data is significantly stale

Caches
------
data/ireland_data/cache_roi_pois.geojson   — POI points (ITM) with poi_weight column.
    Re-used across runs; delete to force re-extraction from PBF.

Needs: Docker (osmctools-roaaads image), local PBF (demographics_config.PBF_PATH).
"""

import glob, os, subprocess, sys
import xml.etree.ElementTree as ET

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import Point

from demographics_config import PROJECTED_CRS, EXCLUDE_AMENITY, POI_WEIGHTS

# ── Paths ────────────────────────────────────────────────────────────────────

_BASE = "data/ireland_data"
WZ_FILE          = f"{_BASE}/Workplace_Zones_ITM/Workplace_Zones_ITM.shp"
SA_BOUNDARY_GLOB = f"{_BASE}/Small_Area_National_Statistical_Boundaries_2022_*.geojson"
SA_WP_OUT        = f"{_BASE}/cache_sa_workplace.csv"
POI_CACHE        = f"{_BASE}/cache_roi_pois.geojson"

OSRM_ROOT        = "/home/matthew/Documents/CodingFun/osrm"
PBF_NAME         = "ireland-and-northern-ireland-latest.osm.pbf"
OSMCTOOLS_IMAGE  = "osmctools-roaaads"

# Intermediate osmctools scratch (prefer reusing edge_index o5m if present)
_WORK_DIR        = os.path.join(OSRM_ROOT, "wz_pois")

SLIVER_AREA_M2   = 100     # intersection pieces smaller than this are dropped

# ── Coordinate transformer (WGS84 → ITM) ─────────────────────────────────────

_to_itm = Transformer.from_crs("EPSG:4326", PROJECTED_CRS, always_xy=True)


# ── osmctools POI extraction ──────────────────────────────────────────────────

def _docker(cmd):
    uidgid = f"{os.getuid()}:{os.getgid()}"
    subprocess.run(
        ["docker", "run", "--rm", "--user", uidgid,
         "-v", f"{OSRM_ROOT}:/osrm",
         OSMCTOOLS_IMAGE, "sh", "-c", cmd],
        check=True,
    )


def _extract_poi_osm():
    """PBF → filtered POI node .osm via Docker.  Returns host path of the .osm."""
    os.makedirs(_WORK_DIR, exist_ok=True)

    # Prefer existing o5m from edge_index pipeline (same PBF, already converted)
    edge_o5m_host   = os.path.join(OSRM_ROOT, "edge_index", "ni.o5m")
    wz_o5m_host     = os.path.join(_WORK_DIR, "ni.o5m")
    pois_osm_host   = os.path.join(_WORK_DIR, "pois.osm")

    if os.path.exists(pois_osm_host):
        print(f"  Reusing existing {pois_osm_host} "
              f"({os.path.getsize(pois_osm_host)/1e6:.0f} MB)")
        return pois_osm_host

    if os.path.exists(edge_o5m_host):
        print(f"  Reusing edge_index o5m …")
        o5m_docker = "/osrm/edge_index/ni.o5m"
    elif os.path.exists(wz_o5m_host):
        print(f"  Reusing existing wz_pois o5m …")
        o5m_docker = "/osrm/wz_pois/ni.o5m"
    else:
        print("  Converting pbf → o5m (osmconvert, streaming) …")
        _docker(
            f"osmconvert /osrm/{PBF_NAME} "
            f"-t=/osrm/wz_pois/_osmconvert_tmp -o=/osrm/wz_pois/ni.o5m"
        )
        o5m_docker = "/osrm/wz_pois/ni.o5m"

    print("  Filtering amenity/shop/office nodes (osmfilter) …")
    _docker(
        f"osmfilter {o5m_docker} -t=/osrm/wz_pois/_osmfilter_tmp "
        f'--keep-nodes="amenity= shop= office=" '
        f"-o=/osrm/wz_pois/pois.osm"
    )
    print(f"  POI osm: {os.path.getsize(pois_osm_host)/1e6:.0f} MB")
    return pois_osm_host


def _parse_poi_nodes(osm_path):
    """Stream <node> elements from the POI osm and return a GeoDataFrame in ITM."""
    records = []
    context = ET.iterparse(osm_path, events=("start", "end"))
    _event, root = next(context)

    for event, elem in context:
        if event != "end":
            continue
        if elem.tag == "node":
            lat_s = elem.get("lat")
            lon_s = elem.get("lon")
            if lat_s is None or lon_s is None:
                root.clear()
                continue
            tags = {t.get("k"): t.get("v") for t in elem.findall("tag")}
            amenity = tags.get("amenity")
            shop    = tags.get("shop")
            office  = tags.get("office")

            if amenity and amenity in EXCLUDE_AMENITY:
                root.clear()
                continue

            # Weight: amenity > shop > office > generic fallback
            if amenity and amenity in POI_WEIGHTS:
                w = POI_WEIGHTS[amenity]
            elif amenity and amenity not in EXCLUDE_AMENITY:
                w = POI_WEIGHTS.get(amenity, 1.0)
            elif shop and shop in POI_WEIGHTS:
                w = POI_WEIGHTS[shop]
            elif shop:
                w = 1.0
            elif office:
                w = POI_WEIGHTS.get(office, 2.0)
            else:
                root.clear()
                continue

            x, y = _to_itm.transform(float(lon_s), float(lat_s))
            records.append((x, y, w))

        if elem.tag in ("node", "way", "relation"):
            root.clear()

    print(f"  {len(records):,} POI nodes parsed")
    xs, ys, ws = zip(*records) if records else ([], [], [])
    return gpd.GeoDataFrame(
        {"poi_weight": list(ws)},
        geometry=[Point(x, y) for x, y in zip(xs, ys)],
        crs=PROJECTED_CRS,
    )


def _load_pois():
    """Load POIs from cache (GeoJSON) or extract fresh from PBF."""
    if os.path.exists(POI_CACHE):
        print(f"Loading POI cache ({POI_CACHE}) …")
        poi_gdf = gpd.read_file(POI_CACHE)
        print(f"  {len(poi_gdf):,} POIs")
        return poi_gdf

    print("Extracting POIs from PBF …")
    osm_path = _extract_poi_osm()
    poi_gdf  = _parse_poi_nodes(osm_path)

    print(f"Saving POI cache → {POI_CACHE}")
    poi_gdf.to_file(POI_CACHE, driver="GeoJSON")
    return poi_gdf


# ── Main apportionment ────────────────────────────────────────────────────────

def main():
    if os.path.exists(SA_WP_OUT):
        print(f"Output already exists: {SA_WP_OUT}")
        print("Delete it to recompute.")
        return

    # ── Load WZ shapefile (already in EPSG:2157) ─────────────────────────────
    print(f"Loading WZ boundaries ({WZ_FILE}) …")
    wz = gpd.read_file(WZ_FILE)[["WORKPLACE", "T1_T", "geometry"]]
    wz = wz.to_crs(PROJECTED_CRS)   # already ITM; this is a no-op but defensive
    wz_total = wz["T1_T"].sum()
    print(f"  {len(wz):,} WZs  |  total T1_T: {wz_total:,.0f}")

    # ── Load SA boundaries (only the code + geometry needed for overlay) ──────
    sa_files = glob.glob(SA_BOUNDARY_GLOB)
    if not sa_files:
        sys.exit(f"ERROR: SA boundary not found matching {SA_BOUNDARY_GLOB}")
    print(f"Loading SA boundaries ({os.path.basename(sa_files[0])}) — large file …")
    sa_full = gpd.read_file(sa_files[0])
    sa = sa_full[["SA_PUB2022", "geometry"]].copy().to_crs(PROJECTED_CRS)
    print(f"  {len(sa):,} SAs")
    del sa_full   # free memory before overlay

    # ── Geometric intersection ────────────────────────────────────────────────
    print("Computing WZ × SA geometric overlay … (may take a few minutes)")
    pieces = gpd.overlay(wz, sa, how="intersection", keep_geom_type=False)
    n_before = len(pieces)
    pieces = pieces[pieces.geometry.area > SLIVER_AREA_M2].copy()
    print(f"  {len(pieces):,} intersection pieces  "
          f"({n_before - len(pieces)} slivers dropped < {SLIVER_AREA_M2} m²)")

    pieces["piece_area"] = pieces.geometry.area

    # ── Load / extract POIs ───────────────────────────────────────────────────
    poi_gdf = _load_pois()

    # ── Spatial join: assign each POI to the intersection piece it falls in ──
    print("Spatial join: POIs → intersection pieces …")
    poi_pts = poi_gdf[["poi_weight", "geometry"]].copy()
    joined  = gpd.sjoin(poi_pts, pieces[["WORKPLACE", "SA_PUB2022", "geometry"]],
                        predicate="within", how="inner")
    piece_poi_weight = (
        joined.groupby(["WORKPLACE", "SA_PUB2022"])["poi_weight"].sum()
    )
    pieces = pieces.join(
        piece_poi_weight, on=["WORKPLACE", "SA_PUB2022"]
    )
    pieces["poi_weight"] = pieces["poi_weight"].fillna(0.0)

    n_poi_pieces = (pieces["poi_weight"] > 0).sum()
    n_area_fallback = len(pieces) - n_poi_pieces
    print(f"  {n_poi_pieces:,} pieces with ≥1 POI  |  "
          f"{n_area_fallback:,} using area-proportional fallback")

    # ── Per-WZ totals ─────────────────────────────────────────────────────────
    wz_poi_total  = pieces.groupby("WORKPLACE")["poi_weight"].sum()
    wz_area_total = pieces.groupby("WORKPLACE")["piece_area"].sum()
    pieces["wz_poi_total"]  = pieces["WORKPLACE"].map(wz_poi_total)
    pieces["wz_area_total"] = pieces["WORKPLACE"].map(wz_area_total)

    # ── Split weight: POI-weighted or area-proportional ───────────────────────
    has_pois = pieces["wz_poi_total"] > 0
    pieces["split_w"] = np.where(
        has_pois,
        pieces["poi_weight"]   / pieces["wz_poi_total"],
        pieces["piece_area"]   / pieces["wz_area_total"],
    )

    # ── Apportion and aggregate ────────────────────────────────────────────────
    pieces["sa_wp_contrib"] = pieces["T1_T"] * pieces["split_w"]
    sa_workplace = (
        pieces.groupby("SA_PUB2022")["sa_wp_contrib"]
        .sum()
        .rename("workplace_pop")
        .reset_index(names="sa_code")
    )

    apportioned_total = sa_workplace["workplace_pop"].sum()
    print(f"\nConservation check:")
    print(f"  WZ total T1_T:           {wz_total:>12,.0f}")
    print(f"  SA apportioned total:     {apportioned_total:>12,.0f}")
    diff = abs(wz_total - apportioned_total)
    if diff > 1.0:
        print(f"  WARNING: discrepancy of {diff:.1f} "
              f"(sliver-dropped pieces lose a small amount — expected if non-zero)")

    # ── Write ─────────────────────────────────────────────────────────────────
    sa_workplace.to_csv(SA_WP_OUT, index=False)
    print(f"\nWrote {SA_WP_OUT}  ({len(sa_workplace):,} SAs)")
    print("Next: python3 simulation/build_census_zones.py")


if __name__ == "__main__":
    main()
