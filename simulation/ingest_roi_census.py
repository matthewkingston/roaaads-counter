"""
Load RoI census data (SA / ED / LEA from CSO) and return standardised
GeoDataFrames consumed by build_census_zones.py.

Public API
----------
load_roi_census() -> (sa_gdf, ed_gdf, lea_gdf)

Standardised columns — same schema as ingest_ni_census.py
    area_code       str   — census code for this area
    parent_code     str   — code of the immediate parent area (NaN where none)
    level           str   — "SA", "ED", or "LEA"
    population      int   — resident population (sa_gdf + aggregates)
    workplace_pop   float — workplace population (see note below)
    geometry        Polygon — projected to PROJECTED_CRS

Workplace note
--------------
No open 2022 per-zone workplace data exists for RoI (2022 POWSCCAR is restricted
microdata).  Workplace population is loaded from a pre-computed cache file produced
by simulation/build_wz_apportionment.py, which intersects 2016 WZ polygons with 2022
SA polygons geometrically and splits each WZ's headcount by POI density.

If the cache is missing, load_roi_census() raises FileNotFoundError telling the user
to run build_wz_apportionment.py first.
"""

import glob, os, sys
import geopandas as gpd
import pandas as pd

from demographics_config import PROJECTED_CRS

# ── File paths (relative to project root) ────────────────────────────────────
_BASE = "data/ireland_data"
SA_BOUNDARY_GLOB = f"{_BASE}/Small_Area_National_Statistical_Boundaries_2022_*.geojson"
SAPS_CSV         = f"{_BASE}/Complete_set_of_Census_2022_SAPs/SAPS_2022_Small_Area_UR_171024.csv"
SA_WP_CACHE      = f"{_BASE}/cache_sa_workplace.csv"


def load_roi_census():
    """
    Load RoI SA/ED/LEA boundaries + CSO 2022 population + pre-computed WZ workplace.
    Returns (sa_gdf, ed_gdf, lea_gdf) with standardised columns.
    """
    # ── SA boundaries ────────────────────────────────────────────────────────
    sa_files = glob.glob(SA_BOUNDARY_GLOB)
    if not sa_files:
        raise FileNotFoundError(f"SA boundary file not found matching: {SA_BOUNDARY_GLOB}")
    sa_path = sa_files[0]
    print(f"RoI: loading SA boundaries ({os.path.basename(sa_path)}) — ~410 MB, takes a moment …")
    sa = gpd.read_file(sa_path).to_crs(PROJECTED_CRS)
    print(f"  {len(sa)} SAs loaded")

    # ── Population (SAPS 2022) ───────────────────────────────────────────────
    # (SAPS T1_1AGETT *column* sums to ~2× the national population because the file
    # carries a "State" aggregate row equal to the sum of all SAs; that row has no
    # SA_PUB2022 match, so the per-SA join below yields the correct 1× population.)
    print("RoI: loading SA population (SAPS 2022) …")
    saps = pd.read_csv(SAPS_CSV, usecols=["GEOGID", "T1_1AGETT"])
    # GEOGID is int64 without leading zero; SA_PUB2022 is a 9-digit string.
    saps["sa_code"] = saps["GEOGID"].astype(str).str.zfill(9)
    pop_lookup = saps.set_index("sa_code")["T1_1AGETT"].to_dict()
    sa["population"] = sa["SA_PUB2022"].map(pop_lookup).fillna(0).astype(int)
    n_matched = sa["SA_PUB2022"].isin(pop_lookup).sum()
    print(f"  {n_matched}/{len(sa)} SAs matched to SAPS; RoI total pop: {sa['population'].sum():,}")

    # ── Workplace (pre-computed WZ→SA apportionment) ──────────────────────────
    print("RoI: loading WZ workplace cache …")
    if not os.path.exists(SA_WP_CACHE):
        raise FileNotFoundError(
            f"Missing {SA_WP_CACHE}\n"
            "Run: python3 simulation/build_wz_apportionment.py"
        )
    wp_df = pd.read_csv(SA_WP_CACHE)
    wp_lookup = wp_df.set_index("sa_code")["workplace_pop"].to_dict()
    sa["workplace_pop"] = sa["SA_PUB2022"].map(wp_lookup).fillna(0)
    print(f"  RoI total workplace_pop: {sa['workplace_pop'].sum():,.0f}")

    # ── Build standardised SA GeoDataFrame ──────────────────────────────────
    sa_gdf = gpd.GeoDataFrame({
        "area_code":        sa["SA_PUB2022"].astype(str).values,
        "parent_code":      sa["ED_ID_STR"].astype(str).values,   # ED code
        "grandparent_code": sa["CSO_LEA"].astype(str).values,     # LEA name
        "level":            "SA",
        "population":       sa["population"].values,
        "workplace_pop":    sa["workplace_pop"].values,
        "geometry":         sa["geometry"].values,
    }, crs=PROJECTED_CRS)

    # ── Derive ED GeoDataFrame by dissolving SAs ─────────────────────────────
    print("RoI: deriving ED boundaries by dissolving SAs …")
    ed_raw = (
        sa_gdf
        .dissolve(
            by="parent_code",
            aggfunc={"population": "sum", "workplace_pop": "sum",
                     "grandparent_code": "first"},
        )
        .reset_index()
    )
    ed_gdf = gpd.GeoDataFrame({
        "area_code":     ed_raw["parent_code"].values,
        "parent_code":   ed_raw["grandparent_code"].values,
        "level":         "ED",
        "population":    ed_raw["population"].values,
        "workplace_pop": ed_raw["workplace_pop"].values,
        "geometry":      ed_raw["geometry"].values,
    }, crs=PROJECTED_CRS)
    print(f"  {len(ed_gdf)} EDs derived")

    # ── Derive LEA GeoDataFrame by dissolving SAs ────────────────────────────
    print("RoI: deriving LEA boundaries by dissolving SAs …")
    lea_raw = (
        sa_gdf
        .dissolve(
            by="grandparent_code",
            aggfunc={"population": "sum", "workplace_pop": "sum"},
        )
        .reset_index()
    )
    lea_gdf = gpd.GeoDataFrame({
        "area_code":   lea_raw["grandparent_code"].values,
        "parent_code": pd.Series([pd.NA] * len(lea_raw), dtype=object),
        "level":       "LEA",
        "population":  lea_raw["population"].values,
        "workplace_pop": lea_raw["workplace_pop"].values,
        "geometry":    lea_raw["geometry"].values,
    }, crs=PROJECTED_CRS)
    print(f"  {len(lea_gdf)} LEAs derived")

    return sa_gdf, ed_gdf, lea_gdf
