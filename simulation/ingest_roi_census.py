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
microdata).  This module uses the 2016 Workplace Zone (WZ) file with a
population-proportional apportionment to SAs:

    SA_workplace = Σ_WZ  WZ_T1T × SA_pop / WZ_pop

where WZ_pop = Σ_{s∈WZ} SA_pop_s.  This is a first-stage approximation;
POI-weighted split is a planned upgrade (recorded in project_roi_extension memory).
"""

import glob, os
import geopandas as gpd
import numpy as np
import pandas as pd

from demographics_config import PROJECTED_CRS

# ── File paths (relative to project root) ────────────────────────────────────
_BASE = "data/ireland_data"
SA_BOUNDARY_GLOB = f"{_BASE}/Small_Area_National_Statistical_Boundaries_2022_*.geojson"
SAPS_CSV         = f"{_BASE}/Complete_set_of_Census_2022_SAPs/SAPS_2022_Small_Area_UR_171024.csv"
WZ_SAPS_FILE     = f"{_BASE}/Workplace_zones_-_SAPS_2016.xlsx"
WZ_LOOKUP_FILE   = f"{_BASE}/Look-up_table_-_Small_Area_to_Workplace_Zone.xlsx"


def load_roi_census():
    """
    Load RoI SA/ED/LEA boundaries + CSO 2022 population + 2016 WZ workplace.
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
    print("RoI: loading SA population (SAPS 2022) …")
    saps = pd.read_csv(SAPS_CSV, usecols=["GEOGID", "T1_1AGETT"])
    # GEOGID is int64 without leading zero; SA_PUB2022 is a 9-digit string.
    saps["sa_code"] = saps["GEOGID"].astype(str).str.zfill(9)
    pop_lookup = saps.set_index("sa_code")["T1_1AGETT"].to_dict()
    sa["population"] = sa["SA_PUB2022"].map(pop_lookup).fillna(0).astype(int)
    n_unmatched = sa["SA_PUB2022"].isin(pop_lookup).sum()
    print(f"  {n_unmatched}/{len(sa)} SAs matched to SAPS; RoI total pop: {sa['population'].sum():,}")

    # ── Workplace (2016 WZ → SA, population-proportional) ────────────────────
    print("RoI: loading WZ workplace (2016) and apportioning to SAs …")
    wz_saps   = pd.read_excel(WZ_SAPS_FILE, usecols=["WORKPLACE_ZONE", "T1_T"])
    wz_t1t    = wz_saps.set_index("WORKPLACE_ZONE")["T1_T"].to_dict()

    # header=3: row 3 is the actual column header (County, ..., WZ1, WZ2, ...)
    lookup    = pd.read_excel(WZ_LOOKUP_FILE, sheet_name="SA to WPZ", header=3)
    wz_cols   = [c for c in lookup.columns if str(c).startswith("WZ")]
    if not wz_cols:
        raise ValueError("No WZ* columns found in SA→WPZ lookup — check header row.")

    # Melt to long (sa_code, wz_code) pairs, drop NaN (SA has <N WZs)
    long = (
        lookup[["Small Area"] + wz_cols]
        .melt(id_vars="Small Area", value_name="wz_code")
        .dropna(subset=["wz_code"])[["Small Area", "wz_code"]]
        .copy()
    )
    long.columns = ["sa_code", "wz_code"]
    long["sa_code"] = long["sa_code"].astype(str).str.zfill(9)

    # Attach SA population and WZ headcount
    long["sa_pop"] = long["sa_code"].map(pop_lookup).fillna(0)
    long["wz_t1t"] = long["wz_code"].map(wz_t1t).fillna(0)

    # WZ_pop = sum of SA_pop for all SAs in each WZ (the denominator for the split)
    wz_pop = long.groupby("wz_code")["sa_pop"].sum().rename("wz_pop")
    long   = long.join(wz_pop, on="wz_code")
    # Guard zero-population WZs: uniform split across member SAs
    zero_wz = long["wz_pop"] == 0
    if zero_wz.any():
        wz_sa_count = long.groupby("wz_code")["sa_code"].transform("count")
        long.loc[zero_wz, "sa_pop"] = 1
        long.loc[zero_wz, "wz_pop"] = wz_sa_count[zero_wz]

    long["sa_wk_contrib"] = long["wz_t1t"] * long["sa_pop"] / long["wz_pop"]
    sa_wp = long.groupby("sa_code")["sa_wk_contrib"].sum()
    sa["workplace_pop"] = sa["SA_PUB2022"].map(sa_wp).fillna(0)

    wz_total = sum(wz_t1t.values())
    print(f"  WZ total headcount (raw 2016): {wz_total:,}")
    print(f"  SA apportioned total: {sa['workplace_pop'].sum():,.0f}")

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
    # After dissolve: "parent_code" column = ED code (was the groupby key)
    # "grandparent_code" column = LEA name (first value from child SAs — EDs don't straddle LEAs)
    ed_gdf = gpd.GeoDataFrame({
        "area_code":   ed_raw["parent_code"].values,       # ED code
        "parent_code": ed_raw["grandparent_code"].values,  # LEA name
        "level":       "ED",
        "population":  ed_raw["population"].values,
        "workplace_pop": ed_raw["workplace_pop"].values,
        "geometry":    ed_raw["geometry"].values,
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
        "area_code":   lea_raw["grandparent_code"].values,  # LEA name
        "parent_code": pd.Series([pd.NA] * len(lea_raw), dtype=object),
        "level":       "LEA",
        "population":  lea_raw["population"].values,
        "workplace_pop": lea_raw["workplace_pop"].values,
        "geometry":    lea_raw["geometry"].values,
    }, crs=PROJECTED_CRS)
    print(f"  {len(lea_gdf)} LEAs derived")

    return sa_gdf, ed_gdf, lea_gdf
