"""
Load NI census data (DZ / SDZ / DEA) from NISRA sources and return standardised
GeoDataFrames consumed by build_census_zones.py.

Public API
----------
load_ni_census() -> (dz_gdf, sdz_gdf, dea_gdf)

Standardised columns on every returned GeoDataFrame
    area_code       str  — census code for this area
    parent_code     str  — code of the immediate parent area (NaN where none)
    level           str  — "DZ", "SDZ", or "DEA"
    population      int  — resident population (dz_gdf only; 0 elsewhere)
    workplace_pop   float— workplace population (dz_gdf only; 0 elsewhere)
    geometry        Polygon — projected to PROJECTED_CRS
"""

import os, sys, urllib.request
from io import StringIO

import geopandas as gpd
import pandas as pd

from demographics_config import PROJECTED_CRS

# ── File paths (relative to project root) ────────────────────────────────────
DZ_BOUNDARY_FILE  = "simulation/dz2021/DZ2021.geojson"
SDZ_BOUNDARY_FILE = "simulation/sdz2021/SDZ2021.geojson"
DEA_BOUNDARY_FILE = "simulation/dea2021/DEA2021.geojson"
WORKPLACE_FILE    = "data/census-2021-apwp001.xlsx"
POPULATION_CACHE  = "data/cache_nisra_population.csv"

POPULATION_API = (
    "https://ws-data.nisra.gov.uk/public/api.restful/"
    "PxStat.Data.Cube_API.ReadDataset/MYE01T011/CSV/1.0/en/"
)


def _find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load_ni_census():
    """
    Load NI DZ/SDZ/DEA boundaries + NISRA population + workplace.
    Returns (dz_gdf, sdz_gdf, dea_gdf) with standardised columns.
    """
    missing = [f for f in [DZ_BOUNDARY_FILE, SDZ_BOUNDARY_FILE,
                            DEA_BOUNDARY_FILE, WORKPLACE_FILE]
               if not os.path.exists(f)]
    if missing:
        print("ERROR: Missing NI data files:")
        for f in missing:
            print(f"  {f}")
        sys.exit(1)

    # ── Population (NISRA API, cached) ──────────────────────────────────────
    print("NI: loading DZ population …")
    if os.path.exists(POPULATION_CACHE):
        pop_df = pd.read_csv(POPULATION_CACHE)
    else:
        print("  Fetching from NISRA API …")
        req = urllib.request.Request(POPULATION_API, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            raw = r.read().decode("utf-8-sig")
        pop_df = pd.read_csv(StringIO(raw))
        pop_df.to_csv(POPULATION_CACHE, index=False)

    pop_df = pop_df[
        (pop_df["TLIST(A1)"] == 2021) &
        (pop_df["DZ2021"].str.startswith("N20"))
    ][["DZ2021", "VALUE"]].rename(columns={"DZ2021": "DZ2021_cd", "VALUE": "population"})
    pop_lookup = pop_df.set_index("DZ2021_cd")["population"].to_dict()
    print(f"  {len(pop_lookup)} DZs, NI total pop: {sum(pop_lookup.values()):,}")

    # ── Workplace (NISRA Excel) ──────────────────────────────────────────────
    print("NI: loading DZ workplace population …")
    wp_df = pd.read_excel(WORKPLACE_FILE, sheet_name="DZ", header=5)
    wp_df = wp_df[["Geography Code", "Workplace population"]].rename(
        columns={"Geography Code": "DZ2021_cd", "Workplace population": "workplace_pop"}
    )
    wp_df = wp_df[wp_df["DZ2021_cd"].astype(str).str.startswith("N20")].copy()
    wp_df["workplace_pop"] = pd.to_numeric(wp_df["workplace_pop"], errors="coerce").fillna(0)
    wp_lookup = wp_df.set_index("DZ2021_cd")["workplace_pop"].to_dict()
    print(f"  {len(wp_lookup)} DZs with workplace data, NI total: {int(sum(wp_lookup.values())):,}")

    # ── Boundary files ───────────────────────────────────────────────────────
    print("NI: loading DZ boundaries …")
    dz = gpd.read_file(DZ_BOUNDARY_FILE).to_crs(PROJECTED_CRS)
    dz["population"]    = dz["DZ2021_cd"].map(pop_lookup).fillna(0)
    dz["workplace_pop"] = dz["DZ2021_cd"].map(wp_lookup).fillna(0)
    print(f"  {len(dz)} DZs loaded")

    print("NI: loading SDZ boundaries …")
    sdz = gpd.read_file(SDZ_BOUNDARY_FILE).to_crs(PROJECTED_CRS)
    print(f"  {len(sdz)} SDZs loaded")

    print("NI: loading DEA boundaries …")
    dea = gpd.read_file(DEA_BOUNDARY_FILE).to_crs(PROJECTED_CRS)
    print(f"  {len(dea)} DEAs loaded")

    # ── Detect code columns ──────────────────────────────────────────────────
    dz_sdz_col   = _find_col(dz,  ["SDZ2021_cd", "SDZ_cd", "SDZ2021", "SDZCODE"])
    sdz_code_col = _find_col(sdz, ["SDZ2021_cd", "SDZ_cd", "SDZ2021", "SDZCODE"])
    sdz_dea_col  = _find_col(sdz, ["DEA2021_cd", "DEA_cd", "DEA2021", "DEACODE"])
    dea_code_col = _find_col(dea, ["DEA2021_cd", "DEA_cd", "DEA2021", "DEACODE"])

    if sdz_code_col is None:
        cands = [c for c in sdz.columns if "code" in c.lower() or "cd" in c.lower()]
        sdz_code_col = cands[0] if cands else sdz.columns[0]
        print(f"  SDZ code column inferred: '{sdz_code_col}'")
    if dea_code_col is None:
        cands = [c for c in dea.columns if "code" in c.lower() or "cd" in c.lower()]
        dea_code_col = cands[0] if cands else dea.columns[0]
        print(f"  DEA code column inferred: '{dea_code_col}'")

    # ── DZ → SDZ parent lookup ───────────────────────────────────────────────
    if dz_sdz_col:
        dz_parent = dz[dz_sdz_col].astype(str)
    else:
        print("  DZ→SDZ: spatial join (no parent-code column in DZ file)")
        cents = dz[["DZ2021_cd", "geometry"]].copy()
        cents["geometry"] = cents.geometry.centroid
        joined = gpd.sjoin(cents, sdz[[sdz_code_col, "geometry"]],
                           how="left", predicate="within")
        dz_parent = joined.set_index("DZ2021_cd")[sdz_code_col].reindex(
            dz["DZ2021_cd"]).values

    # ── SDZ → DEA parent lookup ──────────────────────────────────────────────
    if sdz_dea_col:
        sdz_parent = sdz[sdz_dea_col].astype(str)
    else:
        print("  SDZ→DEA: spatial join (no parent-code column in SDZ file)")
        cents2 = sdz[[sdz_code_col, "geometry"]].copy()
        cents2["geometry"] = cents2.geometry.centroid
        joined2 = gpd.sjoin(cents2, dea[[dea_code_col, "geometry"]],
                            how="left", predicate="within")
        sdz_parent = joined2.set_index(sdz_code_col)[dea_code_col].reindex(
            sdz[sdz_code_col]).values

    # ── Standardise ──────────────────────────────────────────────────────────
    dz_gdf = gpd.GeoDataFrame({
        "area_code":     dz["DZ2021_cd"].astype(str).values,
        "parent_code":   pd.array(dz_parent, dtype=object),
        "level":         "DZ",
        "population":    dz["population"].values,
        "workplace_pop": dz["workplace_pop"].values,
        "geometry":      dz["geometry"].values,
    }, crs=PROJECTED_CRS)

    sdz_gdf = gpd.GeoDataFrame({
        "area_code":   sdz[sdz_code_col].astype(str).values,
        "parent_code": pd.array(sdz_parent, dtype=object),
        "level":       "SDZ",
        "geometry":    sdz["geometry"].values,
    }, crs=PROJECTED_CRS)

    dea_gdf = gpd.GeoDataFrame({
        "area_code": dea[dea_code_col].astype(str).values,
        "level":     "DEA",
        "geometry":  dea["geometry"].values,
    }, crs=PROJECTED_CRS)

    return dz_gdf, sdz_gdf, dea_gdf
