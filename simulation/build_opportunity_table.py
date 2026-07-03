"""Island-wide opportunity table — one row per census small area (NI DZ + RoI SA).

Writes `data/island_opportunity_table.csv`: for every small area on the island, its producer
and attractor masses (population, commute producers/attractor, per-level school producers/demand,
retail parking spaces) plus a WGS84 centroid. This is the per-area aggregation that
`build_census_zones.py` already does for the CENTRE's external zones, but run **island-wide over
every small area with no core/external classification** — each small area is the unit, no hierarchy,
no CENTRE. It is the frozen input feed for `analysis/build_n_of_t.py` (the national n(t) sampler).

Uses the SAME estimators as `build_census_zones.py` (so per-area values aggregate to the external-node
values in `data/census_zones.json`):
  * population, workplace_pop  ← ingest_ni_census.load_ni_census / ingest_roi_census.load_roi_census
  * commute_producers          ← census_supply.load_supply (car-driver resident commuters)
  * school_producers_<level>   ← census_school_producers.load_school_producers (per-level)
  * commute_attractor          ← census_attractor.load_attractor (car-commute jobs)
  * retail_spaces              ← parking_demand.parking_spaces over the island parking cache, sjoin to
                                 area polygons + workplace-derived fallback (replicates build_census_zones)
  * school_demand_<level>      ← per-POI enrolment (school_attractor.add_level_enrolments) sjoin to areas

The retail workplace-fallback ratio is computed over SMALL AREAS here (vs external nodes in
build_census_zones): direct-parking areas match census_zones.json when aggregated; the fallback for
no-parking areas differs slightly (different aggregation level) — a verification-only column, not used
by the sampler (which reads the parking POIs directly).

Run: python3 simulation/build_opportunity_table.py   (from the repo root, like build_census_zones.py)
"""
import os

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj

from demographics_config import PROJECTED_CRS, PARKING_ISLAND_CACHE, SCHOOL_ISLAND_CACHE
from parking_demand import parking_spaces
from school_attractor import add_level_enrolments, LEVEL_ENROL_COLS
import census_supply
import census_school_producers
import census_attractor
from ingest_ni_census import load_ni_census
from ingest_roi_census import load_roi_census

OUTPUT_FILE = "data/island_opportunity_table.csv"
LEVELS = ("primary", "postprimary", "tertiary")

# Final column order for the CSV.
COLUMNS = [
    "area_code", "centroid_lat", "centroid_lon",
    "population", "commute_attractor", "retail_spaces",
    "school_demand_primary", "school_demand_postprimary", "school_demand_tertiary",
    "commute_producers",
    "school_producers_primary", "school_producers_postprimary", "school_producers_tertiary",
]


def load_areas():
    """All island small areas (NI DZ + RoI SA) with population/workplace_pop/geometry in ITM."""
    dz_gdf, _, _ = load_ni_census()
    sa_gdf, _, _ = load_roi_census()
    areas = pd.concat([dz_gdf, sa_gdf], ignore_index=True)
    areas = gpd.GeoDataFrame(areas, geometry="geometry", crs=PROJECTED_CRS)
    areas["area_code"] = areas["area_code"].astype(str)
    print(f"Loaded {len(areas):,} small areas (NI DZ + RoI SA)")
    return areas


def add_census_masses(areas):
    """Map commute/school producers + commute attractor onto areas by area_code."""
    supply = census_supply.load_supply()
    areas["commute_producers"] = areas["area_code"].map(
        lambda c: supply.get(c, {}).get("commute", 0.0))

    school_prod = census_school_producers.load_school_producers()
    for lvl in LEVELS:
        areas["school_producers_" + lvl] = areas["area_code"].map(
            lambda c, _l=lvl: school_prod.get(c, {}).get(_l, 0.0))

    attractor = census_attractor.load_attractor()
    areas["commute_attractor"] = areas["area_code"].map(lambda c: attractor.get(c, 0.0))

    n_sup = int((~areas["area_code"].isin(supply)).sum())
    n_att = int((~areas["area_code"].isin(attractor)).sum())
    print(f"Producers: {len(areas) - n_sup}/{len(areas)} matched to census_supply"
          + (f" ({n_sup} unmatched → 0)" if n_sup else ""))
    print(f"Commute attractor: {len(areas) - n_att}/{len(areas)} matched"
          + (f" ({n_att} unmatched → 0)" if n_att else ""))
    return areas


def add_retail_spaces(areas, zones):
    """Estimated retail parking spaces per area (parking cache sjoin + workplace fallback)."""
    if not os.path.exists(PARKING_ISLAND_CACHE):
        raise SystemExit(f"ERROR: {PARKING_ISLAND_CACHE} not found (run build_parking.py)")
    park = gpd.read_file(PARKING_ISLAND_CACHE).to_crs(PROJECTED_CRS)
    park = park[park.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    tag_cols = [c for c in park.columns if c != "geometry"]
    park["spaces"] = [
        parking_spaces({c: r[c] for c in tag_cols}, g.area)
        for r, g in zip(park.to_dict("records"), park.geometry)
    ]
    park = park[park["spaces"] > 0].copy()
    park["geometry"] = park.geometry.centroid                  # assign lots by centroid

    joined = gpd.sjoin(park[["spaces", "geometry"]], zones, how="inner", predicate="within")
    spaces_by_area = joined.groupby("area_code")["spaces"].sum()
    areas["retail_spaces_raw"] = areas["area_code"].map(spaces_by_area).fillna(0.0)

    # Workplace-derived fallback for areas with no mapped parking (median spaces/workplace over
    # matched areas) — same mechanism as build_census_zones (computed here at small-area level).
    matched = areas[(areas["retail_spaces_raw"] > 0) & (areas["workplace_pop"] > 0)]
    fb_ratio = float(np.median(matched["retail_spaces_raw"] / matched["workplace_pop"])) \
        if len(matched) else 0.0
    fb = (areas["retail_spaces_raw"] <= 0)
    areas["retail_spaces"] = np.where(
        fb, areas["workplace_pop"] * fb_ratio, areas["retail_spaces_raw"]).round(1)

    print(f"Retail: {len(park):,} lots, {spaces_by_area.sum():,.0f} spaces into "
          f"{spaces_by_area.gt(0).sum():,} areas; median spaces/workplace = {fb_ratio:.3f}; "
          f"{int(fb.sum()):,} areas used the workplace fallback")
    return areas


def add_school_demand(areas, zones):
    """Per-level school enrolment (attractor) per area via the level-tagged school cache."""
    if not os.path.exists(SCHOOL_ISLAND_CACHE):
        raise SystemExit(f"ERROR: {SCHOOL_ISLAND_CACHE} not found (run build_admin_schools.py)")
    sch = gpd.read_file(SCHOOL_ISLAND_CACHE).to_crs(PROJECTED_CRS)
    sch["geometry"] = sch.geometry.centroid
    sch = add_level_enrolments(sch)                            # + enrol_primary/postprimary/tertiary
    cols = list(LEVEL_ENROL_COLS)
    joined = gpd.sjoin(sch[[*cols, "geometry"]], zones, how="inner", predicate="within")
    by_area = joined.groupby("area_code")[cols].sum()
    for c in LEVEL_ENROL_COLS:                                 # enrol_<level> → school_demand_<level>
        lvl = c.split("_", 1)[1]
        areas["school_demand_" + lvl] = areas["area_code"].map(by_area[c]).fillna(0.0).round(1)
    tot = sum(float(areas["school_demand_" + lvl].sum()) for lvl in LEVELS)
    print(f"School: {len(sch):,} POIs, {tot:,.0f} pupils into "
          f"{by_area.index.nunique():,} areas")
    return areas


def add_centroids(areas):
    """Per-area polygon centroid in WGS84 (lat/lon)."""
    to_wgs = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)
    cent = areas.geometry.centroid
    lon, lat = to_wgs.transform(cent.x.values, cent.y.values)
    areas["centroid_lat"] = np.round(lat, 6)
    areas["centroid_lon"] = np.round(lon, 6)
    return areas


def main():
    areas = load_areas()
    zones = gpd.GeoDataFrame({"area_code": areas["area_code"].values},
                             geometry=areas.geometry.values, crs=PROJECTED_CRS)
    areas = add_census_masses(areas)
    areas = add_retail_spaces(areas, zones)
    areas = add_school_demand(areas, zones)
    areas = add_centroids(areas)

    areas["population"] = areas["population"].round().astype(int)
    for col in ("commute_attractor", "commute_producers",
                "school_producers_primary", "school_producers_postprimary",
                "school_producers_tertiary"):
        areas[col] = areas[col].astype(float).round(3)

    out = areas[COLUMNS].copy()
    out.to_csv(OUTPUT_FILE, index=False)

    # ── Verification summary ─────────────────────────────────────────────────
    print(f"\nSaved {OUTPUT_FILE}  ({len(out):,} rows)")
    print(f"  Σ population        = {out['population'].sum():,}")
    print(f"  Σ commute_producers = {out['commute_producers'].sum():,.0f}")
    print(f"  Σ commute_attractor = {out['commute_attractor'].sum():,.0f}")
    print(f"  Σ retail_spaces     = {out['retail_spaces'].sum():,.0f}")
    for lvl in LEVELS:
        print(f"  Σ school [{lvl:<11}] producers={out['school_producers_'+lvl].sum():>12,.0f}"
              f"  demand={out['school_demand_'+lvl].sum():>12,.0f}")


if __name__ == "__main__":
    main()
