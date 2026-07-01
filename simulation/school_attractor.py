"""Per-POI school ATTRACTOR enrolment split by level: primary / post-primary / tertiary.

Phase-2 attractor-side counterpart to census_school_producers.py. Reads the unified, level-tagged
school cache (SCHOOL_ISLAND_CACHE, from build_admin_schools: admin school-age + curated OSM
tertiary) and adds three enrolment columns — `enrol_primary` / `enrol_postprimary` /
`enrol_tertiary` — so build_census_zones (external zones, aggregated per zone) and
build_demographics (internal core, POI-distributed to road nodes) build the three attractor node
layers from one source.

Level mapping: cache `level` primary→primary, post_primary→post-primary, tertiary→tertiary.
**Special schools are mixed-age**, so each special POI's enrolment is split into primary/post-primary
by the **NI primary:post-primary enrolment ratio** (the same ratio census_school_producers folds
special by — user's decision), contributing to *both* layers at that school's location. Tertiary
gets no special. The three columns sum to the original `enrolment` (special is split, not dropped),
so Σ level layers ≡ the old single `node_school_demand`.
"""
import geopandas as gpd

from demographics_config import SCHOOL_ISLAND_CACHE, PROJECTED_CRS

LEVEL_ENROL_COLS = ("enrol_primary", "enrol_postprimary", "enrol_tertiary")


def ni_primary_ratio(gdf):
    """NI primary / (NI primary + NI post-primary) enrolment — the special-school split ratio.
    Computed from the cache so it matches census_school_producers' special fold exactly."""
    ni = gdf[gdf["jurisdiction"] == "NI"]
    p  = float(ni.loc[ni["level"] == "primary", "enrolment"].sum())
    pp = float(ni.loc[ni["level"] == "post_primary", "enrolment"].sum())
    return p / (p + pp) if (p + pp) > 0 else 0.5


def add_level_enrolments(gdf):
    """Return a copy of the school-cache GeoDataFrame with the three level-enrolment columns added."""
    r_prim = ni_primary_ratio(gdf)
    e = gdf["enrolment"].astype(float)
    lvl = gdf["level"]
    is_special = (lvl == "special")
    gdf = gdf.copy()
    gdf["enrol_primary"]     = e.where(lvl == "primary", 0.0) + e.where(is_special, 0.0) * r_prim
    gdf["enrol_postprimary"] = e.where(lvl == "post_primary", 0.0) + e.where(is_special, 0.0) * (1.0 - r_prim)
    gdf["enrol_tertiary"]    = e.where(lvl == "tertiary", 0.0)
    return gdf


def load_level_enrolments():
    """Load the school cache in PROJECTED_CRS with the three level-enrolment columns added."""
    return add_level_enrolments(gpd.read_file(SCHOOL_ISLAND_CACHE).to_crs(PROJECTED_CRS))


if __name__ == "__main__":
    g = load_level_enrolments()
    r = ni_primary_ratio(g)
    print(f"{len(g)} school POIs; NI primary:post-primary special-split ratio = {r:.4f}")
    tot = float(g["enrolment"].sum())
    print(f"{'level':>14} {'enrolment':>12}")
    for c in LEVEL_ENROL_COLS:
        print(f"{c:>14} {g[c].sum():>12,.0f}")
    lvlsum = sum(float(g[c].sum()) for c in LEVEL_ENROL_COLS)
    print(f"  Σ level layers = {lvlsum:,.0f}  vs single 'enrolment' total = {tot:,.0f}  "
          f"(diff {lvlsum - tot:+.1f})")
    # by jurisdiction, for the NI:RoI sanity check
    for j in ("NI", "RoI"):
        gj = g[g["jurisdiction"] == j]
        print(f"  {j:>4}: " + "  ".join(f"{c.split('_')[1]}={gj[c].sum():,.0f}" for c in LEVEL_ENROL_COLS))
