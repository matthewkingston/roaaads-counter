"""Per-small-area car-commute *attractor* (jobs reached by car), harmonised NI (DZ) + RoI (SA).

`load_attractor() -> {area_code: car_commute_jobs}` for every NI Data Zone and RoI Small Area,
computed **once, island-wide** (jurisdiction handled internally) and consumed by both
build_census_zones.py (external zones — aggregated) and build_demographics.py (internal core —
POI-distributed). It is the attractor counterpart of the car-commute producer in
census_supply.py; together they make the commute gravity component car-specific (the model
assigns car flow). Keyed purely by `area_code`, so a core anywhere on the island works.

Definitions
-----------
NI  car_jobs[DZ] = apwp001 workplace total[DZ] (place-of-work jobs, incl. WFH)
        × car_share[parent SDZ]
    car_share[SDZ] = (Driving a car or van + Motorcycle/scooter + Taxi)
                     ───────────────────────────────────────────────────
                     Σ(all apwp035 travel-method columns, INCLUDING WFH)
    The DZ workplace mode-split is published only at SDZ level (apwp035), so the DZ magnitude is
    disaggregated by the parent-SDZ car-share. WFH is removed in the WFH-inclusive denominator
    (numerator already excludes it). This is exact, not an approximation, because empirically
    apwp001_total[SDZ] ≡ Σ(all apwp035 method columns)[SDZ] (verified across all 850 SDZs to
    NISRA disclosure rounding) — the DZ total and the ratio denominator share one universe, so
    Σ_DZ car_jobs[DZ] reproduces the SDZ car-commuter count exactly. DZ→SDZ comes from the
    DZ2021.geojson attribute table (DZ2021_cd / SDZ2021_cd — a plain attribute read, no geometry).

RoI car_jobs[SA] = `commute_car` column of cache_sa_workplace.csv (WZ daytime drivers
    T2_M5+T2_M6+T2_M8 apportioned WZ→SA, × national work-driver share 0.9588 — see
    build_wz_apportionment.py). Fails loud if the column is absent (cache predates this change).
"""
import csv
import glob

import geopandas as gpd
import pandas as pd

APWP001_FILE    = "data/census-2021-apwp001.xlsx"      # DZ workplace totals (sheet "DZ")
APWP035_FILE    = "data/ni_census/census-2021-apwp035.xlsx"  # SDZ travel-to-work mode split
DZ_BOUNDARY     = "simulation/dz2021/DZ2021.geojson"   # DZ→SDZ parent map (attributes only)
ROI_WP_CACHE    = "data/ireland_data/cache_sa_workplace.csv"

# apwp035 method columns selected as car-driver modes (the model assigns car flow).
_NI_CAR_METHODS = ("Driving a car or van", "Motorcycle, scooter or moped", "Taxi")


def _dz_to_sdz():
    """DZ code → parent SDZ code, read from the DZ boundary attribute table (no geometry op)."""
    df = gpd.read_file(DZ_BOUNDARY, ignore_geometry=True, columns=["DZ2021_cd", "SDZ2021_cd"])
    return {str(d): str(s) for d, s in zip(df["DZ2021_cd"], df["SDZ2021_cd"])}


def _ni_attractor():
    """NI per-DZ car-commute jobs = apwp001 DZ total × parent-SDZ car-share (apwp035)."""
    # DZ workplace totals (place-of-work jobs, incl. WFH).
    wp = pd.read_excel(APWP001_FILE, sheet_name="DZ", header=5)
    dz_total = {
        str(c): float(v)
        for c, v in zip(wp["Geography Code"], pd.to_numeric(wp["Workplace population"],
                                                            errors="coerce").fillna(0.0))
        if str(c).startswith("N20")
    }

    # SDZ car-share = car methods / Σ(all methods incl. WFH).
    sdz = pd.read_excel(APWP035_FILE, sheet_name="SDZ", header=5)
    method_cols = [c for c in sdz.columns if c not in ("Geography", "Geography Code")]
    car_share = {}
    for _, row in sdz.iterrows():
        code = str(row["Geography Code"])
        if not code.startswith("N21"):
            continue
        denom = sum(float(row[c] or 0) for c in method_cols)
        car   = sum(float(row[c] or 0) for c in _NI_CAR_METHODS)
        car_share[code] = (car / denom) if denom > 0 else 0.0

    dz_sdz = _dz_to_sdz()
    out = {}
    for dz, total in dz_total.items():
        sdz_code = dz_sdz.get(dz)
        if sdz_code is None:
            raise KeyError(f"DZ {dz} has no parent SDZ in {DZ_BOUNDARY}")
        out[dz] = total * car_share.get(sdz_code, 0.0)
    return out


def _roi_attractor():
    """RoI per-SA car-commute jobs from the `commute_car` column of the WZ apportionment cache."""
    hits = glob.glob(ROI_WP_CACHE)
    if not hits:
        raise FileNotFoundError(
            f"{ROI_WP_CACHE} missing — run simulation/build_wz_apportionment.py first.")
    out = {}
    with open(ROI_WP_CACHE, newline="") as f:
        r = csv.DictReader(f)
        if "commute_car" not in (r.fieldnames or []):
            raise KeyError(
                f"{ROI_WP_CACHE} has no 'commute_car' column — re-run "
                f"simulation/build_wz_apportionment.py (delete the cache first) to add it.")
        for row in r:
            out[str(row["sa_code"])] = float(row["commute_car"] or 0.0)
    return out


def load_attractor():
    """{area_code: car_commute_jobs} for all NI DZ + RoI SA (island-wide, computed once)."""
    out = {}
    out.update(_roi_attractor())
    out.update(_ni_attractor())
    return out


if __name__ == "__main__":
    ni  = _ni_attractor()
    roi = _roi_attractor()
    nc  = sum(ni.values()); rc = sum(roi.values())
    print(f"NI DZs:  {len(ni):,}  car-commute jobs = {nc:,.0f}")
    print(f"RoI SAs: {len(roi):,}  car-commute jobs = {rc:,.0f}")
    print(f"TOTAL    car-commute jobs = {nc + rc:,.0f}")
    # Reconcile NI against the apwp035 NI-sheet car total (independent NISRA marginal).
    sdz = pd.read_excel(APWP035_FILE, sheet_name="NI", header=5)
    row = sdz.iloc[0]
    ni_car_marginal = sum(float(row[c] or 0) for c in _NI_CAR_METHODS)
    print(f"NI reconcile: Σ DZ car-jobs {nc:,.0f}  vs apwp035 NI car total "
          f"{ni_car_marginal:,.0f}  (diff {nc - ni_car_marginal:+,.0f})")
