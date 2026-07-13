"""Derive the island car-driver mobilisation level + spatial-dispersion width.

Writes analysis/mobilisation.json — the single island-wide mobilisation factor and
its common-mode uncertainty, consumed by the K-normalisation prior:

  * ``island_level``  — car-driver-ish journeys per person per day, island-wide.
  * ``m_island``      — ``island_level / Σρ_NTS`` — the single global multiplier that
                        rescales the England-NTS generation rates ρ_c to the island
                        level (model.load_generation_rates applies it, preserving the
                        NTS purpose split), so each K_c anchors cleanly at 1.
  * ``sigma_mob``     — the common-mode (all-components-together) prior width for Λ_K:
                        the spatial dispersion of car mobilisation (a place could sit
                        this far from the island anchor).

WHY A SINGLE ISLAND NUMBER (not a NI/RoI split)
-----------------------------------------------
Car-driver mobilisation is driven by settlement density / car-dependence — continuous
covariates that vary *within* each jurisdiction and don't jump at the border.  A hard
NI/RoI rescale would impose a discontinuity that doesn't exist and mis-assign rural-NI /
urban-RoI areas, and would not be CENTRE-portable.  So the level is one island constant
and the genuine spatial spread becomes ``sigma_mob`` (the loose leash the local counts
then refine).  See the plan / agent memory on model portability.

LEVEL SOURCES (locked; journey ≈ car-movement, level × NTS_share_c)
------------------------------------------------------------------
* NI  — TSNI Table 3.2, 2023, journeys/person/year (all ages): Car Driver + Motorcycle
        + Other private (vans/lorries) + Taxi = 455  →  /365 per day.
* RoI — CSO NTA, 2019: Σ_(ageband×sex) NTA04 journeys/person/day × RoI census population
        (SAPS T1_1AGE* single-year cols) × 0.688 car-driver-ish fraction (NTA11 2019:
        car-driver 64.9 + taxi/hackney 0.9 + lorry/motorcycle 2.8 + van 0.2) / RoI total
        population (all ages).  NTA04 covers 18+ only; the under-18 driver gap is
        negligible (only 17-yos drive) and the all-ages denominator matches TSNI.
* Island — population-weighted pool of the two per-capita levels (NI + RoI census pops).

Mixed vintage (NI 2023 / RoI 2019) is a deliberate best-available-per-jurisdiction call
(NI 2023 aligns with the 2023/24 NTS shares + 2023 counts; RoI 2019 is its last pre-COVID
wave).  The mode sets match the NTS {car/van driver, motorcycle, taxi} basis; NI/RoI fold
in a small lorry share NTS mode 3 excludes (negligible).

SIGMA_MOB (spatial dispersion — the common-mode width)
------------------------------------------------------
Two derived components, combined in quadrature:
  * between-region  — pop-weighted SD of the NI vs RoI per-capita levels (relative).
  * within-region   — pop-weighted CV of a census car-ownership intensity across small
                      areas (RoI SAPS T15_1 cars/household; NI person-weighted car
                      availability), a PROXY for the within-jurisdiction spread of car
                      *use*.  (Units differ NI vs RoI — person- vs household-weighted —
                      so this is approximate; documented.)
The recipe (how to weight between vs within, ownership-as-use-proxy) is a modelling
choice — the components are all written out so it can be revisited without a re-run.

Usage:  python3 analysis/derive_mobilisation.py
Needs the TSNI ODS, data/ireland_nts/ (NTA04, NTA11), the RoI SAPS CSV, NI NISRA
population + car-availability CSVs, and analysis/generation_rates.json (for Σρ_NTS).
"""

import glob
import json
import math
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT_FILE = "analysis/mobilisation.json"
GEN_RATES_FILE = "analysis/generation_rates.json"

# ── NI level (TSNI Table 3.2, 2023) ──────────────────────────────────────────
NI_JOURNEYS_PER_YEAR = 434 + 2 + 13 + 6      # Car Driver + Motorcycle + Other private + Taxi
NI_LEVEL_MODES = "Car Driver + Motorcycle + Other private (van/lorry) + Taxi (TSNI Table_3_2, 2023)"
NISRA_POP_CACHE = "data/cache_nisra_population.csv"
NI_CAR_GLOB = "data/ni_census/*hh_car_van_tc5_pers*.csv"

# ── RoI level (CSO NTA, 2019) ────────────────────────────────────────────────
NTA04_GLOB = "data/ireland_nts/NTA04*.csv"
ROI_YEAR = "2019"
ROI_CAR_FRACTION = 0.688                      # NTA11 2019: car-driver+taxi/hackney+lorry/motorcycle+van
SAPS_CSV = ("data/ireland_data/Complete_set_of_Census_2022_SAPs/"
            "SAPS_2022_Small_Area_UR_171024.csv")
STATE_ROW_POP = 100_000                       # drop the CSO "State" aggregate row (pop ~5M)

# NTA04 age band -> SAPS T1_1AGE* single-year/5-year column stems (append sex M/F)
NTA04_BANDS = {
    "18 - 24 years":      ["T1_1AGE18", "T1_1AGE19", "T1_1AGE20_24"],
    "25 - 34 years":      ["T1_1AGE25_29", "T1_1AGE30_34"],
    "35 - 44 years":      ["T1_1AGE35_39", "T1_1AGE40_44"],
    "45 - 54 years":      ["T1_1AGE45_49", "T1_1AGE50_54"],
    "55 - 64 years":      ["T1_1AGE55_59", "T1_1AGE60_64"],
    "65 - 74 years":      ["T1_1AGE65_69", "T1_1AGE70_74"],
    "75 years and over":  ["T1_1AGE75_79", "T1_1AGE80_84", "T1_1AGEGE_85"],
}
SEXES = {"Male": "M", "Female": "F"}


def _one(pattern):
    m = sorted(glob.glob(pattern))
    if not m:
        sys.exit(f"ERROR: no file matching {pattern}")
    return m[0]


def _wmean_sd(vals, wts):
    """Population-weighted mean, SD, and CV (SD/mean)."""
    W = sum(wts)
    mean = sum(v * w for v, w in zip(vals, wts)) / W
    var = sum(w * (v - mean) ** 2 for v, w in zip(vals, wts)) / W
    sd = math.sqrt(var)
    return mean, sd, (sd / mean if mean else 0.0)


def roi_level():
    """RoI car-driver-ish journeys per capita per day (all ages) + total population."""
    nta = {}
    for r in pd.read_csv(_one(NTA04_GLOB)).to_dict("records"):
        if str(r["Year"]) == ROI_YEAR:
            nta[(r["Age Group"], r["Sex"])] = float(r["VALUE"])

    pop_band = {}          # (band, sex) -> RoI population
    roi_pop = 0.0
    saps = pd.read_csv(SAPS_CSV)
    saps = saps[saps["T1_1AGETT"].astype(float) <= STATE_ROW_POP]   # drop State aggregate row
    roi_pop = float(saps["T1_1AGETT"].astype(float).sum())
    for band, stems in NTA04_BANDS.items():
        for sex, sfx in SEXES.items():
            cols = [s + sfx for s in stems]
            pop_band[(band, sex)] = float(saps[cols].astype(float).sum().sum())

    journeys_18plus = sum(nta[(b, s)] * pop_band[(b, s)]
                          for b in NTA04_BANDS for s in SEXES)      # all-mode 18+ journeys/day
    car_driver_journeys = journeys_18plus * ROI_CAR_FRACTION
    level = car_driver_journeys / roi_pop
    return level, roi_pop, {"journeys_18plus_per_day": journeys_18plus,
                            "car_fraction": ROI_CAR_FRACTION, "year": ROI_YEAR}


def ni_level():
    """NI car-driver journeys per capita per day (TSNI direct) + total population."""
    level = NI_JOURNEYS_PER_YEAR / 365.0
    pop = pd.read_csv(NISRA_POP_CACHE)
    pop = pop[(pop["TLIST(A1)"] == 2021) & (pop["DZ2021"].astype(str).str.startswith("N20"))]
    ni_pop = float(pop["VALUE"].sum())
    return level, ni_pop, {"journeys_per_year": NI_JOURNEYS_PER_YEAR, "modes": NI_LEVEL_MODES}


def _roi_ownership_cv():
    """Pop-weighted CV of RoI cars/household (SAPS T15_1_*) across Small Areas."""
    vals, wts = [], []
    for row in pd.read_csv(SAPS_CSV).to_dict("records"):
        if float(row["T1_1AGETT"]) > STATE_ROW_POP:
            continue
        denom = float(row["T15_1_TC"]) - float(row["T15_1_NSC"])      # exclude not-stated
        if denom <= 0:
            continue
        cars = (1 * float(row["T15_1_1C"]) + 2 * float(row["T15_1_2C"])
                + 3 * float(row["T15_1_3C"]) + 4 * float(row["T15_1_GE4C"]))  # GE4 counted as 4
        vals.append(cars / denom)
        wts.append(float(row["T1_1AGETT"]))
    return _wmean_sd(vals, wts)


def _ni_ownership_cv():
    """Pop-weighted CV of NI mean cars per person's household across Data Zones.
    NOTE person-weighted (the NISRA file is *_pers), so not unit-identical to RoI's
    household-weighted metric — an approximation, see module docstring."""
    band = {"0": 0, "1": 1, "2": 2, "3": 3, "4": 4, "5": 5}
    dz = {}
    for r in pd.read_csv(_one(NI_CAR_GLOB)).itertuples(index=False):
        code = str(r[2])
        if code == "-8":
            continue
        dz.setdefault(r[0], {})[code] = float(r[4])
    vals, wts = [], []
    for _, b in dz.items():
        tot = sum(b.values())
        if tot <= 0:
            continue
        vals.append(sum(band[k] * v for k, v in b.items()) / tot)
        wts.append(tot)
    return _wmean_sd(vals, wts)


def main():
    print("Deriving island mobilisation level + spatial dispersion …")
    ni_l, ni_pop, ni_meta = ni_level()
    roi_l, roi_pop, roi_meta = roi_level()
    island_level = (ni_l * ni_pop + roi_l * roi_pop) / (ni_pop + roi_pop)

    sigma_nts = sum(json.load(open(GEN_RATES_FILE))["rates"].values())
    m_island = island_level / sigma_nts

    # sigma_mob components
    _, between_sd, _ = _wmean_sd([ni_l, roi_l], [ni_pop, roi_pop])
    between_cv = between_sd / island_level
    roi_mean, _, roi_cv = _roi_ownership_cv()
    ni_mean, _, ni_cv = _ni_ownership_cv()
    within_cv = math.sqrt((roi_pop * roi_cv ** 2 + ni_pop * ni_cv ** 2) / (roi_pop + ni_pop))
    sigma_mob = math.hypot(between_cv, within_cv)   # default recipe; components stored for revision

    out = {
        "_meta": {
            "purpose": "island car-driver mobilisation level + common-mode prior width",
            "ni_source": ni_meta,
            "roi_source": roi_meta,
            "ni_population": ni_pop,
            "roi_population": roi_pop,
            "sigma_nts": sigma_nts,
            "vintage_note": "NI TSNI 2023, RoI CSO 2019 — deliberate best-available-per-jurisdiction",
            "sigma_mob_recipe": "hypot(between_region_cv, within_region_ownership_cv); "
                                "within uses car-ownership CV as a proxy for use-dispersion "
                                "(RoI household-weighted, NI person-weighted — approximate)",
        },
        "ni_level": ni_l,
        "roi_level": roi_l,
        "island_level": island_level,
        "m_island": m_island,
        "sigma_mob": sigma_mob,
        "sigma_mob_components": {
            "between_region_cv": between_cv,
            "within_region_cv": within_cv,
            "within_roi_ownership_cv": roi_cv,
            "within_ni_ownership_cv": ni_cv,
        },
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\n  NI level      {ni_l:.4f}/cap/day   (pop {ni_pop:,.0f})")
    print(f"  RoI level     {roi_l:.4f}/cap/day   (pop {roi_pop:,.0f})")
    print(f"  island level  {island_level:.4f}/cap/day  (RoI/NI ratio {roi_l/ni_l:.3f})")
    print(f"  Σρ NTS        {sigma_nts:.4f}   →   m_island = {m_island:.4f}")
    print(f"\n  sigma_mob = {sigma_mob:.4f}  "
          f"(between {between_cv:.4f} ⊕ within {within_cv:.4f}; "
          f"ownership CV RoI {roi_cv:.3f} / NI {ni_cv:.3f})")
    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
