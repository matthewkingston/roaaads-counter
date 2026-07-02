"""Derive per-component vehicle-driver trip-generation rates from the NTS microdata.

Writes analysis/generation_rates.json — the a-priori daily car-driver trip rate
per person for each gravity component (commute / retail / school / res).  These
pin the model's *generation* (production) magnitudes so each component's tuned
scale K_c should land at ≈ 1.0 (a verification anchor, not a fit knob): a producer
weight in vehicle-driver-trips/day means K_c·p_i = the node's daily trips.

Source
------
The NTS trip-level microdata (UKDS SN 5340), via ``analysis/nts_microdata.py``.
This replaces the pre-aggregated NTS0409a table: the per-person/day rate is built
directly from trip records (England, SurveyYear in YEARS), so the derivation is no
longer bound to that table's fixed purpose × mode breakdown.  Reproduces the
published NTS0409a rate to rounding (All-purposes veh-driver = 371.2 trips/person/yr).

Rate = trips-per-person-per-day, computed the way NTS grosses NTS0409a:
    numerator   = Σ(JJXSC × W5) over vehicle-driver trips        [diary-week trips]
    denominator = Σ W2 over individuals (via their household W2)  [persons]
    per day     = (numerator / denominator) / 7                  [diary week = 7 days]
JJXSC is the NTS trip count (short walks ×7, series-of-calls ×0); W5 the trip weight
(folds in the household weight); W2 the diary-sample person weight.

Vehicle basis
-------------
The model calibrates against on-road *vehicle* counts, so we keep the modes that
put one vehicle on the road per household-recorded trip — MainMode_B04ID in
{Car / van driver, Motorcycle, Taxi / minicab} (VEHICLE_MODE_CODES).  Using the
*driver* mode makes it vehicles by construction (no occupancy correction), and
"Education or escort education" is already the adult escort (a child can't drive).

Purpose → component mapping  (JUDGMENT ALLOCATIONS — candidate error sources)
----------------------------------------------------------------------------
Trip purpose is TripPurpose_B04ID — the 8-category breakdown that matches NTS0409a's
published columns 1:1 (B04_TO_CANON).  The organising principle is the *attractor*
each component offers:
    commute → workplace (jobs)            retail → retail_spaces (PARKING = all
    school  → school places                        commercial / venue)
    res     → population (HOMES)
The allocations below (in purpose_mapping.COMPONENT_PURPOSES) are deliberate
modelling decisions, the first thing to revisit if the fit is scrutinised:

  * Commuting        → commute   (pure home ↔ own workplace).
  * Business         → RETAIL    (commercial premises = parking, not the jobs count).
  * Personal business→ RETAIL    (services / banks / medical = parking destinations).
  * Shopping         → retail.
  * Education/escort → school.
  * Leisure          → SPLIT  LEISURE_RETAIL_FRAC to retail (venue) and the rest to
                                res (visit friends at home → pop↔pop).  The 0.5 split
                                is the single largest assumption here.  (The finer
                                TripPurpose_B01ID codes can data-derive this split —
                                a planned improvement, not applied here.)
  * Other escort, Other → res    (residual discretionary, pop↔pop).

Usage
-----
  python3 analysis/derive_generation_rates.py
Re-run whenever the microdata (data/NTS), the purpose mapping, or the enrolment
split changes.
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from purpose_mapping import COMPONENT_PURPOSES, LEISURE_RETAIL_FRAC, CANONICAL_PURPOSES
import nts_microdata as nts
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulation"))
import census_school_producers   # island_enrolment_by_level() — for the per-level school ρ split

OUT_FILE = "analysis/generation_rates.json"
YEARS = [2023, 2024]
# MainMode_B04ID vehicle-driver modes: Car/van driver, Motorcycle, Taxi/minicab.
VEHICLE_MODE_CODES = [3, 5, 12]
# TripPurpose_B04ID (== NTS0409a's published purpose breakdown) → canonical purpose.
B04_TO_CANON = {1: "commuting", 2: "business", 3: "education_escort", 4: "shopping",
                5: "other_escort", 6: "personal_business", 7: "leisure", 8: "other"}


def _purpose_rates():
    """{canonical_purpose: vehicle-driver trips/person/day} from the microdata."""
    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B04ID",
                                    "W5", "JJXSC"], years=YEARS)
    veh = tr[tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES)].copy()
    veh["canon"] = veh.TripPurpose_B04ID.map(B04_TO_CANON)
    if veh["canon"].isna().any():
        bad = sorted(veh.loc[veh["canon"].isna(), "TripPurpose_B04ID"].unique())
        sys.exit(f"ERROR: unmapped TripPurpose_B04ID codes {bad} — update B04_TO_CANON")
    trips = veh.assign(w=veh.JJXSC * veh.W5).groupby("canon")["w"].sum()   # Σ(JJXSC×W5)/purpose

    # Persons (diary sample): each individual weighted by their household's W2.
    ind = nts.load("individual", columns=["SurveyYear", "HouseholdID"], years=YEARS)
    hh = nts.load("household", columns=["SurveyYear", "HouseholdID", "W2"], years=YEARS)
    persons = ind.merge(hh[["HouseholdID", "W2"]], on="HouseholdID", how="left")["W2"].sum()

    # trips/person/day = (diary-week trips per person) / 7
    return {p: float(trips.get(p, 0.0)) / persons / 7.0 for p in CANONICAL_PURPOSES}


def main():
    print("Deriving generation rates from the NTS microdata (data/NTS) …")
    purpose_rates = _purpose_rates()
    rates = {comp: sum(w * purpose_rates[p] for p, w in terms)
             for comp, terms in COMPONENT_PURPOSES.items()}

    # Split the single "school" component into primary/post-primary/tertiary, each with its OWN
    # generation rate = ρ_school × that level's island enrolment share.  Frozen here as three
    # explicit constants (enrolment from the admin school cache via census_school_producers — the
    # same source the attractor uses), so each level is a fully independent component: later
    # changing one level's data never shifts another's ρ.  The three sum to ρ_school (education
    # generation is conserved), so the component-partition sanity check below still holds.
    _enrol = census_school_producers.island_enrolment_by_level()      # {primary, postprimary, tertiary}
    _e_tot = sum(_enrol.values())
    _school_shares = {lvl: _enrol[lvl] / _e_tot for lvl in ("primary", "postprimary", "tertiary")}
    _rho_school = rates.pop("school")
    for lvl, share in _school_shares.items():
        rates["school_" + lvl] = _rho_school * share

    # Sanity: components should partition all vehicle-driver trips (Σ over canonical purposes).
    total_comp = sum(rates.values())
    allp = sum(purpose_rates.values())
    if abs(total_comp - allp) > 1e-9:
        print(f"  WARNING: component sum {total_comp:.4f} ≠ all purposes {allp:.4f} "
              f"(diff {total_comp - allp:+.4f}/person/day)")

    out = {
        "_meta": {
            "source": "NTS microdata (UKDS SN 5340) via analysis/nts_microdata.py",
            "years": YEARS,
            "vehicle_mode_codes_B04": VEHICLE_MODE_CODES,
            "trip_count_measure": "JJXSC (short walks ×7, series-of-calls ×0)",
            "weights": "trips W5, persons W2",
            "leisure_retail_frac": LEISURE_RETAIL_FRAC,
            "units": "vehicle-driver trips per person per day",
            "judgment_allocations": [
                "Business -> retail (not commute)",
                "Personal business -> retail",
                f"Leisure split {LEISURE_RETAIL_FRAC} retail / {1 - LEISURE_RETAIL_FRAC} res",
            ],
            "school_split_by_enrolment": {lvl: round(s, 4) for lvl, s in _school_shares.items()},
        },
        "rates": rates,
        "purpose_rates": purpose_rates,
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nVehicle-driver generation rates (/person/day, {YEARS} avg, "
          f"MainMode_B04ID {VEHICLE_MODE_CODES}):")
    for comp, r in rates.items():
        print(f"  {comp:16s} {r:.4f}  ({r / total_comp * 100:4.1f}%)")
    print(f"  {'total':16s} {total_comp:.4f}")
    print(f"  school split by enrolment share: "
          + ", ".join(f"{lvl} {s:.3f}" for lvl, s in _school_shares.items()))
    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
