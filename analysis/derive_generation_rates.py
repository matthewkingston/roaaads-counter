"""Derive per-component vehicle-driver trip-generation rates from the NTS microdata.

Writes analysis/generation_rates.json — the a-priori daily car-driver trip rate
per person for each gravity component (commute / retail / school / res).  These
pin the model's *generation* (production) magnitudes so each component's tuned
scale K_c should land at ≈ 1.0 (a verification anchor, not a fit knob): a producer
weight in vehicle-driver-trips/day means K_c·p_i = the node's daily trips.

Source
------
The NTS trip-level microdata (UKDS SN 5340), via ``analysis/nts_microdata.py``.

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
{Car / van driver, Motorcycle, Taxi / minicab} (VEHICLE_MODE_CODES).

Purpose → component  (data-derived, 23-cat)
-------------------------------------------
Component rates use ``purpose_mapping.B01_COMPONENT`` — the 23-category
TripPurpose_B01ID → component mapping (rule: **res iff endpoint is a home**;
otherwise routed by land-use).  This replaces the old 8-category + LEISURE_RETAIL_FRAC
scheme: leisure is split by data (visit-home→res, venues→retail) and escorts are
routed by destination.  Two allocations remain modelling decisions the codes don't
resolve (Business/Other-work, Personal-business → retail).

The **school** total (Education 4 + Escort-education 21) is split into
primary/post-primary/tertiary by island enrolment share — a placeholder pending a
dedicated school-generation review (self-driven Education trips suggest the
enrolment-share split may overstate tertiary; not addressed here).

``purpose_rates`` (canonical 8-cat, via TripPurpose_B04ID) is retained in the output
**only** for the not-yet-migrated temporal derivation
(analysis/derive_component_profiles.py), which still reads it; it is not the
generation anchor.

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
from purpose_mapping import B01_COMPONENT, B01_EXCLUDE, CANONICAL_PURPOSES
import nts_microdata as nts
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulation"))
import census_school_producers   # island_enrolment_by_level() — for the per-level school ρ split

OUT_FILE = "analysis/generation_rates.json"
YEARS = [2023, 2024]
# MainMode_B04ID vehicle-driver modes: Car/van driver, Motorcycle, Taxi/minicab.
VEHICLE_MODE_CODES = [3, 5, 12]
# TripPurpose_B04ID (== NTS0409a's published breakdown) → canonical purpose, for the
# legacy purpose_rates output consumed by the (not-yet-migrated) temporal derivation.
B04_TO_CANON = {1: "commuting", 2: "business", 3: "education_escort", 4: "shopping",
                5: "other_escort", 6: "personal_business", 7: "leisure", 8: "other"}


def _rates():
    """(component_rates, canonical_purpose_rates) in veh-driver trips/person/day."""
    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B01ID",
                                    "TripPurpose_B04ID", "W5", "JJXSC"], years=YEARS)
    veh = tr[tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES)].copy()
    veh["w"] = veh.JJXSC * veh.W5                       # NTS trip count × trip weight

    # Persons (diary sample): each individual weighted by their household's W2.
    ind = nts.load("individual", columns=["SurveyYear", "HouseholdID"], years=YEARS)
    hh = nts.load("household", columns=["SurveyYear", "HouseholdID", "W2"], years=YEARS)
    persons = ind.merge(hh[["HouseholdID", "W2"]], on="HouseholdID", how="left")["W2"].sum()

    # Component rates via the 23-cat B01 mapping.  Any code neither mapped nor in the
    # intentional-exclude set (17 just-walk, sentinels) is a loud error.
    veh["component"] = veh.TripPurpose_B01ID.map(B01_COMPONENT)
    stray = veh[veh.component.isna() & ~veh.TripPurpose_B01ID.isin(B01_EXCLUDE)]
    if len(stray):
        bad = sorted(stray.TripPurpose_B01ID.unique())
        sys.exit(f"ERROR: unmapped TripPurpose_B01ID codes {bad} — update B01_COMPONENT")
    comp = veh.dropna(subset=["component"]).groupby("component")["w"].sum() / persons / 7.0

    # Legacy canonical purpose_rates (B04ID) — for the temporal derivation only.
    veh["canon"] = veh.TripPurpose_B04ID.map(B04_TO_CANON)
    purp = veh.dropna(subset=["canon"]).groupby("canon")["w"].sum() / persons / 7.0
    purpose_rates = {p: float(purp.get(p, 0.0)) for p in CANONICAL_PURPOSES}
    return comp, purpose_rates


def main():
    print("Deriving generation rates from the NTS microdata (data/NTS) …")
    comp, purpose_rates = _rates()
    rates = {c: float(comp.get(c, 0.0)) for c in ("commute", "retail", "res")}

    # Split the school component (Education 4 + Escort-education 21) into
    # primary/post-primary/tertiary by island enrolment share (placeholder — pending a
    # dedicated school-generation review).  Enrolment from the admin school cache via
    # census_school_producers (the same source the attractor uses), so each level is a
    # fully independent component; the three sum to the school total.
    _enrol = census_school_producers.island_enrolment_by_level()      # {primary, postprimary, tertiary}
    _e_tot = sum(_enrol.values())
    _school_shares = {lvl: _enrol[lvl] / _e_tot for lvl in ("primary", "postprimary", "tertiary")}
    _rho_school = float(comp.get("school", 0.0))
    for lvl, share in _school_shares.items():
        rates["school_" + lvl] = _rho_school * share

    # Sanity: components should partition all vehicle-driver trips (Σ over B01 codes,
    # excluding the intentional drops).
    total_comp = sum(rates.values())
    allp = float(comp.sum())
    if abs(total_comp - allp) > 1e-9:
        print(f"  WARNING: component sum {total_comp:.4f} ≠ mapped total {allp:.4f} "
              f"(diff {total_comp - allp:+.4f}/person/day)")

    out = {
        "_meta": {
            "source": "NTS microdata (UKDS SN 5340) via analysis/nts_microdata.py",
            "years": YEARS,
            "vehicle_mode_codes_B04": VEHICLE_MODE_CODES,
            "trip_count_measure": "JJXSC (short walks ×7, series-of-calls ×0)",
            "weights": "trips W5, persons W2",
            "purpose_mapping": "purpose_mapping.B01_COMPONENT (23-cat; res iff endpoint=home)",
            "units": "vehicle-driver trips per person per day",
            "judgment_allocations": [
                "Business / Other-work -> retail (commercial premises, not the jobs count)",
                "Personal business -> retail",
                "Escorts routed by destination (commuting->commute, shopping/business->retail, "
                "education->school, home->res)",
            ],
            "school_split_by_enrolment": {lvl: round(s, 4) for lvl, s in _school_shares.items()},
            "school_split_note": "enrolment-share placeholder pending school-generation review",
            "purpose_rates_note": "canonical 8-cat (B04ID); legacy field for the temporal "
                                  "derivation only, not the generation anchor",
        },
        "rates": rates,
        "purpose_rates": purpose_rates,
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nVehicle-driver generation rates (/person/day, {YEARS} avg, "
          f"MainMode_B04ID {VEHICLE_MODE_CODES}):")
    for comp_name, r in rates.items():
        print(f"  {comp_name:16s} {r:.4f}  ({r / total_comp * 100:4.1f}%)")
    print(f"  {'total':16s} {total_comp:.4f}")
    print(f"  school split by enrolment share: "
          + ", ".join(f"{lvl} {s:.3f}" for lvl, s in _school_shares.items()))
    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
