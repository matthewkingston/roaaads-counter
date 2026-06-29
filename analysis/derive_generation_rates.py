"""Derive per-component vehicle-driver trip-generation rates from NTS0409a.

Writes analysis/generation_rates.json — the a-priori daily car-driver trip rate
per person for each gravity component (commute / retail / school / res).  These
pin the model's *generation* (production) magnitudes so each component's tuned
scale K_c should land at ≈ 1.0 (a verification anchor, not a fit knob): a producer
weight in vehicle-driver-trips/day means K_c·p_i = the node's daily trips.

Source
------
NTS0409a "Average number of trips by purpose and main mode (trips per person per
year)", England (data/nts0409.ods, sheet NTS0409a_trips), averaged over 2023+2024.

Vehicle basis
-------------
The model calibrates against on-road *vehicle* counts, so we sum the modes that
put one vehicle on the road per household-recorded trip:
  Car or van driver + Motorcycle + Taxi or minicab.
Buses are excluded (many passengers per vehicle; the driver is not a household
trip).  Using the *driver* row makes it vehicles by construction (no occupancy
correction) and means "Education or escort education" is already the adult escort
(a child can't drive) — so no all-mode education down-weighting is needed.

Purpose → component mapping  (JUDGMENT ALLOCATIONS — candidate error sources)
----------------------------------------------------------------------------
The organising principle is the *attractor* each component offers:
    commute → workplace (jobs)            retail → retail_spaces (PARKING = all
    school  → school places                        commercial / venue)
    res     → population (HOMES)
NTS0409a cannot sub-split these purposes by car-driver mode, so the allocations
below are deliberate modelling decisions, NOT data lookups.  They are the first
thing to revisit if the fit is scrutinised:

  * Commuting        → commute   (pure home ↔ own workplace).
  * Business         → RETAIL    (not commute): business visits go to commercial
                                  premises captured by the parking proxy, not to
                                  the workplace-jobs count.  Keeps commute pure.
  * Personal business→ RETAIL    : services / banks / medical are commercial,
                                  parking-attracted destinations.
  * Shopping         → retail.
  * Education/escort → school.
  * Leisure          → SPLIT  LEISURE_RETAIL_FRAC to retail (venue: entertainment,
                                  sport, holiday, day trip → parking) and the rest
                                  to res (visit friends at home → pop↔pop).  The
                                  0.5 split is the single largest assumption here —
                                  leisure is the biggest bucket — and is a pure
                                  judgment (no driver-mode sub-split is published).
  * Other escort, Other → res    (residual discretionary, pop↔pop).

Usage
-----
  python3 analysis/derive_generation_rates.py
Re-run whenever data/nts0409.ods changes or the purpose mapping is revised.
"""

import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from purpose_mapping import COMPONENT_PURPOSES, LEISURE_RETAIL_FRAC, CANONICAL_PURPOSES

NTS0409_FILE = "data/nts0409.ods"
OUT_FILE     = "analysis/generation_rates.json"
NTS_YEARS    = [2023, 2024]
VEHICLE_MODES = ["Car or van driver", "Motorcycle", "Taxi or minicab"]

# Canonical purpose key (purpose_mapping.py) → NTS0409a column display name.
# The purpose→component mapping itself lives in purpose_mapping.COMPONENT_PURPOSES
# (shared with derive_component_profiles.py); here we only resolve the canonical
# keys to this table's columns.
NTS0409_COL = {
    "commuting":         "Commuting",
    "business":          "Business",
    "education_escort":  "Education or escort education",
    "shopping":          "Shopping",
    "other_escort":      "Other escort",
    "personal_business": "Personal business",
    "leisure":           "Leisure",
    "other":             "Other",
}


def _resolve_col(df, name):
    """Match an NTS purpose column, tolerating a trailing ' [note N]' suffix."""
    hits = [c for c in df.columns if c == name or c.startswith(name + " [")]
    if len(hits) != 1:
        sys.exit(f"ERROR: column '{name}' matched {len(hits)} headers (expected 1)")
    return hits[0]


def main():
    print(f"Loading {NTS0409_FILE} …")
    df = pd.read_excel(NTS0409_FILE, sheet_name="NTS0409a_trips",
                       header=5, engine="odf")
    df.columns = [str(c).strip() for c in df.columns]
    ycol, mcol = df.columns[0], df.columns[1]
    df[ycol] = pd.to_numeric(df[ycol], errors="coerce")

    sub = df[df[ycol].isin(NTS_YEARS)
             & df[mcol].astype(str).str.strip().isin(VEHICLE_MODES)]
    n_expected = len(NTS_YEARS) * len(VEHICLE_MODES)
    if len(sub) != n_expected:
        sys.exit(f"ERROR: expected {n_expected} (year×mode) rows, got {len(sub)} "
                 f"— check years {NTS_YEARS} and modes {VEHICLE_MODES} exist")

    # Per-person/day rate for an NTS0409a purpose column = Σ(year,mode) trips/yr ÷ n_years ÷ 365.
    def purpose_rate(name):
        col = _resolve_col(df, name)
        return sub[col].astype(float).sum() / len(NTS_YEARS) / 365.0

    # Per-canonical-purpose car-driver rate ρ_p (consumed by derive_component_profiles.py
    # to ρ-weight each component's temporal-shape blend) and per-component rates.
    purpose_rates = {p: purpose_rate(NTS0409_COL[p]) for p in CANONICAL_PURPOSES}
    rates = {comp: sum(w * purpose_rates[p] for p, w in terms)
             for comp, terms in COMPONENT_PURPOSES.items()}

    # Sanity: components should partition All purposes (vehicle modes only).
    total_comp = sum(rates.values())
    allp = purpose_rate("All purposes")
    if abs(total_comp - allp) > 1e-6:
        print(f"  WARNING: component sum {total_comp:.4f} ≠ All purposes {allp:.4f} "
              f"(diff {total_comp - allp:+.4f}/person/day)")

    out = {
        "_meta": {
            "source": NTS0409_FILE,
            "sheet": "NTS0409a_trips",
            "years": NTS_YEARS,
            "vehicle_modes": VEHICLE_MODES,
            "leisure_retail_frac": LEISURE_RETAIL_FRAC,
            "units": "vehicle-driver trips per person per day",
            "judgment_allocations": [
                "Business -> retail (not commute)",
                "Personal business -> retail",
                f"Leisure split {LEISURE_RETAIL_FRAC} retail / {1 - LEISURE_RETAIL_FRAC} res",
            ],
        },
        "rates": rates,
        "purpose_rates": purpose_rates,
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nVehicle-driver generation rates (/person/day, {NTS_YEARS} avg, "
          f"{'+'.join(VEHICLE_MODES)}):")
    for comp, r in rates.items():
        print(f"  {comp:8s} {r:.4f}  ({r / total_comp * 100:4.1f}%)")
    print(f"  {'total':8s} {total_comp:.4f}")
    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
