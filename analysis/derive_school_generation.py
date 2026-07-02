"""Derive per-student school trip-generation rates (escort + self-drive) from microdata.

Writes analysis/school_generation_rates.json — for each school level (primary /
post-primary / tertiary) the **vehicle-driver trips per full-time student per day**
that the school component should generate, plus the pre-school escort magnitude that
is folded into the retail component.  Consumed by derive_generation_rates.py.

WHY THIS EXISTS (per-student, not per-capita)
---------------------------------------------
The rest of generation is per-capita × population.  School is instead anchored
**per student × student counts**: total level-L school trips = ρ_L(per student) ×
n_students_L, with the student counts supplied per jurisdiction (NI/RoI) by
census_school_producers.  This decouples the transferable *behaviour* (car trips per
school child, assumed the same everywhere given school status) from the local
*demography* (how many children of each level there are, which differs by system).
It avoids importing England's child-per-capita ratio onto Ireland.  Downstream,
model.compute_generation_scales applies these rates with k=1 (rate is already
per-student; producer = students).

TWO COMPONENTS PER LEVEL
------------------------
1. ESCORT  (an adult drives the child; TripPurpose_B01ID == 21, vehicle-driver modes)
   The escort trip is the *parent's* record and carries no child-level info, so it is
   attributed to the child via a **household regression**: per household, regress the
   escort-vehicle trip count on the number of primary / secondary / tertiary students
   (+ a pre-school bucket), weighted by W2.  The β's are the marginal escort trips per
   student of each level.  This keeps the vehicle basis (driver trips) AND handles
   ride-sharing for free (siblings sharing a car → escort ÷ kids).
2. SELF-DRIVE  (the student's own vehicle trip to education; TripPurpose_B01ID == 4,
   vehicle-driver modes = car-driver / motorcycle / taxi).  Combined per level from
   the per-age self-drive rate × the level's student age composition.  ~0 for primary/
   secondary; the majority of (small) tertiary generation.

ENGLAND AGE→LEVEL SPLIT (the hard part)
---------------------------------------
NTS gives the student's actual age but NOT their school level; a 16-18-year-old
student may be secondary (sixth form) or tertiary (FE college), which NTS cannot
distinguish (EducN lumps "school or college").  We split by age using DfE's
"Participation in education, training and employment age 16 to 21" (England,
DFE_FILE): per academic age, the share of full-time education participants in
secondary providers (state schools + independent + sixth-form colleges + special) vs
tertiary (general FE + HE).  **Definition (user-set): sixth form = secondary, FE =
tertiary.**  Two frames are reconciled:
  * DfE age is ACADEMIC age (age on 31 Aug = school year); NTS age is ACTUAL age.
    They differ by the September-cutoff half-year offset, so an NTS actual age a maps
    to 50% academic (a-1) + 50% academic a (the same convention used for the Irish
    school-year producers).  Below actual 16 → all secondary; 5-10 → primary.
  * Full-time only (mode_of_study = Full-time; EducN = full-time), to match
    census_school_producers' tertiary = full-time-students definition.

FUDGES / ASSUMPTIONS ON THE RECORD (candidate error sources)
------------------------------------------------------------
* PRE-SCHOOL → RETAIL.  Escort of 3-4-year-olds to nursery is a real vehicle trip but
  the model excludes pre-school from the school component.  Rather than drop it, its
  magnitude (β_preschool × England pre-school-per-capita ≈ 0.006/capita/day, ~0.8% of
  retail) is added to the retail component as a flat magnitude fudge — no pre-school
  producers/attractors are built, so it has retail's spatial distribution (wrong, but
  the volume is tiny and better than dropping it).  England-per-capita basis.
* SELF-DRIVE IS THE LEAST TRANSFERABLE PIECE.  Escort (parent-drives-child) is a
  near-universal pattern.  Self-drive depends on young-adult car access, which varies
  by economy/settlement — and it dominates the (small) tertiary rate.  Irish factors
  push both ways (stronger live-at-home/commuter-student culture and rurality ↑;
  higher young-driver insurance ↓).  Kept at the England rate as an acceptable
  baseline because tertiary is small on both rate and producer count; flagged as the
  first thing to revisit (vs TSNI/Irish data) if the fit ever points at school.
* SELF-DRIVE LEVEL ATTRIBUTION assumes self-drive propensity is age-driven, not
  level-driven (a sixth-form 18yo drives as much as a college 18yo).  Shakier than for
  escort (college students plausibly drive more), so it slightly over-attributes to
  secondary / under-attributes to tertiary — but the amounts are ~1% of secondary.
* AGE-18 SELF-DRIVE is a genuine spike (0.11, newly-licensed + still at home) — firmed
  from the noisy single-year value to the pooled 2013-24 non-COVID mean (FIRM_AGE18);
  only this one number uses the wider pool.
* Vehicle basis = MainMode_B04ID {car/van driver, motorcycle, taxi}; taxi (child taken
  by taxi) is a real on-road vehicle trip, not self-driving, but counted as the
  student's own vehicle trip.  Series-of-calls excluded (JJXSC).
* England NTS behaviour throughout; no Ireland trip data.  Years 2023-24 (except the
  firmed age-18 self-drive).

Usage:  python3 analysis/derive_school_generation.py
Needs data/NTS (nts_microdata) + the DfE participation CSV (DFE_FILE, gitignored).
"""

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nts_microdata as nts

OUT_FILE = "analysis/school_generation_rates.json"
# DfE "Participation in education, training and employment age 16 to 21" (2025 release,
# all-institutions national file; gitignored reference data — keep local copy).
DFE_FILE = ("data/participation-in-education-training-and-employment-age-16-to-21_2025/"
            "data/Participation_insts_and_quals_16_21_allinsts_040626.csv")
DFE_YEAR = "2024"                      # latest complete year (2025 lacks age-18 HE)
YEARS = [2023, 2024]                   # NTS behaviour window
VEH = [3, 5, 12]                       # MainMode_B04ID: car/van driver, motorcycle, taxi
FIRM_AGE18 = 0.1085                    # pooled 2013-24 (ex-COVID) age-18 self-drive/FT-student/day
# TripPurpose_B01ID codes
P_ESCORT_EDU, P_EDUCATION = 21, 4
# DfE provider_type → level (sixth form = secondary, FE = tertiary)
DFE_SECONDARY = ["All state-funded schools", "Independent schools",
                 "Sixth form colleges", "Special schools"]
DFE_TERTIARY = ["General FE tertiary and specialist colleges",
                "Higher education institutions"]
# Age_B01ID coded age → representative actual age
CODE_AGE = {4: 8, 5: 13, 6: 16, 7: 17, 8: 18, 9: 19, 10: 20, 11: 23, 12: 27}


def _england_secondary_share():
    """actual-age → secondary share of FT students, from DfE (academic age) + half-year
    academic→actual reconciliation.  <=14 → 1.0 (all secondary); 5-10 handled as primary."""
    d = pd.read_csv(DFE_FILE, dtype=str)
    d["number"] = pd.to_numeric(d["number"], errors="coerce")
    q = d[(d.participation_measure == "Education") & (d.mode_of_study == "Full-time")
          & (d.qualification_level == "Total") & (d.qualification_sublevel == "Total")
          & (d.sex == "Total") & (d.time_period == DFE_YEAR)]
    acad = {15: 1.0}                                   # academic 15 = compulsory school
    for a in ("16", "17", "18", "19", "20", "21"):
        sec = q[(q.age == a) & (q.provider_type.isin(DFE_SECONDARY))].number.sum()
        ter = q[(q.age == a) & (q.provider_type.isin(DFE_TERTIARY))].number.sum()
        acad[int(a)] = sec / (sec + ter) if (sec + ter) > 0 else 0.0

    def sec_share(actual_age):                         # actual = 50% academic(a-1) + 50% academic a
        if actual_age <= 14:
            return 1.0
        return 0.5 * acad.get(actual_age - 1, 0.0) + 0.5 * acad.get(actual_age, 0.0)
    return sec_share


def _level_shares(sec_share, code, educn):
    """(primary, secondary, tertiary) membership fractions for an individual.
    5-10 → primary; 11-15 → secondary; 16-19 → sec/ter by age split (FT students only);
    20+ → tertiary (FT students only).  Non-FT 16+ and pre-school → (0,0,0)."""
    if code == 4:
        return (1.0, 0.0, 0.0)
    if code == 5:
        return (0.0, 1.0, 0.0)
    if code in (6, 7, 8, 9):
        if educn != 1:
            return (0.0, 0.0, 0.0)
        s = sec_share(CODE_AGE[code])
        return (0.0, s, 1.0 - s)
    if code >= 10:
        return (0.0, 0.0, 1.0) if educn == 1 else (0.0, 0.0, 0.0)
    return (0.0, 0.0, 0.0)


def derive():
    sec_share = _england_secondary_share()
    ind = nts.load("individual", columns=["SurveyYear", "IndividualID", "Age_B01ID",
                                           "EducN_B01ID", "HouseholdID"], years=YEARS)
    hh = nts.load("household", columns=["SurveyYear", "HouseholdID", "W2"], years=YEARS)
    ind = ind.merge(hh[["HouseholdID", "W2"]], on="HouseholdID", how="left")
    sh = ind.apply(lambda r: _level_shares(sec_share, r.Age_B01ID, r.EducN_B01ID),
                   axis=1, result_type="expand")
    sh.columns = ["prim", "sec", "ter"]
    ind = pd.concat([ind, sh], axis=1)
    ind["pre"] = ind.Age_B01ID.isin([2, 3]).astype(float)      # 1-2, 3-4 = pre-school

    tr = nts.load("trip", columns=["SurveyYear", "TripPurpose_B01ID", "MainMode_B04ID",
                                   "IndividualID", "HouseholdID", "W5", "JJXSC"], years=YEARS)

    # ── ESCORT: household regression of escort-vehicle trips on level-student counts ──
    esc = tr[(tr.TripPurpose_B01ID == P_ESCORT_EDU) & (tr.MainMode_B04ID.isin(VEH))]
    hk = (ind.groupby("HouseholdID")[["pre", "prim", "sec", "ter"]].sum()
          .join(hh.set_index("HouseholdID")["W2"]))
    hk["y"] = esc.groupby("HouseholdID")["JJXSC"].sum().reindex(hk.index).fillna(0.0)
    D = hk[hk[["pre", "prim", "sec", "ter"]].sum(axis=1) > 0]
    X = D[["pre", "prim", "sec", "ter"]].values
    sw = np.sqrt(D["W2"].values)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], D["y"].values * sw, rcond=None)
    esc_pre, esc_prim, esc_sec, esc_ter = beta / 7.0            # per student per day

    # ── SELF-DRIVE: per-age own-vehicle education rate, combined per level ──
    sd = tr[(tr.TripPurpose_B01ID == P_EDUCATION) & (tr.MainMode_B04ID.isin(VEH))].copy()
    sd["w"] = sd.JJXSC * sd.W5
    ind["sd"] = ind.IndividualID.map(sd.groupby("IndividualID")["w"].sum()).fillna(0.0)

    def sd_rate(code):
        g = ind[ind.Age_B01ID == code] if code <= 5 else \
            ind[(ind.Age_B01ID == code) & (ind.EducN_B01ID == 1)]
        if len(g) == 0 or g.W2.sum() == 0:
            return 0.0
        return FIRM_AGE18 if code == 8 else g.sd.sum() / g.W2.sum() / 7.0

    def sd_level(idx):
        num = den = 0.0
        for code in range(4, 22):
            g = ind[ind.Age_B01ID == code]
            if len(g) == 0:
                continue
            share = _level_shares(sec_share, code, 1)[idx]
            students = (g[g.EducN_B01ID.isin([1, 2])].W2.sum() if code >= 6 else g.W2.sum())
            num += sd_rate(code) * students * share
            den += students * share
        return num / den if den else 0.0

    sd_prim, sd_sec, sd_ter = sd_level(0), sd_level(1), sd_level(2)

    # ── pre-school escort → retail magnitude (England per-capita) ──
    pre_percapita = ind[ind.Age_B01ID.isin([2, 3])].W2.sum() / ind.W2.sum()

    rates = {
        "primary":     esc_prim + sd_prim,
        "postprimary": esc_sec + sd_sec,
        "tertiary":    esc_ter + sd_ter,
    }
    detail = {
        "escort":     {"primary": esc_prim, "postprimary": esc_sec, "tertiary": esc_ter},
        "self_drive": {"primary": sd_prim, "postprimary": sd_sec, "tertiary": sd_ter},
    }
    preschool_escort_percapita = esc_pre * pre_percapita
    return rates, detail, esc_pre, preschool_escort_percapita


def main():
    print("Deriving per-student school generation from NTS microdata + DfE England split …")
    rates, detail, esc_pre, pre_pc = derive()
    out = {
        "_meta": {
            "units": "vehicle-driver trips per full-time student per day (escort + self-drive)",
            "nts_years": YEARS,
            "dfe_source": DFE_FILE,
            "dfe_year": DFE_YEAR,
            "age18_self_drive_firmed": FIRM_AGE18,
            "detail": detail,
            "preschool_escort_per_preschooler": esc_pre,
            "note": "sixth form=secondary, FE=tertiary; FT students; half-year academic→actual "
                    "alignment; tertiary self-drive = least transferable (see module docstring)",
        },
        "rates_per_student": rates,
        "preschool_escort_retail_percapita": pre_pc,
    }
    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)
    print("\nPer-student school generation (escort + self-drive, trips/FT-student/day):")
    for lvl in ("primary", "postprimary", "tertiary"):
        print(f"  {lvl:12s} escort {detail['escort'][lvl]:.4f} + self-drive "
              f"{detail['self_drive'][lvl]:.4f} = {rates[lvl]:.4f}")
    print(f"\n  pre-school escort → retail: {pre_pc:.5f}/capita/day "
          f"(β_preschool {esc_pre:.4f}/preschooler)")
    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    main()
