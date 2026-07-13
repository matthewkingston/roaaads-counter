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
otherwise routed by land-use): leisure is split by endpoint (visit-home→res,
venues→retail) and escorts are routed by destination.  Two allocations remain
modelling decisions the codes don't resolve (Business/Other-work, Personal-business →
retail).

**Every rate is per-capita** — it encodes the island TOTAL journeys of its type (rate ×
pop), and the producer/attractor layer only distributes them spatially.  School behaviour
is derived **per full-time student** (escort + self-drive; analysis/derive_school_generation.py)
because that is the transferable quantity, then converted here to per-capita by ×(island
students_L / population) using the node-weight totals the model's k_students divides by
(so ρ_school/k_students recovers the per-student rate exactly).  The per-capita education
codes (4+21) from the B01 mapping are NOT used for school.  Retail additionally absorbs
the pre-school escort magnitude (already per-capita) as a documented fudge (no pre-school
producers exist).  This makes the school per-capita rates all-Ireland-specific — re-run if
the census school producers / population change.

Usage
-----
  python3 analysis/derive_school_generation.py   # first — writes school_generation_rates.json
  python3 analysis/derive_generation_rates.py
Re-run whenever the microdata (data/NTS), the purpose mapping, or the school rates change.
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from purpose_mapping import B01_COMPONENT, B01_EXCLUDE
import nts_microdata as nts

OUT_FILE = "analysis/generation_rates.json"
# Canonical 6-component order for the bootstrap stat_cov (matches SCHOOL_LEVELS + the
# non-school components; the K-prior in tune_assignment.py keys off these names).
SCHOOL_LEVELS = ("primary", "postprimary", "tertiary")
COMP_ORDER = ("res", "commute", "retail",
              "school_primary", "school_postprimary", "school_tertiary")
BOOT_SEED = 20260713          # fixed so the derived stat_cov is reproducible
# Per-student school rates (escort + self-drive) + the pre-school escort magnitude
# routed to retail, from analysis/derive_school_generation.py.
SCHOOL_RATES_FILE = "analysis/school_generation_rates.json"
# Node weights (same island totals the model's k_students uses) — to convert the
# per-student school behaviour to a per-capita rate that cancels exactly with k_students.
NODE_WEIGHTS_FILE = "simulation/node_weights_reduced.json"
YEARS = [2023, 2024]
# MainMode_B04ID vehicle-driver modes: Car/van driver, Motorcycle, Taxi/minicab.
VEHICLE_MODE_CODES = [3, 5, 12]


def _load_gen_frames():
    """Load the trip + person frames once (individual pre-merged with its household W2),
    so a cluster bootstrap can recompute rates without re-reading the microdata."""
    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B01ID",
                                    "HouseholdID", "W5", "JJXSC"], years=YEARS)
    ind = nts.load("individual", columns=["SurveyYear", "HouseholdID"], years=YEARS)
    hh = nts.load("household", columns=["SurveyYear", "HouseholdID", "W2"], years=YEARS)
    ind = ind.merge(hh[["HouseholdID", "W2"]], on="HouseholdID", how="left")
    return tr, ind


def _rates_from(frames, hh_mult=None):
    """Per-component veh-driver trips/person/day (23-cat B01 mapping) from preloaded
    frames, optionally under a household-resample multiplicity (dict HouseholdID→count).
    hh_mult=None ⇒ the point estimate (all weights ×1) — bit-identical to the un-refactored
    path.  A resample scales every trip weight and every person weight by the household's
    multiplicity (numerator and denominator alike), i.e. the household counted that often."""
    tr, ind = frames
    veh = tr[tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES)].copy()
    w = veh.JJXSC * veh.W5                             # NTS trip count × trip weight
    person_w = ind["W2"]                               # each individual ← household W2
    if hh_mult is not None:
        w = w * veh["HouseholdID"].map(hh_mult).fillna(0.0)
        person_w = person_w * ind["HouseholdID"].map(hh_mult).fillna(0.0)
    veh["w"] = w
    persons = float(person_w.sum())

    # Component rates via the 23-cat B01 mapping.  Any code neither mapped nor in the
    # intentional-exclude set (17 just-walk, sentinels) is a loud error.
    veh["component"] = veh.TripPurpose_B01ID.map(B01_COMPONENT)
    stray = veh[veh.component.isna() & ~veh.TripPurpose_B01ID.isin(B01_EXCLUDE)]
    if len(stray):
        bad = sorted(stray.TripPurpose_B01ID.unique())
        sys.exit(f"ERROR: unmapped TripPurpose_B01ID codes {bad} — update B01_COMPONENT")
    return veh.dropna(subset=["component"]).groupby("component")["w"].sum() / persons / 7.0


def _rates():
    """Point-estimate per-component rates (the B01 mapping)."""
    return _rates_from(_load_gen_frames())


def _students_pop():
    """Per-level island students/population ratio (census, held FIXED in the bootstrap —
    it is not survey-sampled)."""
    nw = json.load(open(NODE_WEIGHTS_FILE))
    pop = sum(nw.get("node_population", {}).values())
    return {lvl: sum(nw.get(f"node_school_producers_{lvl}", {}).values()) / pop
            for lvl in SCHOOL_LEVELS}


def _six_vector(comp, srates, pre_pc, students_pop):
    """Assemble the six per-capita component rates from the two derivations, exactly as
    main() builds generation_rates.json (retail absorbs the pre-school escort fudge;
    school per-capita = per-student × island students/pop)."""
    v = {
        "res":     float(comp.get("res", 0.0)),
        "commute": float(comp.get("commute", 0.0)),
        "retail":  float(comp.get("retail", 0.0)) + float(pre_pc),
    }
    for lvl in SCHOOL_LEVELS:
        v[f"school_{lvl}"] = float(srates[lvl]) * students_pop[lvl]
    return v


def bootstrap_stat_cov(n_boot, seed=BOOT_SEED):
    """PSU-cluster bootstrap → the 6×6 covariance of the log generation rates.

    Resamples PSUs with replacement *within SurveyYear* (region is degenerate in the
    2023/24 NTS subset), giving a household multiplicity that is threaded through BOTH
    the non-school B01 rates (_rates_from) AND the school escort+self-drive derivation
    (derive_school_generation._compute), so the six rates carry survey sampling variance
    and their true cross-correlations (shared resampled households).  Reported in
    LOG/relative space because K_c is a multiplier on ρ_c (Var(K_c) ≈ Var(log ρ̂_c)); the
    school students/pop ratio is held fixed (census).

    Returns (list(COMP_ORDER), cov 6×6 ndarray).  Needs the DfE participation CSV
    (school _prepare) + node weights, same as the full derivation."""
    import derive_school_generation as dsg
    gen_frames = _load_gen_frames()
    school_prepared = dsg._prepare()
    students_pop = _students_pop()

    meta = nts.load("household", columns=["SurveyYear", "HouseholdID", "PSUID"], years=YEARS)
    hids = meta["HouseholdID"].to_numpy()
    years = meta["SurveyYear"].to_numpy()
    _, hid_psu_idx = np.unique(meta["PSUID"].to_numpy(), return_inverse=True)
    n_psu = int(hid_psu_idx.max()) + 1
    year_psu_idx = {yr: np.unique(hid_psu_idx[years == yr]) for yr in np.unique(years)}

    rng = np.random.default_rng(seed)
    logs = np.empty((n_boot, len(COMP_ORDER)))
    for b in range(n_boot):
        counts = np.zeros(n_psu)
        for pidx in year_psu_idx.values():                       # stratify by SurveyYear
            np.add.at(counts, rng.choice(pidx, size=len(pidx), replace=True), 1.0)
        hh_mult = dict(zip(hids, counts[hid_psu_idx]))           # per-household multiplicity
        comp = _rates_from(gen_frames, hh_mult)
        srates, _, _, pre_pc = dsg._compute(school_prepared, hh_mult)
        vec = _six_vector(comp, srates, pre_pc, students_pop)
        logs[b] = [math.log(max(vec[c], 1e-30)) for c in COMP_ORDER]
    return list(COMP_ORDER), np.cov(logs, rowvar=False)


def main(n_boot=500, seed=BOOT_SEED):
    print("Deriving generation rates from the NTS microdata (data/NTS) …")
    comp = _rates()
    # commute / res are per-capita from the B01 mapping.  (comp["school"] — the per-capita
    # education codes 4+21 — is NOT used: school is now per-student, see below.)
    rates = {c: float(comp.get(c, 0.0)) for c in ("commute", "res")}

    # School: convert the per-student behaviour (escort + self-drive, from
    # derive_school_generation.py) to a PER-CAPITA rate so every rate is per-capita and
    # encodes an island total — ρ_percapita_L = ρ_perstudent_L × (island students_L / pop).
    # The students/pop ratio uses the SAME node-weight totals the model's k_students
    # divides by, so ρ_school/k_students recovers the per-student rate exactly (the
    # producer layer then does pure spatial distribution).  All-Ireland-specific: re-run
    # if the census school producers / population change.  Retail additionally absorbs the
    # pre-school escort magnitude (already per-capita; a documented fudge — no pre-school
    # producers/attractors exist).
    with open(SCHOOL_RATES_FILE) as f:
        sch = json.load(f)
    nw = json.load(open(NODE_WEIGHTS_FILE))
    pop = sum(nw.get("node_population", {}).values())
    for lvl in ("primary", "postprimary", "tertiary"):
        students = sum(nw.get(f"node_school_producers_{lvl}", {}).values())
        rates["school_" + lvl] = float(sch["rates_per_student"][lvl]) * (students / pop)
    rates["retail"] = (float(comp.get("retail", 0.0))
                       + float(sch["preschool_escort_retail_percapita"]))

    # Data-integrity: the per-capita B01 components (incl. per-capita school) partition all
    # vehicle-driver trips — a check on the B01 mapping, independent of the per-student swap.
    b01_total = float(comp.sum())
    b01_partition = sum(float(comp.get(c, 0.0)) for c in comp.index)
    if abs(b01_total - b01_partition) > 1e-9:
        print(f"  WARNING: B01 partition {b01_partition:.4f} ≠ mapped total {b01_total:.4f}")

    out = {
        "_meta": {
            "source": "NTS microdata (UKDS SN 5340) via analysis/nts_microdata.py",
            "years": YEARS,
            "vehicle_mode_codes_B04": VEHICLE_MODE_CODES,
            "trip_count_measure": "JJXSC (short walks ×7, series-of-calls ×0)",
            "weights": "trips W5, persons W2",
            "purpose_mapping": "purpose_mapping.B01_COMPONENT (23-cat; res iff endpoint=home)",
            "units": "vehicle-driver trips per person per day (per-capita); school_* converted "
                     "from per-student behaviour via island students_L/pop (see derive_school_generation)",
            "judgment_allocations": [
                "Business / Other-work -> retail (commercial premises, not the jobs count)",
                "Personal business -> retail",
                "Escorts routed by destination (commuting->commute, shopping/business->retail, "
                "education->school, home->res)",
            ],
            "school_source": "analysis/derive_school_generation.py (per-student escort + "
                             "self-drive; England age->level split; see that module)",
            "retail_preschool_fudge_percapita": float(sch["preschool_escort_retail_percapita"]),
        },
        "rates": rates,
    }

    # Cluster-bootstrap sampling covariance of the six log-rates (the DERIVED width for the
    # K-normalisation prior's differential/split directions; see analysis/tune_assignment.py).
    if n_boot and n_boot > 0:
        print(f"  bootstrap: PSU-cluster resample × {n_boot} (seed {seed}) …")
        comps, cov = bootstrap_stat_cov(n_boot, seed)
        out["_meta"]["stat_cov"] = {
            "components": comps,
            "space": "log-rate (relative; Var(K_c)≈Var(log ρ̂_c))",
            "cluster": "PSU within SurveyYear (region degenerate in 2023/24)",
            "n_boot": n_boot,
            "seed": seed,
            "cov": cov.tolist(),
        }
        sd = np.sqrt(np.diag(cov))
        print("  per-component sampling CV (√diag of log-cov): "
              + ", ".join(f"{c}={s:.3f}" for c, s in zip(comps, sd)))

    with open(OUT_FILE, "w") as f:
        json.dump(out, f, indent=2)

    print("\nGeneration rates (per person/day, per-capita):")
    print(f"  commute {rates['commute']:.4f}  retail {rates['retail']:.4f}  res {rates['res']:.4f}"
          f"   (retail incl. pre-school escort fudge {sch['preschool_escort_retail_percapita']:.5f})")
    print(f"  school (per-capita; from per-student × island students/pop): "
          f"primary {rates['school_primary']:.4f}  post-primary {rates['school_postprimary']:.4f}  "
          f"tertiary {rates['school_tertiary']:.4f}")
    print(f"\nSaved → {OUT_FILE}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Derive per-component generation rates "
                                             "(+ bootstrap stat_cov for the K-prior).")
    ap.add_argument("--boot", type=int, default=500,
                    help="bootstrap replicates for stat_cov (0 to skip; default 500)")
    ap.add_argument("--seed", type=int, default=BOOT_SEED, help="bootstrap RNG seed")
    args = ap.parse_args()
    main(n_boot=args.boot, seed=args.seed)
