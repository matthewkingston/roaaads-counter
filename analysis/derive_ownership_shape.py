"""INSPECTION ARTIFACT — car-availability -> car-driver-trip band shapes (NTS SN 5340).

Module 1 of the per-area car-ownership work (see .claude/plans/eager-painting-pearl.md):
derive, per gravity component x household car-availability band, the car-driver trip
rate, plus per-cell sample sizes and a diagnostic plot.  This is the transferable donor
curve ("how does car-driver trip generation rise with household car ownership, by
purpose") — the numerator of a future per-area mobilisation multiplier.  DATA-ONLY: it
writes only analysis/ownership_shape.json + reports/ownership_shape.png and touches NO
model artifact (no producer weighting, no census_zones/node_weights change, no sigma_mob
change, no re-tune).

Band variable: household NumCarVan (car-or-van count, plain numeric) -> bands 0/1/2/3+
(BANDS).  Matches the NI census "car or van availability" variable; RoI SAPS T15_1 is
cars-only (a minor definitional gap resolved later in the harmonisation module).

Res / commute / retail  (a pure stratification of the existing per-capita rate):
    rate(c, band) = Sum(JJXSC*W5 over veh-driver trips of component c in band) [numerator]
                  / Sum(W2 persons in band)                                    [denominator]
                  / 7
The same quantity as derive_generation_rates' rho_c but split by band; because it is a
pure partition, the band-population-weighted mean of rate(c, band) recovers rho_c
EXACTLY (the consistency gate below; retail compared to the RAW B01 rate, before the
pre-school escort fudge that derive_generation_rates adds).

School (primary/postprimary/tertiary, full escort-stratified) mirrors
derive_school_generation._compute with a band dimension:
  * ESCORT: run the W2-weighted NON-NEGATIVE household regression (nnls) of escort-
    vehicle-trip count on (pre/prim/sec/ter student counts) WITHIN each band ->
    beta_L(band).  Non-negativity is a no-op where the effect is well identified
    (primary/postprimary) and floors the ill-identified negative tertiary bands (0/1,
    college students are rarely car-escorted) to 0.  Band 0 ~ 0 (no car to drive) is a
    sanity check.  NB a per-band regression re-aggregated by band student-weights only
    APPROXIMATELY recovers the pooled beta (regression is not a simple partition) — so
    the school consistency print is indicative, not exact.
  * SELF-DRIVE: per-age own-vehicle education rate stratified by the student's household
    band, combined per level.  The firmed age-18 self-drive constant (FIRM_AGE18) is
    kept BAND-INVARIANT — the faithful mirror of _compute; a flagged caveat, since it
    hides any ownership effect on young-adult self-drive (revisit if it matters).

Usage:  python3 analysis/derive_ownership_shape.py   [--years 2023,2024]
Needs data/NTS (nts_microdata), the DfE participation CSV (derive_school_generation),
and simulation/node_weights_reduced.json (for the school students/pop consistency ref).
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from purpose_mapping import B01_COMPONENT, B01_EXCLUDE
import nts_microdata as nts
import derive_generation_rates as dgr
import derive_school_generation as dsg

OUT_JSON = "analysis/ownership_shape.json"
OUT_PNG = "reports/ownership_shape.png"
YEARS = [2023, 2024]
VEHICLE_MODE_CODES = dgr.VEHICLE_MODE_CODES         # {3,5,12} car/van driver, motorcycle, taxi
SCHOOL_LEVELS = dgr.SCHOOL_LEVELS                    # ("primary","postprimary","tertiary")

# Car-availability bands: label -> (lo, hi_inclusive); top band is open (3+).
BANDS = [("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3+", 3, 99)]
BAND_LABELS = [b[0] for b in BANDS]


def band_of(n):
    """NumCarVan count -> band label (str), or None if missing."""
    if pd.isna(n):
        return None
    n = int(n)
    for lab, lo, hi in BANDS:
        if lo <= n <= hi:
            return lab
    return None


# ------------------------------------------------------- res / commute / retail

def _bandmap(years):
    """HouseholdID -> band label, from the household NumCarVan (missing dropped)."""
    hh = nts.load("household", columns=["SurveyYear", "HouseholdID", "NumCarVan"],
                  years=years)
    hh["band"] = hh.NumCarVan.map(band_of)
    return hh[["HouseholdID", "band"]]


def b01_band_rates(years):
    """rate(component, band) for res/commute/retail (+ raw B01 retail), with per-cell
    counts.  Returns (rates, trip_n, persons_w, persons_n, point_rates)."""
    bandmap = _bandmap(years)

    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B01ID",
                                   "HouseholdID", "W5", "JJXSC"], years=years)
    ind = nts.load("individual", columns=["SurveyYear", "HouseholdID"], years=years)
    hhw = nts.load("household", columns=["SurveyYear", "HouseholdID", "W2"], years=years)

    # numerator: vehicle-driver trips, tagged with component + band
    veh = tr[tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES)].copy()
    veh["component"] = veh.TripPurpose_B01ID.map(B01_COMPONENT)
    stray = veh[veh.component.isna() & ~veh.TripPurpose_B01ID.isin(B01_EXCLUDE)]
    if len(stray):
        sys.exit(f"ERROR: unmapped TripPurpose_B01ID {sorted(stray.TripPurpose_B01ID.unique())}")
    veh = veh.dropna(subset=["component"]).merge(bandmap, on="HouseholdID", how="left")
    veh = veh.dropna(subset=["band"])
    veh["w"] = veh.JJXSC * veh.W5

    num = (veh.groupby(["component", "band"])["w"].sum()
           .unstack("band").reindex(columns=BAND_LABELS).fillna(0.0))
    trip_n = (veh.groupby(["component", "band"]).size()
              .unstack("band").reindex(columns=BAND_LABELS).fillna(0).astype(int))

    # denominator: persons per band (each individual <- household W2)
    ind = ind.merge(hhw[["HouseholdID", "W2"]], on="HouseholdID", how="left")
    ind = ind.merge(bandmap, on="HouseholdID", how="left").dropna(subset=["band"])
    persons_w = ind.groupby("band")["W2"].sum().reindex(BAND_LABELS).fillna(0.0)
    persons_n = ind.groupby("band").size().reindex(BAND_LABELS).fillna(0).astype(int)

    rates = num.div(persons_w, axis=1) / 7.0            # rate(component, band)
    point_rates = dgr._rates()                           # raw B01 point rates (Series)
    return rates, trip_n, persons_w, persons_n, point_rates


# ------------------------------------------------------------------ school

def school_band_rates(years):
    """Per-level school rate(band) = escort_L(band) + selfdrive_L(band), escort-stratified.
    Returns (rates_by_level{lvl:{band:rate}}, detail, design_counts, students_w)."""
    prepared = dsg._prepare()
    bandmap = _bandmap(years)                            # HouseholdID -> band
    b2 = bandmap.set_index("HouseholdID")["band"]

    ind = prepared["ind"].copy()
    ind["band"] = ind.HouseholdID.map(b2)
    hk = prepared["hk"].copy()
    hk["band"] = hk.index.map(b2)
    sec_share = prepared["sec_share"]

    rates = {lvl: {} for lvl in SCHOOL_LEVELS}
    detail = {"escort": {lvl: {} for lvl in SCHOOL_LEVELS},
              "self_drive": {lvl: {} for lvl in SCHOOL_LEVELS}}
    design = {}                       # band -> #regression households + per-level student mass
    students_w = {lvl: {} for lvl in SCHOOL_LEVELS}

    for lab in BAND_LABELS:
        # -- escort: W2-weighted household regression within this band --
        sub = hk[hk.band == lab]
        Xall = sub[["pre", "prim", "sec", "ter"]].to_numpy()
        keep = Xall.sum(axis=1) > 0
        X = Xall[keep]
        if len(X) >= 4:
            sw = np.sqrt(sub["W2"].to_numpy()[keep])
            # Non-negative escort β — escort trips/student can't be < 0.  A no-op where the
            # plain-lstsq β are already positive (primary/postprimary); floors the ill-
            # identified negative tertiary bands (0/1) to 0 without touching the rest.
            beta, _ = nnls(X * sw[:, None], sub["y"].to_numpy()[keep] * sw)
            esc_pre, esc_prim, esc_sec, esc_ter = beta / 7.0
        else:
            esc_prim = esc_sec = esc_ter = float("nan")
        design[lab] = {"reg_households": int(keep.sum()),
                       "student_mass": {lvl: float(sub[col].mul(sub["W2"]).sum())
                                        for lvl, col in zip(SCHOOL_LEVELS,
                                                            ["prim", "sec", "ter"])}}

        # -- self-drive within this band --
        indb = ind[ind.band == lab]
        age = indb["Age_B01ID"].to_numpy()
        educ = indb["EducN_B01ID"].to_numpy()
        W2 = indb["W2"].to_numpy()
        sdv = indb["sd"].to_numpy()

        def sd_rate(code):
            sel = (age == code) if code <= 5 else ((age == code) & (educ == 1))
            den = W2[sel].sum()
            if den == 0:
                return 0.0
            return dsg.FIRM_AGE18 if code == 8 else sdv[sel].sum() / den / 7.0

        def sd_level(idx):
            num = den = 0.0
            for code in range(4, 22):
                sel = age == code
                if not sel.any():
                    continue
                share = dsg._level_shares(sec_share, code, 1)[idx]
                students = (W2[sel & np.isin(educ, [1, 2])].sum() if code >= 6
                            else W2[sel].sum())
                num += sd_rate(code) * students * share
                den += students * share
            return (num / den) if den else 0.0

        sd_prim, sd_sec, sd_ter = sd_level(0), sd_level(1), sd_level(2)

        esc = {"primary": esc_prim, "postprimary": esc_sec, "tertiary": esc_ter}
        sdl = {"primary": sd_prim, "postprimary": sd_sec, "tertiary": sd_ter}
        for lvl in SCHOOL_LEVELS:
            rates[lvl][lab] = float(esc[lvl] + sdl[lvl])
            detail["escort"][lvl][lab] = float(esc[lvl])
            detail["self_drive"][lvl][lab] = float(sdl[lvl])
            students_w[lvl][lab] = design[lab]["student_mass"][lvl]

    return rates, detail, design, students_w


# ------------------------------------------------------------------ assembly

def _weighted_mean(rate_by_band, wt_by_band):
    num = sum(rate_by_band[b] * wt_by_band[b] for b in BAND_LABELS
              if not (isinstance(rate_by_band[b], float) and np.isnan(rate_by_band[b])))
    den = sum(wt_by_band[b] for b in BAND_LABELS
              if not (isinstance(rate_by_band[b], float) and np.isnan(rate_by_band[b])))
    return (num / den) if den else float("nan")


def _relative(rate_by_band, mean):
    return {b: (rate_by_band[b] / mean if mean else float("nan")) for b in BAND_LABELS}


def plot(rates, out_png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                               # noqa: BLE001
        print(f"  (plot skipped: {e})")
        return
    x = list(range(len(BAND_LABELS)))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    for comp, r in rates.items():
        y = [r["rate"][b] for b in BAND_LABELS]
        ax1.plot(x, y, marker="o", label=comp)
        yr = [r["relative"][b] for b in BAND_LABELS]
        ax2.plot(x, yr, marker="o", label=comp)
    for ax, ttl in ((ax1, "absolute rate (trips/person or /student per day)"),
                    (ax2, "relative to band-weighted mean")):
        ax.set_xticks(x); ax.set_xticklabels(BAND_LABELS)
        ax.set_xlabel("household cars/vans (NumCarVan band)")
        ax.set_title(ttl); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax2.axhline(1.0, color="k", lw=0.8, ls="--")
    fig.suptitle("Car-availability band shapes (NTS SN 5340, %s)" % "+".join(map(str, YEARS)))
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=110)
    print(f"  plot -> {out_png}")


def main(years):
    global YEARS
    YEARS = years
    dgr.YEARS = years
    dsg.YEARS = years
    print(f"Deriving car-availability band shapes (NTS {years}, bands {BAND_LABELS}) …")

    b01_rates, trip_n, persons_w, persons_n, point_rates = b01_band_rates(years)
    sch_rates, sch_detail, sch_design, students_w = school_band_rates(years)

    out = {"_meta": {
        "note": "INSPECTION ARTIFACT (module 1) — car-availability band shapes; no model wiring",
        "source": "NTS microdata (UKDS SN 5340) via analysis/nts_microdata.py",
        "years": years,
        "band_variable": "household NumCarVan (car-or-van count)",
        "bands": BAND_LABELS,
        "vehicle_mode_codes_B04": VEHICLE_MODE_CODES,
        "units": "res/commute/retail = veh-driver trips/person/day; school = trips/FT-student/day",
        "school_escort": "per-band W2-weighted NON-NEGATIVE household regression (nnls); "
                         "indicative, not an exact partition of the pooled beta",
        "school_age18_selfdrive": "FIRM_AGE18 kept band-invariant (faithful mirror of "
                                  "derive_school_generation._compute)",
    }, "components": {}}

    # res / commute / retail
    print("\n== res / commute / retail: rate(component, band) [n trips] ==")
    print("persons/band (W2):  " + "  ".join(f"{b}={persons_w[b]:,.0f}(n={persons_n[b]})"
                                             for b in BAND_LABELS))
    for comp in ("res", "commute", "retail"):
        rbb = {b: float(b01_rates.loc[comp, b]) for b in BAND_LABELS}
        mean = _weighted_mean(rbb, persons_w)
        ref = float(point_rates.get(comp, float("nan")))
        out["components"][comp] = {
            "basis": "per-capita", "rate": rbb, "relative": _relative(rbb, mean),
            "band_weighted_mean": mean, "point_rate_ref": ref,
            "n_trips": {b: int(trip_n.loc[comp, b]) for b in BAND_LABELS},
        }
        # Pure partition ⇒ wmean == rho up to the (negligible) missing-NumCarVan drop;
        # flag only a real divergence well above that (a genuine tagging/denominator bug).
        rel = abs(mean - ref) / ref if ref else 0.0
        flag = "  <-- MISMATCH" if rel > 5e-3 else ""
        print(f"  {comp:8s} " + "  ".join(f"{b}:{rbb[b]:.4f}[{trip_n.loc[comp,b]}]"
                                          for b in BAND_LABELS)
              + f"   wmean={mean:.4f} vs rho={ref:.4f} (Δ={rel*100:.2f}%, "
                f"missing-band drop){flag}")

    # school
    print("\n== school (escort-stratified): rate(level, band) ==")
    dsg_rates, dsg_detail, _, _ = dsg.derive()          # pooled per-student ref
    for lvl in SCHOOL_LEVELS:
        rbb = sch_rates[lvl]
        wt = students_w[lvl]
        mean = _weighted_mean(rbb, wt)
        ref = float(dsg_rates[lvl])
        out["components"]["school_" + lvl] = {
            "basis": "per-student", "rate": rbb, "relative": _relative(rbb, mean),
            "band_weighted_mean": mean, "point_rate_ref": ref,
            "escort": sch_detail["escort"][lvl], "self_drive": sch_detail["self_drive"][lvl],
            "reg_households": {b: sch_design[b]["reg_households"] for b in BAND_LABELS},
            "student_mass": {b: students_w[lvl][b] for b in BAND_LABELS},
        }
        print(f"  {lvl:11s} " + "  ".join(f"{b}:{rbb[b]:.4f}" for b in BAND_LABELS)
              + f"   wmean~{mean:.4f} vs pooled={ref:.4f}  (indicative)")
    print("  reg households/band: " + "  ".join(f"{b}={sch_design[b]['reg_households']}"
                                                for b in BAND_LABELS))

    out["_meta"]["persons_per_band_W2"] = {b: float(persons_w[b]) for b in BAND_LABELS}
    out["_meta"]["persons_per_band_n"] = {b: int(persons_n[b]) for b in BAND_LABELS}

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved -> {OUT_JSON}")
    plot(out["components"], OUT_PNG)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Inspect NTS car-availability band shapes.")
    ap.add_argument("--years", default="2023,2024",
                    help="comma-separated NTS survey years (default 2023,2024)")
    args = ap.parse_args()
    main([int(y) for y in args.years.split(",")])
