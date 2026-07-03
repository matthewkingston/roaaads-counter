"""Derive per-component temporal SHAPE profiles from the NTS trip-level microdata.

Under generation pinning (analysis/derive_generation_rates.py) each component's
absolute magnitude AND the inter-component split are set by data (ρ_c × producers,
K_c ≈ 1).  The temporal profiles f_c therefore carry ONLY each component's relative
time-of-day-and-week SHAPE — they do not partition a shared aggregate.

This is the **car-specific**, microdata-derived replacement for the old NTS0502a/0504b
(all-mode aggregate-table) derivation.  The profile is measured DIRECTLY as the joint
distribution of vehicle-driver trip START times: for each component (the shared
purpose_mapping.B01_COMPONENT scheme, same as the generation rates) we count car-driver
trips (MainMode_B04ID {3,5,12}, weight JJXSC×W5 — the NTS trip-count measure and trip
weight per the Data Extract User Guide) by (day_type, hour) and normalise.  This retires,
in one step, the old PROXY_0502 merged-purpose proxies, the separable V(dow)×H(hour)
assumption, the all-mode weekend fallback, the Bayes-flip, and the ρ-weighting/purpose_rates
dependency (real trip counts already weight by volume).  Calendar day-of-week comes from
TravelWeekDay_B01ID on the Day table (1=Mon..7=Sun), joined to trips via DayID.

Bins and day_type
-----------------
The model collapses day-of-week to day_type — weekday (Mon–Fri, averaged), Saturday,
Sunday — so the estimand is the 72 (day_type, hour) slots per component (weekday pools
Mon–Fri, so those bins are well sampled; Sat/Sun are single-day).  We write the weekday
profile identically into all five weekday rows of the 168-row CSV (the model averages
them; the official observations are per day_type), which also maximises weekday bin stats.

Adaptive per-bin year pooling (MIN_BIN_N)
-----------------------------------------
The temporal profile is a PINNED input with no propagated uncertainty, so a bin's
sampling noise propagates straight into the χ² at that slot — every bin needs adequate
statistics, even near-zero night hours.  But pooling years costs a little year-drift, so
we pool **only where needed**: each (day_type, hour) bin expands its year window through
the ex-COVID tiers (2023/24 → +2018/19 → +2013–17) until its unweighted count reaches
MIN_BIN_N, then stops.  Well-sampled daytime/weekday bins stay on 2023/24 (zero drift);
only thin weekend-night bins pool wider.  This is valid because the shape is near-stationary
across the window (weekday-peak drift ≈ 1 pp), so a bin's fraction is the same estimand at
any window — the **normalisation that makes windows comparable is dividing each bin's pooled
weighted count by its number of years** (a per-year rate).  A handful of intrinsically thin
res weekend-night bins cannot reach MIN_BIN_N even on the full window (car trips at 3 am are
genuinely rare); the residual noise there is removed by a **light within-night smoothing**
(NIGHT_HOURS / NIGHT_SMOOTH_WIN): a small moving average over the flat pre-dawn trough
(00:00–05:59), per day_type, preserving each day_type's night-block total so the day/night
balance is unchanged — it touches only those hours, never a daytime peak, and brings every
bin's effective noise to ≲15% (at the observation σ floor).

Trips that cannot be placed or weighted are **dropped explicitly** (NaN start-hour /
day-of-week, or NaN W5 — ~0.04% lack a trip weight), rather than being silently skipped by a
weighted sum while still inflating the unweighted pooling count.

Normalisation to Σ=7
--------------------
With per-year rates A[dt,h] (weekday = Mon–Fri combined), the CSV columns are
  weekday-row(h) = 7·(A[0,h]/5)/T,  Sat(h) = 7·A[1,h]/T,  Sun(h) = 7·A[2,h]/T,
  T = Σ_h (A[0,h] + A[1,h] + A[2,h]),
so each component column sums to 7 over the 168 rows (W_c = 1) — the existing convention
(magnitude/split live in generation).

School (primary / post-primary / tertiary)
------------------------------------------
School is a **separate path** (`_school_level_grids`).  A school car trip is "when does the
escorting car leave" — a COUNT, so escort **ride-sharing matters** here (unlike the driveshare
*share*, where it cancels), and it needs the escort VEHICLE departure times, not the child's
passenger trip.  Per level the car-departure timing = ESCORT + SELF-DRIVE:
  - ESCORT: individual escort trips carry no child level, so — exactly as the school GENERATION
    rates do — a **household regression** attributes them: regress each household's escort-veh
    trip count on its preschool/primary/secondary/tertiary student counts (W2-weighted).  Run
    **per (day_type,hour) bin**, β_L(dt,h) is the per-student escort timing for level L,
    ride-share-correct (y counts cars) and level-resolved (primary-heavy households drive at
    primary bell times, etc.).  Off-peak regression noise (β<0, where school≈0) is clamped to 0.
    Preschool stays a regressor to avoid contaminating the others but its β is unused (→ retail).
  - SELF-DRIVE: the student's own purpose-4 veh trip, age→level (same `_level_shares`, incl. the
    DfE 16–18 split), mostly tertiary.
The two normalised shapes are blended by the generation escort/self-drive split
(`school_generation_rates.json`): primary ≈ post-primary ≈ escort double-peak (h8 / h15), post-
primary's AM peak a touch earlier, tertiary spread/college-like (self-drive-led, no sharp PM peak).
**Uniform ex-COVID pool** (not the adaptive per-bin path — the regression runs on the whole set;
the signal is thin + weekday-peak-concentrated; bell-time drift is ~nil, and 2yr-vs-9yr is stable
for primary/post-primary and only tertiary genuinely needs the pooling).  Same night smoothing
(harmless — school night ≈ 0).

`mean_fraction` / `std_fraction` (the aggregate columns used elsewhere, e.g. ingest_counts) are
preserved.  This completes the temporal migration — all six columns are microdata-derived, fully
retiring nts0502/0504 and the vestigial `purpose_rates`.

Usage:  python3 analysis/derive_component_profiles.py
Needs data/NTS (nts_microdata) + the DfE participation CSV (via derive_school_generation).
Re-run when the NTS microdata, the purpose mapping, or the school generation split change.
"""

import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nts_microdata as nts
from purpose_mapping import B01_COMPONENT

FRACS_FILE = "analysis/hourly_fractions.csv"
VEHICLE_MODE_CODES = [3, 5, 12]          # car/van driver, motorcycle, taxi
COMPONENTS = ["res", "commute", "retail"]   # adaptive-pooled path; school is derived separately below
MIN_BIN_N = 100                          # per-(day_type,hour) unweighted target (~10% noise)
# Expanding ex-COVID year tiers (England-only since 2012; COVID 2020-22 excluded).
YEAR_TIERS = [[2023, 2024], [2018, 2019], [2017, 2016, 2015, 2014, 2013]]
POOL_YEARS = [y for tier in YEAR_TIERS for y in tier]
DAY_NAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
# Light within-night smoothing: even pooled to the full window, the deepest pre-dawn hours
# stay thin on single-day (Sat/Sun) day-types (car trips at 3 am are genuinely rare).  That
# trough is flat and near-zero, so a small moving average over these hours borrows strength
# across adjacent hours with negligible bias — and it touches ONLY these hours, never the
# daytime peaks.  Each day-type's night-block total is preserved (only the within-night
# distribution is smoothed), so the day/night balance is unchanged.
NIGHT_HOURS = list(range(0, 6))          # 00:00–05:59 pre-dawn trough
NIGHT_SMOOTH_WIN = 3                     # moving-average width (±1); flat trough ⇒ low bias

# ── School temporal (separate path) ──────────────────────────────────────────
# Per-level car-departure timing = ESCORT (attributed to levels by the SAME household regression
# as the school generation rates, run per (day_type,hour) bin — β_L(dt,h) is per-student escort
# timing for level L, ride-share-correct because y counts cars, and level-resolved even though
# escort trips carry no child level) + SELF-DRIVE (purpose-4 veh, age→level), blended by the
# generation escort/self-drive split.  Uniform ex-COVID pool (the regression runs on the whole set;
# the signal is thin + weekday-peak-concentrated; school bell-time drift over 2013-24 is ~nil).
SCHOOL_LEVELS = ["primary", "postprimary", "tertiary"]
SCHOOL_LEVEL_KEY = {"primary": "prim", "postprimary": "sec", "tertiary": "ter"}  # -> _level_shares col
SCHOOL_POOL_YEARS = POOL_YEARS                       # uniform 9yr ex-COVID
SCHOOL_GEN_RATES = "analysis/school_generation_rates.json"


def _load():
    """Vehicle-driver trips with component, day_type (0=wkday,1=Sat,2=Sun), hour, weight."""
    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B01ID",
                                   "TripStartHours", "DayID", "W5", "JJXSC"], years=POOL_YEARS)
    day = nts.load("day", columns=["SurveyYear", "DayID", "TravelWeekDay_B01ID"], years=POOL_YEARS)
    tr = tr.merge(day[["DayID", "TravelWeekDay_B01ID"]], on="DayID", how="left")
    veh = tr[tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES)].copy()
    n0 = len(veh)
    # Explicit drop of trips we cannot place or weight (NaN start-hour / day-of-week, or NaN
    # weight — ~0.04% have a missing W5).  Loud, not a silent groupby.sum skip.
    veh = veh[veh.TripStartHours.notna() & veh.TravelWeekDay_B01ID.notna()
              & veh.W5.notna() & veh.JJXSC.notna()].copy()
    if len(veh) < n0:
        print(f"  dropped {n0 - len(veh):,} of {n0:,} vehicle trips with NaN "
              f"start-hour / day-of-week / weight")
    tr = veh
    dow = tr.TravelWeekDay_B01ID.astype(int) - 1                 # 0=Mon .. 6=Sun
    tr["dt"] = np.where(dow < 5, 0, np.where(dow == 5, 1, 2))     # day_type
    tr["h"] = tr.TripStartHours.astype(int)
    tr["w"] = tr.JJXSC * tr.W5
    tr["component"] = tr.TripPurpose_B01ID.map(B01_COMPONENT)
    return tr


def _bin_rate(dt, h, ncount_by_year, w_by_year):
    """Per-year weighted rate for one (dt,h) bin via adaptive tier expansion.

    Returns (rate, n_years, unweighted_n).  Expands the year window tier by tier until the
    cumulative unweighted count reaches MIN_BIN_N (or tiers are exhausted), then divides the
    cumulative weighted count by the number of years accumulated — the per-year normalisation
    that makes bins comparable across the windows they needed.
    """
    years, n, wsum = [], 0, 0.0
    for tier in YEAR_TIERS:
        for y in tier:
            years.append(y)
            n += ncount_by_year.get((y, dt, h), 0)
            wsum += w_by_year.get((y, dt, h), 0.0)
        if n >= MIN_BIN_N:
            break
    return wsum / len(years), len(years), n


def _smooth_night(grid):
    """Light within-night moving average over NIGHT_HOURS, per day_type, preserving each
    day_type's night-block total.  Only the flat pre-dawn trough is touched; daytime hours
    (incl. every peak) are returned unchanged."""
    out = grid.copy()
    half = NIGHT_SMOOTH_WIN // 2
    lo, hi = NIGHT_HOURS[0], NIGHT_HOURS[-1]
    for dt in range(grid.shape[0]):
        block = grid[dt, lo:hi + 1]
        sm = np.array([block[max(0, i - half):i + half + 1].mean() for i in range(len(block))])
        s = sm.sum()
        if s > 0:
            sm *= block.sum() / s                # preserve the night-block total
        out[dt, lo:hi + 1] = sm
    return out


def _write_grid(rows, col, grid):
    """Write a (3,24) day_type×hour rate grid into `col` of the 168 rows, normalised so the
    column sums to 7: weekday = grid[0]/5 (Mon–Fri pooled) into each of the 5 weekday rows,
    Sat = grid[1], Sun = grid[2], all ×7/T with T = grid.sum()."""
    T = grid.sum()
    if T <= 0:
        sys.exit(f"ERROR: zero total rate for {col}")
    wkday = 7.0 * (grid[0] / 5.0) / T
    sat = 7.0 * grid[1] / T
    sun = 7.0 * grid[2] / T
    for r in rows:
        dow = int(r["day_of_week"]); h = int(r["hour"].split(":")[0])
        v = wkday[h] if dow < 5 else (sat[h] if dow == 5 else sun[h])
        r[col] = f"{v:.10f}"


def _school_level_grids():
    """{level: (3,24) blended escort+self-drive car-departure rate grid} for the school columns,
    plus the (escort_share, sd_share) used.  See the SCHOOL constants block for the method."""
    import json
    import pandas as pd
    from derive_school_generation import (_england_secondary_share, _level_shares,
                                           P_ESCORT_EDU, P_EDUCATION, VEH)

    det = json.load(open(SCHOOL_GEN_RATES))["_meta"]["detail"]
    esc_share, sd_share = {}, {}
    for lvl in SCHOOL_LEVELS:
        e = det["escort"][lvl]; s = det["self_drive"][lvl]
        esc_share[lvl] = e / (e + s); sd_share[lvl] = s / (e + s)

    yrs = SCHOOL_POOL_YEARS
    sec_share = _england_secondary_share()
    ind = nts.load("individual", columns=["SurveyYear", "IndividualID", "HouseholdID",
                                          "Age_B01ID", "EducN_B01ID"], years=yrs)
    sh = ind.apply(lambda r: _level_shares(sec_share, r.Age_B01ID, r.EducN_B01ID),
                   axis=1, result_type="expand")
    sh.columns = ["prim", "sec", "ter"]
    ind = pd.concat([ind, sh], axis=1)
    ind["pre"] = ind.Age_B01ID.isin([2, 3]).astype(float)          # pre-school control (unused col)
    hh = nts.load("household", columns=["SurveyYear", "HouseholdID", "W2"], years=yrs)
    hk = (ind.groupby("HouseholdID")[["pre", "prim", "sec", "ter"]].sum()
          .join(hh.set_index("HouseholdID")["W2"]))

    tr = nts.load("trip", columns=["SurveyYear", "TripPurpose_B01ID", "MainMode_B04ID",
                                   "IndividualID", "HouseholdID", "TripStartHours", "DayID",
                                   "W5", "JJXSC"], years=yrs)
    day = nts.load("day", columns=["SurveyYear", "DayID", "TravelWeekDay_B01ID"], years=yrs)
    tr = tr.merge(day[["DayID", "TravelWeekDay_B01ID"]], on="DayID", how="left")
    tr = tr[tr.MainMode_B04ID.isin(VEH) & tr.TripStartHours.notna()
            & tr.TravelWeekDay_B01ID.notna() & tr.JJXSC.notna() & tr.W5.notna()].copy()
    dow = tr.TravelWeekDay_B01ID.astype(int) - 1
    tr["dt"] = np.where(dow < 5, 0, np.where(dow == 5, 1, 2))
    tr["h"] = tr.TripStartHours.astype(int)
    tr["bin"] = tr.dt * 24 + tr.h

    # escort: per-bin household regression (same X/weighting as the generation escort β, per bin)
    D = hk[(hk[["pre", "prim", "sec", "ter"]].sum(axis=1) > 0) & hk.W2.notna()]
    sw = np.sqrt(D["W2"].values)
    Xw = D[["pre", "prim", "sec", "ter"]].values * sw[:, None]
    esc = tr[tr.TripPurpose_B01ID == P_ESCORT_EDU]
    P = (esc.pivot_table(index="HouseholdID", columns="bin", values="JJXSC",
                         aggfunc="sum", fill_value=0.0)
         .reindex(index=D.index, columns=range(72), fill_value=0.0))
    B, *_ = np.linalg.lstsq(Xw, P.values * sw[:, None], rcond=None)          # (4, 72)
    beta = {lvl: np.maximum(B[i], 0.0).reshape(3, 24)                        # clamp off-peak noise ≥0
            for lvl, i in [("primary", 1), ("postprimary", 2), ("tertiary", 3)]}

    # self-drive: student's own purpose-4 veh trip, age→level (fractional level share × trip weight)
    sd = tr[tr.TripPurpose_B01ID == P_EDUCATION].copy()
    sd["w"] = sd.JJXSC * sd.W5
    sdi = ind.set_index("IndividualID")
    sdgrid = {}
    for lvl in SCHOOL_LEVELS:
        sd["sw"] = sd.IndividualID.map(sdi[SCHOOL_LEVEL_KEY[lvl]]).fillna(0.0) * sd.w
        sdgrid[lvl] = (sd.groupby(["dt", "h"]).sw.sum()
                       .reindex([(d, h) for d in range(3) for h in range(24)], fill_value=0.0)
                       .values.reshape(3, 24))

    # blend the two normalised shapes by the generation escort/sd split, then smooth the night
    grids = {}
    for lvl in SCHOOL_LEVELS:
        eg, sg = beta[lvl], sdgrid[lvl]
        eg = eg / eg.sum() if eg.sum() > 0 else eg
        sg = sg / sg.sum() if sg.sum() > 0 else sg
        grids[lvl] = _smooth_night(esc_share[lvl] * eg + sd_share[lvl] * sg)
    return grids, esc_share, sd_share


def main():
    print(f"Deriving car-specific temporal profiles from NTS microdata "
          f"(pool tiers {YEAR_TIERS}, MIN_BIN_N={MIN_BIN_N}) …")
    tr = _load()
    # Pre-aggregate per (year, dt, hour): unweighted n and Σw, per component.
    A = {}          # component -> (3,24) per-year rate grid
    diag = {}       # component -> (windows list, final n list)
    for comp in COMPONENTS:
        sub = tr[tr.component == comp]
        g = sub.groupby(["SurveyYear", "dt", "h"])
        ncount = g.size().to_dict()
        wsum = g.w.sum().to_dict()
        grid = np.zeros((3, 24))
        wins, ns = [], []
        for dt in range(3):
            for h in range(24):
                rate, nyr, nn = _bin_rate(dt, h, ncount, wsum)
                grid[dt, h] = rate
                wins.append(nyr)
                ns.append(nn)
        A[comp] = _smooth_night(grid)                # light pre-dawn smoothing
        diag[comp] = (np.array(wins), np.array(ns))

    # ── School temporal (separate path: regression escort + age→level self-drive) ──
    print("Deriving school temporal (per-bin household regression escort + age→level self-drive) …")
    school_grids, esc_share, sd_share = _school_level_grids()

    # ── Read the CSV, rewrite all six shape columns, preserve the rest (mean_fraction etc.) ──
    with open(FRACS_FILE, newline="") as f:
        rows = list(csv.DictReader(f))
        fieldnames = list(rows[0].keys())
    out_cols = ([f"mean_fraction_{c}" for c in COMPONENTS]
                + [f"mean_fraction_school_{lvl}" for lvl in SCHOOL_LEVELS])
    for col in out_cols:
        if col not in fieldnames:
            sys.exit(f"ERROR: column {col} missing from {FRACS_FILE}")

    for comp in COMPONENTS:
        _write_grid(rows, f"mean_fraction_{comp}", A[comp])
    for lvl in SCHOOL_LEVELS:
        _write_grid(rows, f"mean_fraction_school_{lvl}", school_grids[lvl])

    with open(FRACS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote res/commute/retail + school primary/postprimary/tertiary shape columns → {FRACS_FILE}")

    # ── Diagnostics + verification ───────────────────────────────────────────
    print("\nAdaptive pooling — bins (of 72) per window size, and final unweighted n:")
    print(f"  {'comp':8s} {'2yr':>4} {'4yr':>4} {'9yr':>4} | {'min_n':>6} {'median_n':>9} "
          f"{'#<MIN':>6}  daytime(wkday 7-19) min_n")
    for comp in COMPONENTS:
        wins, ns = diag[comp]
        hist = {2: int((wins == 2).sum()), 4: int((wins == 4).sum()), 9: int((wins == 9).sum())}
        # daytime weekday bins = dt 0 (first 24), hours 7..19
        dt0 = ns[:24]
        print(f"  {comp:8s} {hist[2]:>4} {hist[4]:>4} {hist[9]:>4} | {ns.min():>6} "
              f"{int(np.median(ns)):>9} {int((ns < MIN_BIN_N).sum()):>6}  {dt0[7:20].min():>6}")

    hrs = [6, 7, 8, 9, 12, 14, 15, 16, 17, 18]
    print("\nWeekday hourly shape % (sanity):")
    print("  component        " + " ".join(f"{h:4d}" for h in hrs))
    for comp in COMPONENTS:
        g = A[comp][0]; g = g / g.sum()
        print(f"  {comp:15s}" + " ".join(f"{100*g[h]:4.1f}" for h in hrs))
    for lvl in SCHOOL_LEVELS:
        g = school_grids[lvl][0]; g = g / g.sum()
        print(f"  school_{lvl:8s}" + " ".join(f"{100*g[h]:4.1f}" for h in hrs)
              + f"   [escort {esc_share[lvl]:.2f} / self-drive {sd_share[lvl]:.2f}]")

    print("\nVerification (each shape column Σ over 168 rows → 7):")
    ok = True
    for col in out_cols:
        s = sum(float(r[col]) for r in rows)
        neg = any(float(r[col]) < 0 for r in rows)
        flag = "" if abs(s - 7) < 1e-6 and not neg else "  <-- FAIL"
        ok = ok and not flag
        print(f"  {col:32s} Σ_168 = {s:.6f}{'  (has negative!)' if neg else ''}{flag}")
    print("  ✓ all checks pass" if ok else "  ✗ CHECK FAILED")


if __name__ == "__main__":
    main()
