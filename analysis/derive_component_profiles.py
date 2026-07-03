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
res weekend-night bins cannot reach MIN_BIN_N even on the full window; they use the widest
window (their best available).

Normalisation to Σ=7
--------------------
With per-year rates A[dt,h] (weekday = Mon–Fri combined), the CSV columns are
  weekday-row(h) = 7·(A[0,h]/5)/T,  Sat(h) = 7·A[1,h]/T,  Sun(h) = 7·A[2,h]/T,
  T = Σ_h (A[0,h] + A[1,h] + A[2,h]),
so each component column sums to 7 over the 168 rows (W_c = 1) — the existing convention
(magnitude/split live in generation).

Scope / staging
---------------
This derives res / commute / retail.  The three school columns are **carried forward
untouched** (they stay on the previous method until the school temporal piece, which must
handle escort ride-sharing — a school car trip is "when does the escorting car leave", a
COUNT, so ride-sharing matters, unlike the driveshare share).  mean_fraction / std_fraction
(the aggregate columns used elsewhere, e.g. ingest_counts) are also preserved.

Usage:  python3 analysis/derive_component_profiles.py
Re-run when the NTS microdata or the purpose mapping change.
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
COMPONENTS = ["res", "commute", "retail"]   # school carried forward (pending its own piece)
MIN_BIN_N = 100                          # per-(day_type,hour) unweighted target (~10% noise)
# Expanding ex-COVID year tiers (England-only since 2012; COVID 2020-22 excluded).
YEAR_TIERS = [[2023, 2024], [2018, 2019], [2017, 2016, 2015, 2014, 2013]]
POOL_YEARS = [y for tier in YEAR_TIERS for y in tier]
DAY_NAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def _load():
    """Vehicle-driver trips with component, day_type (0=wkday,1=Sat,2=Sun), hour, weight."""
    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B01ID",
                                   "TripStartHours", "DayID", "W5", "JJXSC"], years=POOL_YEARS)
    day = nts.load("day", columns=["SurveyYear", "DayID", "TravelWeekDay_B01ID"], years=POOL_YEARS)
    tr = tr.merge(day[["DayID", "TravelWeekDay_B01ID"]], on="DayID", how="left")
    tr = tr[tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES)
            & tr.TripStartHours.notna() & tr.TravelWeekDay_B01ID.notna()].copy()
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
        A[comp] = grid
        diag[comp] = (np.array(wins), np.array(ns))

    # ── Read the CSV, rewrite res/commute/retail, carry everything else forward ──
    with open(FRACS_FILE, newline="") as f:
        rows = list(csv.DictReader(f))
        fieldnames = rows[0].keys()
    for comp in COMPONENTS:
        if f"mean_fraction_{comp}" not in fieldnames:
            sys.exit(f"ERROR: column mean_fraction_{comp} missing from {FRACS_FILE}")
    for lvl in ("primary", "postprimary", "tertiary"):
        if f"mean_fraction_school_{lvl}" not in fieldnames:
            sys.exit(f"ERROR: school column missing from {FRACS_FILE} — school columns are "
                     f"carried forward and must already exist (see module docstring / staging)")

    # Normalise each component to Σ_168 = 7 and write the per-dow rows.
    for comp in COMPONENTS:
        grid = A[comp]
        T = grid.sum()                                   # Σ over the 3 day-types × 24 h
        if T <= 0:
            sys.exit(f"ERROR: zero total rate for {comp}")
        wkday_row = 7.0 * (grid[0] / 5.0) / T             # one weekday (grid[0] pools Mon–Fri)
        sat_row = 7.0 * grid[1] / T
        sun_row = 7.0 * grid[2] / T
        col = f"mean_fraction_{comp}"
        for r in rows:
            dow = int(r["day_of_week"]); h = int(r["hour"].split(":")[0])
            v = wkday_row[h] if dow < 5 else (sat_row[h] if dow == 5 else sun_row[h])
            r[col] = f"{v:.10f}"

    with open(FRACS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote res/commute/retail shape columns → {FRACS_FILE} "
          f"(school columns carried forward)")

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

    print("\nWeekday hourly shape %, h6-9,12,16-18 (sanity):")
    print("  comp     " + "  ".join(f"{h:4d}" for h in [6, 7, 8, 9, 12, 16, 17, 18]))
    for comp in COMPONENTS:
        g = A[comp][0]; g = g / g.sum()
        print(f"  {comp:8s}" + "  ".join(f"{100*g[h]:4.1f}" for h in [6, 7, 8, 9, 12, 16, 17, 18]))

    print("\nVerification (each component column Σ over 168 rows → 7):")
    ok = True
    for comp in COMPONENTS:
        col = f"mean_fraction_{comp}"
        s = sum(float(r[col]) for r in rows)
        neg = any(float(r[col]) < 0 for r in rows)
        flag = "" if abs(s - 7) < 1e-6 and not neg else "  <-- FAIL"
        ok = ok and not flag
        print(f"  {comp:8s} Σ_168 = {s:.6f}{'  (has negative!)' if neg else ''}{flag}")
    print("  ✓ all checks pass" if ok else "  ✗ CHECK FAILED")


if __name__ == "__main__":
    main()
