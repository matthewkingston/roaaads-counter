"""
Derive residential and business temporal profile priors for the two-component
gravity model from NTS trip-purpose-by-start-time data.

Source
------
DfT National Travel Survey Table NTS0502a, "Trip start time by trip purpose
(Monday to Friday only): England, 2002 onwards".
Downloaded from:
  https://www.gov.uk/government/statistical-data-sets/nts05-trips
File: data/nts0502.ods

What the table gives
--------------------
Each row is one (year-range, start-hour) combination.  The percentage columns
give P(purpose | start hour, Mon–Fri), i.e. for all trips that BEGIN in that
hour on a weekday, what fraction belong to each trip purpose.  Rows sum to 100
under the "All purposes" column.

Purpose classification
----------------------
The gravity model's business demand nodes represent workplaces, retail POIs,
and car parks (public and private).  Trips attracted to these nodes are:

  Business ("biz"):
    - Commuting              — direct work trips (largest component)
    - Employer's business    — trips made during/for work (deliveries, etc.)
    - Education              — to/from school or university; schools are
                               tagged in OSM and receive a bonus POI weight
    - Escort education       — taking someone to school (school-run); timed
                               identically to education trips and attracted to
                               the same nodes
    - Shopping               — retail trips; shops are prominent in OSM POIs
                               and shopping car parks dominate the car-park area

  Residential ("res"):
    - Other work, other escort and personal business  (banks, GP, etc.)
    - Visiting friends, entertainment and sport
    - Holiday, day trip and other

  Rationale for "escort education" and "education" as business:
    Schools appear as POIs in the network and are given a bonus weighting in
    build_demographics.py.  School-run trips (escort education) have the same
    temporal signature and destination type as education trips themselves.

Year choice
-----------
"2023 to 2024" — the most recent rolling average in the file.  The individual
years 2020 and 2021 are excluded because COVID reduced commuting anomalously;
2022 is the first post-pandemic year and "2023 to 2024" is more stable.

Weekend handling
----------------
NTS0502a covers Monday to Friday only.  No hour-by-hour purpose split is
available for Saturday or Sunday in this file.

We apply a FLAT business share to all weekend hours, derived as the
unweighted mean of the weekday NTS business shares over hours h10–h14:

  p_biz_weekend = mean( biz_share_weekday(h) for h in 10..14 )

Rationale: on weekends commuting is low and school runs are absent.  The
dominant "business" activity is shopping.  Hours h10–h14 on weekdays are
dominated by shopping and have minimal school-run contamination (education +
escort_education ≈ 2–4% each, vs 26% + 22% at h15 and 16% + 12% at h14).
  - h15 is excluded: it carries the 3 pm school pickup (NTS shows education
    26.2% + escort_edu 22.0% = 48 pp at h15), absent at weekend.
  - h14 is retained despite some school contamination (~16 pp) because the
    non-school business share there (commuting + employer + shopping ≈ 34%)
    still represents plausible weekend afternoon activity.
The volume-weighted weekday mean (~52%) was rejected as it is dominated by
the h07–h08 school run, absent at weekends.
Without dedicated weekend NTS trip-purpose data this remains a documented
approximation; the tuner can move the weekend f_biz away from the prior.

The temporal SHAPE for weekends comes from the existing aggregate mean_fraction
profile (which already varies by hour), so the profiles still rise and fall
through the day — only the overall scaling is estimated from weekdays.

Constraint
----------
The tuner requires mean_fraction_res + mean_fraction_biz = 2 × mean_fraction
at every (day_of_week, hour) slot (the aggregate coupling uses 2×mean_fraction
as the target sum).  This is enforced exactly:

  mean_fraction_biz = 2 × mean_fraction × biz_share
  mean_fraction_res = 2 × mean_fraction × (1 − biz_share)

Neither component necessarily sums to 1 across hours, but that is expected and
correct: the global scale K absorbs the normalisation difference.

Usage
-----
  python3 analysis/derive_component_profiles.py

Overwrites the mean_fraction_res and mean_fraction_biz columns in
analysis/hourly_fractions.csv.  All other columns are preserved unchanged.
Re-run whenever the NTS source file changes or the purpose classification is
revised.
"""

import csv, sys
import pandas as pd
import numpy as np

NTS_FILE   = "data/nts0502.ods"
FRACS_FILE = "analysis/hourly_fractions.csv"
NTS_YEAR   = "2023 to 2024"

# ── Load NTS0502a ─────────────────────────────────────────────────────────────

print(f"Loading {NTS_FILE} …")
raw = pd.read_excel(NTS_FILE, sheet_name="NTS0502a_start_time_by_purpose",
                    header=None, engine="odf")

# Row 5 is the column header; data starts at row 6.
# Columns (by position):
#  0  Year
#  1  Start time     (e.g. "0700 to 0759")
#  2  Commuting (%)
#  3  Business (%)                          ← employer's business
#  4  Education (%)
#  5  Escort education (%)
#  6  Shopping (%)
#  7  Other work, other escort and personal business (%)
#  8  Visiting friends, entertainment and sport (%)
#  9  Holiday, day trip and other (%)
# 10  All purposes (%)                      ← should equal 100
# 11  Unweighted sample size

BIZ_COLS = [2, 3, 4, 5, 6]   # commuting, employer_biz, education, escort_edu, shopping
ALL_COL  = 10

data = raw.iloc[6:].copy()
data.columns = range(data.shape[1])

# Filter to chosen year, drop the summary "All day" row (where col 1 == "All day")
year_rows = data[data[0] == NTS_YEAR].copy()
year_rows = year_rows[year_rows[1] != "All day"].reset_index(drop=True)

if len(year_rows) != 24:
    print(f"ERROR: expected 24 hourly rows for '{NTS_YEAR}', got {len(year_rows)}")
    sys.exit(1)

# Parse hour index (0–23) from start-time string "0700 to 0759"
year_rows["hour"] = year_rows[1].str[:2].astype(int)
year_rows = year_rows.sort_values("hour").reset_index(drop=True)

# Coerce numeric (handles any "[low]" strings → 0)
for c in BIZ_COLS + [ALL_COL]:
    year_rows[c] = pd.to_numeric(year_rows[c], errors="coerce").fillna(0.0)

year_rows["biz_share"] = year_rows[BIZ_COLS].sum(axis=1) / year_rows[ALL_COL]

# Dictionary hour → weekday business share
biz_share_by_hour = dict(zip(year_rows["hour"], year_rows["biz_share"]))

print(f"\nNTS {NTS_YEAR} weekday business share by hour")
print(f"  (business = commuting + employer_biz + education + escort_edu + shopping)")
print(f"  {'Hour':>4}  {'biz%':>6}  {'res%':>6}")
for h in range(24):
    bs = biz_share_by_hour[h]
    print(f"  h{h:02d}    {100*bs:5.1f}%  {100*(1-bs):5.1f}%")

# ── Load hourly fractions ─────────────────────────────────────────────────────

print(f"\nLoading {FRACS_FILE} …")
rows = []
with open(FRACS_FILE, newline="") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

# ── Weekend business share: mean over core shopping-window hours h10–h16 ──────
# Weekday h10–h16 has minimal commute/school-run contamination and is the best
# NTS-grounded proxy for weekend business activity (primarily shopping).
WEEKEND_HOURS = range(10, 15)   # 10:00–14:59; h15 excluded (school pickup)
p_biz_weekend = sum(biz_share_by_hour[h] for h in WEEKEND_HOURS) / len(WEEKEND_HOURS)

print(f"\nWeekend flat business share (mean of h10–h16 weekday NTS shares): {p_biz_weekend:.4f}")
print(f"  Hours used: " + ", ".join(f"h{h:02d}={100*biz_share_by_hour[h]:.1f}%" for h in WEEKEND_HOURS))
print(f"  Residual (residential): {1 - p_biz_weekend:.4f}")

# ── Update component columns ──────────────────────────────────────────────────

for row in rows:
    dow = int(row["day_of_week"])
    h   = int(row["hour"].split(":")[0])
    mfa = float(row["mean_fraction"])

    bs = biz_share_by_hour[h] if dow <= 4 else p_biz_weekend

    # Enforce mean_fraction_biz + mean_fraction_res = 2 × mean_fraction exactly.
    row["mean_fraction_biz"] = f"{2 * mfa * bs:.10f}"
    row["mean_fraction_res"] = f"{2 * mfa * (1 - bs):.10f}"

# ── Write back, preserving column order ──────────────────────────────────────

# Keep original columns, replacing or appending the component columns at the end.
base_cols = [c for c in fieldnames
             if c not in ("mean_fraction_res", "mean_fraction_biz")]
out_cols = base_cols + ["mean_fraction_res", "mean_fraction_biz"]

with open(FRACS_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=out_cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print(f"\nWrote updated component columns → {FRACS_FILE}")

# ── Sanity checks ─────────────────────────────────────────────────────────────

print("\nSanity checks:")

# 1. Sum constraint per row
errors = 0
for row in rows:
    mfa = float(row["mean_fraction"])
    mfr = float(row["mean_fraction_res"])
    mfb = float(row["mean_fraction_biz"])
    diff = abs(mfr + mfb - 2 * mfa)
    if diff > 1e-9:
        print(f"  FAIL: dow={row['day_of_week']} h={row['hour']} "
              f"res+biz={mfr+mfb:.10f} vs 2×mfa={2*mfa:.10f} (Δ={diff:.2e})")
        errors += 1
if errors == 0:
    print("  ✓  res + biz = 2 × agg for all 168 rows")

# 2. All values non-negative
neg = [(r["day_of_week"], r["hour"]) for r in rows
       if float(r["mean_fraction_res"]) < 0 or float(r["mean_fraction_biz"]) < 0]
if neg:
    print(f"  FAIL: negative fractions at {neg}")
else:
    print("  ✓  all fractions non-negative")

# 3. Show weekday biz profile (a useful spot-check)
print("\n  Weekday biz fraction profile (mfb / 2*mfa = biz_share):")
wkday_rows = sorted([r for r in rows if int(r["day_of_week"]) == 0],
                    key=lambda r: int(r["hour"].split(":")[0]))
for r in wkday_rows:
    h   = int(r["hour"].split(":")[0])
    mfa = float(r["mean_fraction"])
    mfb = float(r["mean_fraction_biz"])
    bs  = mfb / (2 * mfa) if mfa > 0 else 0
    bar = "█" * int(bs * 40)
    print(f"  Mon h{h:02d}  {100*bs:4.1f}%  {bar}")
