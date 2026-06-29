"""
Derive residential, business, and school temporal profile priors for the
three-component gravity model from NTS trip-purpose data.

Sources
-------
Weekdays — DfT NTS Table NTS0502a, "Trip start time by trip purpose
  (Monday to Friday only): England, 2002 onwards".
  Gives P(purpose | start hour, Mon–Fri) as percentages summing to 100.
  File: data/nts0502.ods

Weekends — DfT NTS Table NTS0504b, "Average number of trips by day of the
  week and purpose (trips per person per year): England, 2002 onwards".
  Gives absolute trip rates per person per year by day and purpose.
  Used here for Saturday and Sunday rows only.
  File: data/nts0504.ods

Both files downloaded from:
  https://www.gov.uk/government/statistical-data-sets/nts05-trips

Purpose classification
----------------------
The gravity model splits work-/shopping-related travel into two independent
components — commute (attractor = workplace jobs) and retail (attractor = retail
parking spaces) — each with its own temporal profile.  Trips are classified:

  Commute ("commute"):
    - Commuting              — direct work trips
    - Employer's business    — trips made during/for work (deliveries, etc.)

  Retail ("retail"):
    - Shopping               — trips to retail; shops dominate OSM POIs and
                               car parks

  School ("school"):
    - Education              — to/from school or university; DOWNWEIGHTED
                               × EDU_FACTOR (see below)
    - Escort education       — taking someone to school; adult car driver

  Residential ("res"):
    - Other escort and personal business  (GP, bank, etc.)
    - Visiting friends, entertainment and sport
    - Holiday, day trip and other

Education downweighting (EDU_FACTOR = 1/5)
------------------------------------------
NTS records ALL trips regardless of mode.  The "Education" category is
predominantly children walking, cycling, or being bussed — not car trips.
The "Escort education" category captures the adult car driver making the
school-run trip. The child's own journey ("Education") has far lower
per-person car trip generation than commuting.

We therefore downweight the Education contribution by EDU_FACTOR = 1/5
in the school numerator, and reduce the denominator by the same amount
so that shares of all remaining purposes are proportionally increased.

  commute_numerator = commuting + employer_business
  retail_numerator  = shopping
  school_numerator  = education*EDU_FACTOR + escort_education
  adjusted_total    = all_purposes - education*(1 - EDU_FACTOR)
  commute_share     = commute_numerator / adjusted_total
  retail_share      = retail_numerator  / adjusted_total
  school_share      = school_numerator  / adjusted_total
  res_share         = 1 - commute_share - retail_share - school_share

For NTS0504b weekends, "Just walk" trips are additionally removed from
the denominator before this reweighting, because the aggregate traffic
count profiles we are matching contain no pedestrian trips.

Year choice
-----------
"2023 to 2024" — the most recent rolling average in both files.  Individual
years 2020–2021 are avoided as COVID anomalies; 2022 is the first
post-pandemic year; the 2-year rolling average gives a more stable estimate.

Weekday approach (NTS0502a)
---------------------------
For each hour h, commute_share(h) and retail_share(h) are computed from the
2023–2024 per-hour purpose percentages with the education downweighting applied.
This gives distinct shares for every hour of the day.

Weekend approach (NTS0504b)
---------------------------
NTS0502a covers Monday–Friday only.  NTS0504b provides absolute trip rates
per person per year for each day of the week, including Saturday and Sunday.
We extract the Saturday and Sunday rows from the 2023–2024 rolling average,
remove "Just walk", apply the education downweighting, and compute flat
commute/retail shares for each weekend day.  The temporal SHAPE within each
weekend day is inherited from the aggregate mean_fraction profile (which
already captures the Saturday vs. Sunday difference in volume pattern);
only the overall component split is estimated from NTS0504b.

Constraint
----------
The tuner requires
  mean_fraction_res + mean_fraction_commute + mean_fraction_retail
    + mean_fraction_school = mean_fraction
at every (day_of_week, hour) slot.  This is enforced exactly:

  mean_fraction_commute = mean_fraction × commute_share
  mean_fraction_retail  = mean_fraction × retail_share
  mean_fraction_school  = mean_fraction × school_share
  mean_fraction_res     = mean_fraction × (1 − commute − retail − school shares)

Each component is directly interpretable as the fraction of AADT attributable
to that traffic type at that slot.  The global scale K absorbs normalisation.

Usage
-----
  python3 analysis/derive_component_profiles.py

Overwrites the mean_fraction_res, mean_fraction_commute, mean_fraction_retail,
and mean_fraction_school columns in analysis/hourly_fractions.csv.  All other
columns are preserved.
Re-run whenever the NTS files change or the purpose classification is revised.
"""

import csv, sys
import pandas as pd

NTS0502_FILE = "data/nts0502.ods"
NTS0504_FILE = "data/nts0504.ods"
FRACS_FILE   = "analysis/hourly_fractions.csv"
NTS_YEAR     = "2023 to 2024"

EDU_FACTOR   = 1 / 5   # education trips have ~1/5 the car trip generation of commuting

# ── Load NTS0502a (weekday hourly purpose split) ──────────────────────────────

print(f"Loading {NTS0502_FILE} …")
raw502 = pd.read_excel(NTS0502_FILE, sheet_name="NTS0502a_start_time_by_purpose",
                       header=None, engine="odf")

# Row 5 = header; data from row 6.
# Column positions (0-indexed):
#  0  Year       1  Start time
#  2  Commuting  3  Business (employer's)  4  Education  5  Escort education
#  6  Shopping   7  Other…  8  Visiting…  9  Holiday…
# 10  All purposes (%)     11  Sample size
C_COMM, C_BIZ_E, C_EDU, C_ESCORT, C_SHOP = 2, 3, 4, 5, 6
C_ALL = 10

data502 = raw502.iloc[6:].copy()
data502.columns = range(data502.shape[1])

yr502 = data502[data502[0] == NTS_YEAR].copy()
yr502 = yr502[yr502[1] != "All day"].reset_index(drop=True)
if len(yr502) != 24:
    print(f"ERROR: expected 24 rows for '{NTS_YEAR}' in NTS0502a, got {len(yr502)}")
    sys.exit(1)

yr502["hour"] = yr502[1].str[:2].astype(int)
yr502 = yr502.sort_values("hour").reset_index(drop=True)

for c in [C_COMM, C_BIZ_E, C_EDU, C_ESCORT, C_SHOP, C_ALL]:
    yr502[c] = pd.to_numeric(yr502[c], errors="coerce").fillna(0.0)

# Education downweighting: remove (1-EDU_FACTOR) × education from denominator
yr502["adj_total"]     = yr502[C_ALL] - yr502[C_EDU] * (1 - EDU_FACTOR)
yr502["commute_num"]   = yr502[C_COMM] + yr502[C_BIZ_E]
yr502["retail_num"]    = yr502[C_SHOP]
yr502["school_num"]    = yr502[C_EDU] * EDU_FACTOR + yr502[C_ESCORT]
yr502["commute_share"] = yr502["commute_num"] / yr502["adj_total"]
yr502["retail_share"]  = yr502["retail_num"]  / yr502["adj_total"]
yr502["school_share"]  = yr502["school_num"]  / yr502["adj_total"]

commute_share_weekday = dict(zip(yr502["hour"], yr502["commute_share"]))
retail_share_weekday  = dict(zip(yr502["hour"], yr502["retail_share"]))
school_share_weekday  = dict(zip(yr502["hour"], yr502["school_share"]))

print(f"\nNTS {NTS_YEAR} weekday component shares by hour  (education × {EDU_FACTOR})")
print(f"  {'Hour':>4}  {'com%':>6}  {'ret%':>6}  {'sch%':>6}  {'res%':>6}")
for h in range(24):
    cs = commute_share_weekday[h]
    ts = retail_share_weekday[h]
    ss = school_share_weekday[h]
    print(f"  h{h:02d}    {100*cs:5.1f}%  {100*ts:5.1f}%  {100*ss:5.1f}%  {100*(1-cs-ts-ss):5.1f}%")

# ── Load NTS0504b (day-of-week purpose rates, trips/person/year) ──────────────

print(f"\nLoading {NTS0504_FILE} …")
raw504 = pd.read_excel(NTS0504_FILE, sheet_name="NTS0504b_day_purpose",
                       header=None, engine="odf")

# Row 5 = header; data from row 6.
# Column positions (0-indexed):
#  0  Year          1  Day of the week
#  2  Commuting     3  Business (employer's)   4  Education   5  Escort education
#  6  Shopping      7  Other escort            8  Personal business
#  9  Visit friends at private home           10  Visit friends elsewhere
# 11  Sport/entertainment                     12  Holiday or day trip
# 13  Just walk                               14  Other
# 15  All purposes                            16  Sample size
C504_COMM, C504_BIZ_E, C504_EDU, C504_ESCORT, C504_SHOP = 2, 3, 4, 5, 6
C504_WALK = 13
C504_ALL  = 15

data504 = raw504.iloc[6:].copy()
data504.columns = range(data504.shape[1])

yr504 = data504[data504[0] == NTS_YEAR].copy()
for c in [C504_COMM, C504_BIZ_E, C504_EDU, C504_ESCORT,
          C504_SHOP, C504_WALK, C504_ALL]:
    yr504[c] = pd.to_numeric(yr504[c], errors="coerce").fillna(0.0)

commute_share_weekend = {}
retail_share_weekend  = {}
school_share_weekend  = {}
print(f"\nNTS {NTS_YEAR} weekend component shares from NTS0504b  (education × {EDU_FACTOR}, walk removed)")
print(f"  {'Day':<10}  {'com_trips':>9}  {'ret_trips':>9}  {'sch_trips':>9}  {'adj_total':>9}"
      f"  {'com%':>6}  {'ret%':>6}  {'sch%':>6}  {'res%':>6}")

for day in ["Saturday", "Sunday"]:
    row = yr504[yr504[1] == day]
    if row.empty:
        print(f"ERROR: '{day}' not found in NTS0504b for year '{NTS_YEAR}'")
        sys.exit(1)
    row = row.iloc[0]

    comm    = float(row[C504_COMM])
    biz_e   = float(row[C504_BIZ_E])
    edu     = float(row[C504_EDU])
    escort  = float(row[C504_ESCORT])
    shop    = float(row[C504_SHOP])
    walk    = float(row[C504_WALK])
    all_p   = float(row[C504_ALL])

    # Remove walking trips (not in vehicle counts), then downweight education
    no_walk     = all_p - walk
    adj_total   = no_walk - edu * (1 - EDU_FACTOR)
    commute_num = comm + biz_e
    retail_num  = shop
    school_num  = edu * EDU_FACTOR + escort

    cs = commute_num / adj_total
    ts = retail_num  / adj_total
    ss = school_num  / adj_total
    commute_share_weekend[day] = cs
    retail_share_weekend[day]  = ts
    school_share_weekend[day]  = ss

    print(f"  {day:<10}  {commute_num:9.3f}  {retail_num:9.3f}  {school_num:9.3f}  {adj_total:9.3f}"
          f"  {100*cs:5.1f}%  {100*ts:5.1f}%  {100*ss:5.1f}%  {100*(1-cs-ts-ss):5.1f}%")

# ── Load hourly fractions CSV ─────────────────────────────────────────────────

print(f"\nLoading {FRACS_FILE} …")
rows = []
with open(FRACS_FILE, newline="") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        rows.append(row)

# ── Compute and write component columns ──────────────────────────────────────

DOW_TO_DAY = {5: "Saturday", 6: "Sunday"}

for row in rows:
    dow = int(row["day_of_week"])
    h   = int(row["hour"].split(":")[0])
    mfa = float(row["mean_fraction"])

    if dow <= 4:
        cs = commute_share_weekday[h]
        ts = retail_share_weekday[h]
        ss = school_share_weekday[h]
    else:
        cs = commute_share_weekend[DOW_TO_DAY[dow]]
        ts = retail_share_weekend[DOW_TO_DAY[dow]]
        ss = school_share_weekend[DOW_TO_DAY[dow]]

    row["mean_fraction_commute"] = f"{mfa * cs:.10f}"
    row["mean_fraction_retail"]  = f"{mfa * ts:.10f}"
    row["mean_fraction_school"]  = f"{mfa * ss:.10f}"
    row["mean_fraction_res"]     = f"{mfa * (1 - cs - ts - ss):.10f}"

_comp_cols = ("mean_fraction_res", "mean_fraction_commute",
              "mean_fraction_retail", "mean_fraction_school",
              "mean_fraction_biz")   # drop legacy biz column if present
base_cols = [c for c in fieldnames if c not in _comp_cols]
out_cols = base_cols + ["mean_fraction_res", "mean_fraction_commute",
                        "mean_fraction_retail", "mean_fraction_school"]

with open(FRACS_FILE, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=out_cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print(f"\nWrote updated component columns → {FRACS_FILE}")

# ── Sanity checks ─────────────────────────────────────────────────────────────

print("\nSanity checks:")

errors = 0
for row in rows:
    mfa = float(row["mean_fraction"])
    mfr = float(row["mean_fraction_res"])
    mfc = float(row["mean_fraction_commute"])
    mft = float(row["mean_fraction_retail"])
    mfs = float(row["mean_fraction_school"])
    diff = abs(mfr + mfc + mft + mfs - mfa)
    if diff > 1e-9:
        print(f"  FAIL: dow={row['day_of_week']} h={row['hour']} "
              f"res+com+ret+sch={mfr+mfc+mft+mfs:.10f} vs mfa={mfa:.10f} (Δ={diff:.2e})")
        errors += 1
if errors == 0:
    print("  ✓  res + commute + retail + school = agg for all 168 rows")

neg = [(r["day_of_week"], r["hour"]) for r in rows
       if (float(r["mean_fraction_res"]) < 0 or float(r["mean_fraction_commute"]) < 0
           or float(r["mean_fraction_retail"]) < 0 or float(r["mean_fraction_school"]) < 0)]
if neg:
    print(f"  FAIL: negative fractions at {neg}")
else:
    print("  ✓  all fractions non-negative")

print("\n  Weekday component shares by hour (Mon shown; other weekdays identical):")
wkday_rows = sorted([r for r in rows if int(r["day_of_week"]) == 0],
                    key=lambda r: int(r["hour"].split(":")[0]))
for r in wkday_rows:
    h   = int(r["hour"].split(":")[0])
    mfa = float(r["mean_fraction"])
    mfc = float(r["mean_fraction_commute"])
    mft = float(r["mean_fraction_retail"])
    mfs = float(r["mean_fraction_school"])
    cs  = mfc / mfa if mfa > 0 else 0
    ts  = mft / mfa if mfa > 0 else 0
    ss  = mfs / mfa if mfa > 0 else 0
    bar_c = "█" * int(cs * 30)
    bar_t = "▓" * int(ts * 30)
    bar_s = "░" * int(ss * 30)
    print(f"  h{h:02d}  com={100*cs:4.1f}%  ret={100*ts:4.1f}%  sch={100*ss:4.1f}%  {bar_c}{bar_t}{bar_s}")

print("\n  Weekend flat shares applied:")
for dow, day in DOW_TO_DAY.items():
    cs = commute_share_weekend[day]
    ts = retail_share_weekend[day]
    ss = school_share_weekend[day]
    print(f"  {day}: com={100*cs:.1f}%  ret={100*ts:.1f}%  sch={100*ss:.1f}%  "
          f"res={100*(1-cs-ts-ss):.1f}%")
