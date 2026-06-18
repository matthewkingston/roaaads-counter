"""
Derive residential and business temporal profile priors for the two-component
gravity model from NTS trip-purpose data.

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
The gravity model's business demand nodes represent workplaces, retail POIs,
and car parks (public and private).  Trips attracted to these nodes are:

  Business ("biz"):
    - Commuting              — direct work trips
    - Employer's business    — trips made during/for work (deliveries, etc.)
    - Education              — to/from school or university; schools are
                               tagged in OSM and receive a bonus POI weight;
                               DOWNWEIGHTED × EDU_FACTOR (see below)
    - Escort education       — taking someone to school; same destination type
    - Shopping               — trips to retail; shops dominate OSM POIs and
                               car parks

  Residential ("res"):
    - Other escort and personal business  (GP, bank, etc.)
    - Visiting friends, entertainment and sport
    - Holiday, day trip and other

Education downweighting (EDU_FACTOR = 1/5)
------------------------------------------
NTS records ALL trips regardless of mode.  The "Education" category is
predominantly children walking, cycling, or being bussed — not car trips.
The "Escort education" category already captures the adult car driver making
the school-run trip; the standalone "Education" trip is the child's own
journey, which has far lower per-person car trip generation than commuting.

We therefore downweight the Education contribution by EDU_FACTOR = 1/5:
the denominator is also reduced by the same amount so that the shares of
all remaining purposes are proportionally increased (consistent reweighting).

  biz_numerator   = commuting + business + education*EDU_FACTOR
                    + escort_education + shopping
  adjusted_total  = all_purposes - education*(1 - EDU_FACTOR)
  biz_share       = biz_numerator / adjusted_total

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
For each hour h, biz_share(h) is computed from the 2023–2024 per-hour
purpose percentages with the education downweighting applied.  This gives
a distinct share for every hour of the day.

Weekend approach (NTS0504b)
---------------------------
NTS0502a covers Monday–Friday only.  NTS0504b provides absolute trip rates
per person per year for each day of the week, including Saturday and Sunday.
We extract the Saturday and Sunday rows from the 2023–2024 rolling average,
remove "Just walk", apply the education downweighting, and compute a single
flat business share for each weekend day.  The temporal SHAPE within each
weekend day is inherited from the aggregate mean_fraction profile (which
already captures the Saturday vs. Sunday difference in volume pattern);
only the overall business / residential split is estimated from NTS0504b.

Constraint
----------
The tuner requires mean_fraction_res + mean_fraction_biz = 2 × mean_fraction
at every (day_of_week, hour) slot.  This is enforced exactly:

  mean_fraction_biz = 2 × mean_fraction × biz_share
  mean_fraction_res = 2 × mean_fraction × (1 − biz_share)

Neither component needs to sum to 1 across hours; the global scale K absorbs
the normalisation.

Usage
-----
  python3 analysis/derive_component_profiles.py

Overwrites the mean_fraction_res and mean_fraction_biz columns in
analysis/hourly_fractions.csv.  All other columns are preserved.
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
yr502["adj_total"] = yr502[C_ALL] - yr502[C_EDU] * (1 - EDU_FACTOR)
yr502["biz_num"]   = (yr502[C_COMM] + yr502[C_BIZ_E]
                      + yr502[C_EDU] * EDU_FACTOR
                      + yr502[C_ESCORT] + yr502[C_SHOP])
yr502["biz_share"] = yr502["biz_num"] / yr502["adj_total"]

biz_share_weekday = dict(zip(yr502["hour"], yr502["biz_share"]))

print(f"\nNTS {NTS_YEAR} weekday business share by hour  (education × {EDU_FACTOR})")
print(f"  {'Hour':>4}  {'biz%':>6}  {'res%':>6}")
for h in range(24):
    bs = biz_share_weekday[h]
    print(f"  h{h:02d}    {100*bs:5.1f}%  {100*(1-bs):5.1f}%")

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

biz_share_weekend = {}
print(f"\nNTS {NTS_YEAR} weekend business share from NTS0504b  (education × {EDU_FACTOR}, walk removed)")
print(f"  {'Day':<10}  {'biz_trips':>9}  {'adj_total':>9}  {'biz%':>6}  {'res%':>6}")

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
    no_walk    = all_p - walk
    adj_total  = no_walk - edu * (1 - EDU_FACTOR)
    biz_num    = comm + biz_e + edu * EDU_FACTOR + escort + shop

    bs = biz_num / adj_total
    biz_share_weekend[day] = bs

    print(f"  {day:<10}  {biz_num:9.3f}  {adj_total:9.3f}  {100*bs:5.1f}%  {100*(1-bs):5.1f}%")

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
        bs = biz_share_weekday[h]
    else:
        bs = biz_share_weekend[DOW_TO_DAY[dow]]

    row["mean_fraction_biz"] = f"{2 * mfa * bs:.10f}"
    row["mean_fraction_res"] = f"{2 * mfa * (1 - bs):.10f}"

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

neg = [(r["day_of_week"], r["hour"]) for r in rows
       if float(r["mean_fraction_res"]) < 0 or float(r["mean_fraction_biz"]) < 0]
if neg:
    print(f"  FAIL: negative fractions at {neg}")
else:
    print("  ✓  all fractions non-negative")

print("\n  Weekday biz share by hour (Mon shown; other weekdays identical):")
wkday_rows = sorted([r for r in rows if int(r["day_of_week"]) == 0],
                    key=lambda r: int(r["hour"].split(":")[0]))
for r in wkday_rows:
    h   = int(r["hour"].split(":")[0])
    mfa = float(r["mean_fraction"])
    mfb = float(r["mean_fraction_biz"])
    bs  = mfb / (2 * mfa) if mfa > 0 else 0
    bar = "█" * int(bs * 40)
    print(f"  h{h:02d}  {100*bs:4.1f}%  {bar}")

print("\n  Weekend flat shares applied:")
for dow, day in DOW_TO_DAY.items():
    bs = biz_share_weekend[day]
    print(f"  {day}: biz={100*bs:.1f}%  res={100*(1-bs):.1f}%")
