"""
Compute AADT-normalised hourly traffic fractions from the 2023 NI count dataset,
broken down by day of the week.

For each site the divisor is the site's own weekly average daily traffic
(weekly total / 7), so fractions carry day-of-week volume information.
A fraction f[D, H] means: "hour H on day D carries fraction f of AADT."

Usage:
  python3 analysis/hourly_fractions.py

Conversion formula:
  AADT = observed_count / (mean_fraction[D, H] * (duration_s / 3600))
"""

import numpy as np
import pandas as pd

ODS_FILE = "data/2023-northern-ireland-traffic-count-data-in-ods-format.ods"
OUT_CSV  = "analysis/hourly_fractions.csv"

DAY_NAMES = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]

# Row/column positions confirmed by inspection of the ODS layout
DATA_ROW_START = 15   # first hour (00:00), 0-indexed
DATA_ROW_END   = 39   # exclusive; rows 15-38 = 24 hours
DAY_COLS       = [1, 2, 3, 4, 5, 6, 7]   # Mon=1 … Sun=7

print(f"Reading {ODS_FILE} …")
xl = pd.ExcelFile(ODS_FILE, engine="odf")

# fractions_by_day[d] = list of 24-element arrays, one per valid site for day d
fractions_by_day = [[] for _ in range(7)]

for name in xl.sheet_names:
    df = pd.read_excel(xl, sheet_name=name, header=None)

    if df.shape[0] < DATA_ROW_END or df.shape[1] <= max(DAY_COLS):
        continue

    block = df.iloc[DATA_ROW_START:DATA_ROW_END, DAY_COLS]
    block = block.apply(pd.to_numeric, errors="coerce")

    if block.shape[0] != 24:
        continue

    # Weekly average daily traffic for this site (divisor for AADT normalisation)
    weekly_total = block.values.astype(float).sum()
    if weekly_total <= 0 or np.isnan(weekly_total):
        continue
    aadt = weekly_total / 7.0

    for day_idx, col in enumerate(block.columns):
        counts = block[col].values.astype(float)
        if np.isnan(counts).any():
            continue
        fractions_by_day[day_idx].append(counts / aadt)

hours = [f"{h:02d}:00" for h in range(24)]

rows = []
for day_idx, day_name in enumerate(DAY_NAMES):
    arrs = fractions_by_day[day_idx]
    N = len(arrs)
    if N == 0:
        continue
    arr    = np.array(arrs)   # shape: (N_sites, 24)
    mean_f = arr.mean(axis=0)
    std_f  = arr.std(axis=0, ddof=1) if N > 1 else np.zeros(24)

    day_total = mean_f.sum()   # should equal day_volume / AADT for this day
    print(f"\n{day_name}  ({N} sites, day/AADT ratio={day_total:.3f})")
    print(f"{'Hour':<7}  {'Mean %AADT':>10}  {'Std %AADT':>9}")
    print("-" * 32)
    for h, (m, s) in enumerate(zip(mean_f, std_f)):
        print(f"{hours[h]:<7}  {m*100:10.3f}%  {s*100:9.3f}%")

    for h in range(24):
        rows.append({
            "day_of_week":   day_idx,
            "day_name":      day_name,
            "hour":          hours[h],
            "mean_fraction": mean_f[h],
            "std_fraction":  std_f[h],
        })

out = pd.DataFrame(rows, columns=["day_of_week", "day_name", "hour", "mean_fraction", "std_fraction"])
out.to_csv(OUT_CSV, index=False)
print(f"\nSaved {len(out)} rows → {OUT_CSV}")
