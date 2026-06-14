"""
Compute the typical hourly traffic distribution from the 2023 NI count dataset.

For each site and each day of the week, calculates the fraction of daily
traffic occurring each hour. Reports mean and sample std across all valid
(site, day) pairs, giving one figure ± uncertainty per hour.

Usage:
  python3 analysis/hourly_fractions.py
"""

import numpy as np
import pandas as pd

ODS_FILE = "data/2023-northern-ireland-traffic-count-data-in-ods-format.ods"
OUT_CSV  = "analysis/hourly_fractions.csv"

# Row/column positions confirmed by inspection of the ODS layout
DATA_ROW_START = 15   # first hour (00:00), 0-indexed
DATA_ROW_END   = 39   # exclusive; rows 15-38 = 24 hours
DAY_COLS       = [1, 2, 3, 4, 5, 6, 7]   # Mon=1 … Sun=7

print(f"Reading {ODS_FILE} …")
xl = pd.ExcelFile(ODS_FILE, engine="odf")

all_fractions = []   # list of 24-element arrays, one per valid (site, day)

for name in xl.sheet_names:
    df = pd.read_excel(xl, sheet_name=name, header=None)

    if df.shape[0] < DATA_ROW_END or df.shape[1] <= max(DAY_COLS):
        continue

    block = df.iloc[DATA_ROW_START:DATA_ROW_END, DAY_COLS]
    block = block.apply(pd.to_numeric, errors="coerce")

    if block.shape[0] != 24:
        continue

    for col in block.columns:
        counts = block[col].values.astype(float)
        total  = counts.sum()
        if total <= 0 or np.isnan(total):
            continue
        all_fractions.append(counts / total)

arr = np.array(all_fractions)   # shape: (N, 24)
N   = arr.shape[0]

mean_f = arr.mean(axis=0)
std_f  = arr.std(axis=0, ddof=1)

hours = [f"{h:02d}:00" for h in range(24)]

print(f"\nResults from {N} valid site-day samples ({N // 7} sites × 7 days equiv.)\n")
print(f"{'Hour':<7}  {'Mean %':>8}  {'Std %':>7}")
print("-" * 28)
for h, (m, s) in enumerate(zip(mean_f, std_f)):
    print(f"{hours[h]:<7}  {m*100:8.3f}%  {s*100:7.3f}%")

print(f"\nSum of mean fractions: {mean_f.sum():.6f}  (should be 1.0)")

out = pd.DataFrame({
    "hour":          hours,
    "mean_fraction": mean_f,
    "std_fraction":  std_f,
})
out.to_csv(OUT_CSV, index=False)
print(f"\nSaved to {OUT_CSV}")
