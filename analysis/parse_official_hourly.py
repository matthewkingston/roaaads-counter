"""
Parse per-hour vehicle count data for official AADT sites 444, 507, 508 from
the Northern Ireland 2023 traffic count ODS file.

Writes data/official_hourly.json with one observation per (site, day_type, hour),
in the same format consumed by the two-component tuner:
  count  — average vehicles/hour for that slot
  sigma  — uncertainty (between-weekday std, 10% floor; 15% floor for Sat/Sun)
  T_s    — 3600 (observation is exactly 1 hour)

These replace the three single-AADT observations previously used in COUNT_SITES.
"""

import json, math, os, sys
sys.path.insert(0, "simulation")
from model import COUNT_SITES

try:
    from odf.opendocument import load
    from odf.table import Table, TableRow, TableCell
    from odf.text import P
except ImportError:
    print("ERROR: odfpy not installed — run:  pip install odfpy")
    raise SystemExit(1)

ODS_FILE   = "data/2023-northern-ireland-traffic-count-data-in-ods-format.ods"
OUTPUT     = "data/official_hourly.json"

# Map ODS sheet name → COUNT_SITES entry
SITE_MAP = {s["label"].split(",")[0].strip().split()[-1]: s for s in COUNT_SITES}


def cell_text(cell):
    return " ".join(str(p) for p in cell.getElementsByType(P)).strip()


def parse_sheet(sheet):
    """Return list of (hour, mon, tue, wed, thu, fri, sat, sun) for hours 0–23.

    Sheets contain multiple sub-sections (All directions, per-lane, per-direction).
    Only the first 24 hour-rows are taken — these belong to the 'All directions'
    section that always appears first.
    """
    rows = sheet.getElementsByType(TableRow)
    data = []
    seen_hours = set()
    for row in rows:
        if len(data) == 24:
            break
        cells = row.getElementsByType(TableCell)
        vals = [cell_text(c) for c in cells]
        if not vals or ":" not in vals[0] or vals[0].count(":") != 2:
            continue
        try:
            hour = int(vals[0].split(":")[0])
        except ValueError:
            continue
        if hour in seen_hours:
            break   # second section starting — stop
        try:
            counts = [int(vals[i]) if i < len(vals) and vals[i].isdigit() else 0
                      for i in range(1, 8)]
        except ValueError:
            continue
        if len(counts) == 7:
            data.append((hour, *counts))
            seen_hours.add(hour)
    return sorted(data, key=lambda x: x[0])


print(f"Loading {ODS_FILE} …")
doc    = load(ODS_FILE)
sheets = {s.getAttribute("name"): s
          for s in doc.spreadsheet.getElementsByType(Table)}

output   = {}
total_n  = 0
MIN_REL  = 0.10   # 10% relative floor for weekday sigma
MIN_REL_WE = 0.15 # 15% relative floor for Sat/Sun

for site_label, site_cfg in SITE_MAP.items():
    sheet_name = site_label  # "507", "508", "444"
    if sheet_name not in sheets:
        print(f"  WARNING: sheet '{sheet_name}' not found — skipping")
        continue

    rows = parse_sheet(sheets[sheet_name])
    if len(rows) != 24:
        print(f"  WARNING: site {sheet_name}: expected 24 hours, got {len(rows)}")

    observations = []
    for tup in rows:
        hour, mon, tue, wed, thu, fri, sat, sun = tup

        # Weekday (day_type 0): average of Mon–Fri, sigma from between-day variance.
        # Three floors: between-day std, 10% relative, sqrt(count).
        # sqrt(count) matters at overnight hours where between-day variance is near
        # zero (all days similar) and 10% of a small count is also small, giving an
        # unrealistically tight sigma and inflated z-scores.
        wd_vals = [mon, tue, wed, thu, fri]
        wd_mean = sum(wd_vals) / 5.0
        wd_var  = sum((v - wd_mean) ** 2 for v in wd_vals) / 4.0  # sample variance
        wd_sig  = max(math.sqrt(wd_var), MIN_REL * wd_mean, math.sqrt(max(wd_mean, 0.5)))

        # Saturday (day_type 1) and Sunday (day_type 2): single column, Poisson + rel floor
        sat_sig = max(math.sqrt(max(sat, 0.5)), MIN_REL_WE * sat)
        sun_sig = max(math.sqrt(max(sun, 0.5)), MIN_REL_WE * sun)

        for day_type, count, sigma in [
            (0, wd_mean, wd_sig),
            (1, float(sat), sat_sig),
            (2, float(sun), sun_sig),
        ]:
            observations.append({
                "time_slot": [day_type, hour],
                "count":     round(count, 3),
                "sigma":     round(sigma, 3),
                "T_s":       3600,
            })

    output[sheet_name] = {
        "label": site_cfg["label"],
        "node":  site_cfg["node"],
        "links": [[u, v] for u, v in site_cfg["links"]] if site_cfg["links"] else None,
        "observations": observations,
    }
    total_n += len(observations)
    print(f"  site {sheet_name}: {len(observations)} observations")

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
with open(OUTPUT, "w") as f:
    json.dump(output, f, indent=2)

print(f"\nWrote {total_n} observations → {OUTPUT}")
