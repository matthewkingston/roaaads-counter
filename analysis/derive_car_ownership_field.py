"""INSPECTION ARTIFACT — per-area car-availability field (NI DZ + RoI SA, harmonised).

Module 2 of the per-area car-ownership work (see .claude/plans/eager-painting-pearl.md):
build the island-wide per-small-area distribution of PERSONS over household car-
availability bands (0/1/2/3+, matching the M1 band shapes), and inspect its spatial
structure + the NI-vs-RoI gap.  DATA-ONLY: writes analysis/car_ownership_field.json +
reports/car_ownership_field.png; touches NO model artifact.

Basis = persons-by-band (the M1 rates are trips-per-person by the person's household
band), keyed per area_code (NI Data-Zone code, RoI Small-Area code).

NI (DZ)  — data/ni_census/…hh_car_van_tc5_pers…csv is already PERSONS by household
    car-or-van availability (0..5+).  Collapse 3/4/5 -> 3+; drop code -8 ("no code
    required", ~1.4%, communal/non-household).

RoI (SA) — SAPS 2022 T15_1 is HOUSEHOLDS by cars (NC/1C/2C/3C/GE4C/NSC).  There is no
    persons-by-car census cross-tab, so convert households->persons with a persons-per-
    household-by-band profile.  DECISION (see plan): use the NTS actual-household profile
    (computed here from SN 5340), which an ecological (Goodman) regression on RoI's own
    marginals validated on the identifiable bands (0/1/3+ agree; RoI marginals cannot
    reliably pin band 2, so the actual-household NTS value is preferred there).  Collapse
    3C+GE4C -> 3+; drop NSC (not-stated, ~7.1%) — proportional, so person-shares are
    unaffected.

Definitional gaps on the record: NI counts car-OR-VAN, RoI counts cars only (T15_1);
RoI code alignment to the model's padded SA id is an M3 concern (keyed here by GEOGID).

Usage:  python3 analysis/derive_car_ownership_field.py
Needs data/NTS (NTS profile), the RoI SAPS Small-Area CSV, and the NI persons CSV.
"""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nts_microdata as nts

OUT_JSON = "analysis/car_ownership_field.json"
OUT_PNG = "reports/car_ownership_field.png"
ROI_SAPS = ("data/ireland_data/Complete_set_of_Census_2022_SAPs/"
            "SAPS_2022_Small_Area_UR_171024.csv")
NI_CAR_GLOB = "data/ni_census/*hh_car_van_tc5_pers*.csv"
BANDS = ["0", "1", "2", "3+"]
NTS_YEARS = [2023, 2024]
# Band midpoints for a scalar summary = person-weighted mean HOUSEHOLD car-availability
BAND_CARS = {"0": 0.0, "1": 1.0, "2": 2.0, "3+": 3.5}


def nts_persons_per_hh():
    """persons-per-household by car band (0/1/2/3+), W2-weighted, from NTS SN 5340 —
    the households->persons conversion profile for RoI (validated by the ecological
    regression on RoI marginals; see module docstring)."""
    hh = nts.load("household", columns=["SurveyYear", "NumCarVan", "HHoldNumPeople", "W2"],
                  years=NTS_YEARS).dropna(subset=["NumCarVan", "HHoldNumPeople", "W2"])
    hh["band"] = hh.NumCarVan.map(lambda n: str(int(min(n, 3))) if int(n) < 3 else "3+")
    hh.loc[hh.NumCarVan >= 3, "band"] = "3+"
    prof = {}
    for b in BANDS:
        sub = hh[hh.band == b]
        prof[b] = float((sub.HHoldNumPeople * sub.W2).sum() / sub.W2.sum())
    return prof


def ni_persons_by_band():
    """{DZ_code: {band: persons}} from the NI persons-by-car-availability CSV."""
    path = sorted(glob.glob(NI_CAR_GLOB))[0]
    df = pd.read_csv(path)
    code_col, avail_col, cnt_col = df.columns[0], df.columns[2], df.columns[4]
    out = {}
    for _, r in df.iterrows():
        c = str(r[avail_col]).strip()
        if c == "-8":                       # "no code required" — drop (non-household)
            continue
        band = c if c in ("0", "1", "2") else "3+"   # 3/4/5 -> 3+
        dz = str(r[code_col]).strip()
        out.setdefault(dz, {b: 0.0 for b in BANDS})
        out[dz][band] += float(r[cnt_col])
    return out


def roi_persons_by_band(prof):
    """{SA_code: {band: persons}} from RoI SAPS households-by-car × NTS persons/hh profile."""
    cols = ["GEOGID", "T1_1AGETT", "T15_1_NC", "T15_1_1C", "T15_1_2C",
            "T15_1_3C", "T15_1_GE4C"]           # NSC dropped (proportional)
    d = pd.read_csv(ROI_SAPS, usecols=cols)
    d = d[pd.to_numeric(d.T1_1AGETT, errors="coerce") <= 100000]   # drop CSO 'State' row
    for c in cols[2:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=cols[2:])
    out = {}
    for _, r in d.iterrows():
        hh_band = {"0": r.T15_1_NC, "1": r.T15_1_1C, "2": r.T15_1_2C,
                   "3+": r.T15_1_3C + r.T15_1_GE4C}
        sa = str(r.GEOGID).strip()
        out[sa] = {b: float(hh_band[b]) * prof[b] for b in BANDS}
    return out


def _shares(pb):
    tot = sum(pb.values())
    return ({b: pb[b] / tot for b in BANDS}, tot) if tot > 0 else (None, 0.0)


def _hh_car_index(share):
    return sum(BAND_CARS[b] * share[b] for b in BANDS)


def build():
    prof = nts_persons_per_hh()
    print("NTS persons/household by band (RoI conversion profile): "
          + "  ".join(f"{b}={prof[b]:.2f}" for b in BANDS))
    ni = ni_persons_by_band()
    roi = roi_persons_by_band(prof)
    print(f"areas: NI DZ={len(ni)}  RoI SA={len(roi)}")

    field = {}
    for code, pb in list(ni.items()):
        share, tot = _shares(pb)
        if share:
            field[code] = {"j": "NI", "persons_total": tot, "share": share,
                           "hh_car_index": _hh_car_index(share)}
    for code, pb in list(roi.items()):
        share, tot = _shares(pb)
        if share:
            field[code] = {"j": "RoI", "persons_total": tot, "share": share,
                           "hh_car_index": _hh_car_index(share)}
    return prof, field


def inspect(field):
    """Print the NI-vs-RoI aggregate band split + per-area index distribution."""
    for j in ("NI", "RoI"):
        areas = {k: v for k, v in field.items() if v["j"] == j}
        totp = sum(v["persons_total"] for v in areas.values())
        agg = {b: sum(v["persons_total"] * v["share"][b] for v in areas.values()) / totp
               for b in BANDS}
        idx = np.array([v["hh_car_index"] for v in areas.values()])
        wts = np.array([v["persons_total"] for v in areas.values()])
        wmean = float((idx * wts).sum() / wts.sum())
        print(f"\n{j}: {len(areas)} areas, {totp:,.0f} persons")
        print("  aggregate person-share by band: "
              + "  ".join(f"{b}={agg[b]*100:4.1f}%" for b in BANDS))
        print(f"  0-car person share = {agg['0']*100:.1f}%   "
              f"person-wtd mean HOUSEHOLD car-availability = {wmean:.3f} (NOT cars/person)")
        print(f"  per-area index: p10={np.percentile(idx,10):.2f} "
              f"median={np.median(idx):.2f} p90={np.percentile(idx,90):.2f}")


def plot(field, out_png):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:                                       # noqa: BLE001
        print(f"  (plot skipped: {e})")
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    x = np.arange(len(BANDS))
    for j, col in (("NI", "tab:blue"), ("RoI", "tab:green")):
        areas = {k: v for k, v in field.items() if v["j"] == j}
        totp = sum(v["persons_total"] for v in areas.values())
        agg = [sum(v["persons_total"] * v["share"][b] for v in areas.values()) / totp
               for b in BANDS]
        ax1.bar(x + (0.2 if j == "RoI" else -0.2), agg, width=0.4, label=j, color=col)
        idx = [v["hh_car_index"] for v in areas.values()]
        ax2.hist(idx, bins=40, density=True, alpha=0.5, label=j, color=col)
    ax1.set_xticks(x); ax1.set_xticklabels(BANDS); ax1.set_xlabel("cars/vans band")
    ax1.set_ylabel("person share"); ax1.set_title("aggregate person-share by band")
    ax1.legend()
    ax2.set_xlabel("per-area mean household car-availability (index)"); ax2.set_ylabel("density")
    ax2.set_title("per-area ownership index distribution"); ax2.legend()
    fig.suptitle("Per-area car-availability field (NI DZ + RoI SA)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=110)
    print(f"  plot -> {out_png}")


def main():
    print("Building per-area car-availability field (NI DZ + RoI SA) …")
    prof, field = build()
    inspect(field)
    out = {"_meta": {
        "note": "INSPECTION ARTIFACT (module 2) — per-area persons-by-band car availability",
        "bands": BANDS,
        "basis": "persons by household car-availability band; person-share per area_code",
        "ni_source": "NI 2021 hh_car_van_tc5_pers (persons; car-or-van; -8 dropped)",
        "roi_source": "RoI SAPS 2022 T15_1 households × NTS persons/hh profile (NSC dropped)",
        "roi_conversion_profile_persons_per_hh": prof,
        "caveats": ["NI car-or-van vs RoI cars-only", "RoI keyed by GEOGID (M3 pads to model SA id)"],
    }, "areas": field}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f)
    print(f"\nSaved -> {OUT_JSON}  ({len(field)} areas)")
    plot(field, OUT_PNG)


if __name__ == "__main__":
    main()
