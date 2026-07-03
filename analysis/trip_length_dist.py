"""Per-component car-driver trip-length distributions (TLDs) from the NTS microdata.

Writes ``analysis/trip_length_distributions.json`` — for each gravity component
(res / commute / retail) the empirical distribution of **car-driver** trip lengths
(miles), as a normalised histogram on a fixed bin grid, with per-bin effective-n so a
downstream fit can be uncertainty-weighted.  This is the *numerator* of the intended
empirical kernel ``f = TLD / n(t)`` (memory: project_tld_retail_mixture,
project_n_eng_source_geometry): the TLD ÷ a national opportunity geometry n(t), then a
double-exp fit, recovers the willingness decay.  **This module builds only the
distributions** — not the ÷n(t) divide, not the kernel fit, not the three school levels
(deferred, as with driveshare/temporal; the deployed design uses one shared tau_school
with no per-level school TLD).

Mode / axis / range (settled with the user)
--------------------------------------------
* **Car-driver basis:** MainMode_B04ID in {3,5,12} (car/van driver, motorcycle, taxi) —
  the numerator of f = TLD_car/n, same basis as generation and driveshare.
* **Distance axis, miles** (TripDisIncSW).  The reliable NTS quantity; reconcile to the
  model's time axis later at the ÷n(t) step via analysis/equiv_miles.py.
* **Full range, no cap.**  The fat tail carries the willingness long-scale; unlike
  driveshare (rise only, capped at 25 mi) the whole range is kept.

Body vs tail
------------
* **Body (<50 mi): the trip table alone** — Sum(JJXSC x W5) reproduces the published TLD
  bit-for-bit across all lengths (memory: project_nts_ldj_tld_weighting).
* **Tail (>=50 mi): LDJ-boosted shape.**  The long-distance-journey table (`ldj`) adds a
  *contemporaneous* extra week of >=50 mi journeys (interview-recalled prior week,
  TripID='NA'), roughly doubling the tail sample with no decade drift — preferred over
  year-pooling.  Handled correctly (a NAIVE splice double-counts the tail):
  - **Purpose:** LDJPurpose_B01ID is byte-identical to TripPurpose_B01ID -> B01_COMPONENT
    directly, no crosswalk.
  - **Mode:** LDJMode_B01ID has no driver/passenger split (4 Car, 6 Van, 5 Motorcycle,
    16 Taxi, 17 Minicab).  Restrict to those vehicle modes (excludes rail/coach/air, whose
    long-haul lengths would pollute the car tail); a driver-fraction correction on Car+Van
    is measured from the trip table's own >=50 mi rows.
  - **Period exposure:** the LDJ >=50 block is **total-pinned** to the trip table's
    (unbiased, 1-week) >=50 mi car-driver mass per component, so LDJ contributes only the
    within-tail *shape* and the ~/2 period factor is absorbed exactly (with any per-year
    calibration imperfection).  The raw ratio and the tranche-(a) consistency check are
    reported as diagnostics.
  - **Sentinels:** TripID / LDJDistance store missing as the string 'NA' — use
    nts_microdata.na_mask (numeric coercion also turns 'NA' -> NaN).
  LDJ reaches only >=50 mi, so the 25-50 mi band still rests on 2023/24 (a coarse bin);
  year-pooling there is the fallback if it prints thin.

Re-derive with:  python3 analysis/trip_length_dist.py --build
(prints a per-bin table per component + writes the JSON and reports/trip_length_dist.png)
"""

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nts_microdata as nts
from purpose_mapping import B01_COMPONENT
from driveshare import _load_simple, SIMPLE_COMPONENTS   # reuse the car-driver body loader

OUT_JSON = "analysis/trip_length_distributions.json"
PLOT_PATH = "reports/trip_length_dist.png"
YEARS = [2023, 2024]

# Fine at short range (resolve the peak/rise), coarse in the tail (guarantee counts).
BIN_EDGES = np.array([0, .25, .5, .75, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 25, 35,
                      50, 75, 100, 150, 300.])
LDJ_MIN_MI = 50.0                      # LDJ covers >=50 mi journeys only
# LDJMode_B01ID vehicle codes: 4 Car, 5 Motorcycle, 6 Van/lorry, 16 Taxi, 17 Minicab.
LDJ_VEHICLE_MODES = [4, 5, 6, 16, 17]
LDJ_CARVAN_MODES = [4, 6]              # driver/passenger not split -> driver-fraction correction


def _mids(edges):
    return np.array([np.sqrt(edges[i] * edges[i + 1]) if edges[i] > 0 else edges[1] / 2
                     for i in range(len(edges) - 1)])


def _bin_stats(dist, w, edges):
    """Per-bin weighted mass Sum(w) and effective-n = (Sum w)^2 / Sum(w^2)."""
    idx = pd.cut(np.asarray(dist, float), bins=edges, right=False, labels=False)
    df = pd.DataFrame({"b": idx, "w": np.asarray(w, float)}).dropna(subset=["b"])
    g = df.groupby("b")
    n = len(edges) - 1
    mass = g.w.sum().reindex(range(n)).fillna(0.0).values
    sw2 = g.w.apply(lambda s: float((s ** 2).sum())).reindex(range(n)).fillna(0.0).values
    neff = np.zeros_like(mass)
    np.divide(mass ** 2, sw2, out=neff, where=sw2 > 0)      # effective-n; 0 where empty
    return mass, neff


def _load_body():
    """Car-driver trip records + the >=50 mi driver fraction (single trip-table stream).

    ``driveshare._load_simple`` keeps all modes (it fits a Bernoulli), so the passenger
    rows needed for the driver fraction are present; y==1 selects the car-driver TLD trips.
    """
    tr = _load_simple()                                   # dist, w=JJXSC*W5, y, component, MainMode
    tail = tr[tr.dist >= LDJ_MIN_MI]
    drv = tail[tail.MainMode_B04ID == 3].w.sum()          # car/van driver
    pas = tail[tail.MainMode_B04ID == 4].w.sum()          # car/van passenger
    driver_frac = float(drv / (drv + pas)) if (drv + pas) > 0 else 1.0
    body = tr[(tr.y == 1) & tr.component.isin(SIMPLE_COMPONENTS)].copy()
    return body, driver_frac


def _load_ldj(driver_frac):
    """LDJ >=50 mi vehicle journeys, component-mapped, with the car-driver weight ``wv``."""
    ldj = nts.load("ldj", columns=["SurveyYear", "LDJPurpose_B01ID", "LDJMode_B01ID",
                                    "LDJDistance", "TripID", "W4"], years=YEARS)
    ldj["dist"] = pd.to_numeric(ldj.LDJDistance, errors="coerce")   # 'NA' string -> NaN
    ldj["component"] = ldj.LDJPurpose_B01ID.map(B01_COMPONENT)
    ldj["is_diary"] = ~nts.na_mask(ldj.TripID)                      # tranche (a): also in trip table
    veh = ldj[ldj.LDJMode_B01ID.isin(LDJ_VEHICLE_MODES)
              & ldj.dist.notna() & (ldj.dist >= LDJ_MIN_MI)
              & ldj.W4.notna() & (ldj.W4 > 0)
              & ldj.component.isin(SIMPLE_COMPONENTS)].copy()
    # Car+Van carry drivers+passengers; scale to car-driver.  (Motorcycle/taxi/minicab as-is.)
    veh["wv"] = veh.W4 * np.where(veh.LDJMode_B01ID.isin(LDJ_CARVAN_MODES), driver_frac, 1.0)
    return veh


def build():
    print(f"Building per-component car-driver TLDs — NTS microdata {YEARS}, "
          f"body Sum(JJXSC x W5), tail LDJ-boosted >=50 mi\n")
    body, driver_frac = _load_body()
    veh = _load_ldj(driver_frac)
    edges = BIN_EDGES
    mids = _mids(edges)
    tail0 = int(np.searchsorted(edges, LDJ_MIN_MI))       # first bin whose left edge is >=50
    print(f"  driver fraction (>=50 mi car+van, trip table): {driver_frac:.3f}")

    # --- period-exposure diagnostics (per year) -------------------------------------
    print("  LDJ period/consistency diagnostics (per year):")
    print(f"    {'year':>6} {'body>=50 Sw':>12} {'ldj a+b Swv':>12} {'rho~0.5':>8} "
          f"{'ldj(a) Swv':>11} {'trancheA~1':>11}")
    for y in YEARS:
        b50 = body[(body.SurveyYear == y) & (body.dist >= LDJ_MIN_MI)].w.sum()
        vy = veh[veh.SurveyYear == y]
        ab = vy.wv.sum()
        a = vy[vy.is_diary].wv.sum()
        rho = b50 / ab if ab > 0 else float("nan")        # ~0.5: a+b is ~2 weeks
        chk = a / b50 if b50 > 0 else float("nan")        # ~1: tranche (a) reproduces trip table
        print(f"    {y:>6} {b50:>12,.0f} {ab:>12,.0f} {rho:>8.3f} {a:>11,.0f} {chk:>11.3f}")
    print()

    out = {"_meta": {
        "source": "NTS microdata (UKDS SN 5340) via analysis/nts_microdata.py",
        "years": YEARS,
        "mode_basis": "car-driver MainMode_B04ID {3,5,12}; LDJ vehicle LDJMode_B01ID {4,5,6,16,17}",
        "weight": "body JJXSC x W5; tail LDJ W4 x driver_fraction, total-pinned to trip-table >=50 mass",
        "axis": "trip length, miles (TripDisIncSW / LDJDistance)",
        "driver_fraction_ge50": driver_frac,
        "bin_edges_miles": edges.tolist(),
        "tail_source": f"LDJ >=50 mi shape, total-pinned per component (first tail bin index {tail0})",
    }, "bin_edges_miles": edges.tolist(), "bin_mid_miles": mids.tolist(), "components": {}}

    widths = np.diff(edges)
    for comp in SIMPLE_COMPONENTS:
        b = body[body.component == comp]
        mass_tt, neff_tt = _bin_stats(b.dist, b.w, edges)          # trip-table-only (all lengths)

        v = veh[veh.component == comp]
        ldj_mass, ldj_neff = _bin_stats(v.dist, v.wv, edges)       # LDJ, populated only in tail
        # Total-pin the LDJ >=50 block to the trip table's >=50 mass, borrow its shape.
        tt_tail_total = mass_tt[tail0:].sum()
        ldj_tail_total = ldj_mass[tail0:].sum()
        mass_lb = mass_tt.copy()
        neff_lb = neff_tt.copy()
        if ldj_tail_total > 0:
            scale = tt_tail_total / ldj_tail_total
            mass_lb[tail0:] = ldj_mass[tail0:] * scale
            neff_lb[tail0:] = ldj_neff[tail0:]                     # effective-n unaffected by scaling

        def _norm(mass):
            tot = mass.sum()
            share = mass / tot if tot > 0 else mass
            return share, share / widths                           # share (Sum=1), density (per mile)

        share_tt, dens_tt = _norm(mass_tt)
        share_lb, dens_lb = _norm(mass_lb)

        out["components"][comp] = {
            "n_trips_body": int(len(b)), "sum_w_body": float(b.w.sum()),
            "n_ldj_tail": int(len(v)),
            "triptable_only": {"mass": mass_tt.tolist(), "share": share_tt.tolist(),
                               "density": dens_tt.tolist(), "eff_n": neff_tt.tolist()},
            "ldj_boosted": {"mass": mass_lb.tolist(), "share": share_lb.tolist(),
                            "density": dens_lb.tolist(), "eff_n": neff_lb.tolist()},
        }

        # --- printed per-bin table (LDJ-boosted; trip-table tail eff-n in brackets) ----
        print(f"=== {comp}  (body n {len(b):,}, Sum_w {b.w.sum():,.0f}; LDJ tail rows {len(v):,}) ===")
        print(f"  {'lo':>6} {'hi':>6} {'share%':>7} {'dens/mi':>9} {'eff_n':>9} "
              f"{'[tt eff_n]':>11}")
        for i in range(len(edges) - 1):
            tag = " <LDJ" if i >= tail0 else ""
            print(f"  {edges[i]:>6.2f} {edges[i+1]:>6.2f} {100*share_lb[i]:>7.2f} "
                  f"{dens_lb[i]:>9.4f} {neff_lb[i]:>9.0f} {neff_tt[i]:>11.0f}{tag}")
        print(f"  >=50 mi share: trip-table-only {100*share_tt[tail0:].sum():.2f}%  "
              f"LDJ-boosted {100*share_lb[tail0:].sum():.2f}%\n")

    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved -> {OUT_JSON}")
    _plot(out, edges, mids, tail0)


def _plot(out, edges, mids, tail0):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(plot skipped: {e})")
        return
    colors = {"commute": "tab:blue", "retail": "tab:green", "res": "tab:red"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    for comp in SIMPLE_COMPONENTS:
        c = colors[comp]
        d_tt = np.array(out["components"][comp]["triptable_only"]["density"])
        d_lb = np.array(out["components"][comp]["ldj_boosted"]["density"])
        for ax in axes:
            ax.plot(mids, d_lb, "-", color=c, lw=2, label=f"{comp} (LDJ-boosted tail)")
            ax.plot(mids, d_tt, "--", color=c, lw=1, alpha=0.7,
                    label=f"{comp} (trip-table only)")
    axes[0].set_xlim(0, 25)
    axes[0].set_title("body (0-25 mi, linear)", fontsize=10)
    axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[1].axvline(edges[tail0], color="grey", ls=":", lw=1)
    axes[1].set_title("full range (log-log; dotted = LDJ tail start 50 mi)", fontsize=10)
    for ax in axes:
        ax.set_xlabel("trip length (miles)")
        ax.set_ylabel("density (share per mile)")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(alpha=0.3)
    fig.suptitle("Per-component car-driver trip-length distributions (NTS microdata)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(PLOT_PATH), exist_ok=True)
    fig.savefig(PLOT_PATH, dpi=120)
    print(f"saved {PLOT_PATH}")


if __name__ == "__main__":
    if "--build" in sys.argv:
        build()
    else:
        print(__doc__.split("\n\n")[0])
        print("\nRe-derive with: python3 analysis/trip_length_dist.py --build")
