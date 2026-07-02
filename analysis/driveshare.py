"""Vehicle-driver mode share as a function of trip length (miles), per component.

The gravity kernel's short-range rise is *mode substitution*: short trips are
walked, not driven, so car demand is suppressed at short cost.  `driveshare(d, c)`
is that factor for component `c`; the kernel is
    f(c) = driveshare(equiv_miles(c), component) * exp(-c / tau_c)
(see simulation/model.py `_modesub_kernel`; `equiv_miles` maps OSRM seconds→miles).

**Per-component curves.** The mode-substitution rise is *not* identical across trip
purposes (a short commute is walked far more readily than a short residential run),
so each component gets its own fixed driveshare curve, derived from the NTS trip-level
microdata (UKDS SN 5340) via the shared `purpose_mapping.B01_COMPONENT` mapping — the
same purpose→component scheme as the generation rates.  This adds **zero tuned
parameters** (driveshare stays empirical/data-fixed; only tau_c is tuned) and de-confounds
tau_c (a purpose whose real short-range walk-substitution differs from a shared curve no
longer pushes that difference into its willingness decay).

Derivation (`--fit`): survey-weighted **binomial MLE** of the share form directly on the
trip records — each trip a Bernoulli `veh-driver?` outcome at its actual `TripDisIncSW`
length, frequency-weighted by `JJXSC×W5` (short walks ×7, series-of-calls ×0), vehicle-driver
modes = car/van driver + motorcycle + taxi (MainMode_B04ID {3,5,12}, same set as generation),
years 2023/24.  Two corrections that the old aggregate-table derivation needed are now
structural, not assumptions:

  1. JUST-WALK removed exactly at the record level (TripPurpose_B01ID == 17 dropped),
     replacing the old "allocate the all-modes just-walk total ∝ walk-length profile".
  2. Continuous length pins the sub-1-mile rise directly (evaluated at each trip's real
     length), replacing the coarse-band + origin-pinning band-integral trick.  The form's
     K>0 gives driveshare(0)=0 (⇒ f(0)=0) structurally.

The England long-distance rail decline (past ~25 mi) is handled by **capping the fit at
25 mi**: the monotone saturating form cannot represent that decline, those trips are killed
by the willingness exp(-c/tau) and carry no weight in the model, and (D0,K) are stable to
the cap.  This is the NI/RoI transfer stance (sparser rail ⇒ car share holds at its
plateau), equivalent to the old plateau-clip.

Form: driveshare(d, c) = PLATEAU_c * (1 - exp(-(d/D0_c)^K_c)).  PLATEAU_c is interpretive
only — a constant factor that cancels *within each component* in the production constraint,
so only the rise shape (D0_c, K_c) is load-bearing.

School is **not yet component-fitted** — its aggregate purpose (education 4 + escort 21)
driveshare is confounded by the child's own non-driver trip (purpose 4 ≈ 0–13% driver, all
passenger/walk), so the car school-run curve must come from the escort/self-drive trips, a
separate derivation (see project-nts-microdata-gains).  Until then `school` carries the shared
placeholder curve, and `component=None` returns that same legacy all-purpose blend so the
not-yet-per-component-wired model runs unchanged.

The module constants below are the *single* authoritative source of truth; the imported
`driveshare` does no file I/O.  Re-derive with `--fit` after the NTS data change (prints
refreshed constants to paste back, and draws the plot).
"""

import numpy as np

# --- Fitted per-component constants (NTS microdata 2023/24; weighted binomial MLE, cap 25 mi)
# Re-derive with: python3 analysis/driveshare.py --fit
#   component: (PLATEAU, D0 [miles], K)
CURVES = {
    "commute": (0.6940, 1.2867, 0.9888),   # slow rise, high plateau
    "retail":  (0.5492, 0.8706, 1.6052),   # fast/steep rise, low plateau
    "res":     (0.6148, 0.7409, 1.3967),   # earliest rise
    # School is not component-fitted yet (aggregate is confounded by the child's own
    # non-driver trip); shared placeholder until the escort-based school derivation.
    "school":  (0.5779, 1.2573, 1.2610),
}

# Shared/default curve (component=None): the legacy all-purpose blend, kept until model.py
# is wired to pass a component through _modesub_kernel.  Equals the old single fit.
PLATEAU, D0, K = 0.5779, 1.2573, 1.2610


def driveshare(d_miles, component=None):
    """Vehicle-driver share for a trip of length d_miles in `component`.

    Float or numpy array.  `component` is one of CURVES (commute/retail/res/school); the
    default `None` uses the shared legacy curve (for callers not yet passing a component).
    """
    pl, d0, k = CURVES[component] if component is not None else (PLATEAU, D0, K)
    return pl * (1.0 - np.exp(-(np.asarray(d_miles, dtype=float) / d0) ** k))


# --- Re-fit / diagnostics / plot (offline only; not on the import path) ---------

YEARS = [2023, 2024]                 # pool COVID-excluded extra years only for thin cells
VEHICLE_MODE_CODES = [3, 5, 12]      # Car/van driver, Motorcycle, Taxi/minicab
JUST_WALK_B01 = 17                   # TripPurpose_B01ID for "Just walk" (dropped exactly)
FIT_CAP_MILES = 25.0                 # exclude the England rail-decline tail from the fit
FIT_COMPONENTS = ["commute", "retail", "res"]   # school derived separately (escort-based)


def _load_trips():
    """Trip records with dist/weight/veh/component columns (just-walk dropped)."""
    import os
    import sys
    import pandas as pd
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nts_microdata as nts
    from purpose_mapping import B01_COMPONENT

    tr = nts.load("trip", columns=["SurveyYear", "MainMode_B04ID", "TripPurpose_B01ID",
                                   "TripDisIncSW", "W5", "JJXSC"], years=YEARS)
    tr["dist"] = pd.to_numeric(tr.TripDisIncSW, errors="coerce")
    tr["w"] = tr.JJXSC * tr.W5
    tr["y"] = tr.MainMode_B04ID.isin(VEHICLE_MODE_CODES).astype(float)
    tr = tr[(tr.TripPurpose_B01ID != JUST_WALK_B01)
            & tr.dist.notna() & (tr.dist > 0) & (tr.w > 0)].copy()
    tr["component"] = tr.TripPurpose_B01ID.map(B01_COMPONENT)
    return tr


def _fit_component(d, y, w, cap=FIT_CAP_MILES):
    """Weighted binomial MLE of PLATEAU*(1-exp(-(d/D0)^K)), trips with d<=cap."""
    from scipy.optimize import minimize
    m = d <= cap
    d, y, w = d[m].values, y[m].values, w[m].values

    def nll(p):
        pl, d0, k = p
        pr = np.clip(pl * (1.0 - np.exp(-(d / d0) ** k)), 1e-9, 1 - 1e-9)
        return -np.sum(w * (y * np.log(pr) + (1 - y) * np.log(1 - pr)))

    r = minimize(nll, [0.58, 1.26, 1.26], method="L-BFGS-B",
                 bounds=[(0.2, 0.95), (0.1, 30), (0.5, 6)])
    return tuple(r.x)


def _fit(plot_path="reports/driveshare.png"):
    import os
    import numpy as np
    import pandas as pd
    tr = _load_trips()

    # fine bins for the empirical scatter (plot + printed table)
    edges = np.array([0, .25, .5, .75, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 25, 50, 300.])
    mids = np.array([np.sqrt(edges[i] * edges[i + 1]) if edges[i] > 0 else edges[1] / 2
                     for i in range(len(edges) - 1)])

    def emp(sub):
        b = pd.cut(sub.dist, bins=edges, right=False, labels=False)
        num = sub[sub.y > 0].groupby(b).w.sum()
        den = sub.groupby(b).w.sum()
        return (num / den).reindex(range(len(edges) - 1)), den.reindex(range(len(edges) - 1))

    fits = {}
    print(f"NTS microdata {YEARS}, vehicle modes {VEHICLE_MODE_CODES}, "
          f"weight JJXSC×W5, just-walk removed, fit cap {FIT_CAP_MILES:g} mi\n")
    for comp in FIT_COMPONENTS:
        sub = tr[tr.component == comp]
        pl, d0, k = _fit_component(sub.dist, sub.y, sub.w)
        fits[comp] = (pl, d0, k)
        e, den = emp(sub)
        print(f"=== {comp}  (Σw {sub.w.sum():,.0f}, n {len(sub):,})  "
              f"PLATEAU={pl:.4f} D0={d0:.4f} K={k:.4f} ===")
        print(f"  {'mid_mi':>7} {'emp%':>6} {'fit%':>6} {'Σw':>12}")
        for i in range(len(edges) - 1):
            fv = 100 * pl * (1 - np.exp(-(mids[i] / d0) ** k))
            ev = 100 * e[i] if pd.notna(e[i]) else float("nan")
            dv = den[i] if pd.notna(den[i]) else 0.0
            print(f"  {mids[i]:7.2f} {ev:6.1f} {fv:6.1f} {dv:12,.0f}")
        print()

    print("Paste back into CURVES:")
    for comp in FIT_COMPONENTS:
        pl, d0, k = fits[comp]
        print(f'    "{comp}":{"":>{9-len(comp)}}({pl:.4f}, {d0:.4f}, {k:.4f}),')

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        dd = np.linspace(0.01, 25, 800)
        colors = {"commute": "tab:blue", "retail": "tab:green", "res": "tab:red"}
        fig, ax = plt.subplots(figsize=(9, 5.2))
        for comp in FIT_COMPONENTS:
            pl, d0, k = fits[comp]
            sub = tr[tr.component == comp]
            e, _ = emp(sub)
            c = colors[comp]
            ax.plot(dd, pl * (1 - np.exp(-(dd / d0) ** k)), "-", color=c, lw=2,
                    label=f"{comp}: {pl:.2f}(1-exp(-(d/{d0:.2f})^{k:.2f}))")
            ax.scatter(mids, e, color=c, s=22, zorder=5, alpha=0.8)
        ax.plot(dd, PLATEAU * (1 - np.exp(-(dd / D0) ** K)), "k--", lw=1.3, alpha=0.7,
                label=f"shared (legacy): {PLATEAU:.2f}(1-exp(-(d/{D0:.2f})^{K:.2f}))")
        ax.set_xlim(0, 25)
        ax.set_ylim(0, 0.8)
        ax.set_xlabel("trip length (miles)")
        ax.set_ylabel("vehicle-driver share")
        ax.set_title("Per-component driveshare (NTS microdata 2023/24, points = binned empirical)")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        os.makedirs(os.path.dirname(plot_path), exist_ok=True)
        fig.savefig(plot_path, dpi=120)
        print(f"\nsaved {plot_path}")
    except Exception as e:
        print(f"\n(plot skipped: {e})")


if __name__ == "__main__":
    import sys
    if "--fit" in sys.argv:
        _fit()
    else:
        print("driveshare(d_miles, component) module. Re-derive with: "
              "python3 analysis/driveshare.py --fit\n")
        for comp in ["commute", "retail", "res"]:
            vals = "  ".join(f"{d}mi={float(driveshare(d, comp)):.3f}"
                             for d in [0.25, 0.5, 1, 2, 5])
            print(f"  {comp:8s}: {vals}")
