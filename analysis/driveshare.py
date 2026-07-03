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

Non-school derivation (`--fit`): survey-weighted **binomial MLE** of the share form
directly on the trip records — each trip a Bernoulli `veh-driver?` outcome at its actual
`TripDisIncSW` length, frequency-weighted by `JJXSC×W5` (short walks ×7, series-of-calls ×0),
vehicle-driver modes = car/van driver + motorcycle + taxi (MainMode_B04ID {3,5,12}, same set
as generation), years 2023/24.  Two corrections the old aggregate-table derivation needed are
now structural, not assumptions:

  1. JUST-WALK removed exactly at the record level (TripPurpose_B01ID == 17 dropped),
     replacing the old "allocate the all-modes just-walk total ∝ walk-length profile".
  2. Continuous length pins the sub-1-mile rise directly (evaluated at each trip's real
     length), replacing the coarse-band + origin-pinning band-integral trick.  The form's
     K>0 gives driveshare(0)=0 (⇒ f(0)=0) structurally.

**School is per-level (primary / post-primary / tertiary)** and measured differently — as a
per-trip **mode share of the child's own school journey**, because the school producer (the
student) is usually *not* the driver.  For each level, the curve is the share of that level's
**Education** trips (TripPurpose_B01ID == 4) made **by car** — child as car driver *or*
passenger, motorcycle or taxi (MainMode_B04ID {3,4,5,12}) — out of all modes, as a function of
distance.  Using the child's own trip (not the parent's escort record, which carries no
child-level tag) makes each trip level-taggable **by the child's age**, reusing the *exact*
generation age→level machinery (`derive_school_generation._level_shares` /
`_england_secondary_share`: 5–10 → primary, 11–15 → post-primary, 16–18 split secondary/
tertiary by the DfE England full-time shares, 19+ FT → tertiary).  The by-car numerator
captures escort (child = passenger, the same physical journey as the parent's escort trip, at
the home→school distance) and self-drive (child = driver) in one curve.  Ride-sharing does not
bias a *share* (siblings sit at the same distance in the same category), so — unlike the
generation magnitude — no ride-share correction is needed here.  Tertiary is thin on 2023/24, so
it is pooled over ex-COVID years (2013–19 + 23–24); primary/post-primary use 2023/24.  The
school run is non-monotonic (walk → car → **bus** at range), which the saturating form cannot
represent; the fit is capped at 25 mi and the long-distance bus decline is left to the
willingness `exp(-c/tau_school)` (which kills long school trips) + the per-component plateau-cancel.

The England long-distance rail decline (past ~25 mi) is likewise handled by **capping the fit
at 25 mi**: the monotone saturating form cannot represent that decline, those trips are killed by
the willingness exp(-c/tau) and carry no weight in the model, and (D0,K) are stable to the cap.
This is the NI/RoI transfer stance (sparser rail ⇒ car share holds at its plateau).

Form: driveshare(d, c) = PLATEAU_c * (1 - exp(-(d/D0_c)^K_c)).  PLATEAU_c is interpretive only —
a constant factor that cancels *within each component* in the production constraint, so only the
rise shape (D0_c, K_c) is load-bearing.

The module constants below are the *single* authoritative source of truth; the imported
`driveshare` does no file I/O.  Re-derive with `--fit` after the NTS data change (prints
refreshed constants to paste back, and draws the plot).
"""

import numpy as np

# --- Fitted per-component constants (NTS microdata; weighted binomial MLE, cap 25 mi)
# Re-derive with: python3 analysis/driveshare.py --fit
#   component: (PLATEAU, D0 [miles], K)
CURVES = {
    "commute": (0.6940, 1.2867, 0.9888),   # slow rise, high plateau
    "retail":  (0.5492, 0.8706, 1.6052),   # fast/steep rise, low plateau
    "res":     (0.6148, 0.7409, 1.3967),   # earliest rise
    # School by-car share of the child's education trip, per level (see module docstring).
    # tertiary pooled over ex-COVID years (2013-19 + 23-24) for sample size.
    "school_primary":     (0.8474, 0.9126, 2.1546),   # high plateau, driven once past walking
    "school_postprimary": (0.4273, 0.9896, 2.3069),   # lower — more bus/cycle/walk independence
    "school_tertiary":    (0.3199, 1.6968, 1.9983),   # lowest, rises latest (drive at range)
}

# Shared/default curve (component=None): the legacy all-purpose blend, kept until model.py is
# wired to pass a component through _modesub_kernel (then removed — no back-compat retained).
PLATEAU, D0, K = 0.5779, 1.2573, 1.2610


def driveshare(d_miles, component=None):
    """Vehicle-driver / by-car share for a trip of length d_miles in `component`.

    Float or numpy array.  `component` is one of CURVES (commute/retail/res and the three
    school_* levels); the default `None` uses the shared legacy curve (for callers not yet
    passing a component).
    """
    pl, d0, k = CURVES[component] if component is not None else (PLATEAU, D0, K)
    return pl * (1.0 - np.exp(-(np.asarray(d_miles, dtype=float) / d0) ** k))


# --- Re-fit / diagnostics / plot (offline only; not on the import path) ---------

YEARS = [2023, 2024]                  # non-school + school primary/post-primary window
VEHICLE_MODE_CODES = [3, 5, 12]       # Car/van driver, Motorcycle, Taxi/minicab (driver basis)
JUST_WALK_B01 = 17                    # TripPurpose_B01ID "Just walk" (dropped exactly)
FIT_CAP_MILES = 25.0                  # exclude the England rail / school-bus decline tail
SIMPLE_COMPONENTS = ["commute", "retail", "res"]   # driver-share, purpose→component mapping

EDUCATION_B01 = 4                     # TripPurpose_B01ID "Education" (the child's own trip)
SCHOOL_BYCAR_MODES = [3, 4, 5, 12]    # child in a car: driver OR passenger, motorcycle, taxi
SCHOOL_POOL_YEARS = [2013, 2014, 2015, 2016, 2017, 2018, 2019, 2023, 2024]   # ex-COVID
# (level-share key, CURVES key, year window) — tertiary pooled, prim/sec recent
SCHOOL_LEVELS = [("prim", "school_primary", YEARS),
                 ("sec", "school_postprimary", YEARS),
                 ("ter", "school_tertiary", SCHOOL_POOL_YEARS)]


def _load_simple():
    """Non-school trip records with dist/weight/veh/component columns (just-walk dropped)."""
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


def _load_school(years):
    """Education trips with by-car outcome + per-trip prim/sec/ter level-share columns.

    Reuses the generation age→level machinery so the level definition (DfE 16–18 split, FT
    filter, age bands) is literally the same trips the school generation rates use.
    """
    import os
    import sys
    import pandas as pd
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import nts_microdata as nts
    from derive_school_generation import _england_secondary_share, _level_shares

    sec_share = _england_secondary_share()
    ind = nts.load("individual", columns=["SurveyYear", "IndividualID", "Age_B01ID",
                                          "EducN_B01ID"], years=years)
    sh = ind.apply(lambda r: _level_shares(sec_share, r.Age_B01ID, r.EducN_B01ID),
                   axis=1, result_type="expand")
    sh.columns = ["prim", "sec", "ter"]
    ind = pd.concat([ind[["IndividualID"]], sh], axis=1).set_index("IndividualID")

    tr = nts.load("trip", columns=["SurveyYear", "TripPurpose_B01ID", "MainMode_B04ID",
                                   "IndividualID", "TripDisIncSW", "W5", "JJXSC"], years=years)
    tr = tr[tr.TripPurpose_B01ID == EDUCATION_B01].copy()
    tr["dist"] = pd.to_numeric(tr.TripDisIncSW, errors="coerce")
    tr["w"] = tr.JJXSC * tr.W5
    tr["y"] = tr.MainMode_B04ID.isin(SCHOOL_BYCAR_MODES).astype(float)
    tr = tr[tr.dist.notna() & (tr.dist > 0) & (tr.w > 0)]
    for lv in ("prim", "sec", "ter"):
        tr[lv] = tr.IndividualID.map(ind[lv]).fillna(0.0)
    return tr


def _fit_curve(d, y, w, cap=FIT_CAP_MILES):
    """Weighted binomial MLE of PLATEAU*(1-exp(-(d/D0)^K)) over records with d<=cap, w>0."""
    from scipy.optimize import minimize
    m = (d <= cap) & (w > 0)
    d, y, w = d[m].values, y[m].values, w[m].values

    def nll(p):
        pl, d0, k = p
        pr = np.clip(pl * (1.0 - np.exp(-(d / d0) ** k)), 1e-9, 1 - 1e-9)
        return -np.sum(w * (y * np.log(pr) + (1 - y) * np.log(1 - pr)))

    r = minimize(nll, [0.58, 1.0, 1.3], method="L-BFGS-B",
                 bounds=[(0.2, 0.98), (0.1, 30), (0.5, 6)])
    return tuple(r.x)


def _emp(d, y, w, edges):
    """Binned weighted empirical share + per-bin Σw (for the printed table + plot points)."""
    import pandas as pd
    b = pd.cut(d, bins=edges, right=False, labels=False)
    num = pd.Series(w * y).groupby(b).sum()
    den = pd.Series(w).groupby(b).sum()
    share = (num / den).reindex(range(len(edges) - 1))
    return share, den.reindex(range(len(edges) - 1)).fillna(0.0)


def _fit(plot_path="reports/driveshare.png"):
    import os
    import pandas as pd

    edges = np.array([0, .25, .5, .75, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 25, 50, 300.])
    mids = np.array([np.sqrt(edges[i] * edges[i + 1]) if edges[i] > 0 else edges[1] / 2
                     for i in range(len(edges) - 1)])

    def report(title, d, y, w, pl, d0, k):
        e, den = _emp(d, y, w, edges)
        print(f"=== {title}  PLATEAU={pl:.4f} D0={d0:.4f} K={k:.4f} ===")
        print(f"  {'mid_mi':>7} {'emp%':>6} {'fit%':>6} {'Σw':>12}")
        for i in range(len(edges) - 1):
            fv = 100 * pl * (1 - np.exp(-(mids[i] / d0) ** k))
            ev = 100 * e[i] if pd.notna(e[i]) else float("nan")
            print(f"  {mids[i]:7.2f} {ev:6.1f} {fv:6.1f} {den[i]:12,.0f}")
        print()

    fits, emps = {}, {}

    # non-school: driver-share by the B01 component mapping
    print(f"NON-SCHOOL — NTS microdata {YEARS}, veh modes {VEHICLE_MODE_CODES}, "
          f"weight JJXSC×W5, just-walk removed, cap {FIT_CAP_MILES:g} mi\n")
    simple = _load_simple()
    for comp in SIMPLE_COMPONENTS:
        sub = simple[simple.component == comp]
        pl, d0, k = _fit_curve(sub.dist, sub.y, sub.w)
        fits[comp] = (pl, d0, k)
        emps[comp] = _emp(sub.dist, sub.y, sub.w, edges)[0]
        report(f"{comp}  (Σw {sub.w.sum():,.0f}, n {len(sub):,})", sub.dist, sub.y, sub.w, pl, d0, k)

    # school: by-car share of the child's education trip, per level
    print("SCHOOL — by-car share of the child's Education trip, level via generation "
          "age→level machinery (tertiary pooled ex-COVID)\n")
    cache = {}
    for lv_key, comp, years in SCHOOL_LEVELS:
        yk = tuple(years)
        if yk not in cache:
            cache[yk] = _load_school(years)
        tr = cache[yk]
        lw = tr[lv_key] * tr.w
        pl, d0, k = _fit_curve(tr.dist, tr.y, lw)
        fits[comp] = (pl, d0, k)
        emps[comp] = _emp(tr.dist, tr.y, lw, edges)[0]
        report(f"{comp}  (years {years[0]}-{years[-1]}, Σw {lw.sum():,.0f})",
               tr.dist, tr.y, lw, pl, d0, k)

    print("Paste back into CURVES:")
    for comp in SIMPLE_COMPONENTS + [c for _, c, _ in SCHOOL_LEVELS]:
        pl, d0, k = fits[comp]
        print(f'    "{comp}":{" " * max(1, 20 - len(comp))}({pl:.4f}, {d0:.4f}, {k:.4f}),')

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        dd = np.linspace(0.01, 25, 800)
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
        groups = [("non-school", SIMPLE_COMPONENTS,
                   {"commute": "tab:blue", "retail": "tab:green", "res": "tab:red"}),
                  ("school (by-car share of child's education trip)",
                   [c for _, c, _ in SCHOOL_LEVELS],
                   {"school_primary": "tab:purple", "school_postprimary": "tab:orange",
                    "school_tertiary": "tab:brown"})]
        for ax, (title, comps, colors) in zip(axes, groups):
            for comp in comps:
                pl, d0, k = fits[comp]
                c = colors[comp]
                ax.plot(dd, pl * (1 - np.exp(-(dd / d0) ** k)), "-", color=c, lw=2,
                        label=f"{comp}: {pl:.2f}(1-exp(-(d/{d0:.2f})^{k:.2f}))")
                ax.scatter(mids, emps[comp], color=c, s=20, zorder=5, alpha=0.8)
            if title.startswith("non"):
                ax.plot(dd, PLATEAU * (1 - np.exp(-(dd / D0) ** K)), "k--", lw=1.2, alpha=0.6,
                        label="shared (legacy)")
            ax.set_xlim(0, 25)
            ax.set_ylim(0, 1.0)
            ax.set_xlabel("trip length (miles)")
            ax.set_title(title, fontsize=10)
            ax.legend(fontsize=7.5, loc="lower right")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("driveshare (by-car share)")
        fig.suptitle("Per-component driveshare (NTS microdata; points = binned empirical)")
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
        for comp in CURVES:
            vals = "  ".join(f"{d}mi={float(driveshare(d, comp)):.3f}"
                             for d in [0.25, 0.5, 1, 2, 5])
            print(f"  {comp:20s}: {vals}")
