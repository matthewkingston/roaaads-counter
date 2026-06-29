"""Vehicle-driver mode share as a function of trip length (miles).

The gravity kernel's short-range rise is *mode substitution*: short trips are
walked, not driven, so car demand is suppressed at short cost — a mode/network
property, shared across trip purposes and empirical, NOT a per-purpose "desire".
`driveshare(d_miles)` is that factor; the kernel is
    f(c) = driveshare(equiv_miles(c)) * exp(-c / tau_c)
(see simulation/model.py `_modesub_kernel`; `equiv_miles` maps OSRM seconds→miles).

Derived from DfT NTS0308a (trips by trip length × main mode, England, 2023/24),
vehicle-driver modes = car/van driver + motorcycle + taxi/minicab (same set as the
generation rates), as a share of all destination journeys, with three corrections:

  1. JUST-WALK removed — NTS "Just walk" is a recreational purpose with no
     destination, not a point-to-point journey, so it is excluded from the
     mode-split denominator (consistent with generation, which also drops it).
     There is no purpose×length cross-tab on disk, so the just-walk total
     (NTS0409a, all-modes) is allocated across length bands ∝ the walk-mode length
     profile.  *** This allocation is the single biggest assumption here — if
     just-walk skews longer than utility walking, the short-band lift is slightly
     overstated. ***
  2. PLATEAU-CLIPPED — the England decline past ~50 mi (rail/coach substitution) is
     clipped flat: NI+RoI rail is far sparser, so car share holds at its plateau.
     (All long-range suppression then comes from the willingness exp(-c/tau).)
  3. ORIGIN-PINNED — fit so driveshare(0)=0 (⇒ f(0)=0), with the <1 mi band treated
     as an integral (band-average matching) so the *pointwise* value goes to ~0 at
     the origin rather than sitting at the band mean (otherwise a 40 m car trip
     would get spurious appeal, since the willingness exp is ~1 there).

Form: driveshare(d) = PLATEAU * (1 - exp(-(d/D0)^K)).  PLATEAU is interpretive only
— it is a shared constant factor that cancels in the production constraint.

The module constants below are the *single* authoritative source of truth; the
imported `driveshare` does no file I/O.  Re-derive with `--fit` after the NTS files
change (prints refreshed constants to paste back, and draws the plot).
"""

import numpy as np

# --- Fitted constants (NTS0308a 2023/24; just-walk removed, plateau-clipped) -----
# Re-derive with: python3 analysis/driveshare.py --fit
PLATEAU = 0.5779    # vehicle-driver share ceiling (cancels in the production constraint)
D0      = 1.2573    # scale (miles)
K       = 1.2610    # shape (>1 ⇒ gentle origin start)


def driveshare(d_miles):
    """Vehicle-driver share for a trip of length d_miles. Float or numpy array."""
    return PLATEAU * (1.0 - np.exp(-(np.asarray(d_miles, dtype=float) / D0) ** K))


# --- Re-fit / diagnostics / plot (offline only; not on the import path) ---------

NTS0308_FILE = "data/nts0308.ods"
NTS0409_FILE = "data/nts0409.ods"
NTS_YEARS    = [2023, 2024]
VEH_MODES    = ["Car or van driver", "Motorcycle", "Taxi or minicab"]
BANDS = ["Under 1 mile", "1 to under 2 miles", "2 to under 5 miles", "5 to under 10 miles",
         "10 to under 25 miles", "25 to under 50 miles", "50 to under 100 miles",
         "100 miles and over"]
EDGES = np.array([0, 1, 2, 5, 10, 25, 50, 100, 300.])   # last band capped at 300 mi


def _band_data():
    """(veh, total, walk) trips/person/yr per length band + just-walk total."""
    import pandas as pd
    df = pd.read_excel(NTS0308_FILE, sheet_name="NTS0308a_trips", header=5, engine="odf")
    df.columns = [str(c).strip() for c in df.columns]
    yc, mc = df.columns[0], df.columns[1]
    df[yc] = pd.to_numeric(df[yc], errors="coerce")

    def col(name):
        return [c for c in df.columns if c == name or c.startswith(name)][0]

    def msum(prefix, band):
        s = df[df[yc].isin(NTS_YEARS) & df[mc].astype(str).str.strip().str.startswith(prefix)]
        return s[col(band)].astype(float).sum() / len(NTS_YEARS)

    veh   = np.array([sum(msum(v, b) for v in VEH_MODES) for b in BANDS])
    total = np.array([msum("All modes", b) for b in BANDS])
    walk  = np.array([msum("Walk", b) for b in BANDS])      # 'Walk [notes 2, 3]'

    # just-walk total (all modes) from NTS0409a
    d9 = pd.read_excel(NTS0409_FILE, sheet_name="NTS0409a_trips", header=5, engine="odf")
    d9.columns = [str(c).strip() for c in d9.columns]
    y9, m9 = d9.columns[0], d9.columns[1]
    d9[y9] = pd.to_numeric(d9[y9], errors="coerce")
    jwcol = [c for c in d9.columns if c.startswith("Just walk")][0]
    am = d9[d9[y9].isin(NTS_YEARS) & (d9[m9].astype(str).str.strip() == "All modes")]
    just_walk = am[jwcol].astype(float).sum() / len(NTS_YEARS)
    return veh, total, walk, just_walk


def _fit(plot_path="reports/driveshare.png"):
    from scipy.optimize import least_squares
    veh, total, walk, just_walk = _band_data()
    # (1) remove just-walk, allocated ∝ walk-mode length profile
    jw_band = just_walk * walk / walk.sum()
    total_c = total - jw_band
    ds_raw  = veh / total
    ds_corr = veh / total_c
    ds_clip = np.maximum.accumulate(ds_corr)             # (2) plateau-clip

    # (3) origin-pinned band-AVERAGE least-squares fit of PLATEAU*(1-exp(-(d/D0)^K))
    def ds_of(p, d):
        pl, d0, k = p
        return pl * (1.0 - np.exp(-(d / d0) ** k))

    def band_avg(p, lo, hi, n=400):
        xs = np.linspace(lo, hi, n)
        return np.trapz(ds_of(p, xs), xs) / (hi - lo)

    def resid(p):
        return np.array([band_avg(p, EDGES[i], EDGES[i + 1]) for i in range(len(BANDS))]) - ds_clip

    sol = least_squares(resid, [0.58, 2.0, 1.5], bounds=([0.2, 0.1, 0.5], [0.9, 30, 6]))
    pl, d0, k = sol.x

    print(f"NTS0308a {NTS_YEARS}, vehicle modes = {'+'.join(VEH_MODES)}")
    print(f"just-walk removed: {just_walk:.1f}/yr (allocated ∝ walk-mode length profile)\n")
    print(f"{'mid_mi':>7} {'raw%':>6} {'corr%':>6} {'clip%':>6} {'fit%':>6}")
    for i, b in enumerate(BANDS):
        mid = np.sqrt(EDGES[i] * EDGES[i + 1]) if EDGES[i] > 0 else 0.5
        print(f"{mid:7.1f} {100*ds_raw[i]:6.1f} {100*ds_corr[i]:6.1f} {100*ds_clip[i]:6.1f} "
              f"{100*band_avg(sol.x, EDGES[i], EDGES[i+1]):6.1f}")
    print(f"\nPaste back into the module constants:")
    print(f"  PLATEAU = {pl:.4f}\n  D0      = {d0:.4f}\n  K       = {k:.4f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        dd = np.linspace(0.01, 30, 1000)
        mids = [np.sqrt(EDGES[i] * EDGES[i + 1]) if EDGES[i] > 0 else 0.5 for i in range(len(BANDS))]
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(dd, ds_of(sol.x, dd), "b-", lw=2,
                label=f"fit {pl:.2f}*(1-exp(-(d/{d0:.2f})^{k:.2f}))")
        ax.scatter(mids, ds_corr, c="green", zorder=5, label="band (just-walk removed)")
        ax.scatter(mids, ds_clip, c="red", marker="x", zorder=6, label="band (plateau-clipped)")
        ax.set_xlim(0, 30); ax.set_ylim(0, 0.7)
        ax.set_xlabel("trip length (miles)"); ax.set_ylabel("vehicle-driver share")
        ax.set_title("driveshare(distance): NTS0308, just-walk removed, plateau-clipped, origin-pinned")
        ax.legend(fontsize=8); ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig(plot_path, dpi=110)
        print(f"\nsaved {plot_path}")
    except Exception as e:
        print(f"\n(plot skipped: {e})")


if __name__ == "__main__":
    import sys
    if "--fit" in sys.argv:
        _fit()
    else:
        print("driveshare(d_miles) module. Re-derive constants with: python3 analysis/driveshare.py --fit")
        for d in [0.05, 0.25, 0.5, 1, 2, 3, 5, 10]:
            print(f"  driveshare({d:5.2f} mi) = {float(driveshare(d)):.4f}")
