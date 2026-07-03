"""Anchor the gravity willingness kernel from the national TLD / n_Ire(t) divide.

The kernel is `f(c) = driveshare(equiv_miles(c), comp) · W(c)`, c = OSRM seconds
(simulation/model.py `_modesub_kernel`; W is currently the tuned `exp(-c/TAU_c)`).
The intended empirical anchor is `TLD(c) ∝ f(c)·n(c)` ⇒ dividing the national car
trip-length distribution by the national opportunity-geometry n_Ire(t), and dividing
out the fixed empirical driveshare rise, leaves the **willingness** decay:

    W(t) ∝ [TLD(t) / n(t)] / driveshare(equiv_miles(t), comp)

This module computes W(t) per component and fits it two ways — a single exponential
(the anchored `tau_c`, a drop-in for the current kernel) and a double exponential
(short scale fixed, long-scale + weight free — tests whether the tail wants a fatter
form).  **Artifact only** — writes analysis/kernel_fit.json + reports/kernel_fit.png
and prints a table; it does NOT touch model.py or re-tune (the `f≠TLD` caveat: the
divide is a default SHAPE to anchor, refined later by local counts).

Inputs (both on main):
  * analysis/trip_length_distributions.json — 6 car TLDs, MILES bins (share/density/eff_n).
  * data/national_n_of_t.json — 6 n_Ire(t), SECONDS bins (P·A mass per 30 s bin, unnormalised).

Axis reconciliation (the crux).  n(t) is native seconds; the TLD is miles.  We work on
the seconds grid (model + n axis) and map the TLD across analysis/equiv_miles.py WITH the
density Jacobian dd/dt: `TLD_dens_s(t) = TLD_dens_mi(equiv_miles(t)) · d(equiv_miles)/dt`.
driveshare is an exact fixed curve, so dividing it out only rescales — it preserves the
TLD's relative per-bin uncertainty (eff_n), which is the dominant noise (n(t) is sampled
from ~1M+ pairs/purpose, so its per-bin noise is second-order).

School: one shared willingness is the deployed design, so the three levels are each
divided/fit independently and their `tau_L` compared (the shared-τ test); a mass-weighted
shared `tau_school` is reported.

Re-derive:  python3 analysis/fit_kernel.py --build
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from equiv_miles import equiv_miles, _C0, _C1, _C2          # noqa: F401 (constants for Jacobian)
from driveshare import driveshare, CURVES

TLD_FILE = "analysis/trip_length_distributions.json"
NT_FILE = "data/national_n_of_t.json"
TUNED_FILE = "simulation/tuned_params.json"
OUT_JSON = "analysis/kernel_fit.json"
PLOT_PATH = "reports/kernel_fit.png"

NONSCHOOL = ["res", "commute", "retail"]
SCHOOL = ["school_primary", "school_postprimary", "school_tertiary"]

# Fit domain: willingness is a mid/long-range property — drop the short-range rise
# (driveshare owns it; f/driveshare is a noisy 0/0 there) and thin TLD bins.
DRIVESHARE_MIN_FRAC = 0.85    # fit only where driveshare has ~plateaued (pure willingness;
                              # below this the rise leaks into f/driveshare)
EFFN_MIN = 30.0               # TLD mapped-bin effective-n floor
T_REF_S = 600.0              # reference time for plot normalisation (~10 min)


def equiv_miles_deriv(t):
    """dd/dt for d = equiv_miles(t) = exp(C0 + C1·ln t + C2·(ln t)^2)."""
    lt = np.log(t)
    return equiv_miles(t) * (_C1 + 2.0 * _C2 * lt) / t


def _load_inputs():
    tld = json.load(open(TLD_FILE))
    nt = json.load(open(NT_FILE))
    tuned = json.load(open(TUNED_FILE)) if os.path.exists(TUNED_FILE) else {}
    return tld, nt, tuned


def _nt_curve(nt, purpose):
    """Finite-bin seconds centres + per-second density for one purpose's n(t)."""
    edges = np.array(nt["bin_edges_s"], float)
    vals = np.array(nt["n_of_t"][purpose]["n_of_t"], float)      # 481 masses
    fin = np.isfinite(edges[1:])                                 # drop the [14400,inf) overflow
    tc = 0.5 * (edges[:-1][fin] + edges[1:][fin])
    width = np.diff(edges)[fin]
    return tc, vals[fin] / width                                 # centres (s), mass per second


def _tld_interp(tld, comp):
    """Return (density_at_miles(d), effn_at_miles(d)) callables for one component."""
    dist = tld["components"][comp]["distribution"]
    mids = np.array(tld["bin_mid_miles"], float)
    edges = np.array(tld["bin_edges_miles"], float)
    dens = np.array(dist["density"], float)
    effn = np.array(dist["eff_n"], float)
    good = dens > 0
    lm, ld = np.log(mids[good]), np.log(dens[good])              # log-log interpolation

    def dens_at(d):
        return np.exp(np.interp(np.log(d), lm, ld))

    def effn_at(d):
        idx = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, len(effn) - 1)
        return effn[idx]

    return dens_at, effn_at


def willingness(tld, nt, comp):
    """W(t) on the seconds grid for `comp`, plus per-point weight and fit-domain mask."""
    tc, n_dens_s = _nt_curve(nt, comp)
    dens_at, effn_at = _tld_interp(tld, comp)
    d = equiv_miles(tc)                                          # seconds -> miles
    jac = equiv_miles_deriv(tc)                                  # dd/dt
    tld_dens_s = dens_at(d) * jac                                # per-mile -> per-second (Jacobian)
    ds = driveshare(d, comp)
    plateau = CURVES[comp][0]
    with np.errstate(divide="ignore", invalid="ignore"):
        f = tld_dens_s / n_dens_s                               # whole kernel (shape, unnormalised)
        W = f / ds                                             # willingness (driveshare divided out)
    effn = effn_at(d)
    dom = (n_dens_s > 0) & (tld_dens_s > 0) & np.isfinite(W) & (W > 0) \
        & (ds >= DRIVESHARE_MIN_FRAC * plateau) & (effn >= EFFN_MIN)
    return tc, d, W, effn, dom


def fit_single(t, W, wts):
    """Weighted log-linear fit W ∝ exp(-t/tau); returns (tau, logA)."""
    slope, intercept = np.polyfit(t, np.log(W), 1, w=np.sqrt(wts))
    return -1.0 / slope, intercept


def fit_double(t, W, wts):
    """Weighted fit W = A·(w·exp(-t/tau_s) + (1-w)·exp(-t/tau_l)); BOTH scales free.

    Returns (A, w, tau_s, tau_l) ordered tau_s <= tau_l (w = weight on the fast head).
    """
    from scipy.optimize import least_squares
    logW = np.log(W)
    sw = np.sqrt(wts)

    def resid(p):
        logA, wgt, ts, tl = p
        model = np.exp(logA) * (wgt * np.exp(-t / ts) + (1 - wgt) * np.exp(-t / tl))
        return (np.log(np.clip(model, 1e-300, None)) - logW) * sw

    A0 = np.exp(np.average(logW, weights=wts))
    r = least_squares(resid, [np.log(A0), 0.6, 300.0, 4000.0],
                      bounds=([-np.inf, 0.0, 30.0, 100.0], [np.inf, 1.0, 1e5, 1e6]))
    logA, wgt, ts, tl = r.x
    if ts > tl:                                    # order fast head then fat tail
        ts, tl, wgt = tl, ts, 1.0 - wgt
    return float(np.exp(logA)), float(wgt), float(ts), float(tl)


def _wrms(t, W, wts, model):
    """Weighted RMS of log residual for a fitted model callable."""
    r = np.log(model(t)) - np.log(W)
    return float(np.sqrt(np.average(r ** 2, weights=wts)))


def _fit_component(tld, nt, comp):
    tc, d, W, effn, dom = willingness(tld, nt, comp)
    t, Wd, wts = tc[dom], W[dom], effn[dom]
    tau, logA = fit_single(t, Wd, wts)
    A, wgt, tau_s, tau_l = fit_double(t, Wd, wts)
    single = lambda tt: np.exp(logA) * np.exp(-tt / tau)                       # noqa: E731
    double = lambda tt: A * (wgt * np.exp(-tt / tau_s) + (1 - wgt) * np.exp(-tt / tau_l))  # noqa: E731
    return {
        "tau_single_s": float(tau),
        "double": {"w": wgt, "tau_s_s": float(tau_s), "tau_l_s": float(tau_l), "A": A},
        "wrms_single": _wrms(t, Wd, wts, single),
        "wrms_double": _wrms(t, Wd, wts, double),
        "fit_domain_s": [float(t.min()), float(t.max())],
        "n_points": int(dom.sum()),
        # samples (for the JSON + plot); full grid W with the domain mask
        "_t": tc, "_W": W, "_effn": effn, "_dom": dom, "_single": single, "_double": double,
    }


def build():
    tld, nt, tuned = _load_inputs()
    print("Anchoring willingness from TLD / n_Ire(t) — seconds axis, driveshare divided out\n")
    results = {}
    for comp in NONSCHOOL + SCHOOL:
        results[comp] = _fit_component(tld, nt, comp)

    # --- report table ---
    tuned_map = {"res": "TAU_res", "commute": "TAU_commute", "retail": "TAU_retail"}
    print(f"{'component':20s} {'tau_single':>11} {'tuned_tau':>10} {'double: w':>10} "
          f"{'tau_s':>7} {'tau_l':>8} {'wrms 1->2':>14}")
    for comp in NONSCHOOL:
        r = results[comp]
        tt = tuned.get(tuned_map[comp], float("nan"))
        dbl = r["double"]
        print(f"{comp:20s} {r['tau_single_s']:>9.0f}s {tt:>9.0f}s {dbl['w']:>10.2f} "
              f"{dbl['tau_s_s']:>6.0f}s {dbl['tau_l_s']:>7.0f}s "
              f"{r['wrms_single']:>6.3f}->{r['wrms_double']:<6.3f}")
    print("  (tau in seconds; double-exp = fast head tau_s + fat tail tau_l, both free; "
          "wrms = weighted log-resid, single->double)\n")

    # --- school: shared-tau test ---
    print("SCHOOL — shared-tau test (three levels fit independently):")
    tsch = tuned.get("TAU_school", float("nan"))
    taus, masses = [], []
    for comp in SCHOOL:
        r = results[comp]
        dbl = r["double"]
        m = float(nt["n_of_t"][comp]["total_mass"])
        taus.append(r["tau_single_s"]); masses.append(m)
        print(f"  {comp:20s} tau_single {r['tau_single_s']:>6.0f}s   double(w={dbl['w']:.2f}, "
              f"tau_l={dbl['tau_l_s']:.0f}s)   wrms {r['wrms_single']:.3f}->{r['wrms_double']:.3f}")
    tau_shared = float(np.average(taus, weights=masses))
    spread = f"{min(taus):.0f}-{max(taus):.0f}s"
    print(f"  -> per-level tau spread {spread}; mass-weighted shared tau_school ~ {tau_shared:.0f}s "
          f"(current tuned {tsch:.0f}s)\n")

    _write_json(results, nt, tuned, tau_shared)
    _plot(results)


def _write_json(results, nt, tuned, tau_shared):
    out = {"_meta": {
        "purpose": "empirical willingness anchor W(t)=[TLD/n]/driveshare (artifact; not wired)",
        "tld_file": TLD_FILE, "n_of_t_file": NT_FILE,
        "axis": "OSRM seconds; TLD mapped miles->seconds via equiv_miles + Jacobian dd/dt",
        "n_of_t_v1": {k: nt["_meta"].get(k) for k in ("unconstrained", "leg", "mode")},
        "fit_domain_rule": f"driveshare>= {DRIVESHARE_MIN_FRAC}*plateau and TLD eff_n>= {EFFN_MIN}",
        "school_shared_tau_s": tau_shared,
        "finding": "willingness is two-scale: a robust fast head (tau_s ~7-13 min) + a heavier-"
                   "than-single-exp tail. The FORM change (single-exp -> double-exp) is warranted "
                   "regardless of how much the tail value is trusted.",
        "head_trust": "tau_s (fast head) is trustworthy — short-range, well-sampled, insensitive "
                      "to n(t) far-field limitations.",
        "tail_caveat": "tau_l is QUALITATIVE-ONLY. n(t) is v1 UNCONSTRAINED (no 1/D_i — the very "
                       "production constraint that suppresses far trips; dividing by unconstrained "
                       "n mis-attributes that suppression to willingness and inflates the tail), "
                       "plus finite-island truncation and the uncorrected n_Eng/n_Ire geometry "
                       "ratio. Building the constrained (1/D_i) n(t) is what would firm it up.",
        "school_shared_tau_note": "per-level tau_single 235/338/1143 s. Primary vs post-primary "
                                  "(similar geometry) are roughly consistent; tertiary is "
                                  "confounded by its distinctive big-city (university) clustering "
                                  "already visible in its n(t) — so this is weak, not clean, "
                                  "evidence against the shared-tau_school design.",
        "tuned_reference": {k: tuned.get(k) for k in
                            ("TAU_res", "TAU_commute", "TAU_retail", "TAU_school")},
        "tuned_reference_note": "STALE — these tuned tau predate the mode-sub re-tune (K far from "
                                "1); not a fair comparison, shown for reference only.",
    }, "components": {}}
    for comp, r in results.items():
        dom = r["_dom"]
        out["components"][comp] = {
            "tau_single_s": r["tau_single_s"], "double": r["double"],
            "wrms_single": r["wrms_single"], "wrms_double": r["wrms_double"],
            "fit_domain_s": r["fit_domain_s"], "n_points": r["n_points"],
            "W_samples": {"t_s": r["_t"][dom].tolist(), "W": r["_W"][dom].tolist(),
                          "eff_n": r["_effn"][dom].tolist()},
        }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved -> {OUT_JSON}")


def _norm_ref(t, W):
    """Scale W so W(T_REF_S)=1 (log-interp), for readable overlays."""
    ref = np.exp(np.interp(np.log(T_REF_S), np.log(t), np.log(W)))
    return W / ref if ref > 0 else W


def _plot(results):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(plot skipped: {e})")
        return
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    panels = NONSCHOOL + ["school"]
    for ax, comp in zip(axes.flat, panels):
        if comp != "school":
            r = results[comp]; dom = r["_dom"]
            t, W = r["_t"][dom], _norm_ref(r["_t"][dom], r["_W"][dom])
            wts = r["_effn"][dom]
            ref = np.exp(np.interp(np.log(T_REF_S), np.log(r["_t"][dom]), np.log(r["_W"][dom])))
            ax.scatter(t, W, s=3 + 30 * wts / wts.max(), color="tab:gray", alpha=0.5, label="W(t) empirical")
            tt = np.linspace(t.min(), t.max(), 300)
            ax.plot(tt, r["_single"](tt) / ref, "-", color="tab:blue", lw=2,
                    label=f"single exp (τ={r['tau_single_s']:.0f}s)")
            ax.plot(tt, r["_double"](tt) / ref, "--", color="tab:red", lw=1.8,
                    label=f"double (τs={r['double']['tau_s_s']:.0f}s, τl={r['double']['tau_l_s']:.0f}s)")
            ax.set_title(comp, fontsize=11)
        else:
            colors = {"school_primary": "tab:purple", "school_postprimary": "tab:orange",
                      "school_tertiary": "tab:brown"}
            for lvl, c in colors.items():
                r = results[lvl]; dom = r["_dom"]
                t, Wraw = r["_t"][dom], r["_W"][dom]
                ref = np.exp(np.interp(np.log(T_REF_S), np.log(t), np.log(Wraw)))   # W @600s
                ax.scatter(t, Wraw / ref, s=6, color=c, alpha=0.4)
                tt = np.linspace(t.min(), t.max(), 300)
                ax.plot(tt, r["_single"](tt) / ref, "-", color=c, lw=1.8,
                        label=f"{lvl.split('_')[1]} (τ={r['tau_single_s']:.0f}s)")
            ax.set_title("school levels (shared-τ test)", fontsize=11)
        ax.set_yscale("log")
        ax.set_xlabel("travel time c (seconds)")
        ax.set_ylabel("willingness W(c)  (norm. @600s)")
        ax.legend(fontsize=7.5, loc="upper right")
        ax.grid(alpha=0.3)
    fig.suptitle("Empirical willingness W(c) = [TLD/n_Ire] / driveshare  (anchor artifact)")
    fig.tight_layout()
    os.makedirs(os.path.dirname(PLOT_PATH), exist_ok=True)
    fig.savefig(PLOT_PATH, dpi=120)
    print(f"saved {PLOT_PATH}")


if __name__ == "__main__":
    if "--build" in sys.argv:
        build()
    else:
        print(__doc__.split("\n\n")[0])
        print("\nRe-derive with: python3 analysis/fit_kernel.py --build")
