"""Journey-time -> equivalent road-distance (miles), as a cheap closed form.

The gravity model works in least-time seconds, but surveyed trip-length
distributions (NTS/TSNI) are published in miles. `equiv_miles(t_seconds)` maps a
journey *time* to the road *distance* a typical journey of that duration covers,
so a seconds-based kernel can be compared against mile-based TLDs.

Average speed is strongly trip-length-dependent (~17 mph for short in-town trips
rising to ~49 mph on long inter-urban legs), so this is NOT a constant
`miles = v * t`. The relationship is fit from the Google routing cache (the same
free-flow `TRAFFIC_UNAWARE` times used for the OSRM profile calibration), Google
best route per OD: `g_dur` (seconds) vs `g_dist` (metres -> miles).

Chosen form: log-log quadratic
    ln(miles) = C0 + C1 * ln(t) + C2 * (ln t)^2
i.e. equiv_miles(t) = exp(C0 + C1*ln t + C2*(ln t)^2).

One log + one exp -> cheap enough for the tuning/eval hot loop, and it vectorises
over a numpy array of OD times. Monotone increasing across the whole realistic
domain (d ln(mi)/d ln(t) = C1 + 2*C2*ln t stays > 0), so no clamp is needed.

The module constants below are the authoritative source of truth: the imported
`equiv_miles` does no file I/O. Re-derive them with `--fit` after the Google cache
changes (it prints the refreshed constants to paste back in).
"""

import numpy as np

# --- Fitted coefficients (Google best-route, 989 points, span 119-8141 s) -------
# Re-derive with: python3 analysis/equiv_miles.py --fit
# log-log quadratic: ln(miles) = C0 + C1*ln(t) + C2*(ln t)^2
_C0 = -5.6726
_C1 = 0.9591
_C2 = 0.02319

# Simpler 2-param power-law alternative: miles = _POW_A * t**_POW_B
# (kept for reference / quick analytic use; mild short-trip bias vs the quadratic).
_POW_A = 1.04e-03
_POW_B = 1.2955

MILES_PER_M = 1609.344


def equiv_miles(t_seconds):
    """Equivalent road distance (miles) for a journey time in seconds.

    Accepts a Python float or a numpy array (returns the same shape).
    """
    lt = np.log(t_seconds)
    return np.exp(_C0 + _C1 * lt + _C2 * lt * lt)


def equiv_miles_pow(t_seconds):
    """Power-law variant: miles = _POW_A * t**_POW_B. Float or numpy array."""
    return _POW_A * np.power(t_seconds, _POW_B)


# --- Re-fit / diagnostics / plot (offline only; not on the import path) ---------

def _load_google_points(results_path):
    """Google best-route (t_seconds, miles) points from results.jsonl."""
    import json

    pts = []
    with open(results_path) as fh:
        for line in fh:
            r = json.loads(line)
            best = r.get("google_best_dur")
            if not best:
                continue
            for rt in r.get("routes", []):
                if (
                    rt.get("valid")
                    and rt.get("g_dur")
                    and rt.get("g_dist")
                    and abs(rt["g_dur"] - best) < 1e-6
                ):
                    pts.append((rt["g_dur"], rt["g_dist"] / MILES_PER_M))
                    break
    return pts


def _fit(results_path, plot_path, json_path):
    import datetime
    import json

    pts = _load_google_points(results_path)
    t = np.array([p[0] for p in pts])
    miles = np.array([p[1] for p in pts])
    lt = np.log(t)
    ly = np.log(miles)

    # log-log quadratic (numpy returns highest-degree coeff first)
    c2, c1, c0 = np.polyfit(lt, ly, 2)
    # power law (degree-1 in log-log)
    b, lna = np.polyfit(lt, ly, 1)
    a = float(np.exp(lna))

    def quad(x):
        L = np.log(x)
        return np.exp(c0 + c1 * L + c2 * L * L)

    def powf(x):
        return a * np.power(x, b)

    def diagnostics(name, pred):
        lr = np.abs(np.log(pred(t) / miles))
        med = float(np.median(lr))
        bands = [(0, 300), (300, 600), (600, 1200), (1200, 2400), (2400, 4800), (4800, 9000)]
        cells = []
        for lo, hi in bands:
            m = (t >= lo) & (t < hi)
            if m.any():
                cells.append(f"{lo}-{hi}:{np.median(pred(t[m]) / miles[m]):.2f}")
        print(f"  {name:10s} median|log-ratio|={med:.4f} ({np.exp(med) - 1:+.1%})")
        print("    band pred/actual: " + "  ".join(cells))
        return med

    print(f"n Google best-route points: {len(pts)}  span {t.min():.0f}-{t.max():.0f} s\n")
    print("log-log quadratic: ln(miles) = C0 + C1*ln t + C2*(ln t)^2")
    print(f"  _C0 = {c0:.4f}\n  _C1 = {c1:.4f}\n  _C2 = {c2:.5f}")
    print(f"power law:         miles = {a:.3e} * t**{b:.4f}  (_POW_A, _POW_B)\n")
    med_quad = diagnostics("quadratic", quad)
    diagnostics("power", powf)
    speeds = ", ".join(
        f"{float(quad(np.array([s]))[0]) / (s / 3600):.0f}"
        for s in (120, 300, 600, 1200, 2400, 4800, 8000)
    )
    print(f"\n  implied mph @120/300/600/1200/2400/4800/8000 s: {speeds}")

    # provenance record
    import os

    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, "w") as fh:
        json.dump(
            {
                "form": "ln(miles) = C0 + C1*ln(t_s) + C2*(ln t_s)^2",
                "C0": float(c0),
                "C1": float(c1),
                "C2": float(c2),
                "pow_a": a,
                "pow_b": float(b),
                "n_points": len(pts),
                "t_span_s": [float(t.min()), float(t.max())],
                "median_abs_log_ratio": med_quad,
                "source": os.path.relpath(results_path),
                "basis": "Google best route, TRAFFIC_UNAWARE free-flow",
                "fitted": datetime.date.today().isoformat(),
            },
            fh,
            indent=2,
        )
    print(f"\nwrote {json_path}")

    # plot
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        grid = np.logspace(np.log10(t.min()), np.log10(t.max()), 200)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(t, miles, s=8, alpha=0.3, color="#446", label=f"Google best route (n={len(pts)})")
        ax.plot(grid, quad(grid), color="#c33", lw=2.2, label="log-log quadratic (fit)")
        ax.plot(grid, powf(grid), color="#393", lw=1.4, ls="--", label="power law (ref)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("journey time (s)")
        ax.set_ylabel("road distance (miles)")
        ax.set_title("equiv_miles(t_seconds) — Google free-flow")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(plot_path, dpi=120)
        print(f"wrote {plot_path}")
    except ImportError:
        print("matplotlib not available — skipping plot")


def _main():
    import argparse
    import os

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fit", action="store_true", help="re-derive coefficients + plot from the cache")
    ap.add_argument("--results", default=os.path.join(repo, "data/google_cache/results.jsonl"))
    ap.add_argument("--plot", default=os.path.join(repo, "reports/equiv_miles.png"))
    ap.add_argument("--json", default=os.path.join(repo, "reports/equiv_miles_fit.json"))
    args = ap.parse_args()

    if args.fit:
        _fit(args.results, args.plot, args.json)
    else:
        for s in (120, 300, 600, 1200, 2400, 4800, 8000):
            print(f"  equiv_miles({s:5d} s) = {float(equiv_miles(s)):6.2f} miles")


if __name__ == "__main__":
    _main()
