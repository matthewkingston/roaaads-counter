"""
External-focused profile tuner.

Fits the per-(class x speed-band) speed **factors** to minimise the offline
log-ratio time error vs Google over the skeleton cache, weighted toward the
external corridors (X2B / B2X / X2X) that drive flow into the core. INT is
down-weighted to zero by default: verify_profile showed the offline turn model
under-counts in-town junctions, so INT is not yet a trustworthy tuning target.

Only **bucket factors** are tuned here; the global turn params are held at their
defaults. On the external legs turn time is a small fraction (~3-4%) and INT is
excluded, so the identifiable signal is the road-class x speed-band factors.
(Turn-cost tuning waits until the in-town turn model is improved.)

With turn params and base speeds fixed, predicted time is *linear* in the factor
vector, so this is a smooth low-dim fit:

    predict_r = const_r + sum_b  factor_b * (len_rb * 3.6 / base_speed_b)

where const_r = fixed turn time + edge time of non-tuned buckets. Only buckets
with enough (weighted) matched distance are tuned; the rest stay at 1.0.

  python3 analysis/tune_profile.py                       # default external weights
  python3 analysis/tune_profile.py --min-km 100 --lam 0.02
  python3 analysis/tune_profile.py --leg-weights X2B=1,B2X=1,X2X=1,INT=0.2
  python3 analysis/tune_profile.py --init simulation/tuned_profile.json
"""

import argparse
import collections
import datetime
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "simulation"))
import profile_spec as ps          # noqa: E402
import skeleton_model as sm         # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
SKELETONS_FILE = os.path.join(CACHE_DIR, "skeletons.jsonl")
BASE_SPEEDS_FILE = os.path.join(CACHE_DIR, "base_speeds.json")
OUT_SPEC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "simulation", "tuned_profile.json")
HISTORY = os.path.join(CACHE_DIR, "profile_tuning_history.jsonl")

DEFAULT_LEG_WEIGHTS = {"X2B": 1.0, "B2X": 1.0, "X2X": 1.0, "INT": 0.0}
FACTOR_LO, FACTOR_HI = 0.2, 5.0


def parse_leg_weights(s):
    if not s:
        return dict(DEFAULT_LEG_WEIGHTS)
    out = dict(DEFAULT_LEG_WEIGHTS)
    for part in s.split(","):
        k, v = part.split("=")
        out[k.strip()] = float(v)
    return out


def build_matrices(skels, leg_w, min_km):
    """Return (active_buckets, Ls, const, g, w, legs, rows).

    Ls[r,b] = len_rb * 3.6 / base_speed_b (active buckets only);
    const_r = edge time of non-active buckets (factor 1) + unbucketed + turn time
    (default turn params, factor-independent); predict_r = const_r + Ls[r] @ factor.
    """
    rows = [s for s in skels
            if s.get("valid") and s.get("g_dur", 0) > 0 and leg_w.get(s["leg_type"], 0) > 0]

    # active buckets = enough weighted matched distance to be identifiable
    cov = collections.defaultdict(float)
    for s in rows:
        wl = leg_w[s["leg_type"]]
        for k, m in s["length_by_bucket"].items():
            cov[k] += wl * m
    active = sorted(k for k, m in cov.items() if m / 1000.0 >= min_km)
    aidx = {k: i for i, k in enumerate(active)}
    nB = len(active)

    R = len(rows)
    Ls = np.zeros((R, nB))
    const = np.zeros(R)
    g = np.zeros(R)
    w = np.zeros(R)
    legs = []
    default_spec = ps.ProfileSpec()       # default turn params, factors 1.0
    for r, s in enumerate(rows):
        g[r] = s["g_dur"]
        w[r] = leg_w[s["leg_type"]]
        legs.append(s["leg_type"])
        c = 0.0
        for k, m in s["length_by_bucket"].items():
            cls, band = k.split("|", 1)
            coef = m * 3.6 / ps.base_speed_for(cls, band)
            if k in aidx:
                Ls[r, aidx[k]] += coef
            else:
                c += coef
        ub = s.get("unbucketed_m", 0.0)
        if ub:
            c += ub * 3.6 / ps.STOCK_SPEED_KMH["unclassified"]
        c += sm.predict_components(s, default_spec)[1]   # turn time (factor-free)
        const[r] = c
    return active, Ls, const, g, w, np.array(legs), rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skeletons", default=SKELETONS_FILE)
    ap.add_argument("--leg-weights", help="e.g. X2B=1,B2X=1,X2X=1,INT=0")
    ap.add_argument("--min-km", type=float, default=50.0,
                    help="min weighted matched km to tune a bucket (else fixed 1.0)")
    ap.add_argument("--lam", type=float, default=0.01,
                    help="L2 regularisation of ln(factor) toward 0 (factor 1.0)")
    ap.add_argument("--init", help="ProfileSpec JSON to start from (default: factors 1.0)")
    ap.add_argument("--out", default=OUT_SPEC)
    ap.add_argument("--note", default="")
    ap.add_argument("--dry-run", action="store_true", help="don't write spec/history")
    args = ap.parse_args()

    leg_w = parse_leg_weights(args.leg_weights)
    n_emp = ps.load_empirical_base_speeds(BASE_SPEEDS_FILE)
    print(f"Base speeds: {'empirical (%d buckets)' % n_emp if n_emp else 'analytical'}")
    skels = sm.load_skeletons(args.skeletons)

    active, Ls, const, g, w, legs, rows = build_matrices(skels, leg_w, args.min_km)
    lng = np.log(g)
    wsum = w.sum()
    print(f"Tuning {len(active)} buckets over {len(rows)} routes "
          f"(weights {leg_w}); lambda={args.lam}")

    lo, hi = math.log(FACTOR_LO), math.log(FACTOR_HI)

    # initial factor vector
    theta0 = np.zeros(len(active))
    if args.init and os.path.exists(args.init):
        init = ps.ProfileSpec.load(args.init)
        for i, k in enumerate(active):
            f = init.factors.get(k, 1.0)
            theta0[i] = np.clip(math.log(f), lo, hi)

    def loss(theta):
        pred = const + Ls @ np.exp(theta)
        resid = np.log(pred) - lng
        data = float((w * resid * resid).sum() / wsum)
        reg = args.lam * float((theta * theta).sum())
        return data + reg

    l0 = loss(theta0)
    res = minimize(loss, theta0, method="L-BFGS-B",
                   bounds=[(lo, hi)] * len(active),
                   options={"maxiter": 2000, "ftol": 1e-10})
    theta = res.x
    factors = {k: float(round(math.exp(theta[i]), 4)) for i, k in enumerate(active)}

    # report: data-only loss (no reg) before/after + per-leg medians via evaluate
    def data_loss(th):
        pred = const + Ls @ np.exp(th)
        resid = np.log(pred) - lng
        return float((w * resid * resid).sum() / wsum)

    spec = ps.ProfileSpec(factors={k: v for k, v in factors.items() if abs(v - 1.0) > 1e-3})
    print(f"\nweighted data loss: {data_loss(theta0):.4f} -> {data_loss(theta):.4f} "
          f"(reg+data {l0:.4f} -> {res.fun:.4f})")

    # per-leg median ratio before/after (all legs, for context; INT shown but not fit)
    base_spec = ps.ProfileSpec()
    for label, sp in (("before", base_spec), ("after", spec)):
        ev = sm.evaluate(skels, sp)
        legs_str = "  ".join(f"{lg}={st['median_ratio']:.3f}"
                             for lg, st in sorted(ev["leg_stats"].items()))
        print(f"  {label:6} median ratio: {legs_str}")

    # biggest factor moves by weighted coverage
    cov = collections.defaultdict(float)
    for s in rows:
        for k, m in s["length_by_bucket"].items():
            cov[k] += leg_w[s["leg_type"]] * m
    print("\ntop factor moves (by coverage):")
    print(f"  {'bucket':24} {'factor':>7} {'km':>8}")
    for k in sorted(factors, key=lambda x: -cov.get(x, 0))[:15]:
        print(f"  {k:24} {factors[k]:>7.3f} {cov.get(k,0)/1000:>8.0f}")

    if args.dry_run:
        print("\n[dry-run] not writing spec/history")
        return
    spec.save(args.out)
    print(f"\nWrote {os.path.abspath(args.out)} ({len(spec.factors)} non-unit factors)")
    with open(HISTORY, "a") as f:
        f.write(json.dumps({
            "ts": datetime.datetime.now().isoformat(timespec="seconds"),
            "note": args.note, "leg_weights": leg_w, "min_km": args.min_km,
            "lam": args.lam, "n_routes": len(rows), "n_buckets": len(active),
            "data_loss_before": data_loss(theta0), "data_loss_after": data_loss(theta),
            "factors": factors,
        }) + "\n")
    print(f"Logged to {HISTORY}")
    print("\nNext: compile_profile.py --spec %s -> rebuild :5000 -> verify_profile.py" % args.out)


if __name__ == "__main__":
    main()
