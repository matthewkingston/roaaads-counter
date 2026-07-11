"""Seed the tuner's willingness anchor + its uncertainty from the constrained kernel fit.

Reads the per-component constrained double-exp willingness fit (`analysis/kernel_fit_constrained.json`,
from `analysis/iterate_kernel.py`) and writes into `simulation/tuner_config.json`:
  * `gravity_ref` — the 18 flat willingness keys `{comp}_taus/_taul/_w`, the anchor / start point
    for `tune_assignment.py`;
  * `gravity_stat_cov` — the per-component 3×3 statistical covariance `Cov_stat` of the anchor, in
    the tuner's internal coords (log τs, log τl, logit w), from the fit's `double_iterated_cov`
    (Jacobian for converged components, robust trace-spread for the weakly-identified school tails).

The tuner combines `Cov_stat` with the `anchor_floor` (two scalars — a head-tight / tail-loose
epistemic floor) into the per-component prior precision block `Λ_c = (Cov_stat + Cov_epi)⁻¹`; the
old scalar/per-key `gravity_lambda` regularizer and the hard `gravity_fixed` pins are RETIRED and
stripped here.  `anchor_floor` is source config: seeded to defaults only if absent (user edits
preserved), then adjustable by hand — one edit + re-tune, no re-sync.

The willingness kernel is `f(c)=driveshare(equiv_miles(c),comp)·[w·exp(−c/τs)+(1−w)·exp(−c/τl)]`;
only the shape `{w, τs, τl}` is carried (the fit's amplitude `A` is absorbed by K in the
production constraint).

Usage:  python3 simulation/sync_kernel_anchor.py [--head-sigma 0.182] [--tail-sigma 0.693]
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import WILLINGNESS_COMPONENTS, willingness_keys

ANCHOR_FILE = "analysis/kernel_fit_constrained.json"   # constrained (1/D_i) iterated double-exp
TUNER_CONFIG = "simulation/tuner_config.json"
DEFAULT_HEAD_SIGMA = 0.182     # ±20% on τs  (internal log-coord 1σ) → λ_head ≈ 30
DEFAULT_TAIL_SIGMA = 0.693     # ×÷2 on τl and on tail-mass (1−w)    → λ_tail ≈ 2.1


def _read_anchor(path):
    """Return ({component: (w, τs, τl)}, {component: 3×3 Cov_stat}) from the constrained kernel-fit
    JSON's per-component `double_iterated` (anchor) + `double_iterated_cov` (internal-coord σ_stat)
    blocks.  Adjust this one function if the anchor source/format changes.  NB: school_postprimary/
    tertiary carry `tail_weakly_identified=True` (thin long-tail TLD) — their `double_iterated_cov`
    is the robust trace-spread (loose τl/w), so the tuner's derived prior lets them run."""
    d = json.load(open(path))
    comps = d["components"]
    anchor, cov = {}, {}
    for c in WILLINGNESS_COMPONENTS:
        b = comps[c]["double_iterated"]
        anchor[c] = (float(b["w"]), float(b["tau_s_s"]), float(b["tau_l_s"]))
        cov[c] = comps[c]["double_iterated_cov"]["cov"]      # 3×3, coords (log τs, log τl, logit w)
    return anchor, cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--head-sigma", dest="head_sigma", type=float, default=DEFAULT_HEAD_SIGMA,
                    help="anchor_floor head σ (τs) seeded only if anchor_floor is absent")
    ap.add_argument("--tail-sigma", dest="tail_sigma", type=float, default=DEFAULT_TAIL_SIGMA,
                    help="anchor_floor tail σ (τl, w) seeded only if anchor_floor is absent")
    args = ap.parse_args()

    anchor, cov = _read_anchor(ANCHOR_FILE)

    # Flat gravity_ref (natural units) + per-component internal-coord statistical covariance.
    ref = {}
    for c in WILLINGNESS_COMPONENTS:
        w, tau_s, tau_l = anchor[c]
        ref[f"{c}_taus"] = round(tau_s, 4)
        ref[f"{c}_taul"] = round(tau_l, 4)
        ref[f"{c}_w"]    = round(w, 6)

    cfg = json.load(open(TUNER_CONFIG))
    cfg["gravity_ref"] = ref
    cfg["gravity_stat_cov"] = {c: cov[c] for c in WILLINGNESS_COMPONENTS}
    # anchor_floor = the two epistemic-trust knobs (source config); seed defaults only if absent so
    # hand edits survive a re-sync.  Retire the old regularizers.
    if "anchor_floor" not in cfg:
        cfg["anchor_floor"] = {"head_sigma": args.head_sigma, "tail_sigma": args.tail_sigma}
        print(f"  anchor_floor seeded: head_sigma={args.head_sigma} tail_sigma={args.tail_sigma}")
    else:
        print(f"  anchor_floor: kept existing {cfg['anchor_floor']}")
    for dead in ("gravity_lambda", "gravity_fixed"):
        if cfg.pop(dead, None) is not None:
            print(f"  stripped retired key: {dead}")

    with open(TUNER_CONFIG, "w") as f:
        json.dump(cfg, f, indent=2)

    print(f"Seeded gravity_ref + gravity_stat_cov from {ANCHOR_FILE} → {TUNER_CONFIG}")
    for c in WILLINGNESS_COMPONENTS:
        w, tau_s, tau_l = anchor[c]
        sd = [cov[c][i][i] ** 0.5 for i in range(3)]         # σ per internal coord (τs, τl, w)
        print(f"  {c:20s} w={w:.3f}  τs={tau_s:.0f}s  τl={tau_l:.0f}s   "
              f"σ_stat(logτs,logτl,logitw)=({sd[0]:.2f},{sd[1]:.2f},{sd[2]:.2f})")


if __name__ == "__main__":
    main()
