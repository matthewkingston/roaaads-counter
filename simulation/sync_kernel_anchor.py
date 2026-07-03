"""Seed the tuner's willingness anchor from the TLD ÷ n_Ire(t) kernel fit.

Reads the per-component double-exp willingness fit (`analysis/kernel_fit.json`, from
`analysis/fit_kernel.py`) and writes the 18 flat willingness keys into
`simulation/tuner_config.json`'s `gravity_ref` — the anchor / start point + L2-pull target
for `tune_assignment.py`.  Also (re)initialises `gravity_lambda` to a light uniform value on
those 18 keys if it isn't already willingness-keyed, and strips the dead single-exp
`TAU_*`/`THETA`/`TAU_school` keys.

The willingness kernel is `f(c)=driveshare(equiv_miles(c),comp)·[w·exp(−c/τs)+(1−w)·exp(−c/τl)]`;
only the shape `{w, τs, τl}` is carried (the fit's amplitude `A` is absorbed by K in the
production constraint).

**PATCH POINT for the constrained n_Ire iteration:** when that agent delivers its improved
double-exp params, repoint `ANCHOR_FILE` and adjust `_read_anchor` to its format — everything
downstream keys off the 18 `{comp}_taus/_taul/_w` names in `model.willingness_keys()`.

Usage:  python3 simulation/sync_kernel_anchor.py [--lambda 0.2]
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import WILLINGNESS_COMPONENTS, willingness_keys

ANCHOR_FILE = "analysis/kernel_fit.json"       # ← repoint to the constrained-n_Ire output later
TUNER_CONFIG = "simulation/tuner_config.json"
DEFAULT_LAMBDA = 0.2                            # light anchor-reg (per willingness param, internal coords)


def _read_anchor(path):
    """Return {component: (w, τs, τl)} from the kernel-fit JSON's per-component `double` block.
    Adjust this one function when the anchor source/format changes (the n_Ire hand-off)."""
    d = json.load(open(path))
    comps = d["components"]
    out = {}
    for c in WILLINGNESS_COMPONENTS:
        b = comps[c]["double"]
        out[c] = (float(b["w"]), float(b["tau_s_s"]), float(b["tau_l_s"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lambda", dest="lam", type=float, default=DEFAULT_LAMBDA,
                    help="anchor-reg strength written to gravity_lambda if not already willingness-keyed")
    args = ap.parse_args()

    anchor = _read_anchor(ANCHOR_FILE)
    keys = willingness_keys()

    # Flat gravity_ref (natural units).
    ref = {}
    for c in WILLINGNESS_COMPONENTS:
        w, tau_s, tau_l = anchor[c]
        ref[f"{c}_taus"] = round(tau_s, 4)
        ref[f"{c}_taul"] = round(tau_l, 4)
        ref[f"{c}_w"]    = round(w, 6)

    cfg = json.load(open(TUNER_CONFIG))
    cfg["gravity_ref"] = ref
    # (Re)seed gravity_lambda on the willingness keys only if it isn't already (preserve a
    # hand-tuned per-key lambda across re-syncs; migrate away from the dead TAU_* keys otherwise).
    lam = cfg.get("gravity_lambda")
    if not (isinstance(lam, dict) and all(k in lam for k in keys)):
        cfg["gravity_lambda"] = {k: args.lam for k in keys}
        print(f"  gravity_lambda (re)initialised to {args.lam} on all 18 willingness keys")
    else:
        # keep existing per-key values but drop any stale non-willingness keys
        cfg["gravity_lambda"] = {k: lam[k] for k in keys}
        print("  gravity_lambda: kept existing per-key values (stale keys dropped)")

    with open(TUNER_CONFIG, "w") as f:
        json.dump(cfg, f, indent=2)

    print(f"Seeded gravity_ref from {ANCHOR_FILE} → {TUNER_CONFIG}")
    for c in WILLINGNESS_COMPONENTS:
        w, tau_s, tau_l = anchor[c]
        print(f"  {c:20s} w={w:.3f}  τs={tau_s:.0f}s  τl={tau_l:.0f}s")


if __name__ == "__main__":
    main()
