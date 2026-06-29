"""
Reset gravity parameters in tuned_params.json to the reference values in
tuner_config.json.  External zone params and the temporal profile (slot_fracs_*)
are preserved unchanged.

Gravity shape params — every key in tuner_config gravity_ref, i.e. P, ALPHA,
BETA, P_commute, ALPHA_commute, P_retail, ALPHA_retail, THETA, P_school,
ALPHA_school — are reset to the ref.  The four component scales
K_res/K_commute/K_retail/K_sch are reset to 1.0: with generation pinned
(model.compute_generation_scales puts producer weights in vehicle-driver
trips/day) each K_c is a ≈1 verification anchor, and the tuner's convex
solve_scales recomputes them on the first step regardless.

Legacy 3-component biz keys (K, K_biz, W_BIZ, W_SCHOOL, P_biz, ALPHA_biz) and
dead MU/SIGMA are stripped so the result is a clean 4-component param file.

Usage:
  python3 simulation/reset_gravity_params.py
"""

import json, os

TUNER_CONFIG = "simulation/tuner_config.json"
TUNED_PARAMS = "simulation/tuned_params.json"

with open(TUNER_CONFIG) as f:
    config = json.load(f)

grav_ref     = config["gravity_ref"]
GRAVITY_KEYS = tuple(grav_ref)                       # all shape params in the ref
SCALE_KEYS   = ("K_res", "K_commute", "K_retail", "K_sch")
STALE_KEYS   = ("K", "K_biz", "W_BIZ", "W_SCHOOL", "P_biz", "ALPHA_biz",
                "ALPHA", "ALPHA_commute", "ALPHA_retail", "ALPHA_school",
                "MU", "SIGMA")

existing = {}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        existing = json.load(f)


def _fmt(v):
    return f"{v:.6g}" if isinstance(v, float) else (str(v) if v is not None else "—")


print("Resetting gravity params (tuner_config gravity_ref → tuned_params.json):\n")
print(f"  {'param':<14}  {'before':>12}  {'after':>12}")
for key in SCALE_KEYS + GRAVITY_KEYS:
    after = 1.0 if key in SCALE_KEYS else grav_ref.get(key)
    print(f"  {key:<14}  {_fmt(existing.get(key)):>12}  {_fmt(after):>12}")
for key in STALE_KEYS:
    if key in existing:
        print(f"  {key:<14}  {_fmt(existing.get(key)):>12}  {'(removed)':>12}")

existing.update({k: grav_ref[k] for k in GRAVITY_KEYS})
for k in SCALE_KEYS:
    existing[k] = 1.0
for k in STALE_KEYS:
    existing.pop(k, None)

with open(TUNED_PARAMS, "w") as f:
    json.dump(existing, f, indent=2)

print(f"\nSaved → {TUNED_PARAMS}")
