"""
Reset gravity parameters in tuned_params.json to the reference values in
tuner_config.json.  External zone params are preserved unchanged.

K is reset to 1.0 — the tuner recomputes it analytically on the first step.

Usage:
  python3 simulation/reset_gravity_params.py
"""

import json, os

TUNER_CONFIG = "simulation/tuner_config.json"
TUNED_PARAMS = "simulation/tuned_params.json"

with open(TUNER_CONFIG) as f:
    config = json.load(f)

grav_ref = config["gravity_ref"]

existing = {}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        existing = json.load(f)

GRAVITY_KEYS = ("W_BIZ", "MU", "SIGMA", "ALPHA")

print("Resetting gravity params (tuner_config gravity_ref → tuned_params.json):\n")
print(f"  {'param':<8}  {'before':>12}  {'after':>12}")
for key in ("K",) + GRAVITY_KEYS:
    before = existing.get(key, "—")
    after  = 1.0 if key == "K" else grav_ref[key]
    before_str = f"{before:.6g}" if isinstance(before, float) else str(before)
    print(f"  {key:<8}  {before_str:>12}  {after:>12.6g}")

existing.update({k: grav_ref[k] for k in GRAVITY_KEYS})
existing["K"] = 1.0

with open(TUNED_PARAMS, "w") as f:
    json.dump(existing, f, indent=2)

print(f"\nSaved → {TUNED_PARAMS}")
