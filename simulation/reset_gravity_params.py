"""
Reset gravity parameters in tuned_params.json to the reference values in
tuner_config.json.  External zone params and the temporal profile (slot_fracs_*)
are preserved unchanged.

Gravity shape params — every key in tuner_config gravity_ref, i.e. the 18 flat
willingness keys `<comp>_taus/_taul/_w` (the double-exp {τs, τl, w} for each of the 6
components incl. the 3 independent school levels; seed them from the TLD/n_Ire divide
with `sync_kernel_anchor.py` first) — are reset to the ref.  The six component scales
K_res/K_commute/K_retail/K_primary/K_postprimary/K_tertiary are reset to 1.0: with
generation pinned (model.compute_generation_scales puts producer weights in
vehicle-driver trips/day) each K_c is a ≈1 verification anchor, and the tuner's convex
solve_scales recomputes them on the first step regardless.  `kernel` is set to
"modesub_double".

Legacy keys — the single-exp willingness (TAU_res/TAU_commute/TAU_retail/TAU_school,
THETA), the 3-component biz keys (K, K_biz, W_BIZ, …), the pre-split single school
scale/shape (K_sch, slot_fracs_school) and dead MU/SIGMA — are stripped so the result
is a clean 6-component double-exp param file.  The `doubly_constrained` list (which
components are ALSO attraction-constrained via Furness) is carried across from
tuner_config.json if present.

Usage:
  python3 simulation/sync_kernel_anchor.py     # seed gravity_ref (18 keys) first
  python3 simulation/reset_gravity_params.py
"""

import json, os

TUNER_CONFIG = "simulation/tuner_config.json"
TUNED_PARAMS = "simulation/tuned_params.json"

with open(TUNER_CONFIG) as f:
    config = json.load(f)

grav_ref     = config["gravity_ref"]
GRAVITY_KEYS = tuple(grav_ref)                       # all shape params in the ref
SCALE_KEYS   = ("K_res", "K_commute", "K_retail",
                "K_primary", "K_postprimary", "K_tertiary")
STALE_KEYS   = ("K", "K_biz", "K_sch", "W_BIZ", "W_SCHOOL", "P_biz", "ALPHA_biz",
                "ALPHA", "ALPHA_commute", "ALPHA_retail", "ALPHA_school",
                "P", "BETA", "P_commute", "BETA_commute", "P_retail", "BETA_retail",
                "P_school", "BETA_school", "MU", "SIGMA", "TAU", "slot_fracs_school",
                # single-exp willingness (superseded by the <comp>_taus/_taul/_w keys):
                "TAU_res", "TAU_commute", "TAU_retail", "TAU_school", "THETA")

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
existing["kernel"] = "modesub_double"

# Doubly-constrained (Furness) components — carried from tuner_config.json (a structural
# model choice: which components are ALSO attraction-constrained, Σ_i T_ij ∝ attractor_j).
# Absent ⇒ singly (production) constrained everywhere.
_dblc = config.get("doubly_constrained")
if _dblc:
    existing["doubly_constrained"] = list(_dblc)
    existing["furness_max_sweeps"] = config.get("furness_max_sweeps", 12)
    print(f"  {'doubly_constrained':<14}  {'':>12}  {sorted(_dblc)}")
    print(f"  {'furness_max_sweeps':<14}  {'':>12}  {existing['furness_max_sweeps']}")
else:
    existing.pop("doubly_constrained", None)
    existing.pop("furness_max_sweeps", None)

with open(TUNED_PARAMS, "w") as f:
    json.dump(existing, f, indent=2)

print(f"\nSaved → {TUNED_PARAMS}")
