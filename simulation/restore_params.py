"""
Restore tuned_params.json from a specific entry in tuning_history.jsonl.

All params stored in the history entry's 'params' dict are written directly to
tuned_params.json (existing keys are overwritten; unrelated keys are preserved).
For gravity-only runs the external zone params are absent from history, so any
existing external zone values in tuned_params.json are retained unchanged.

Usage:
  python3 simulation/restore_params.py --list          # show all runs
  python3 simulation/restore_params.py <id>            # restore by ID (partial OK)
"""

import json, os, sys

HISTORY_FILE = "simulation/tuning_history.jsonl"
TUNED_PARAMS = "simulation/tuned_params.json"

# ── Load history ──────────────────────────────────────────────────────────────

entries = []
if os.path.exists(HISTORY_FILE):
    with open(HISTORY_FILE) as f:
        for line in f:
            line = line.strip()
            if line:
                entries.append(json.loads(line))

# ── --list ────────────────────────────────────────────────────────────────────

if "--list" in sys.argv or len(sys.argv) < 2:
    if not entries:
        print("No runs in history yet.")
        sys.exit(0)
    print(f"  {'#':>3}  {'ID':<10}  {'Date':<10}  {'Stage':<8}  {'N':>4}  {'χ²/N':>6}  {'git':<8}  Note")
    for i, e in enumerate(entries, 1):
        date = e.get("timestamp", "")[:10]
        note = e.get("note", "")
        print(f"  {i:>3}  {e['id']:<10}  {date:<10}  {e['stage']:<8}"
              f"  {e['n_obs']:>4}  {e['chi2_per_n']:>6.3f}  {e['git_hash']:<8}  {note}")
    sys.exit(0)

# ── Resolve ID ────────────────────────────────────────────────────────────────

arg = sys.argv[1]
matches = [e for e in entries if e["id"].startswith(arg)]

if len(matches) == 0:
    print(f"No history entry matches '{arg}'.")
    print("Run with --list to see available IDs.")
    sys.exit(1)

if len(matches) > 1:
    print(f"Ambiguous: '{arg}' matches {len(matches)} entries:")
    for e in matches:
        print(f"  {e['id']}  {e['timestamp'][:10]}  {e['stage']}  χ²/N={e['chi2_per_n']:.3f}")
    print("Provide more characters.")
    sys.exit(1)

entry = matches[0]

if entry["params"].get("kernel") != "rational":
    print("WARNING: this history entry uses the old lognormal kernel. "
          "Restoring it will produce incorrect results with the current "
          "rational-kernel pipeline. Aborting.")
    sys.exit(1)

# ── Restore ───────────────────────────────────────────────────────────────────

existing = {}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        existing = json.load(f)

params = entry["params"]

print(f"Restoring from run {entry['id']}  "
      f"({entry['stage']}, {entry['timestamp'][:10]}, χ²/N={entry['chi2_per_n']:.4f})")
if entry.get("note"):
    print(f"  note: {entry['note']}")
hp = entry.get("tuner_hyperparams", {})
if hp:
    print(f"  hyperparams: " +
          "  ".join(f"{k}={v}" for k, v in hp.items()))
print()

# Scalar gravity params: show before/after
SCALAR_KEYS = [k for k in ("K", "K_res", "K_biz", "K_sch", "W_BIZ", "P", "ALPHA", "BETA",
                            "P_biz", "ALPHA_biz", "THETA",
                            "W_SCHOOL", "P_school", "ALPHA_school")
               if k in params]
if SCALAR_KEYS:
    print(f"  {'param':<8}  {'before':>14}  {'after':>14}")
    for key in SCALAR_KEYS:
        before = existing.get(key, "—")
        after  = params[key]
        before_s = f"{before:.6g}" if isinstance(before, (int, float)) else str(before)
        print(f"  {key:<8}  {before_s:>14}  {after:>14.6g}")

# Report what else is being written
n_slot_r = len(params.get("slot_fracs_res",    {}))
n_slot_b = len(params.get("slot_fracs_biz",    {}))
n_slot_s = len(params.get("slot_fracs_school", {}))
if n_slot_r or n_slot_b or n_slot_s:
    _sch_str = f" + {n_slot_s} school" if n_slot_s else ""
    print(f"\n  Slot fracs: {n_slot_r} res + {n_slot_b} biz{_sch_str} slots")
ext_keys = [k for k in params if k.startswith("external_")]
if ext_keys:
    print(f"  External zone params: {', '.join(ext_keys)}")

existing.update(params)
for _stale in ("MU", "SIGMA", "slot_fracs"):
    existing.pop(_stale, None)
existing["chi2"]       = entry["chi2"]
existing["chi2_per_n"] = entry["chi2_per_n"]
existing["n_obs"]      = entry["n_obs"]
existing["stage"]      = entry["stage"]

with open(TUNED_PARAMS, "w") as f:
    json.dump(existing, f, indent=2)

print(f"\nSaved → {TUNED_PARAMS}")
