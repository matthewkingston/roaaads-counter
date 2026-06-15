"""
Restore tuned_params.json from a specific entry in tuning_history.jsonl.

For gravity-only history entries, only the 5 gravity params (K, W_BIZ, MU,
SIGMA, ALPHA) are written; existing external zone entries in tuned_params.json
are left unchanged.  For full-stage entries, all tuned params are restored.

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
        date  = e.get("timestamp", "")[:10]
        note  = e.get("note", "")
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

# ── Restore ───────────────────────────────────────────────────────────────────

existing = {}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        existing = json.load(f)

params = entry["params"]

GRAVITY_KEYS = ("K", "W_BIZ", "MU", "SIGMA", "ALPHA")
ALL_KEYS = list(params.keys())

print(f"Restoring from run {entry['id']}  "
      f"({entry['stage']}, {entry['timestamp'][:10]}, χ²/N={entry['chi2_per_n']:.3f})")
if entry.get("note"):
    print(f"  note: {entry['note']}")
print()

# Print before/after for scalar gravity params
print(f"  {'param':<8}  {'before':>12}  {'after':>12}")
for key in GRAVITY_KEYS:
    before = existing.get(key, "—")
    after  = params[key]
    before_str = f"{before:.6g}" if isinstance(before, float) else str(before)
    print(f"  {key:<8}  {before_str:>12}  {after:>12.6g}")

non_scalar = [k for k in ALL_KEYS if k not in GRAVITY_KEYS]
if non_scalar:
    print(f"\n  External zone params restored: {', '.join(non_scalar)}")

existing.update(params)
existing["chi2"]       = entry["chi2"]
existing["chi2_per_n"] = entry["chi2_per_n"]
existing["n_obs"]      = entry["n_obs"]
existing["stage"]      = entry["stage"]

with open(TUNED_PARAMS, "w") as f:
    json.dump(existing, f, indent=2)

print(f"\nSaved → {TUNED_PARAMS}")
