"""
Aggregate per-session AADT estimates into per-link combined estimates.

Reads data/counts_processed.json and writes data/link_aadt.json.
Always regenerates from scratch — output is a pure function of the input,
so adding new sessions just requires re-running this script.

Combination method: inverse-variance weighted mean (optimal for independent
estimates with known uncertainties).

Usage:
  python3 analysis/aggregate_counts.py
"""

import json, math
from collections import defaultdict

PROCESSED_FILE = "data/counts_processed.json"
LINK_AADT_FILE = "data/link_aadt.json"
MIN_DURATION_S = 60


def link_key(u, v):
    return f"{u},{v}"


with open(PROCESSED_FILE) as f:
    processed = json.load(f)

# ── Collect observations per directed link ────────────────────────────────────

# buckets[(u, v)] = [{"session_id", "direction", "aadt", "aadt_uncertainty"}, ...]
buckets = defaultdict(list)

skipped_duration = 0
skipped_null = 0

for sid, rec in processed["sessions"].items():
    if (rec.get("duration_s") or 0) < MIN_DURATION_S:
        skipped_duration += 1
        continue

    for direction, link_field, aadt_field, unc_field, count_field in (
        ("with",    "matched_link_with",    "with_aadt",    "with_aadt_uncertainty",    "with_count"),
        ("against", "matched_link_against", "against_aadt", "against_aadt_uncertainty", "against_count"),
    ):
        link  = rec.get(link_field)
        aadt  = rec.get(aadt_field)
        unc   = rec.get(unc_field)
        count = rec.get(count_field, None)

        if link is None or aadt is None or unc is None:
            skipped_null += 1
            continue

        u, v = link[0], link[1]
        buckets[(u, v)].append({
            "session_id":       sid,
            "direction":        direction,
            "aadt":             aadt,
            "aadt_uncertainty": unc,
            "time_slot":        rec.get("time_slot"),
            "frac_rel_std":     rec.get("frac_rel_std"),
            "n_eff":            round(count + 0.5, 1) if count is not None else None,
            "duration_s":       rec.get("duration_s"),
        })

# ── Combine each bucket ───────────────────────────────────────────────────────

links_out = {}

for (u, v), obs in sorted(buckets.items()):
    weights = [1.0 / (o["aadt_uncertainty"] ** 2) for o in obs]
    sum_w   = sum(weights)
    aadt    = round(sum(w * o["aadt"] for w, o in zip(weights, obs)) / sum_w)
    sigma   = round(1.0 / math.sqrt(sum_w))

    links_out[link_key(u, v)] = {
        "aadt":             aadt,
        "aadt_uncertainty": sigma,
        "n_observations":   len(obs),
        "observations":     obs,
    }

# ── Write output ──────────────────────────────────────────────────────────────

with open(LINK_AADT_FILE, "w") as f:
    json.dump({"links": links_out}, f, indent=2)

max_obs = max(v["n_observations"] for v in links_out.values()) if links_out else 0
print(f"{len(links_out)} directed link(s)  |  "
      f"{skipped_duration} session(s) skipped (<{MIN_DURATION_S}s)  |  "
      f"{skipped_null} direction(s) skipped (null)  |  "
      f"max {max_obs} obs on one link")
print(f"Saved → {LINK_AADT_FILE}")
