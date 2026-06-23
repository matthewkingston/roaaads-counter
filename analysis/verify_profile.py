"""
Fidelity gate: does the offline model agree with REAL OSRM built from the same
profile?

After compile_profile.py writes car_roaaads.lua and the deployed OSRM (:5000) is
rebuilt from it, this re-matches a validation subset of the Google cache through
that real instance and compares the real /match duration to the offline
predict_duration(skeleton, spec). If they agree, the fast benchmark can be
trusted; if not, the offline model (or the skeleton) is unfaithful and the tuned
profile must NOT be adopted until fixed.

Gate (per the plan):
  (i)  median (predict/real) ratio in [0.98, 1.02] for every leg type, and
  (ii) p90 absolute residual below --gate-resid (default 0.05) — i.e. the
       model's unexplained error is smaller than the loss change tuning intends
       to detect, so the tuner is not fitting model noise.

Read-only: OSRM /match (local, free). No Google calls.

  python3 analysis/verify_profile.py --spec simulation/tuned_profile.json --limit 250
"""

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "simulation"))
import profile_spec as ps           # noqa: E402
import skeleton_model as sm         # noqa: E402
from google_routing_common import (  # noqa: E402
    CONF_MIN, decode_polyline, downsample_by_distance, osrm_match)

REPO_ROOT = "/home/matthew/Documents/CodingFun/roaaads"
CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
RAW_DIR = os.path.join(CACHE_DIR, "raw")
SKELETONS_FILE = os.path.join(CACHE_DIR, "skeletons.jsonl")


def _median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2]) if n else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", help="ProfileSpec JSON the deployed OSRM was built from")
    ap.add_argument("--legacy-factors", action="store_true",
                    help="reference spec = HIGHWAY_COST_FACTOR (deployed legacy profile)")
    ap.add_argument("--osrm-url", default="http://localhost:5000")
    ap.add_argument("--skeletons", default=SKELETONS_FILE)
    ap.add_argument("--limit", type=int, default=250,
                    help="validation subset size (every Nth valid route)")
    ap.add_argument("--gate-resid", type=float, default=0.05)
    args = ap.parse_args()

    if args.legacy_factors:
        spec = sm.legacy_spec_from_highway_cost_factor()
    elif args.spec:
        spec = ps.ProfileSpec.load(args.spec)
    else:
        sys.exit("ERROR: pass --spec or --legacy-factors.")

    skels = [s for s in sm.load_skeletons(args.skeletons) if s.get("valid")]
    by_key = {(s["od_id"], s["route_idx"]): s for s in skels}

    # Held-out subset: stride so it spans all leg types/bands.
    stride = max(1, len(skels) // args.limit)
    subset = skels[::stride][:args.limit]
    print(f"Validating {len(subset)} routes against real OSRM at {args.osrm_url} …")

    rows = []
    for s in subset:
        raw_p = os.path.join(RAW_DIR, f"{s['od_id']}.json")
        if not os.path.exists(raw_p):
            continue
        route = json.load(open(raw_p))["routes"][s["route_idx"]]
        coords = downsample_by_distance(
            decode_polyline(route["polyline"]["encodedPolyline"]))
        m = osrm_match(args.osrm_url, coords)
        if m is None:
            continue
        real_dur, _, conf, _ = m
        if conf < CONF_MIN or real_dur <= 0:
            continue
        pred = sm.predict_duration(s, spec)
        rows.append({"leg_type": s["leg_type"], "ratio": pred / real_dur,
                     "resid": abs(math.log(pred / real_dur))})

    if not rows:
        sys.exit("No comparable routes — is the deployed OSRM rebuilt and serving?")

    by_leg = {}
    for r in rows:
        by_leg.setdefault(r["leg_type"], []).append(r)

    print(f"\n{'leg':6} {'n':>5} {'median ratio':>13} {'p90 |resid|':>12}")
    gate_ok = True
    for leg, rs in sorted(by_leg.items()):
        ratios = [x["ratio"] for x in rs]
        resids = sorted(x["resid"] for x in rs)
        med = _median(ratios)
        p90 = resids[min(len(resids) - 1, int(0.9 * len(resids)))]
        leg_ok = (0.98 <= med <= 1.02) and (p90 < args.gate_resid)
        gate_ok = gate_ok and leg_ok
        print(f"{leg:6} {len(rs):>5} {med:>13.3f} {p90:>12.3f}  {'ok' if leg_ok else 'FAIL'}")

    overall_med = _median([r["ratio"] for r in rows])
    overall_p90 = sorted(r["resid"] for r in rows)[min(len(rows) - 1, int(0.9 * len(rows)))]
    print(f"\nOverall: n={len(rows)}  median ratio={overall_med:.3f}  "
          f"p90 |resid|={overall_p90:.3f}")
    print(f"\nGATE: {'PASS — offline model tracks real OSRM' if gate_ok else 'FAIL — do NOT trust the fast loop; fix the model/skeleton'}")
    sys.exit(0 if gate_ok else 2)


if __name__ == "__main__":
    main()
