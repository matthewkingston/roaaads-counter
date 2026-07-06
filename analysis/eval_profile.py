"""
Benchmark entrypoint: score a candidate OSRM profile against the Google cache,
offline and instantly.

Reads the skeleton cache (simulation/build_skeleton_index.py) and a ProfileSpec
JSON (default = all factors 1.0 = stock OSRM), prints the aggregate log-ratio
loss, the predicted/Google ratio distribution, a per-leg-type and per-cell
breakdown, the per-bucket coverage (which factors are well-determined), and the
turn-time fraction. No OSRM, no Docker, no spend.

  python3 analysis/eval_profile.py                          # stock profile (factors=1)
  python3 analysis/eval_profile.py --spec simulation/tuned_profile.json
  python3 analysis/eval_profile.py --legacy-factors         # reproduce deployed car_roaaads.lua
  python3 analysis/eval_profile.py --all-routes             # include conf<0.5 / off-coverage

The --legacy-factors run is the model-faithfulness sanity check: results.jsonl
was matched under the deployed HIGHWAY_COST_FACTOR profile, so the predicted
median ratio here should track the te_matched the runner reported.
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "simulation"))
import profile_spec as ps           # noqa: E402
import skeleton_model as sm         # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKELETONS_FILE = os.path.join(REPO_ROOT, "data", "google_cache", "skeletons.jsonl")
BASE_SPEEDS_FILE = os.path.join(REPO_ROOT, "data", "google_cache", "base_speeds.json")


def _ratio_pctiles(rows):
    rs = sorted(r["ratio"] for r in rows)
    if not rs:
        return {}
    def pct(p):
        return rs[min(len(rs) - 1, int(p * len(rs)))]
    return {"p10": pct(0.10), "p25": pct(0.25), "p50": pct(0.50),
            "p75": pct(0.75), "p90": pct(0.90)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", help="ProfileSpec JSON (default: all factors 1.0)")
    ap.add_argument("--skeletons", default=SKELETONS_FILE)
    ap.add_argument("--legacy-factors", action="store_true",
                    help="use HIGHWAY_COST_FACTOR per class (deployed-profile reference)")
    ap.add_argument("--all-routes", action="store_true",
                    help="include invalid routes (low conf / off-coverage)")
    ap.add_argument("--coverage-buckets", type=int, default=20,
                    help="how many top buckets to print in the coverage table")
    ap.add_argument("--no-empirical", action="store_true",
                    help="ignore base_speeds.json; use analytical base speeds")
    args = ap.parse_args()

    if not os.path.exists(args.skeletons):
        sys.exit(f"ERROR: {args.skeletons} not found — run build_skeleton_index.py first.")
    skels = sm.load_skeletons(args.skeletons)

    n_emp = 0 if args.no_empirical else ps.load_empirical_base_speeds(BASE_SPEEDS_FILE)
    print(f"Base speeds: {'empirical (%d buckets)' % n_emp if n_emp else 'analytical'}")

    if args.legacy_factors:
        spec = sm.legacy_spec_from_highway_cost_factor()
        label = "legacy HIGHWAY_COST_FACTOR"
    elif args.spec:
        spec = ps.ProfileSpec.load(args.spec)
        label = args.spec
    else:
        spec = ps.ProfileSpec.default()
        label = "stock (factors=1.0)"

    res = sm.evaluate(skels, spec, valid_only=not args.all_routes)
    n_total = len(skels)
    print(f"Profile: {label}")
    print(f"Skeletons: {n_total} routes ({res['n_valid']} scored"
          f"{'' if args.all_routes else ', valid only'})")
    print(f"\nAggregate log-ratio loss  : {res['loss']:.4f}")
    print(f"Median predicted/Google   : {res['median_ratio']:.3f}  "
          f"(>1 = model slower than Google, <1 = model faster)")
    pc = _ratio_pctiles(res["rows"])
    if pc:
        print(f"Ratio pctiles p10..p90    : "
              + "  ".join(f"{k}={v:.3f}" for k, v in pc.items()))

    print("\nPer leg type:")
    print(f"  {'leg':6} {'n':>5} {'med ratio':>10} {'loss':>8} {'turn%':>7}")
    for leg, st in res["leg_stats"].items():
        print(f"  {leg:6} {st['n']:>5} {st['median_ratio']:>10.3f} "
              f"{st['loss']:>8.4f} {100*st['turn_frac']:>6.1f}%")

    print("\nPer (leg x band) median ratio:")
    for (leg, band), st in res["breakdown"].items():
        print(f"  {leg:6} band{band}  n={st['n']:>4}  "
              f"med={st['median_ratio']:.3f}  |lnr|={st['mean_abs_lr']:.3f}")

    cov, unbk = sm.bucket_coverage(skels, valid_only=not args.all_routes)
    total_m = sum(cov.values()) + unbk
    print(f"\nBucket coverage (matched km; total {total_m/1000:.0f} km, "
          f"unbucketed {unbk/1000:.1f} km):")
    for key, m in sorted(cov.items(), key=lambda kv: -kv[1])[:args.coverage_buckets]:
        fac = spec.factors.get(key, 1.0)
        print(f"  {key:24} {m/1000:>8.1f} km   factor={fac:.3f}")
    thin = sum(1 for m in cov.values() if m < 1000)
    print(f"  ({len(cov)} buckets populated; {thin} below 1 km — weakly determined)")


if __name__ == "__main__":
    main()
