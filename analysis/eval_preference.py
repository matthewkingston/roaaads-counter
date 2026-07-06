"""
Offline route-preference benchmark.

Scores a tuned_preference.json against the skeleton cache: what fraction of
violation pairs (Google chose a slower route) are correctly ranked by the
preference multipliers?  Also checks for regressions on concordant pairs.

No OSRM calls, no spend.

  python3 analysis/eval_preference.py                                 # default pref file
  python3 analysis/eval_preference.py --pref simulation/tuned_preference.json
  python3 analysis/eval_preference.py --unit                          # score p_c=1 (baseline)
"""

import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "simulation"))
import profile_spec as ps          # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
SKELETONS_FILE = os.path.join(CACHE_DIR, "skeletons.jsonl")
BASE_SPEEDS_FILE = os.path.join(CACHE_DIR, "base_speeds.json")
DEFAULT_PREF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "simulation", "tuned_preference.json")

SPLIT_CLASSES = {"trunk", "primary", "secondary", "tertiary"}
MAIN_CLASSES = [
    "motorway",
    "trunk", "trunk_rural",
    "primary", "primary_rural",
    "secondary", "secondary_rural",
    "tertiary", "tertiary_rural",
    "unclassified", "residential", "living_street", "service",
]
LINK_PARENT = {
    "motorway_link": "motorway",
    "trunk_link": "trunk",
    "primary_link": "primary",
    "secondary_link": "secondary",
    "tertiary_link": "tertiary",
}
OSRM_SPEED_KMH = {
    "motorway": 90, "motorway_link": 45,
    "trunk": 85, "trunk_link": 40,
    "primary": 65, "primary_link": 30,
    "secondary": 55, "secondary_link": 25,
    "tertiary": 40, "tertiary_link": 20,
    "unclassified": 25, "residential": 25,
    "living_street": 10, "service": 15,
}
URBAN_KMH = 48.3


def load_skeletons(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def _pref_key(raw_cls, cls_key, band):
    if cls_key not in SPLIT_CLASSES:
        return cls_key
    try:
        eff_kmh = int(band) * 1.60934
    except (ValueError, TypeError):
        eff_kmh = OSRM_SPEED_KMH.get(raw_cls, OSRM_SPEED_KMH.get(cls_key, 25))
    return cls_key + ("_rural" if eff_kmh > URBAN_KMH else "")


def route_class_times(skel, base_speeds):
    times = collections.defaultdict(float)
    for bucket, metres in skel.get("length_by_bucket", {}).items():
        raw_cls, band = bucket.split("|", 1)
        cls_key = LINK_PARENT.get(raw_cls, raw_cls)
        if cls_key not in SPLIT_CLASSES and cls_key not in MAIN_CLASSES:
            cls_key = "unclassified"
        key = _pref_key(raw_cls, cls_key, band)
        if key not in MAIN_CLASSES:
            key = cls_key
        spd_info = base_speeds.get(bucket)
        if isinstance(spd_info, dict):
            spd = spd_info.get("harmonic_kmh") or 0.0
        else:
            spd = spd_info or 0.0
        if spd <= 0:
            cls_part, band_part = bucket.split("|", 1)
            spd = ps.base_speed_for(cls_part, band_part)
        times[key] += metres / (spd / 3.6)
    return dict(times)


def route_cost(class_times, pref):
    """Preference-weighted routing cost = sum_c p_c * time_c."""
    return sum(pref.get(cls, 1.0) * t for cls, t in class_times.items())


def _median(xs):
    if not xs:
        return None
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skeletons", default=SKELETONS_FILE)
    ap.add_argument("--pref", default=DEFAULT_PREF,
                    help="tuned_preference.json (default: simulation/tuned_preference.json)")
    ap.add_argument("--unit", action="store_true",
                    help="score with p_c=1 for all classes (baseline)")
    args = ap.parse_args()

    ps.load_empirical_base_speeds(BASE_SPEEDS_FILE)
    with open(BASE_SPEEDS_FILE) as f:
        base_speeds = json.load(f)

    if args.unit:
        pref = {}
        print("Preference: unit (p_c=1 baseline)")
    else:
        if not os.path.exists(args.pref):
            sys.exit(f"ERROR: {args.pref} not found. Run tune_preference.py first.")
        with open(args.pref) as f:
            pref_data = json.load(f)
        pref = pref_data["preference"]
        meta = pref_data.get("meta", {})
        print(f"Preference: {args.pref}")
        if meta:
            print(f"  Tuned {meta.get('ts','')}  pairs={meta.get('n_pairs')}  "
                  f"lam={meta.get('lam')}  log_margin={meta.get('log_margin')}")
        print(f"  Non-unit classes: "
              + ", ".join(f"{c}={v:.3f}" for c, v in sorted(pref.items())
                          if abs(v - 1.0) > 1e-4))
    print()

    skels = load_skeletons(args.skeletons)
    by_od = collections.defaultdict(dict)
    for s in skels:
        by_od[s["od_id"]][s["route_idx"]] = s

    # Score all multi-route ODs.
    # Categorise by (g_dur violation) × (offline violation) to correctly separate
    # true preference violations from timing-model divergences.
    # Unit offline_dur = sum of class times at p=1 (proxy for OSRM routing cost).
    true_violations = []     # g_dur AND offline both say rk is faster; Google chose r0
    timing_diverge  = []     # g_dur says rk faster but offline already routes r0
    concordant_ok   = []     # g_dur says r0 faster, offline agrees
    concordant_bad  = []     # g_dur says r0 faster but offline routes rk (timing error)
    leg_viol = collections.defaultdict(lambda: {"total": 0, "resolved": 0})

    for od, routes in by_od.items():
        if 0 not in routes or len(routes) < 2:
            continue
        r0 = routes[0]
        alts = {k: routes[k] for k in routes if k > 0 and routes[k].get("g_dur", 0) > 0}
        if not alts or not r0.get("g_dur", 0):
            continue
        if not r0.get("valid"):
            continue

        fastest_k = min(alts, key=lambda k: alts[k]["g_dur"])
        rk = alts[fastest_k]
        if not rk.get("valid"):
            continue

        r0_times = route_class_times(r0, base_speeds)
        rk_times = route_class_times(rk, base_speeds)
        offline_r0 = sum(r0_times.values())
        offline_rk = sum(rk_times.values())

        c0 = route_cost(r0_times, pref)
        ck = route_cost(rk_times, pref)
        correctly_ranked = c0 < ck

        g_violation = r0["g_dur"] > rk["g_dur"]
        offline_violation = offline_r0 > offline_rk   # OSRM would pick rk

        entry = {
            "od_id": od,
            "leg_type": r0["leg_type"],
            "g_margin_s": r0["g_dur"] - rk["g_dur"],
            "offline_margin_s": offline_r0 - offline_rk,
            "cost_margin_s": c0 - ck,
            "correctly_ranked": correctly_ranked,
        }

        if g_violation and offline_violation:
            true_violations.append(entry)
            leg_viol[r0["leg_type"]]["total"] += 1
            if correctly_ranked:
                leg_viol[r0["leg_type"]]["resolved"] += 1
        elif g_violation:
            timing_diverge.append(entry)   # offline already routes r0 correctly
        elif offline_violation:
            concordant_bad.append(entry)   # timing model error (not a preference issue)
        else:
            concordant_ok.append(entry)

    # --- Report true violations (the calibration target) ---
    n_true = len(true_violations)
    n_resolved = sum(1 for e in true_violations if e["correctly_ranked"])
    g_margins = sorted(e["g_margin_s"] for e in true_violations)
    print(f"=== True preference violations (g_dur AND offline both route rk; Google chose r0) ===")
    print(f"  N = {n_true}  (these are the calibration target; baseline = 0/{n_true} correct)")
    print(f"  Correctly ranked after preference: "
          f"{n_resolved}/{n_true} "
          f"({100*n_resolved/n_true if n_true else 0:.1f}%)")
    if g_margins:
        print(f"  G-duration margin: median={_median(g_margins):.0f}s  "
              f"mean={sum(g_margins)/len(g_margins):.0f}s  max={g_margins[-1]:.0f}s")
    cost_gaps = sorted(e["cost_margin_s"] for e in true_violations)
    print(f"  Cost gap after pref (r0-rk, <0=resolved): "
          f"median={_median(cost_gaps):.1f}s  mean={sum(cost_gaps)/len(cost_gaps):.1f}s")
    print()
    print("  Per leg type:")
    for leg in sorted(leg_viol):
        s = leg_viol[leg]
        pct = 100 * s["resolved"] / s["total"] if s["total"] else 0
        print(f"    {leg}: {s['resolved']}/{s['total']} ({pct:.0f}%)")

    # --- Report timing divergences (not violations; OSRM already routes r0) ---
    n_td = len(timing_diverge)
    td_still_ok = sum(1 for e in timing_diverge if e["correctly_ranked"])
    print()
    print(f"=== Timing divergences (g_dur says rk faster; offline already routes r0) ===")
    print(f"  N = {n_td}  (OSRM already correct on timing; no preference fix needed)")
    print(f"  Still correctly ranked after preference: {td_still_ok}/{n_td} "
          f"({100*td_still_ok/n_td if n_td else 0:.1f}%)")
    if n_td - td_still_ok:
        print(f"  WARNING: {n_td - td_still_ok} flipped wrong by preference factors")

    # --- Report concordant pairs ---
    n_ok = len(concordant_ok)
    n_bad = len(concordant_bad)
    bad_still_wrong = sum(1 for e in concordant_bad if not e["correctly_ranked"])
    ok_flipped = sum(1 for e in concordant_ok if not e["correctly_ranked"])
    print()
    print(f"=== Concordant pairs (Google chose fastest by g_dur) ===")
    print(f"  Offline also correct: {n_ok}  |  Offline wrong (timing error): {n_bad}")
    if ok_flipped:
        print(f"  WARNING: preference flipped {ok_flipped}/{n_ok} previously-correct concordant pairs")
    if n_bad:
        print(f"  Timing errors still wrong after preference: {bad_still_wrong}/{n_bad} "
              f"(expected — preference shouldn't compensate for timing errors)")


if __name__ == "__main__":
    main()
