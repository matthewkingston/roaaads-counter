"""
Offline route-preference calibration.

Fits per-highway-class preference multipliers p_c such that OSRM would route
through Google's preferred route rather than its (sometimes faster) alternatives.

Uses a scale-invariant log-ratio ranking loss:

    hinge_i = max(0, log(cost(r0_i) / cost(rk_i)) + log_margin)²

where cost(r) = p @ times_r  (dot product of class multipliers and class times).

log(cost(r0)/cost(rk)) < 0 means r0 is preferred.  The log form is scale-
invariant: scaling all p_c by a constant k cancels in the ratio, so the
optimizer must find *relative* class preferences rather than collapsing all p_c
toward the lower bound.

The p_c multipliers are applied to forward_rate only (not forward_speed), so
route choice changes but reported durations are unchanged.  p_c < 1 makes a
class preferred; p_c > 1 makes it avoided.

  python3 analysis/tune_preference.py            # fit + write tuned_preference.json
  python3 analysis/tune_preference.py --dry-run  # report pair counts, don't write
  python3 analysis/tune_preference.py --lam 0.02 --log-margin 0.03
"""

import argparse
import collections
import datetime
import json
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "simulation"))
import profile_spec as ps          # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE_DIR = os.path.join(REPO_ROOT, "data", "google_cache")
SKELETONS_FILE = os.path.join(CACHE_DIR, "skeletons.jsonl")
BASE_SPEEDS_FILE = os.path.join(CACHE_DIR, "base_speeds.json")
OUT_PREF = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "simulation", "tuned_preference.json")
HISTORY = os.path.join(CACHE_DIR, "preference_tuning_history.jsonl")

# Main highway classes to fit.
# Split classes get two entries: base (urban, ≤30mph) and _rural (>30mph).
# Urban/rural is determined by the tagged maxspeed; for untagged roads the
# OSRM class default speed is used as a fallback (same logic OSRM applies
# when no maxspeed tag is present: forward_speed = speed_profile[class]).
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
# OSRM speed_profile defaults (km/h) — used when no maxspeed tag is present.
# Sourced from deployed car_roaaads.lua speed_profile table.
OSRM_SPEED_KMH = {
    "motorway": 90, "motorway_link": 45,
    "trunk": 85, "trunk_link": 40,
    "primary": 65, "primary_link": 30,
    "secondary": 55, "secondary_link": 25,
    "tertiary": 40, "tertiary_link": 20,
    "unclassified": 25, "residential": 25,
    "living_street": 10, "service": 15,
}
URBAN_KMH = 48.3   # 30 mph — effective speed above this → rural suffix

PREF_LO, PREF_HI = 0.33, 3.0


def load_skeletons(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def _pref_key(raw_cls, cls_key, band):
    """Preference lookup key: 'class' for urban/non-split, 'class_rural' for rural split."""
    if cls_key not in SPLIT_CLASSES:
        return cls_key
    # Determine effective speed: tagged maxspeed for known bands, OSRM default for untagged/other.
    try:
        eff_kmh = int(band) * 1.60934   # e.g. "60" → 96.5 km/h
    except (ValueError, TypeError):
        # untagged or other: fall back to OSRM class default (using raw_cls for links)
        eff_kmh = OSRM_SPEED_KMH.get(raw_cls, OSRM_SPEED_KMH.get(cls_key, 25))
    return cls_key + ("_rural" if eff_kmh > URBAN_KMH else "")


def route_class_times(skel, base_speeds):
    """Return dict {pref_key: seconds} from a skeleton.

    pref_key is 'class' for urban/non-split-class segments and 'class_rural'
    for segments on split classes (trunk/primary/secondary/tertiary) where the
    tagged maxspeed (or OSRM class default for untagged roads) exceeds 30 mph.
    """
    times = collections.defaultdict(float)
    for bucket, metres in skel.get("length_by_bucket", {}).items():
        raw_cls, band = bucket.split("|", 1)
        cls_key = LINK_PARENT.get(raw_cls, raw_cls)
        if cls_key not in SPLIT_CLASSES and cls_key not in MAIN_CLASSES:
            cls_key = "unclassified"   # fallback for road/other/etc.
        key = _pref_key(raw_cls, cls_key, band)
        if key not in MAIN_CLASSES:
            key = cls_key   # shouldn't happen but guard anyway
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


def offline_dur(class_times):
    """Total offline route duration = sum of class times."""
    return sum(class_times.values())


def build_pairs(skels, base_speeds):
    """
    Return violation pairs for preference calibration.

    A violation is a case where:
      1. Google chose r0 over a faster alternative rk by Google's own timing
         (g_dur(r0) > g_dur(rk)), AND
      2. The offline/OSRM model also says rk is faster (offline_dur(r0) > offline_dur(rk))
         — so OSRM as a pure time-minimiser would also route via rk.

    Cases where offline_dur(r0) < offline_dur(rk) despite g_dur(r0) > g_dur(rk) are
    OSRM/Google timing divergences: OSRM already routes via r0 on timing alone, so
    no preference factor is needed.  Concordant pairs where offline disagrees are
    timing model errors; preference should not compensate for those.
    """
    by_od = collections.defaultdict(dict)
    for s in skels:
        by_od[s["od_id"]][s["route_idx"]] = s

    pairs = []
    g_violation_osrm_ok = 0    # Google violation but OSRM already picks r0
    concordant = 0
    invalid_skip = 0
    for od, routes in by_od.items():
        if 0 not in routes or len(routes) < 2:
            continue
        r0 = routes[0]
        alts = {k: routes[k] for k in routes if k > 0 and routes[k].get("g_dur", 0) > 0}
        if not alts or not r0.get("g_dur", 0):
            continue
        fastest_k = min(alts, key=lambda k: alts[k]["g_dur"])
        rk = alts[fastest_k]

        if r0["g_dur"] <= rk["g_dur"]:
            concordant += 1
            continue

        # Google violation: r0 is slower by g_dur
        if not r0.get("valid") or not rk.get("valid"):
            invalid_skip += 1
            continue

        r0_times = route_class_times(r0, base_speeds)
        rk_times = route_class_times(rk, base_speeds)

        if offline_dur(r0_times) <= offline_dur(rk_times):
            # OSRM/Google timing divergence: offline already routes via r0
            g_violation_osrm_ok += 1
            continue

        # True preference violation: both g_dur and offline_dur say rk is faster,
        # but Google chose r0.
        pairs.append({
            "od_id": od,
            "leg_type": r0["leg_type"],
            "r0_times": r0_times,
            "rk_times": rk_times,
            "g_dur_r0": r0["g_dur"],
            "g_dur_rk": rk["g_dur"],
            "margin_s": r0["g_dur"] - rk["g_dur"],
        })
    return pairs, concordant, g_violation_osrm_ok, invalid_skip


def times_vector(class_times, class_idx):
    """Return array of per-class seconds aligned to class_idx ordering."""
    v = np.zeros(len(class_idx))
    for cls, t in class_times.items():
        if cls in class_idx:
            v[class_idx[cls]] = t
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skeletons", default=SKELETONS_FILE)
    ap.add_argument("--lam", type=float, default=0.01,
                    help="L2 regularisation of ln(p) toward 0 (p=1)")
    ap.add_argument("--log-margin", dest="log_margin", type=float, default=0.05,
                    help="log-ratio hinge margin (default 0.05 ≈ 5%% cost gap)")
    ap.add_argument("--out", default=OUT_PREF)
    ap.add_argument("--note", default="")
    ap.add_argument("--dry-run", action="store_true",
                    help="report counts and dry-run fit but don't write files")
    args = ap.parse_args()

    ps.load_empirical_base_speeds(BASE_SPEEDS_FILE)
    with open(BASE_SPEEDS_FILE) as f:
        base_speeds = json.load(f)

    skels = load_skeletons(args.skeletons)
    print(f"Loaded {len(skels)} skeletons")

    pairs, concordant, g_viol_osrm_ok, invalid_skip = build_pairs(skels, base_speeds)
    n_g_viol = len(pairs) + g_viol_osrm_ok + invalid_skip
    print(f"Multi-route ODs: {concordant + n_g_viol} total")
    print(f"  Concordant (r0 fastest by g_dur):         {concordant}")
    print(f"  Google violations (r0 slower by g_dur):   {n_g_viol}")
    print(f"    OSRM also routes rk (true violations):  {len(pairs) + invalid_skip}")
    print(f"      both valid (calibration set):         {len(pairs)}")
    print(f"      skipped (match invalid):              {invalid_skip}")
    print(f"    OSRM already routes r0 (timing ok):     {g_viol_osrm_ok}")

    if not pairs:
        sys.exit("No usable violation pairs — nothing to fit.")

    margins = sorted(p["margin_s"] for p in pairs)
    print(f"Violation margin (s): "
          f"median={margins[len(margins)//2]:.0f}  "
          f"mean={sum(margins)/len(margins):.0f}  "
          f"max={margins[-1]:.0f}")
    print()

    if args.dry_run:
        print("[dry-run] Stopping before fit.")
        return

    # Build numpy arrays for optimisation.
    # Use a scale-invariant log-ratio loss:
    #   log(cost(r0) / cost(rk)) = log(p @ t0) - log(p @ tk)
    # Scaling all p_c by k cancels: log(k * p @ t0) - log(k * p @ tk) = unchanged.
    # The optimizer must find *relative* class preferences.
    class_idx = {c: i for i, c in enumerate(MAIN_CLASSES)}
    nC = len(MAIN_CLASSES)

    T0 = np.array([times_vector(p["r0_times"], class_idx) for p in pairs])  # (N, nC)
    Tk = np.array([times_vector(p["rk_times"], class_idx) for p in pairs])  # (N, nC)

    lo = math.log(PREF_LO)
    hi = math.log(PREF_HI)
    theta0 = np.zeros(nC)

    def loss(theta):
        p_vec = np.exp(theta)
        c0 = T0 @ p_vec   # (N,) cost of preferred route
        ck = Tk @ p_vec   # (N,) cost of faster alternative
        # log(c0/ck) < 0 means r0 preferred — hinge fires when log-ratio >= -log_margin
        log_ratio = np.log(np.maximum(c0, 1e-9)) - np.log(np.maximum(ck, 1e-9))
        slack = log_ratio + args.log_margin   # want slack < 0
        hinge = np.maximum(0.0, slack)
        data = float((hinge * hinge).sum()) / len(pairs)
        reg = args.lam * float((theta * theta).sum())
        return data + reg

    def violation_count(theta):
        p_vec = np.exp(theta)
        c0 = T0 @ p_vec
        ck = Tk @ p_vec
        return int((c0 >= ck).sum())

    print(f"Fitting {nC} class preference multipliers over {len(pairs)} violation pairs "
          f"(log_margin={args.log_margin}, lam={args.lam})")
    # By construction all pairs have offline_dur(r0) > offline_dur(rk), so
    # at theta0 (p=1) all are wrong — violation_count == len(pairs).
    assert violation_count(theta0) == len(pairs), \
        "Unexpected: some pairs already resolved at p=1 — check build_pairs filter"
    print(f"Violations correctly ranked before fit: 0/{len(pairs)} (by construction)")

    res = minimize(loss, theta0, method="L-BFGS-B",
                   bounds=[(lo, hi)] * nC,
                   options={"maxiter": 2000, "ftol": 1e-12})
    theta = res.x
    p_fitted = np.exp(theta)

    n_resolved = len(pairs) - violation_count(theta)
    print(f"Violations correctly ranked after fit:  {n_resolved}/{len(pairs)}")
    print(f"Loss: {loss(theta0):.4f} -> {res.fun:.4f}")
    print()
    print("Per-class preference multipliers (p<1 = preferred, p>1 = avoided):")
    print("  (split classes: urban=base key, rural=_rural suffix; ≤30mph vs >30mph)")
    for i, cls in enumerate(MAIN_CLASSES):
        p = p_fitted[i]
        marker = " *" if abs(p - 1.0) > 0.05 else ""
        print(f"  {cls:25s}: {p:.4f}{marker}")

    # Build output: fitted main classes + link classes inheriting parent.
    # Link classes resolve to parent class + same rural/urban sub-bucket in Lua,
    # so we don't need separate link entries — the Lua does parent resolution.
    preference = {}
    for i, cls in enumerate(MAIN_CLASSES):
        p = float(round(p_fitted[i], 4))
        if abs(p - 1.0) > 1e-4:
            preference[cls] = p

    meta = {
        "ts": datetime.datetime.now().isoformat(timespec="seconds"),
        "note": args.note,
        "lam": args.lam,
        "log_margin": args.log_margin,
        "n_pairs": len(pairs),
        "n_resolved_before": 0,
        "n_resolved_after": n_resolved,
        "loss_before": float(loss(theta0)),
        "loss_after": float(res.fun),
        "classes": MAIN_CLASSES,
    }
    out = {"preference": preference, "meta": meta}

    print(f"\nWriting {args.out}")
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Appending history {HISTORY}")
    with open(HISTORY, "a") as f:
        f.write(json.dumps({**meta, "factors": {c: float(round(p_fitted[i], 4))
                                                 for i, c in enumerate(MAIN_CLASSES)}}) + "\n")
    print("Done.")


if __name__ == "__main__":
    main()
