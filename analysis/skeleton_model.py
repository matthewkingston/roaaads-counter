"""
Offline OSRM time model — the fast core of the profile benchmark.

Given a route *skeleton* (simulation/build_skeleton_index.py) and a ProfileSpec
(simulation/profile_spec.py), predict the route's total travel time **without
touching OSRM**, then score the prediction against Google's duration. Because a
bucket factor enters the edge time linearly,

    edge_time = sum_b  factor_b * length_b * 3.6 / base_speed_kmh_b

and the turn time is a cheap closed form, scoring all ~2000 routes is a few
milliseconds — the whole point of the skeleton cache.

The turn model mirrors OSRM car.lua process_turn: a sigmoid penalty applied only
at real junctions (degree > 2) or u-turns, a flat traffic_light_penalty per
matched signal node, and the u_turn_penalty. NI is left-hand driving, so the
turn bias is inverted exactly as the stock profile does.

Pure stdlib. predict_duration() and the car.lua factor block (osrm_lua.
emit_factor_block) apply the *same* per-bucket factor, so the offline score and
real OSRM agree up to base-speed/surface residual — quantified by the verify gate.
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "simulation"))
import profile_spec as ps   # noqa: E402

_UNBUCKETED_BASE_KMH = ps.STOCK_SPEED_KMH["unclassified"]


def load_skeletons(path):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                sys.stderr.write(f"WARNING: skipping malformed line in {path} "
                                 f"(likely a partial write from an interrupted run)\n")
    return out


def _turn_cost(angle, degree, uturn, turn):
    """OSRM-style turn duration (seconds) for one maneuver."""
    if angle is None or not (degree > 2 or uturn):
        return 0.0
    tp = turn["turn_penalty"]
    tb = 1.0 / turn["turn_bias"]          # left-hand driving (NI), as in car.lua
    if angle >= 0:
        pen = tp / (1.0 + math.exp(-((13.0 / tb) * angle / 180.0 - 6.5 * tb)))
    else:
        pen = tp / (1.0 + math.exp(-((13.0 * tb) * (-angle) / 180.0 - 6.5 / tb)))
    if uturn:
        pen += turn["u_turn_penalty"]
    return pen


def predict_components(skel, spec):
    """Return (edge_time_s, turn_time_s) for a skeleton under a spec."""
    edge = 0.0
    for key, length in skel["length_by_bucket"].items():
        cls, band = key.split("|", 1)
        base = ps.base_speed_for(cls, band)
        edge += spec.factor_for(cls, band) * length * 3.6 / base
    unbk = skel.get("unbucketed_m", 0.0)
    if unbk:
        edge += unbk * 3.6 / _UNBUCKETED_BASE_KMH

    turn = 0.0
    for t in skel.get("turns", []):
        turn += _turn_cost(t.get("angle"), t.get("degree", 0),
                           t.get("uturn", False), spec.turn)
    turn += skel.get("n_signals", 0) * spec.turn["traffic_light_penalty"]
    return edge, turn


def predict_duration(skel, spec):
    edge, turn = predict_components(skel, spec)
    return edge + turn


def _median(xs):
    if not xs:
        return None
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def evaluate(skeletons, spec, valid_only=True):
    """Score a spec against the skeleton cache.

    Per-route residual is the squared log time-ratio e = (ln pred - ln g)^2:
    scale-free (a long external leg and a short in-town leg with equal
    proportional error count equally) and symmetric (25% fast == 25% slow).
    Aggregate is the equal-weighted mean over valid routes. Per-cell stats are
    diagnostics, not weights.
    """
    rows = []
    for s in skeletons:
        if valid_only and not s.get("valid"):
            continue
        g = s.get("g_dur", 0.0)
        if g <= 0:
            continue
        edge, turn = predict_components(s, spec)
        p = edge + turn
        if p <= 0:
            continue
        lr = math.log(p) - math.log(g)
        rows.append({"leg_type": s["leg_type"], "len_band": s["len_band"],
                     "g": g, "p": p, "ratio": p / g, "e": lr * lr,
                     "edge": edge, "turn": turn})

    n = len(rows)
    loss = sum(r["e"] for r in rows) / n if n else float("nan")

    # Per (leg_type x len_band) breakdown.
    cells = {}
    for r in rows:
        cells.setdefault((r["leg_type"], r["len_band"]), []).append(r)
    breakdown = {}
    for key, rs in sorted(cells.items()):
        ratios = [x["ratio"] for x in rs]
        breakdown[key] = {
            "n": len(rs),
            "median_ratio": _median(ratios),
            "mean_abs_lr": sum(abs(math.log(x)) for x in ratios) / len(ratios),
            "turn_frac": sum(x["turn"] for x in rs) / sum(x["edge"] + x["turn"] for x in rs),
        }

    # Per-leg_type rollup.
    by_leg = {}
    for r in rows:
        by_leg.setdefault(r["leg_type"], []).append(r)
    leg_stats = {}
    for leg, rs in sorted(by_leg.items()):
        ratios = [x["ratio"] for x in rs]
        leg_stats[leg] = {
            "n": len(rs),
            "median_ratio": _median(ratios),
            "loss": sum(x["e"] for x in rs) / len(rs),
            "turn_frac": sum(x["turn"] for x in rs) / sum(x["edge"] + x["turn"] for x in rs),
        }

    return {"n_valid": n, "loss": loss,
            "median_ratio": _median([r["ratio"] for r in rows]),
            "leg_stats": leg_stats, "breakdown": breakdown, "rows": rows}


def bucket_coverage(skeletons, valid_only=True):
    """Total matched metres per bucket across the cache (factor identifiability)."""
    cov = {}
    unbk = 0.0
    for s in skeletons:
        if valid_only and not s.get("valid"):
            continue
        for key, length in s["length_by_bucket"].items():
            cov[key] = cov.get(key, 0.0) + length
        unbk += s.get("unbucketed_m", 0.0)
    return cov, unbk


def legacy_spec_from_highway_cost_factor():
    """Build the ProfileSpec equivalent to the deployed class-only profile:
    factor = HIGHWAY_COST_FACTOR[class] for every band. Used only as a
    one-off faithfulness reference against results.jsonl (see eval_profile.py).
    """
    from routing_config import HIGHWAY_COST_FACTOR
    factors = {}
    for cls in ps.CLASSES:
        hc = HIGHWAY_COST_FACTOR.get(cls)
        if hc is None or abs(hc - 1.0) < 1e-9:
            continue
        for band in ps.BANDS:
            factors[ps.bucket_key(cls, band)] = hc
    return ps.ProfileSpec(factors=factors)
