"""
OSRM-equivalent internal edge travel-time model — shared by build_paths.py and
reduce_deadends.py (and stamped into the paths-cache signature by model.py).

Replaces the hand-picked, class-only HIGHWAY_COST_FACTOR for *internal* road-
network path-finding. Each edge's routing cost is the same time the deployed,
Google-calibrated OSRM instance assigns it (simulation/tuned_profile.json,
compiled into osrm/car_roaaads.lua by compile_profile.py):

    edge_time_s = factor(class, band) × length_m × 3.6 / base_speed_kmh(class, band)

— identical to analysis/skeleton_model.py's edge term. All bucketisation,
base-speed and factor logic is reused from simulation/profile_spec.py (the single
source of truth shared with the OSRM Lua emitter), so the offline routing model and
real OSRM key on the same (class × band) buckets.

Base speeds are the empirical, OSRM-measured per-bucket speeds from
data/google_cache/base_speeds.json when present (profile_spec.base_speed_for
prefers them), with profile_spec's analytical estimate as fallback.

NOTE: this models edge *impedance* only. The node-based internal Dijkstra does not
apply OSRM's turn/junction/signal penalties (it never did) — internal legs remain
turn-free. Route-*preference* biasing (the old trunk-favouring 0.67) is intentionally
gone; internal routes are now chosen on realistic time alone.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import profile_spec as ps

# Paths are relative to the repo root (scripts are run from there).
TUNED_PROFILE = "simulation/tuned_profile.json"
BASE_SPEEDS   = "data/google_cache/base_speeds.json"

# Synthetic dead-end super-node edges (reduce_deadends.py) encode an OSRM-measured
# intra-region time via osmnx's maxspeed cascade; build_paths.py keeps that osmnx
# time directly (factor 1.0) instead of bucketing them, so they must be recognised.
DEADEND_HIGHWAY = "deadend_collapsed"


def load_profile(profile_path=TUNED_PROFILE, base_speeds_path=BASE_SPEEDS):
    """Load the tuned ProfileSpec + empirical base speeds used by internal routing.

    Fails loud if the tuned profile is missing — using it is the whole point of
    this module (compile/deploy a profile with tune_profile.py + compile_profile.py
    first). A missing base-speeds file is only a warning (analytical fallback works).
    """
    n_base = ps.load_empirical_base_speeds(base_speeds_path)
    if n_base == 0:
        print(f"WARNING: {base_speeds_path} not found — using analytical base speeds "
              f"(less faithful to the deployed OSRM). Run "
              f"build_skeleton_index.py --base-speeds to generate it.", file=sys.stderr)
    if not os.path.exists(profile_path):
        raise SystemExit(
            f"ERROR: tuned profile {profile_path} not found.\n"
            f"Internal routing now uses the Google-calibrated (class × band) profile.\n"
            f"Run: python3 analysis/tune_profile.py  (writes simulation/tuned_profile.json)")
    return ps.ProfileSpec.load(profile_path)


def edge_time_seconds(tags, length_m, spec):
    """OSRM-equivalent travel time (seconds) for one directed edge.

    `tags` is the edge's attribute dict (osmnx/graphml) — only `highway` and the
    maxspeed key set are read, via profile_spec. Robust to `highway`/`maxspeed`
    stored as a list on consolidated edges.
    """
    hw = tags.get("highway")
    if isinstance(hw, list):
        hw = hw[0] if hw else None
    cls  = ps.norm_class(hw)
    band = ps.band_from_tags(tags)          # honours OSRM maxspeed-key precedence
    base = ps.base_speed_for(cls, band)     # km/h; empirical preferred
    return spec.factor_for(cls, band) * float(length_m) * 3.6 / base
