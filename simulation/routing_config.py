"""
Shared routing constants for the probit noise + the legacy OSRM profile generator.

HIGHWAY_COST_FACTOR: multiply the travel time of each road class by this factor
before path-finding.  Values < 1 make that class preferred; > 1 make it avoided.
**No longer used by internal routing** (build_paths.py) or the dead-end reducer
(reduce_deadends.py): those now compute OSRM-equivalent (class × band) edge times
from the Google-calibrated profile via simulation/edge_speed.py. This constant is
retained only for the legacy tooling that still references it — the standalone
OSRM profile generator (build_osrm_profile.py), build_skeleton_index.py's
--base-speeds defactor, and skeleton_model.legacy_spec_from_highway_cost_factor.

PROBIT_CV / PROBIT_LL_SIGMA: probit stochastic route-choice noise (build_paths.py).
Each edge cost is perturbed by exp(eps * w), eps ~ N(0, PROBIT_CV), with a
length-dependent gain w = PROBIT_LL_SIGMA / (PROBIT_LL_SIGMA + PROBIT_CV * cost):
  - short legs (w→1): multiplicative noise, sigma ≈ PROBIT_CV * cost (unchanged);
  - long legs (w→PROBIT_LL_SIGMA/(PROBIT_CV*cost)): noise saturates to an *absolute*
    offset ~ N(0, PROBIT_LL_SIGMA) cost-seconds, so boundary selection on long
    external↔boundary legs is decided by real time differences, not noise.
Crossover at cost ≈ PROBIT_LL_SIGMA / PROBIT_CV.  These are part of the paths-cache
staleness signature (model.paths_cache_signature) — changing them forces a rebuild.
"""

HIGHWAY_COST_FACTOR = {
    "trunk":         0.67,
    "trunk_link":    0.67,
    "primary":       0.67,
    "primary_link":  0.67,
    "secondary":     1.0,
    "tertiary":      1.0,
    "tertiary_link": 1.0,
    "residential":   1.2,
    "unclassified":  1.2,
    "living_street": 1.2,
}

PROBIT_CV       = 0.25     # log-normal noise CV applied to edge costs
PROBIT_LL_SIGMA = 120.0    # long-leg absolute noise sigma (cost-seconds, ≈2 min)
