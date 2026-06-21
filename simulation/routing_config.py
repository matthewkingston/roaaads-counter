"""
Shared routing constants used by both the internal Dijkstra (build_paths.py)
and the OSRM profile generator (build_osrm_profile.py).

HIGHWAY_COST_FACTOR: multiply the travel time of each road class by this factor
before path-finding.  Values < 1 make that class preferred; > 1 make it avoided.
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
