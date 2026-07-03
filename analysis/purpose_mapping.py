"""Shared NTS trip-purpose → gravity-component mapping (single source of truth).

`B01_COMPONENT` maps the 23-category microdata purpose (TripPurpose_B01ID) directly to
a gravity component, imported by both analysis/derive_generation_rates.py (magnitudes)
and analysis/derive_component_profiles.py (temporal shapes) so generation and temporal
use the same split.

Organising rule: a trip is **res iff its endpoint is a home**; otherwise it is routed by
the land-use it serves — workplace→commute, commercial/venue→retail, school→school.  So
leisure is split by endpoint (visit-home→res, all other leisure/venue→retail) and escorts
are routed by destination (escort-commuting→commute, escort-shopping/business→retail,
escort-education→school, escort-home→res).  Two allocations remain modelling decisions the
finer codes do not resolve: Business/Other-work and Personal-business → retail (commercial
premises = parking, not the home-workplace commute count).
"""

# ── TripPurpose_B01ID (23-cat) → component ────────────────────────────────────
# Codes (labels for reference): 1 Commuting, 2 Business, 3 Other work, 4 Education,
# 5 Food shop, 6 Non-food shop, 7/8/9 Personal business (medical/eat-drink/other),
# 10 Visit-home, 11 Eat/drink-friends, 12 Other social, 13 Entertainment, 14 Sport,
# 15 Holiday, 16 Day trip, 18 Other non-escort, 19 Escort-commuting,
# 20 Escort-business, 21 Escort-education, 22 Escort-shopping/pb, 23 Escort-home.
# 17 Just-walk never occurs in car-driver mode; sentinels (-8/-10) excluded.
B01_COMPONENT = {
    1: "commute", 19: "commute",
    2: "retail", 3: "retail", 5: "retail", 6: "retail", 7: "retail", 8: "retail",
    9: "retail", 11: "retail", 12: "retail", 13: "retail", 14: "retail", 15: "retail",
    16: "retail", 20: "retail", 22: "retail",
    4: "school", 21: "school",
    10: "res", 18: "res", 23: "res",
}
# B01 codes intentionally NOT assigned a component (dropped, not an error):
B01_EXCLUDE = {17, -8, -10}   # 17 Just-walk (no car-driver trips), NA/DEAD sentinels
