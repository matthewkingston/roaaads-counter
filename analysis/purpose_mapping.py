"""Shared NTS trip-purpose → gravity-component mapping (single source of truth).

Two mappings live here during the microdata migration:

1. B01_COMPONENT — the **generation-side** mapping (used by
   analysis/derive_generation_rates.py): TripPurpose_B01ID (23-cat microdata) →
   component, directly.  This is the current, data-derived scheme.

2. COMPONENT_PURPOSES / CANONICAL_PURPOSES / LEISURE_RETAIL_FRAC — the **legacy
   canonical** scheme, still used by analysis/derive_component_profiles.py (temporal
   shapes from NTS0502a/0504b) until that derivation is migrated to the microdata.
   Do not delete until then.

Organising rule for B01_COMPONENT: a trip is **res iff its endpoint is a home**;
otherwise it is routed by the land-use it serves — workplace→commute,
commercial/venue→retail, school→school.  This dissolves the old LEISURE_RETAIL_FRAC
0.5 judgment (visit-home→res, all other leisure/venue→retail) and routes escorts by
destination (escort-commuting→commute, escort-shopping/business→retail,
escort-education→school, escort-home→res).  Two allocations remain modelling
decisions the finer codes do not resolve: Business/Other-work and Personal-business
→ retail (commercial premises = parking, not the home-workplace commute count).
"""

# ── Generation-side: TripPurpose_B01ID (23-cat) → component ───────────────────
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

# ── Legacy canonical scheme (temporal derivation only) ────────────────────────
LEISURE_RETAIL_FRAC = 0.5

# component -> [(canonical_purpose, weight)].  Leisure appears in both retail and
# res via LEISURE_RETAIL_FRAC.
COMPONENT_PURPOSES = {
    "commute": [("commuting", 1.0)],
    "retail":  [("shopping", 1.0), ("business", 1.0), ("personal_business", 1.0),
                ("leisure", LEISURE_RETAIL_FRAC)],
    "school":  [("education_escort", 1.0)],
    "res":     [("other_escort", 1.0), ("other", 1.0),
                ("leisure", 1.0 - LEISURE_RETAIL_FRAC)],
}

COMPONENTS = ("res", "commute", "retail", "school")

# Every canonical purpose referenced above (deduped, stable order).
CANONICAL_PURPOSES = ("commuting", "business", "education_escort", "shopping",
                      "other_escort", "personal_business", "leisure", "other")
