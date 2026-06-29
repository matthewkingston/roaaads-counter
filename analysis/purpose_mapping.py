"""Shared NTS trip-purpose → gravity-component mapping (single source of truth).

Imported by BOTH analysis/derive_generation_rates.py (generation magnitudes from
NTS0409a) and analysis/derive_component_profiles.py (temporal shapes from
NTS0502a/0504b), so the purpose→component assignment cannot drift between the two.

JUDGMENT ALLOCATIONS — candidate error sources.  The organising principle is the
attractor each component offers: workplace(jobs)→commute, retail_spaces(PARKING =
all commercial/venue)→retail, school→school, population(HOMES)→res.  NTS cannot
sub-split these purposes by car-driver mode, so each below is a modelling
*decision*, not a data lookup:
  - Business          → retail   (commercial premises = parking, not the jobs count)
  - Personal business → retail   (services/banks/medical = parking)
  - Leisure           → split LEISURE_RETAIL_FRAC retail / rest res
                                 (venue leisure vs visit-friends-at-home)
  - Education/escort  → school   (the car school-run is the adult escort)

Canonical purpose keys (each derive script resolves these to its own NTS table's
columns): commuting, business, education_escort, shopping, other_escort,
personal_business, leisure, other.  "Just walk" is non-vehicle and excluded.
"""

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
