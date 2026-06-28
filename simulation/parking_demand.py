"""Estimate retail parking *spaces* from an OSM parking polygon.

Shared by `build_demographics.py` (internal core nodes) and `build_census_zones.py`
(external census zones) so both jurisdictions and both node types use ONE estimator.
Replaces the old `area/25` (public) / `area/50` (private) "equivalent persons" hack.

Output is an estimated number of parking spaces (a count). The retail demand
component's magnitude is carried by spaces; the tuner's `K_retail` absorbs the global
spaces→trips scale, so the constants (in `demographics_config.py`) only need to be
physically sane and jurisdiction-uniform.

Recipe (see plan_parking_retail / CLAUDE.md):
  1. Exclude `access ∈ {private, no, permit}` (residential/staff parking, not retail).
  2. Decks (`parking ∈ {multi-storey, underground, rooftop}` or `building=parking`):
     trust `capacity=` if present, else `area × levels / 30`. Gate-exempt — a deck's
     capacity legitimately exceeds its footprint estimate.
  3. Else use `capacity=` only if it passes the plausibility gate (implied
     `area/capacity` within [GATE_LO, GATE_HI]); otherwise fall back to an area factor.
  4. Area fallback: `÷13` for on-street (`street_side`/`lane`), `÷30` otherwise
     (surface lots and *untagged* parking default to off-street).
Tiny mis-tagged residential pads self-zero (11 m² ÷ 30 ≈ 0 spaces).
"""

from demographics_config import (
    PARKING_M2_PER_SPACE_OFFSTREET,
    PARKING_M2_PER_SPACE_ONSTREET,
    PARKING_GATE_LO,
    PARKING_GATE_HI,
    PARKING_EXCLUDE_ACCESS,
    PARKING_DECK_TYPES,
    PARKING_ONSTREET_TYPES,
)


def _first(v):
    """OSM tags can arrive as a list (osmnx) or scalar (osmium) — take the first."""
    if isinstance(v, (list, tuple)):
        return v[0] if v else None
    return v


def _to_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def parking_spaces(tags, area_m2):
    """Estimated retail parking spaces for one parking polygon.

    `tags`   : dict-like of OSM tags (values may be scalars or single-element lists).
    `area_m2`: polygon area in m² (projected, e.g. EPSG:2157).
    Returns a float space count; 0.0 for excluded / degenerate features.
    """
    access = _first(tags.get("access"))
    if access in PARKING_EXCLUDE_ACCESS:
        return 0.0
    if area_m2 is None or area_m2 <= 0:
        return 0.0

    ptype    = _first(tags.get("parking"))
    building = _first(tags.get("building"))
    cap      = _to_int(_first(tags.get("capacity")))

    # ── Decks: capacity legitimately exceeds footprint (multi-level) ──────────────
    if ptype in PARKING_DECK_TYPES or building == "parking":
        if cap and cap > 0:
            return float(cap)
        levels = (_to_int(_first(tags.get("building:levels")))
                  or _to_int(_first(tags.get("parking:levels")))
                  or 1)
        return area_m2 * max(levels, 1) / PARKING_M2_PER_SPACE_OFFSTREET

    # ── Surface / on-street / untagged: gated capacity, else area factor ──────────
    if cap and cap > 0:
        implied = area_m2 / cap
        if PARKING_GATE_LO <= implied <= PARKING_GATE_HI:
            return float(cap)
        # gate-failed ⇒ capacity tag is junk; fall through to the area estimate

    divisor = (PARKING_M2_PER_SPACE_ONSTREET if ptype in PARKING_ONSTREET_TYPES
               else PARKING_M2_PER_SPACE_OFFSTREET)
    return area_m2 / divisor
