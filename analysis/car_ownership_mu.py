"""Per-area car-ownership mobilisation multiplier μ — the M1×M2 combiner.

Module 3 of the per-area car-ownership work (see .claude/plans/eager-painting-pearl.md).
Combines the two committed inspection artifacts into a per-area, per-component multiplier
on the producer weights:

    μ_c[area] = Σ_band person_share_band[area] (M2 field)  ·  relrate_c[band] (M1 shape)

where `relrate_c[band]` = M1's *relative* band shape (rate / band-weighted mean, so its
own mean is ~1) and `person_share_band[area]` is M2's per-area person distribution over the
0/1/2/3+ car-availability bands. μ is thus a dimensionless per-area intensity: >1 where car
ownership (and hence car-driver trip generation for that purpose) is above the field mean, <1
below. It is **producer-side** and **redistribute-only** — the level is preserved by
normalising μ_c to producer-weighted mean 1 in the model callers (not here).

Components with a μ (the home-end raw-headcount producer legs): res, retail, and the three
school levels. Commute is deliberately absent — its census producer is already car-driver-
mode-filtered, so it already embeds ownership (μ would double-count).

Pure/read-only: reads the two committed artifacts, no writes. Reused by
build_census_zones.py (external-zone aggregation) and build_demographics.py (internal core).
"""

import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
SHAPE_FILE = os.path.join(_REPO, "analysis", "ownership_shape.json")     # M1
FIELD_FILE = os.path.join(_REPO, "analysis", "car_ownership_field.json")  # M2

# The components that carry a μ (home-end raw-headcount producers). Commute excluded.
MU_COMPONENTS = ("res", "retail", "school_primary", "school_postprimary", "school_tertiary")


def load_relrate(shape_file=SHAPE_FILE):
    """{component: {band: relrate}} — the M1 relative band shapes (rate/band-weighted mean)."""
    with open(shape_file) as f:
        comps = json.load(f)["components"]
    return {c: dict(comps[c]["relative"]) for c in MU_COMPONENTS}


def load_field(field_file=FIELD_FILE):
    """{area_code: {"j","persons_total","share":{band:frac}}} — the M2 per-area field."""
    with open(field_file) as f:
        return json.load(f)["areas"]


def mu_of_share(share, relrate_c):
    """μ for one area+component: Σ_band share_band · relrate_c[band].

    `share` is the area's person-share dict (Σ=1); `relrate_c` the component's band shape.
    Bands are taken from `relrate_c`; a band missing from `share` contributes 0."""
    return sum(float(share.get(b, 0.0)) * float(v) for b, v in relrate_c.items())


def mu_for_area(share, relrate):
    """{component: μ} for one area's share dict, over all MU_COMPONENTS."""
    return {c: mu_of_share(share, relrate[c]) for c in MU_COMPONENTS}


def mu_all_areas(field=None, relrate=None):
    """{area_code: {component: μ}} for every area in the M2 field."""
    field = load_field() if field is None else field
    relrate = load_relrate() if relrate is None else relrate
    return {code: mu_for_area(a["share"], relrate) for code, a in field.items()}
