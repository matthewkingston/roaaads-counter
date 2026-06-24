"""
Single source of truth for the Google-calibrated OSRM time profile.

A profile is a grid of multiplicative speed **factors**, one per
(highway_class x speed_band) bucket, plus a handful of global turn/junction
penalty parameters. In OSRM (Lua) a bucket factor divides the segment speed;
in the offline benchmark (analysis/skeleton_model.py) it multiplies the segment
travel time. The two are exactly equivalent:

    time = length / (base_speed / factor) = factor * length / base_speed

so the offline model never re-implements OSRM's speed cascade — the factor is
the only tunable term on the edge side. `factor = 1.0` means "stock-OSRM base
speed for that bucket". OSRM is currently ~21-26% too fast on the approach
corridors, so fitted factors are expected to land > 1 (factor > 1 => slower).

This module is **pure stdlib** so it can be imported by both the stdlib
analysis tooling and the simulation-side Lua generators. It replaces
simulation/routing_config.py (HIGHWAY_COST_FACTOR) for the calibration work.

Bucketisation (`norm_class`, `parse_band`, `bucket_id`/`bucket_of`) is kept
deliberately simple so the *identical* logic can be emitted as Lua by
simulation/osrm_lua.py — the probe profile and the compiled profile both key
on the same (class, band) pair this module defines.
"""

import json
import os
import re

# ── Road-class axis (full DRIVE_HIGHWAYS classification + fallbacks) ──────────
# Order matters: the index into this list is half of the integer bucket id the
# probe profile encodes. Append-only — never reorder once a skeleton cache is
# built against it. `other` catches any drivable highway tag not listed.
CLASSES = [
    "motorway", "motorway_link", "trunk", "trunk_link",
    "primary", "primary_link", "secondary", "secondary_link",
    "tertiary", "tertiary_link", "unclassified",
    "residential", "living_street", "service", "road", "other",
]

# Stock OSRM car.lua highway base speeds (km/h) — the speed OSRM assigns to an
# *untagged* way of each class. Used as base_speed_b for the "untagged" band.
# `road`/`other` fall back to the unclassified speed.
STOCK_SPEED_KMH = {
    "motorway": 90, "motorway_link": 45,
    "trunk": 85, "trunk_link": 40,
    "primary": 65, "primary_link": 30,
    "secondary": 55, "secondary_link": 25,
    "tertiary": 40, "tertiary_link": 20,
    "unclassified": 25, "residential": 25,
    "living_street": 10, "service": 15,
    "road": 25, "other": 25,
}

# ── Speed-band axis ──────────────────────────────────────────────────────────
# NI maxspeed tags are quantised mph values; everything untagged falls back to
# the class base speed. Order matters (index = second half of the bucket id).
MPH_BANDS = [20, 30, 40, 50, 60, 70]
BANDS = ["untagged"] + [str(v) for v in MPH_BANDS] + ["other"]

# ── maxspeed resolution (mirrors OSRM car.lua exactly) ────────────────────────
# OSRM's WayHandlers.maxspeed consults these keys, first present wins (see
# lib/way_handlers.lua). We resolve the same way so a segment's band matches the
# speed OSRM actually assigned it.
MAXSPEED_KEYS = ("maxspeed:advisory", "maxspeed", "source:maxspeed", "maxspeed:type")

# Symbolic maxspeed values OSRM resolves via profile.maxspeed_table (km/h). Only
# the entries that occur in GB/IE/NI data are listed; the country-prefixed forms
# (gb:/uk:) plus the bare highway_type defaults from maxspeed_table_default.
# Source: /home/matthew/Documents/CodingFun/osrm/car_roaaads.lua.
_SYMBOLIC_MAXSPEED_KMH = {
    "gb:nsl_single": (60 * 1609) / 1000, "uk:nsl_single": (60 * 1609) / 1000,
    "gb:nsl_dual":   (70 * 1609) / 1000, "uk:nsl_dual":   (70 * 1609) / 1000,
    "gb:motorway":   (70 * 1609) / 1000, "uk:motorway":   (70 * 1609) / 1000,
    "none": 140,
}
# maxspeed_table_default — used by OSRM when the value is "XX:highway_type" and
# the full string isn't in the table (e.g. "de:urban").
_MAXSPEED_DEFAULT_KMH = {"urban": 50, "rural": 90, "trunk": 110, "motorway": 130}


def osrm_maxspeed_kmh(value):
    """Resolve a single maxspeed tag string to km/h, mirroring OSRM's
    WayHandlers.parse_maxspeed + Measure.parse_value_speed. Returns 0.0 when OSRM
    would fall back to the class base speed (no usable limit)."""
    if not value:
        return 0.0
    s = str(value).strip()
    m = re.match(r"(\d+)", s)            # OSRM: tonumber(source:match("%d*"))
    if m:
        n = float(m.group(1))
        sl = s.lower()
        if "mph" in sl or "mp/h" in sl:
            n *= MPH_KMH
        return n
    sl = s.lower()
    if sl in _SYMBOLIC_MAXSPEED_KMH:
        return _SYMBOLIC_MAXSPEED_KMH[sl]
    mm = re.match(r"[a-z][a-z]:([a-z]+)", sl)   # OSRM: string.match(source,"%a%a:(%a+)")
    if mm:
        return _MAXSPEED_DEFAULT_KMH.get(mm.group(1), 0.0)
    return 0.0

N_CLASSES = len(CLASSES)
N_BANDS = len(BANDS)
N_BUCKETS = N_CLASSES * N_BANDS

# Probe profile encodes the bucket id as a speed (km/h); offset keeps every
# value strictly positive and plausible so map-matching is not distorted.
PROBE_SPEED_OFFSET = 10

MPH_KMH = 1.609344
SPEED_REDUCTION = 0.8        # OSRM car.lua: tagged maxspeed is used at 0.8x
BAND_SNAP_TOL_MPH = 3        # snap a parsed limit to a standard band within +-3 mph

# ── Global turn / junction penalty defaults (stock OSRM car.lua values) ───────
TURN_DEFAULTS = {
    "turn_penalty": 7.5,
    "traffic_light_penalty": 2.0,
    "u_turn_penalty": 20.0,
    "turn_bias": 1.075,
}


# ── Bucketisation (must stay Lua-reproducible — see osrm_lua.py) ─────────────

def norm_class(highway):
    """Map a raw OSM highway tag to one of CLASSES ('other' if unknown)."""
    h = (highway or "").strip()
    return h if h in CLASSES and h != "other" else "other"


def _band_from_kmh(kmh):
    """km/h posted limit -> band label (snap to nearest NI mph band, else 'other').
    A non-positive limit means OSRM used the class base speed -> 'untagged'."""
    if not kmh or kmh <= 0:
        return "untagged"
    mph = kmh / MPH_KMH
    for b in MPH_BANDS:
        if abs(mph - b) <= BAND_SNAP_TOL_MPH:
            return str(b)
    return "other"


def parse_band(maxspeed_raw):
    """Map a single OSM maxspeed value to a band label in BANDS.

    Resolves numeric ('30 mph', '50') *and* symbolic ('GB:nsl_single', 'none')
    values via OSRM's own maxspeed logic (osrm_maxspeed_kmh), then snaps to the
    nearest standard NI band. No tag / unresolved -> 'untagged' (OSRM uses the
    class base speed). Off-grid numeric -> 'other'. National-speed-limit roads
    (nsl_single -> 60, nsl_dual/motorway -> 70) no longer collapse into the
    'untagged' band, which is the point of the nsl-aware resolution.
    """
    return _band_from_kmh(osrm_maxspeed_kmh(maxspeed_raw))


def band_from_tags(tags):
    """Resolve the speed band from a full OSM tag dict, honouring OSRM's key
    precedence (MAXSPEED_KEYS, first present wins)."""
    for k in MAXSPEED_KEYS:
        v = tags.get(k)
        if v:
            return parse_band(v)
    return "untagged"


def bucket_of(tags):
    """Full OSM tag dict -> (class_label, band_label).

    `tags` is the way's complete tag dict (as cached by build_edge_index.py);
    only `highway` and the maxspeed key set are read here, but the rest is kept
    in the cache for future model versions."""
    return norm_class(tags.get("highway")), band_from_tags(tags)


def bucket_key(cls, band):
    """Canonical string key for the factor grid / JSON."""
    return f"{cls}|{band}"


def bucket_index(cls, band):
    """(class, band) -> integer bucket id used by the probe encoding."""
    return CLASSES.index(cls) * N_BANDS + BANDS.index(band)


def bucket_from_index(idx):
    """Inverse of bucket_index: integer id -> (class_label, band_label)."""
    ci, bi = divmod(int(idx), N_BANDS)
    return CLASSES[ci], BANDS[bi]


def bucket_from_probe_speed(speed_kmh):
    """Decode a probe /match segment speed (km/h) back to (class, band).

    Returns None if the rounded value is outside the valid bucket range (a
    segment whose way the probe did not re-bucket, e.g. ferries)."""
    idx = int(round(speed_kmh)) - PROBE_SPEED_OFFSET
    if idx < 0 or idx >= N_BUCKETS:
        return None
    return bucket_from_index(idx)


# Empirical base-speed override (km/h per bucket_key), measured from real OSRM
# by build_skeleton_index.py --base-speeds. When loaded, it replaces the
# analytical estimate below for any bucket it covers — capturing OSRM's actual
# realised forward_speed (nsl resolution, surface caps, default_speed, the 0.8
# reduction, …) which the analytical formula only approximates.
_EMPIRICAL_BASE = {}


def load_empirical_base_speeds(path):
    """Load a {bucket_key: km/h} override table; absent/None -> analytical only.
    Returns the number of buckets loaded."""
    global _EMPIRICAL_BASE
    if path and os.path.exists(path):
        with open(path) as f:
            _EMPIRICAL_BASE = {k: float(v) for k, v in json.load(f).items()}
    else:
        _EMPIRICAL_BASE = {}
    return len(_EMPIRICAL_BASE)


def base_speed_for(cls, band):
    """Base speed (km/h) OSRM assigns to a bucket *before* any factor.

    Prefers the measured empirical value when loaded; otherwise the analytical
    estimate: untagged/other -> class base table; a tagged mph band ->
    mph*1.609*0.8, mirroring OSRM's WayHandlers.maxspeed (which overrides the
    class base). The analytical form misses nsl/advisory/source maxspeed and
    surface caps — hence the empirical override."""
    if _EMPIRICAL_BASE:
        v = _EMPIRICAL_BASE.get(bucket_key(cls, band))
        if v:
            return v
    if band in ("untagged", "other"):
        return float(STOCK_SPEED_KMH.get(cls, STOCK_SPEED_KMH["other"]))
    return float(band) * MPH_KMH * SPEED_REDUCTION


# ── Profile spec object ──────────────────────────────────────────────────────

class ProfileSpec:
    """A candidate profile: per-bucket factors + global turn params.

    `factors` is a sparse dict {bucket_key: factor}; any bucket absent from it
    defaults to 1.0. `turn` holds the four global penalty parameters.
    """

    def __init__(self, factors=None, turn=None):
        self.factors = dict(factors or {})
        self.turn = dict(TURN_DEFAULTS)
        if turn:
            self.turn.update(turn)

    @classmethod
    def default(cls):
        """All factors 1.0 (stock-OSRM base speeds), stock turn params."""
        return cls()

    def factor_for(self, cls_label, band_label):
        return self.factors.get(bucket_key(cls_label, band_label), 1.0)

    def to_dict(self):
        return {"factors": self.factors, "turn": self.turn}

    @classmethod
    def from_dict(cls, d):
        return cls(factors=d.get("factors"), turn=d.get("turn"))

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, sort_keys=True)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def __repr__(self):
        return (f"ProfileSpec({len(self.factors)} non-unit factors, "
                f"turn={self.turn})")
