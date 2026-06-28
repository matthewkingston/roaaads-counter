"""
Shared configuration for build_demographics.py (node-weight builder) and
build_map.py (interactive map builder).

Pure constants only — no heavy imports — so both scripts agree on the study-area
centre, file paths, OSM tag handling, and map styling without drifting apart.
"""

# ── Study area ──────────────────────────────────────────────────────────────────
# CENTRE lives in zones_config.py (single source); re-exported here so the
# demographics/map scripts can keep importing it from demographics_config.
from zones_config import CENTRE  # noqa: F401  (re-exported)
OUT_DIR = "simulation"

# ── Projected CRS (single source of truth) ────────────────────────────────────
# Irish Transverse Mercator (ITM) covers all of Ireland with uniform accuracy,
# avoiding the zone-boundary distortion that UTM Zone 30N introduces for
# Republic of Ireland towns in the west (Zone 29N territory).
PROJECTED_CRS = "EPSG:2157"

# OSM download radius and DZ selection are bounded by the core polygon
# (data/census_zones.json), not a fixed circle. NETWORK_MARGIN_M sizes the OSM
# POI/building/parking download circle around the core polygon. (The road
# network's own extent is governed by BOUNDARY_BBOX_MARGIN_M below, not this.)
NETWORK_MARGIN_M = 1000

# ── Road-network source (build_network.py) ───────────────────────────────────────
# build_network.py reads the road graph from the local NI .osm.pbf — the same
# snapshot OSRM is built from (see build_osrm_profile.py / build_external_links.py).
# Sharing one OSM snapshot keeps boundary/internal node IDs consistent with OSRM's
# route node IDs. The graph is read for a bounding box of the core polygon buffered
# by BOUNDARY_BBOX_MARGIN_M, which supersedes the old 1 km Overpass download margin
# (it only needs to reach boundary nodes' external neighbours; 5 km is generous).
PBF_PATH = ("/home/matthew/Documents/CodingFun/osrm/"
            "ireland-and-northern-ireland-latest.osm.pbf")
BOUNDARY_BBOX_MARGIN_M = 5000

# ── File paths ──────────────────────────────────────────────────────────────────
POPULATION_API = (
    "https://ws-data.nisra.gov.uk/public/api.restful/"
    "PxStat.Data.Cube_API.ReadDataset/MYE01T011/CSV/1.0/en/"
)
DZ_BOUNDARY_FILE     = "simulation/dz2021/DZ2021.geojson"
GRAPH_PATH           = "simulation/newtownards_consolidated.graphml"
RAW_GRAPH_PATH       = "simulation/newtownards_network.graphml"
WORKPLACE_DATA_FILE  = "data/census-2021-apwp001.xlsx"
POPULATION_CACHE     = "data/cache_nisra_population.csv"
POI_CACHE            = "data/cache_osm_pois.geojson"
BUILDING_CACHE       = "data/cache_osm_buildings.geojson"
PARKING_CACHE        = "data/cache_osm_parking.geojson"          # legacy per-CENTRE Overpass cache (unused)
PARKING_ISLAND_CACHE = "data/cache_osm_parking_island.geojson"   # island-wide parking (build_parking.py)
CENSUS_ZONES_FILE    = "data/census_zones.json"
TUNER_CONFIG_FILE    = "simulation/tuner_config.json"
NODE_WEIGHTS_FILE    = "simulation/node_weights.json"
DEMOGRAPHICS_GEOJSON = "simulation/newtownards_demographics.geojson"
FLOWS_FILE           = "simulation/newtownards_flows.json"
MAP_HTML             = "simulation/newtownards_map.html"

# ── OSM tag handling ────────────────────────────────────────────────────────────
EXCLUDE_AMENITY = {
    "parking", "parking_space", "parking_entrance",
    "vending_machine", "post_box", "waste_basket",
    "bench", "bicycle_parking", "recycling",
    "shelter", "telephone", "grit_bin",
}

# Per-tag trip-generation weights relative to baseline (café/small shop = 1.0).
# Parking layer already handles large retail anchors; weights here add signal
# from institutional employers and high-turnover stops without parking polygons.
# Schools are intentionally excluded — they have their own node_school_demand layer.
POI_WEIGHTS = {
    # amenity tag → weight
    "hospital":        5.0,
    "cinema":          3.0,
    "theatre":         3.0,
    "fuel":            2.0,
    "fast_food":       1.5,
    "place_of_worship": 0.5,
    "atm":             0.5,
    "toilets":         0.25,
    # shop tag → weight
    "supermarket":     1.5,
    # any office tag → 2.0 (applied inline; not listed here as the value covers all subtypes)
}

# Fallback enrollment estimates (pupils/students) for school demand when OSM
# capacity tag is absent. Used as the control total per school POI — units of
# node_school_demand are then pupils, making W_SCHOOL interpretable as a
# trip-production ratio relative to residential population.
# Sources: DE NI school census 2023/24 (primary ~234, post-primary ~854 per school);
# FE college and university figures are rough estimates.
SCHOOL_ENROLL_FALLBACK = {
    "school":           300,   # generic tag — may be primary or secondary; use midpoint
    "secondary_school": 900,
    "college":          2000,
    "university":       3000,
}

# ── Parking → retail-spaces estimator ────────────────────────────────────────────
# Used by simulation/parking_demand.py (shared by build_demographics.py for internal
# core nodes and build_census_zones.py for external census zones) to turn an OSM
# parking polygon into an estimate of retail parking *spaces* (a count, not the old
# "equivalent persons"). The retail demand component's magnitude is then carried by
# spaces; K_retail (tuner) absorbs the spaces→trips scale, so these constants only
# need to be physically sane and jurisdiction-uniform — NOT branched on NI vs RoI.
#
# Validation (this session, island-wide OSM): destination car parks converge to
# ~30 m²/space in BOTH jurisdictions once mis-tagged residential micro-parking is
# excluded (raw NI 24 / RoI 14.7 was a tagging artefact); public on-street bays run
# ~13.9 m²/space (parallel parking, the carriageway is the aisle).
PARKING_M2_PER_SPACE_OFFSTREET = 30.0    # surface lots incl. aisles/landscaping
PARKING_M2_PER_SPACE_ONSTREET  = 13.0    # parallel on-street bays (no aisle)
PARKING_GATE_LO = 8.0     # implied m²/space below this ⇒ capacity= tag implausible
PARKING_GATE_HI = 80.0    # implied m²/space above this ⇒ capacity= tag implausible
PARKING_EXCLUDE_ACCESS = {"private", "no", "permit"}   # residential/staff, not retail
PARKING_DECK_TYPES     = {"multi-storey", "underground", "rooftop"}  # capacity > footprint
PARKING_ONSTREET_TYPES = {"street_side", "lane"}        # denser m²/space than a lot

# ── Map styling ─────────────────────────────────────────────────────────────────
HIGHWAY_STYLE = {
    "trunk":         {"color": "#f5a623", "weight": 4},
    "trunk_link":    {"color": "#f5a623", "weight": 2},
    "primary":       {"color": "#f5d623", "weight": 3},
    "primary_link":  {"color": "#f5d623", "weight": 2},
    "secondary":     {"color": "#a8d08d", "weight": 2},
    "tertiary":      {"color": "#7bafd4", "weight": 2},
    "tertiary_link": {"color": "#7bafd4", "weight": 1},
    "residential":   {"color": "#cccccc", "weight": 1},
    "unclassified":  {"color": "#bbbbbb", "weight": 1},
    "living_street": {"color": "#dddddd", "weight": 1},
}

ROAD_TYPE_LABELS = {
    "trunk":         "Roads · trunk",
    "trunk_link":    "Roads · trunk (links)",
    "primary":       "Roads · primary",
    "primary_link":  "Roads · primary (links)",
    "secondary":     "Roads · secondary",
    "tertiary":      "Roads · tertiary",
    "tertiary_link": "Roads · tertiary (links)",
    "residential":   "Roads · residential",
    "unclassified":  "Roads · unclassified",
    "living_street": "Roads · living street",
}
