"""
Single source of truth for the study-area geographic knobs.

These define where the model is centred and how the census hierarchy is carved
into core / SDZ / DEA zones. build_census_zones.py consumes the radii (the only
script that resolves core membership); build_network.py and demographics_config.py
import CENTRE from here so it is never defined in more than one place.

Edit values here only. CORE_RADIUS / SDZ_ZONE_RADIUS take effect after re-running
build_census_zones.py (which writes the resulting core_polygon into
data/census_zones.json); CENTRE feeds build_census_zones.py and build_network.py.
"""

CENTRE          = (54.5933779, -5.6960935)   # (lat, lon) of Newtownards town centre
CORE_RADIUS     = 3000     # metres — SDZs intersecting this circle → their DZs become core
SDZ_ZONE_RADIUS = 10000    # metres — DEAs intersecting this circle → broken into SDZs
