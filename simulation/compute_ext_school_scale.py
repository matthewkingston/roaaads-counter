"""
Compute the pupils-per-person ratio for the core area.

Reads node_weights.json (internal nodes only — numeric string keys),
sums school demand and population, and prints the ratio to add to
tuner_config.json as "ext_school_per_pop".

This mirrors ext_biz_scale: the value is ad-hoc (core-area OSM schools
applied uniformly to external zones) and is kept explicit in config so
its approximation is visible.

Usage:
    python3 simulation/compute_ext_school_scale.py

Then add the printed value to simulation/tuner_config.json:
    "ext_school_per_pop": <value>
"""

import json, sys, os

WEIGHTS_FILE = "simulation/node_weights.json"

if not os.path.exists(WEIGHTS_FILE):
    print(f"ERROR: {WEIGHTS_FILE} not found. Run build_demographics.py first.")
    sys.exit(1)

with open(WEIGHTS_FILE) as f:
    w = json.load(f)

school = w.get("node_school_demand", {})
pop    = w.get("node_population", {})

# External node IDs are census-area-code strings (e.g. "N21000219").
# Internal node IDs are numeric strings. Filter to internal only.
int_school = {k: v for k, v in school.items() if k[:1].isdigit()}
int_pop    = {k: v for k, v in pop.items()    if k[:1].isdigit()}

total_school = sum(int_school.values())
total_pop    = sum(int_pop.values())

if total_pop == 0:
    print("ERROR: no internal population found in node_weights.json.")
    sys.exit(1)

ratio = total_school / total_pop

print(f"Core area internal nodes: {len(int_pop)}")
print(f"  Total population:    {total_pop:,.0f}")
print(f"  Total school demand: {total_school:,.1f} pupils")
print(f"  Ratio (pupils/person): {ratio:.6f}")
print()
print(f'Add to tuner_config.json:  "ext_school_per_pop": {ratio:.6f}')
