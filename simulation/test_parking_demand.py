"""Unit tests for the parking → retail-spaces estimator (simulation/parking_demand.py).

Run from the simulation/ directory (so zones_config / demographics_config import):
    python3 simulation/test_parking_demand.py
Pure-stdlib, no OSRM/Docker/network. Exits non-zero on any failure.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from parking_demand import parking_spaces
from demographics_config import (
    PARKING_M2_PER_SPACE_OFFSTREET as OFF,
    PARKING_M2_PER_SPACE_ONSTREET as ON,
)

EPS = 1e-6
fails = []


def check(name, got, want, tol=EPS):
    ok = abs(got - want) <= tol
    print(f"  [{'ok' if ok else 'FAIL'}] {name}: got {got:.3f}, want {want:.3f}")
    if not ok:
        fails.append(name)


# 1. Residential micro-parking pad → excluded (access=private)
check("private street_side pad → 0",
      parking_spaces({"parking": "street_side", "access": "private", "capacity": "4"}, 11.0),
      0.0)

# 2. Public on-street with plausible capacity → trust capacity
check("public street_side cap=36 → 36",
      parking_spaces({"parking": "street_side", "access": "yes", "capacity": "36"}, 914.0),
      36.0)

# 3. Multi-storey deck with capacity → trust capacity (gate-exempt)
check("multi-storey cap=400 (tiny footprint) → 400",
      parking_spaces({"parking": "multi-storey", "capacity": "400"}, 1500.0),
      400.0)

# 4. building=parking deck, no capacity, levels=3 → area×levels/30
check("building=parking levels=3 → area×3/30",
      parking_spaces({"building": "parking", "building:levels": "3"}, 1500.0),
      1500.0 * 3 / OFF)

# 5. Untagged surface lot, no capacity → area/30
check("untagged 2000 m² lot → ~67",
      parking_spaces({}, 2000.0),
      2000.0 / OFF)

# 6. parking=surface, no capacity → off-street area/30
check("surface lot 600 m² → 600/30",
      parking_spaces({"parking": "surface"}, 600.0),
      600.0 / OFF)

# 7. On-street, no capacity → area/13 (denser)
check("street_side 130 m², no cap → 130/13",
      parking_spaces({"parking": "street_side"}, 130.0),
      130.0 / ON)

# 8. Gate-failed capacity (implausible) → fall back to area factor, NOT the junk cap
check("cap=120 on 250 m² (2.1 m²/sp, junk) → area/30",
      parking_spaces({"parking": "surface", "capacity": "120"}, 250.0),
      250.0 / OFF)

# 9. access=no and access=permit also excluded
check("access=no → 0", parking_spaces({"access": "no", "capacity": "50"}, 1500.0), 0.0)
check("access=permit → 0", parking_spaces({"access": "permit"}, 1500.0), 0.0)

# 10. Degenerate / robustness
check("zero area → 0", parking_spaces({"parking": "surface"}, 0.0), 0.0)
check("None area → 0", parking_spaces({}, None), 0.0)
check("list-valued tags handled", parking_spaces({"access": ["customers"], "capacity": ["20"]}, 500.0), 20.0)

# 11. Customers access kept (retail) — plausible cap honoured
check("access=customers cap=50 on 1500 m² (30 m²/sp) → 50",
      parking_spaces({"access": "customers", "capacity": "50"}, 1500.0), 50.0)

print()
if fails:
    print(f"FAILED {len(fails)} test(s): {fails}")
    sys.exit(1)
print("All parking_demand tests passed.")
