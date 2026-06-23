"""
Generate a road-class-biased OSRM car profile (car_roaaads.lua).

Pulls the default car.lua from the OSRM Docker image, injects a speed-biasing
block that matches HIGHWAY_COST_FACTOR in simulation/routing_config.py, and
writes the result to the OSRM data directory.

NOTE: this is the legacy, class-only biasing profile (one factor per highway
class, conflating impedance and route preference). The Google-calibrated
successor lives in simulation/compile_profile.py (per-(class x speed-band)
factors + tuned turn costs). Shared Lua/Docker plumbing is in osrm_lua.py.

Run from the repo root:
  python3 simulation/build_osrm_profile.py

Then re-preprocess OSRM (commands printed at end of script).
"""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import osrm_lua
from routing_config import HIGHWAY_COST_FACTOR

OSRM_DATA_DIR = os.path.join(os.path.dirname(__file__), "../../osrm")
OUTPUT_LUA = os.path.join(OSRM_DATA_DIR, "car_roaaads.lua")
PBF_NAME = "ireland-and-northern-ireland-latest.osm.pbf"


def _lua_table(factors):
    return "{" + ", ".join(f"{k}={v}" for k, v in factors.items()) + "}"


INJECTION = (
"\n"
"  -- === Road-class biasing (simulation/routing_config.py HIGHWAY_COST_FACTOR) ===\n"
"  -- Applied after maxspeed capping; dividing forward_speed by a factor < 1\n"
"  -- (trunk, primary) shortens the reported duration and lowers routing cost,\n"
"  -- matching the internal Dijkstra biasing in simulation/build_paths.py.\n"
"  do\n"
f"    local _pref = ({_lua_table(HIGHWAY_COST_FACTOR)})[way:get_value_by_key(\"highway\") or \"\"] or 1.0\n"
"    if _pref ~= 1.0 then\n"
"      if (result.forward_speed  or 0) > 0 then\n"
"        result.forward_speed  = result.forward_speed  / _pref\n"
"        result.forward_rate   = result.forward_speed  / 3.6\n"
"      end\n"
"      if (result.backward_speed or 0) > 0 then\n"
"        result.backward_speed = result.backward_speed / _pref\n"
"        result.backward_rate  = result.backward_speed / 3.6\n"
"      end\n"
"    end\n"
"  end\n"
"  -- ============================================================================\n"
)


def main():
    print(f"Pulling car.lua from Docker image {osrm_lua.OSRM_IMAGE} …")
    base_lua, found_path = osrm_lua.pull_base_lua()
    print(f"  Found at {found_path}  ({len(base_lua)} bytes)")

    patched, strategy = osrm_lua.inject(base_lua, INJECTION)
    print(f"  Injection strategy: {strategy}")

    os.makedirs(os.path.abspath(OSRM_DATA_DIR), exist_ok=True)
    with open(OUTPUT_LUA, "w") as f:
        f.write(patched)
    print(f"\nWrote {os.path.abspath(OUTPUT_LUA)}")

    ok, lib_src = osrm_lua.copy_lib(found_path, OSRM_DATA_DIR)
    print("  lib/ copied." if ok else f"  WARNING: could not copy {lib_src}/ — copy manually.")

    abs_data = os.path.abspath(OSRM_DATA_DIR)
    osrm_base = PBF_NAME.replace(".osm.pbf", "")
    lua_name = os.path.basename(OUTPUT_LUA)
    print(f"""
If the injected block looks correct, re-preprocess OSRM:

  cd {abs_data}
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-extract   -p /data/{lua_name} /data/{PBF_NAME}
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-partition  /data/{osrm_base}.osrm
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-customize  /data/{osrm_base}.osrm
  docker run -t -i -p 5000:5000 -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-routed --algorithm mld /data/{osrm_base}.osrm
""")


if __name__ == "__main__":
    main()
