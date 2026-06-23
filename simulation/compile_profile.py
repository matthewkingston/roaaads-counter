"""
Compile a tuned ProfileSpec into a deployable OSRM car profile.

Reads simulation/tuned_profile.json (or --spec), pulls the stock car.lua, applies
the tuned global turn-penalty values in the setup block, injects the per-bucket
speed-factor block (osrm_lua.emit_factor_block — keyed on the *same* integer
(class x band) bucket id as profile_spec.bucket_index, so it cannot drift from
the offline model), writes car_roaaads.lua, and prints the re-preprocess commands.

  python3 simulation/compile_profile.py                       # uses tuned_profile.json
  python3 simulation/compile_profile.py --spec path/to/spec.json
  python3 simulation/compile_profile.py --out /tmp/car_test.lua

After writing, rebuild the DEPLOYED OSRM (:5000) with the printed commands, then
gate it with analysis/verify_profile.py before adopting downstream.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import osrm_lua
import profile_spec as ps

OSRM_DATA_DIR = "/home/matthew/Documents/CodingFun/osrm"
DEFAULT_SPEC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tuned_profile.json")
PBF_NAME = "ireland-and-northern-ireland-latest.osm.pbf"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", default=DEFAULT_SPEC)
    ap.add_argument("--out", default=os.path.join(OSRM_DATA_DIR, "car_roaaads.lua"))
    args = ap.parse_args()

    if not os.path.exists(args.spec):
        sys.exit(f"ERROR: spec {args.spec} not found. Write/tune a ProfileSpec first "
                 f"(see analysis/eval_profile.py / profile_spec.ProfileSpec).")
    spec = ps.ProfileSpec.load(args.spec)
    n_fac = sum(1 for v in spec.factors.values() if abs(v - 1.0) > 1e-9)
    print(f"Spec: {args.spec}  ({n_fac} non-unit factors, turn={spec.turn})")

    print(f"Pulling stock car.lua from {osrm_lua.OSRM_IMAGE} …")
    base_lua, found_path = osrm_lua.pull_base_lua()
    base_lua = osrm_lua.apply_turn_overrides(base_lua, spec.turn)
    patched, strategy = osrm_lua.inject(base_lua, osrm_lua.emit_factor_block(spec))
    print(f"  Injection strategy: {strategy}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(patched)
    print(f"  Wrote {os.path.abspath(args.out)}")

    if os.path.dirname(os.path.abspath(args.out)) == os.path.abspath(OSRM_DATA_DIR):
        ok, lib_src = osrm_lua.copy_lib(found_path, OSRM_DATA_DIR)
        print("  lib/ copied." if ok else f"  WARNING: copy {lib_src}/ manually.")

    base = PBF_NAME.replace(".osm.pbf", "")
    lua_name = os.path.basename(args.out)
    print(f"""
Rebuild the DEPLOYED OSRM (:5000) with the calibrated profile:

  cd {OSRM_DATA_DIR}
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-extract   -p /data/{lua_name} /data/{PBF_NAME}
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-partition  /data/{base}.osrm
  docker run --rm -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-customize  /data/{base}.osrm
  docker run -t -i -p 5000:5000 -v "$(pwd):/data" {osrm_lua.OSRM_IMAGE} osrm-routed --algorithm mld /data/{base}.osrm

Then gate it:
  python3 analysis/verify_profile.py --spec {args.spec}
""")


if __name__ == "__main__":
    main()
