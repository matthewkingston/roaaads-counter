"""
Generate a road-class-biased OSRM car profile (car_roaaads.lua).

Pulls the default car.lua from the OSRM Docker image, injects a speed-biasing
block that matches HIGHWAY_COST_FACTOR in simulation/routing_config.py, and
writes the result to the OSRM data directory.

Run from the repo root:
  python3 simulation/build_osrm_profile.py

Then re-preprocess OSRM (commands printed at end of script).
"""

import subprocess, sys, os

OSRM_IMAGE     = "osrm/osrm-backend"
OSRM_DATA_DIR  = os.path.join(os.path.dirname(__file__),
                               "../../osrm")          # ../osrm relative to repo root
OUTPUT_LUA     = os.path.join(OSRM_DATA_DIR, "car_roaaads.lua")
PBF_NAME       = "ireland-and-northern-ireland-latest.osm.pbf"

# ── Import factors from routing_config ────────────────────────────────────────

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from routing_config import HIGHWAY_COST_FACTOR

# ── Pull car.lua from Docker ───────────────────────────────────────────────────

LUA_CANDIDATES = [
    "/opt/car.lua",
    "/usr/local/share/osrm/profiles/car.lua",
    "/usr/share/osrm/profiles/car.lua",
]

print(f"Pulling car.lua from Docker image {OSRM_IMAGE} …")
base_lua = None
found_path = None
for candidate in LUA_CANDIDATES:
    result = subprocess.run(
        ["docker", "run", "--rm", OSRM_IMAGE, "cat", candidate],
        capture_output=True, text=True,
    )
    if result.returncode == 0 and result.stdout.strip():
        base_lua = result.stdout
        found_path = candidate
        print(f"  Found at {candidate}  ({len(base_lua)} bytes, "
              f"{base_lua.count(chr(10))} lines)")
        break

if base_lua is None:
    # Last-resort: ask Docker to find it
    r = subprocess.run(
        ["docker", "run", "--rm", OSRM_IMAGE, "find", "/", "-name", "car.lua"],
        capture_output=True, text=True,
    )
    hits = [l.strip() for l in r.stdout.splitlines() if l.strip()]
    if hits:
        r2 = subprocess.run(
            ["docker", "run", "--rm", OSRM_IMAGE, "cat", hits[0]],
            capture_output=True, text=True,
        )
        if r2.returncode == 0 and r2.stdout.strip():
            base_lua = r2.stdout
            found_path = hits[0]
            print(f"  Found (via find) at {found_path}  ({len(base_lua)} bytes)")
    if base_lua is None:
        print("ERROR: Could not locate car.lua in the OSRM Docker image.")
        print("Check that the image is pulled:  docker pull osrm/osrm-backend")
        sys.exit(1)

# ── Build Lua road-preference table from HIGHWAY_COST_FACTOR ──────────────────

def _lua_table(factors: dict) -> str:
    entries = ", ".join(f"{k}={v}" for k, v in factors.items())
    return "{" + entries + "}"

lua_table = _lua_table(HIGHWAY_COST_FACTOR)

INJECTION = (
"\n"
"  -- === Road-class biasing (simulation/routing_config.py HIGHWAY_COST_FACTOR) ===\n"
"  -- Applied after maxspeed capping; dividing forward_speed by a factor < 1\n"
"  -- (trunk, primary) shortens the reported duration and lowers routing cost,\n"
"  -- matching the internal Dijkstra biasing in simulation/build_paths.py.\n"
f"  do\n"
f"    local _pref = ({lua_table})[way:get_value_by_key(\"highway\") or \"\"] or 1.0\n"
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

# ── Find injection point ───────────────────────────────────────────────────────
# We want to inject at the end of the way-processing function, after all speed
# and rate computation is done.  Strategy (tried in order):
#
# 1. After WayHandlers.run() inside process_way: modern profiles delegate all
#    speed/maxspeed/rate work to lib/way_handlers.lua via this call, so
#    injecting after it is always safe.
# 2. After the last result.forward_rate assignment inside process_way/
#    way_function: for profiles that set the rate inline.
# 3. Last resort: just before the final 'return' in the file.

lines = base_lua.splitlines(keepends=True)
inject_idx = None
strategy   = None

# Strategy 1: WayHandlers.run() inside process_way
in_way_fn = False
for i, line in enumerate(lines):
    stripped = line.strip()
    if "function process_way" in stripped or "function way_function" in stripped:
        in_way_fn = True
    if in_way_fn and "WayHandlers.run(" in stripped:
        inject_idx = i
        strategy = "after-WayHandlers.run"
        # don't break — keep scanning so we use the LAST occurrence if needed

if inject_idx is None:
    # Strategy 2: last result.forward_rate inside a function body
    in_way_fn = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if "function process_way" in stripped or "function way_function" in stripped:
            in_way_fn = True
        if not in_way_fn:
            continue
        if stripped.startswith("end") and stripped in ("end", "end\n"):
            in_way_fn = False   # exited function — reset
        if "result.forward_rate" in line and "=" in line:
            inject_idx = i
            strategy = "after-forward_rate-in-function"

if inject_idx is None:
    # Strategy 3: just before the last return
    print("WARNING: could not find a safe injection point inside way function.")
    print("         Falling back to just before the last 'return'.")
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip().startswith("return"):
            inject_idx = i - 1
            strategy = "before-last-return"
            break

if inject_idx is None:
    print("ERROR: No suitable injection point found in car.lua.")
    print("       Please inspect the file manually.")
    sys.exit(1)

print(f"  Injection strategy: {strategy} (after line {inject_idx + 1}: "
      f"{lines[inject_idx].rstrip()!r})")

# ── Inject and write ───────────────────────────────────────────────────────────

patched_lines = (
    lines[:inject_idx + 1]
    + [INJECTION]
    + lines[inject_idx + 1:]
)
patched_lua = "".join(patched_lines)

os.makedirs(os.path.abspath(OSRM_DATA_DIR), exist_ok=True)
with open(OUTPUT_LUA, "w") as f:
    f.write(patched_lua)

out_abs = os.path.abspath(OUTPUT_LUA)
print(f"\nWrote {out_abs}  ({len(patched_lines)} lines)")

# ── Copy lib/ from Docker alongside the profile ───────────────────────────────
# car.lua uses require("lib/...") — the lib/ directory must sit next to the
# profile file when osrm-extract runs.

lib_src = os.path.dirname(found_path) + "/lib"   # e.g. /opt/lib
lib_dst = os.path.join(os.path.abspath(OSRM_DATA_DIR), "lib")
print(f"\nCopying {lib_src} → {lib_dst} …")
copy_result = subprocess.run(
    ["docker", "run", "--rm",
     "-v", f"{os.path.abspath(OSRM_DATA_DIR)}:/data",
     OSRM_IMAGE,
     "sh", "-c", f"cp -r {lib_src} /data/lib"],
    capture_output=True, text=True,
)
if copy_result.returncode == 0:
    print("  lib/ copied.")
else:
    print(f"  WARNING: could not copy lib/ automatically "
          f"(rc={copy_result.returncode}).")
    print(f"  Run manually before osrm-extract:")
    print(f"    docker run --rm -v \"$(pwd):/data\" {OSRM_IMAGE} "
          f"sh -c \"cp -r {lib_src} /data/lib\"")
print()

# ── Show the injected block in context ────────────────────────────────────────

ctx_start = max(0, inject_idx - 2)
ctx_end   = min(len(patched_lines), inject_idx + len(INJECTION.splitlines()) + 6)
print("─" * 70)
print("Injected block in context (inspect before re-processing OSRM):")
print("─" * 70)
for i, line in enumerate(patched_lines[ctx_start:ctx_end], start=ctx_start + 1):
    marker = ">>>" if "Road-class biasing" in line or "_pref" in line else "   "
    print(f"{marker} {i:4d}: {line}", end="" if line.endswith("\n") else "\n")
print("─" * 70)

# ── Print re-preprocessing commands ───────────────────────────────────────────

abs_data = os.path.abspath(OSRM_DATA_DIR)
osrm_base = PBF_NAME.replace(".osm.pbf", "")
lua_name  = os.path.basename(OUTPUT_LUA)

print(f"""
If the injected block looks correct, re-preprocess OSRM:

  cd {abs_data}

  # 1. Extract  (~10–20 min for Ireland+NI)
  docker run --rm -v "$(pwd):/data" {OSRM_IMAGE} \\
    osrm-extract -p /data/{lua_name} /data/{PBF_NAME}

  # 2. Partition  (~2 min)
  docker run --rm -v "$(pwd):/data" {OSRM_IMAGE} \\
    osrm-partition /data/{osrm_base}.osrm

  # 3. Customise  (~2 min)
  docker run --rm -v "$(pwd):/data" {OSRM_IMAGE} \\
    osrm-customize /data/{osrm_base}.osrm

  # 4. Start server
  docker run -t -i -p 5000:5000 -v "$(pwd):/data" {OSRM_IMAGE} \\
    osrm-routed --algorithm mld /data/{osrm_base}.osrm

Then from the repo root:
  python3 simulation/build_external_links.py
  python3 simulation/build_paths.py
  python3 analysis/tune_assignment.py
""")
