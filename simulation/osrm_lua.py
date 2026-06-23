"""
Shared OSRM car.lua plumbing: pull the stock profile from Docker, find the
in-way injection point, and emit the two Lua blocks the calibration tooling
needs:

  * emit_probe_block()      — sets every drivable way's speed to its integer
                              (class x band) bucket id (+offset), so a single
                              OSRM /match reveals the bucket of each matched
                              segment via the `speed` annotation.
  * emit_factor_block(spec) — divides each way's speed by the tuned per-bucket
                              factor (the deployable, calibrated bias).

Both blocks compute the bucket index with the *same* logic as
simulation/profile_spec.py (norm_class / parse_band), generated from its
CLASSES / BANDS / MPH_BANDS constants so the Lua and the Python offline model
can never disagree about which bucket a way falls in.

Used by build_osrm_profile.py (legacy class-only profile), build_skeleton_index.py
(probe), and compile_profile.py (compiled profile).
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import profile_spec as ps

OSRM_IMAGE = "osrm/osrm-backend"

LUA_CANDIDATES = [
    "/opt/car.lua",
    "/usr/local/share/osrm/profiles/car.lua",
    "/usr/share/osrm/profiles/car.lua",
]


# ── Pull stock car.lua + lib/ from the Docker image ──────────────────────────

def pull_base_lua(image=OSRM_IMAGE):
    """Return (lua_text, found_path) for the stock car.lua in the image."""
    for cand in LUA_CANDIDATES:
        r = subprocess.run(["docker", "run", "--rm", image, "cat", cand],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            return r.stdout, cand
    # Last resort: locate it.
    r = subprocess.run(["docker", "run", "--rm", image, "find", "/", "-name", "car.lua"],
                       capture_output=True, text=True)
    hits = [l.strip() for l in r.stdout.splitlines() if l.strip()]
    if hits:
        r2 = subprocess.run(["docker", "run", "--rm", image, "cat", hits[0]],
                            capture_output=True, text=True)
        if r2.returncode == 0 and r2.stdout.strip():
            return r2.stdout, hits[0]
    raise RuntimeError("Could not locate car.lua in the OSRM Docker image "
                       f"({image}); try `docker pull {image}`.")


def copy_lib(found_path, data_dir, image=OSRM_IMAGE):
    """Copy the profile's lib/ (require('lib/...')) next to the output profile."""
    lib_src = os.path.dirname(found_path) + "/lib"
    data_dir = os.path.abspath(data_dir)
    r = subprocess.run(
        ["docker", "run", "--rm", "-v", f"{data_dir}:/data", image,
         "sh", "-c", f"cp -r {lib_src} /data/lib"],
        capture_output=True, text=True)
    return r.returncode == 0, lib_src


# ── Injection point finder (3 fallback strategies) ───────────────────────────

def find_injection_point(base_lua):
    """Return (lines, inject_idx, strategy). Inject *after* lines[inject_idx]."""
    lines = base_lua.splitlines(keepends=True)
    inject_idx = strategy = None

    # 1. after WayHandlers.run() inside the way function (use the last hit)
    in_fn = False
    for i, line in enumerate(lines):
        s = line.strip()
        if "function process_way" in s or "function way_function" in s:
            in_fn = True
        if in_fn and "WayHandlers.run(" in s:
            inject_idx, strategy = i, "after-WayHandlers.run"

    # 2. last result.forward_rate assignment inside the way function
    if inject_idx is None:
        in_fn = False
        for i, line in enumerate(lines):
            s = line.strip()
            if "function process_way" in s or "function way_function" in s:
                in_fn = True
            if not in_fn:
                continue
            if s in ("end", "end\n"):
                in_fn = False
            if "result.forward_rate" in line and "=" in line:
                inject_idx, strategy = i, "after-forward_rate-in-function"

    # 3. just before the final return
    if inject_idx is None:
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].strip().startswith("return"):
                inject_idx, strategy = i - 1, "before-last-return"
                break

    if inject_idx is None:
        raise RuntimeError("No suitable injection point found in car.lua.")
    return lines, inject_idx, strategy


def inject(base_lua, block):
    """Insert `block` at the resolved injection point. Returns (text, strategy)."""
    lines, idx, strategy = find_injection_point(base_lua)
    patched = lines[:idx + 1] + [block] + lines[idx + 1:]
    return "".join(patched), strategy


# ── Lua snippet generators (mirror profile_spec bucketisation) ───────────────

def _class_index_table():
    """Lua table {highway_tag = class_index} for known classes (excl. 'other')."""
    pairs = ", ".join(f"{c}={ps.CLASSES.index(c)}"
                      for c in ps.CLASSES if c != "other")
    return "{" + pairs + "}"


def _bucket_index_lua(indent):
    """Lua that sets local `_bucket` to the integer (class x band) bucket id.

    Mirrors profile_spec.norm_class + parse_band exactly. Assumes `way` is in
    scope. Emits indices, not labels, so callers can key on integers."""
    pad = " " * indent
    other_ci = ps.CLASSES.index("other")
    nb = ps.N_BANDS
    tol = ps.BAND_SNAP_TOL_MPH
    snaps = []
    for k, b in enumerate(ps.MPH_BANDS):          # band index 1..len(MPH_BANDS)
        kw = "if" if k == 0 else "elseif"
        snaps.append(f"{pad}      {kw} math.abs(_mph-{b})<={tol} then _bi={k + 1}")
    snaps.append(f"{pad}      else _bi={ps.BANDS.index('other')} end")
    snap_block = "\n".join(snaps)
    return (
        f"{pad}local _hw = way:get_value_by_key(\"highway\") or \"\"\n"
        f"{pad}local _ci = ({_class_index_table()})[_hw] or {other_ci}\n"
        f"{pad}local _bi = 0\n"
        f"{pad}local _ms = way:get_value_by_key(\"maxspeed\")\n"
        f"{pad}if _ms then\n"
        f"{pad}  local _msl = string.lower(_ms)\n"
        f"{pad}  local _n = tonumber(string.match(_msl, \"%d+%.?%d*\"))\n"
        f"{pad}  if _n then\n"
        f"{pad}    local _mph = _n\n"
        f"{pad}    if not string.find(_msl, \"mph\", 1, true) then _mph = _n / {ps.MPH_KMH} end\n"
        f"{snap_block}\n"
        f"{pad}  end\n"
        f"{pad}end\n"
        f"{pad}local _bucket = _ci * {nb} + _bi\n"
    )


def emit_probe_block():
    """Lua block: set every drivable way's speed to PROBE_SPEED_OFFSET+bucket id."""
    off = ps.PROBE_SPEED_OFFSET
    return (
        "\n"
        "  -- === Bucket-id PROBE (simulation/build_skeleton_index.py) ===\n"
        "  -- Overwrites the speed of each drivable way with its integer\n"
        "  -- (class x speed-band) bucket id (+offset). A single OSRM /match then\n"
        "  -- recovers the bucket of every matched segment from annotation.speed.\n"
        "  do\n"
        + _bucket_index_lua(4) +
        "    local _v = " + str(off) + " + _bucket\n"
        "    if (result.forward_speed  or 0) > 0 then\n"
        "      result.forward_speed = _v; result.forward_rate = _v / 3.6\n"
        "    end\n"
        "    if (result.backward_speed or 0) > 0 then\n"
        "      result.backward_speed = _v; result.backward_rate = _v / 3.6\n"
        "    end\n"
        "  end\n"
        "  -- ============================================================\n"
    )


def _factor_table_lua(spec):
    """Lua table {[bucket_id]=factor} for every non-unit factor in the spec."""
    items = []
    for key, fac in sorted(spec.factors.items()):
        if abs(float(fac) - 1.0) < 1e-9:
            continue
        try:
            cls, band = key.split("|", 1)
            bid = ps.bucket_index(cls, band)
        except (ValueError, KeyError):
            continue
        items.append(f"[{bid}]={fac}")
    return "{" + ", ".join(items) + "}"


def emit_factor_block(spec):
    """Lua block: divide each way's speed by its tuned per-bucket factor."""
    return (
        "\n"
        "  -- === Calibrated per-bucket speed factors (compile_profile.py) ===\n"
        "  -- factor > 1 slows the bucket down (OSRM was too fast); keyed by the\n"
        "  -- same integer (class x band) bucket id as profile_spec.bucket_index.\n"
        "  do\n"
        "    local _FAC = " + _factor_table_lua(spec) + "\n"
        + _bucket_index_lua(4) +
        "    local _f = _FAC[_bucket] or 1.0\n"
        "    if _f ~= 1.0 then\n"
        "      if (result.forward_speed  or 0) > 0 then\n"
        "        result.forward_speed  = result.forward_speed  / _f\n"
        "        result.forward_rate   = result.forward_speed  / 3.6\n"
        "      end\n"
        "      if (result.backward_speed or 0) > 0 then\n"
        "        result.backward_speed = result.backward_speed / _f\n"
        "        result.backward_rate  = result.backward_speed / 3.6\n"
        "      end\n"
        "    end\n"
        "  end\n"
        "  -- ================================================================\n"
    )


def apply_turn_overrides(text, turn):
    """Replace the four global turn-penalty property values in the setup block.

    Only the numeric definitions match `key = <number>`; the `local x =
    profile.x` reads in process_turn have a non-numeric RHS and are untouched.
    """
    import re
    for key in ("turn_penalty", "traffic_light_penalty", "u_turn_penalty", "turn_bias"):
        if key not in turn:
            continue
        val = turn[key]
        text, n = re.subn(rf"(\b{key}\s*=\s*)[\d.]+", rf"\g<1>{val}", text, count=1)
        if n == 0:
            print(f"  WARNING: turn param '{key}' not found in car.lua — not overridden")
    return text
