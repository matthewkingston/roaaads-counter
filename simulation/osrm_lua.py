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


def _lua_table(d):
    """Lua table literal {["k"]=v, ...} from a {str: number} dict."""
    return "{" + ", ".join(f'["{k}"]={v}' for k, v in d.items()) + "}"


def _bucket_index_lua(indent):
    """Lua that sets local `_bucket` to the integer (class x band) bucket id.

    Mirrors profile_spec.norm_class + band_from_tags exactly, including OSRM's
    maxspeed key precedence (maxspeed:advisory > maxspeed > source:maxspeed >
    maxspeed:type) and symbolic/national-speed-limit resolution. Assumes `way` is
    in scope. Emits integer indices, generated from profile_spec constants so the
    Lua bucket and the Python offline bucket can never drift."""
    pad = " " * indent
    other_ci = ps.CLASSES.index("other")
    nb = ps.N_BANDS
    tol = ps.BAND_SNAP_TOL_MPH
    other_bi = ps.BANDS.index("other")
    keys = ", ".join(f'"{k}"' for k in ps.MAXSPEED_KEYS)
    sym = _lua_table({k: round(v, 4) for k, v in ps._SYMBOLIC_MAXSPEED_KMH.items()})
    dft = _lua_table(ps._MAXSPEED_DEFAULT_KMH)
    snaps = []
    for k, b in enumerate(ps.MPH_BANDS):          # band index 1..len(MPH_BANDS)
        kw = "if" if k == 0 else "elseif"
        snaps.append(f"{pad}    {kw} math.abs(_mph-{b})<={tol} then _bi={k + 1}")
    snaps.append(f"{pad}    else _bi={other_bi} end")
    snap_block = "\n".join(snaps)
    return (
        f"{pad}local _hw = way:get_value_by_key(\"highway\") or \"\"\n"
        f"{pad}local _ci = ({_class_index_table()})[_hw] or {other_ci}\n"
        f"{pad}-- resolve maxspeed by OSRM key precedence (first non-empty wins)\n"
        f"{pad}local _ms = nil\n"
        f"{pad}for _, _k in ipairs({{{keys}}}) do\n"
        f"{pad}  if _ms == nil or _ms == \"\" then _ms = way:get_value_by_key(_k) end\n"
        f"{pad}end\n"
        f"{pad}local _kmh = 0\n"
        f"{pad}if _ms and _ms ~= \"\" then\n"
        f"{pad}  local _msl = string.lower(_ms)\n"
        f"{pad}  local _d = string.match(_msl, \"^%s*(%d+)\")\n"
        f"{pad}  if _d then\n"
        f"{pad}    local _n = tonumber(_d)\n"
        f"{pad}    if string.find(_msl,\"mph\",1,true) or string.find(_msl,\"mp/h\",1,true) then\n"
        f"{pad}      _n = _n * {ps.MPH_KMH}\n"
        f"{pad}    end\n"
        f"{pad}    _kmh = _n\n"
        f"{pad}  else\n"
        f"{pad}    local _sym = ({sym})[_msl]\n"
        f"{pad}    if _sym then _kmh = _sym\n"
        f"{pad}    else\n"
        f"{pad}      local _ht = string.match(_msl, \"%a%a:(%a+)\")\n"
        f"{pad}      _kmh = (_ht and ({dft})[_ht]) or 0\n"
        f"{pad}    end\n"
        f"{pad}  end\n"
        f"{pad}end\n"
        f"{pad}local _bi = 0\n"
        f"{pad}if _kmh > 0 then\n"
        f"{pad}  local _mph = _kmh / {ps.MPH_KMH}\n"
        f"{snap_block}\n"
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


def _pref_table_lua(pref_dict):
    """Lua table {key=p} for non-unit preference multipliers.

    Keys are highway class names or 'class_rural' for the rural sub-bucket of
    split classes (trunk/primary/secondary/tertiary).  Non-unit entries only.
    """
    items = [f'["{cls}"]={v}' for cls, v in sorted(pref_dict.items())
             if abs(float(v) - 1.0) > 1e-9]
    return "{" + ", ".join(items) + "}"


def emit_factor_block(spec, pref_dict=None):
    """Lua block: timing factors (class×band) + optional urban/rural preference.

    The two tables are separate within one do…end block (highway tag read once).
    pref_dict maps preference keys → p_c (float); p_c < 1 makes a class
    preferred in routing without changing reported duration.

    Preference keys:
      'class'       — urban (≤30mph effective speed) or non-split class
      'class_rural' — rural (>30mph) for trunk/primary/secondary/tertiary
    For untagged roads the OSRM speed_profile default speed determines urban/rural.
    Link classes are resolved to their parent class for the lookup.

    Timing (forward_speed, forward_rate) and preference (forward_rate only) are
    applied in sequence:
      1. forward_speed /= _f          -- timing correction (duration changes)
         forward_rate  = forward_speed / 3.6
      2. forward_rate /= _p           -- preference bias (duration unchanged)
    """
    pref_block = ""
    if pref_dict:
        # Lua constants for the urban/rural split.
        # OSRM speed_profile defaults (km/h) used when no maxspeed tag (_bi==0 or 7).
        osrm_def = (
            "{trunk=85, primary=65, secondary=55, tertiary=40,"
            " trunk_link=40, primary_link=30, secondary_link=25, tertiary_link=20}"
        )
        # Split classes: get rural sub-bucket key when effective speed > 30 mph (48.3 km/h).
        # Link classes are resolved to their parent before the rural check.
        split_set = "{trunk=true, primary=true, secondary=true, tertiary=true}"
        link_parent = (
            "{motorway_link=\"motorway\", trunk_link=\"trunk\","
            " primary_link=\"primary\", secondary_link=\"secondary\","
            " tertiary_link=\"tertiary\"}"
        )
        # MPH values for band indices 1-6 (20/30/40/50/60/70 mph).
        mph_for_bi = "{[1]=20, [2]=30, [3]=40, [4]=50, [5]=60, [6]=70}"
        pref_block = (
            "    -- Preference bias: forward_rate only (duration unchanged)\n"
            "    -- Split classes (trunk/primary/secondary/tertiary) get separate urban/rural\n"
            "    -- parameters: ≤30mph effective speed → base class key; >30mph → class_rural.\n"
            "    -- For untagged roads (_bi==0) the OSRM speed_profile default is used.\n"
            "    -- Link classes are resolved to their parent for the lookup.\n"
            "    do\n"
            "      local _PREF     = " + _pref_table_lua(pref_dict) + "\n"
            "      local _SPLIT    = " + split_set + "\n"
            "      local _LPAR     = " + link_parent + "\n"
            "      local _OSRMDEF  = " + osrm_def + "\n"
            "      local _MPHBI    = " + mph_for_bi + "\n"
            "      local _hw_base  = _LPAR[_hw] or _hw\n"
            "      local _eff_kmh\n"
            "      if _bi == 0 or _bi == 7 then\n"
            "        _eff_kmh = (_OSRMDEF[_hw] or 25)\n"
            "      else\n"
            "        _eff_kmh = (_MPHBI[_bi] or 25) * 1.60934\n"
            "      end\n"
            "      local _pref_key = (_SPLIT[_hw_base] and _eff_kmh > 48.3)\n"
            "                        and (_hw_base .. \"_rural\") or _hw_base\n"
            "      local _p = _PREF[_pref_key] or _PREF[_hw_base] or 1.0\n"
            "      if _p ~= 1.0 then\n"
            "        if (result.forward_rate  or 0) > 0 then\n"
            "          result.forward_rate  = result.forward_rate  / _p\n"
            "        end\n"
            "        if (result.backward_rate or 0) > 0 then\n"
            "          result.backward_rate = result.backward_rate / _p\n"
            "        end\n"
            "      end\n"
            "    end\n"
        )
    return (
        "\n"
        "  -- === Calibrated speed factors + route preference (compile_profile.py) ===\n"
        "  -- _FAC (class×band): timing corrections — divides forward_speed (affects\n"
        "  --   reported duration and routing weight equally).\n"
        "  -- _PREF (class-only): preference bias — divides forward_rate only\n"
        "  --   (routing weight only; reported duration unchanged).\n"
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
        + pref_block +
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
