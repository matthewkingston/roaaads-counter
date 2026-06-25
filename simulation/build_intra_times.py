#!/usr/bin/env python3
"""Sample intra-zonal OSRM travel times for each external census zone.

Each external zone is a single centroid node in the production-constrained gravity
model, so the k=i diagonal of its per-origin denominator D^c_i (its intra-zonal
trips) is missing — which over-allocates the zone's fixed trip budget to the
observed core (see CLAUDE.md "external intra-zonal self-term" and the memory note
project_production_constrained_gravity).

This script measures that missing diagonal *directly*: for each external zone it
draws M uniform random point-pairs inside the zone's census polygon and routes each
on the local OSRM instance, recording the trip durations.  model.constrained_od_flows
then adds  a^c_i · mean_m F_c(t_im)  to each denominator (denominator-only — no link
flow).  Sampling the real geometry + real OSRM times dissolves the characteristic-
distance constant, the speed assumption, zone-shape irregularity and the urban/rural
speed gap, and using the mean kernel over the sample (E[F], not F(mean)) is correct
for the steep kernel tail.

Inputs:
  data/census_zones.json            — external node list (id, level)
  simulation/sdz2021/SDZ2021.geojson — SDZ polygons (keyed SDZ2021_cd)
  simulation/dz2021/DZ2021.geojson   — DZ polygons  (keyed DZ2021_cd)
  simulation/dea2021/DEA2021.geojson — DEA polygons (keyed FinalR_DEA)
  Local OSRM at OSRM_HOST:OSRM_PORT (same car profile build_external_links.py uses).

Output:
  data/external_intra_times.json  — {"<census_code>": [t1..tM seconds], ...} + _meta.

Run after build_census_zones.py, with OSRM up.  Independent of build_paths.py — the
self-term lives in the model layer, not the paths cache, so re-running this does NOT
require a paths-cache rebuild.  Re-run only when external zones change.
"""
import http.client
import json
import sys
import time

import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from shapely.prepared import prep

CENSUS_ZONES_FILE = "data/census_zones.json"
SDZ_FILE          = "simulation/sdz2021/SDZ2021.geojson"
DZ_FILE           = "simulation/dz2021/DZ2021.geojson"
DEA_FILE          = "simulation/dea2021/DEA2021.geojson"
OUTPUT_FILE       = "data/external_intra_times.json"
OSRM_HOST         = "localhost"
OSRM_PORT         = 5000

M_DEFAULT = 30          # point-pairs sampled per zone
SEED      = 20260625    # deterministic

# geojson property holding the zone code, per level
CODE_COL = {"SDZ": "SDZ2021_cd", "DZ": "DZ2021_cd", "DEA": "FinalR_DEA"}

# ── OSRM helper (mirrors build_external_links.py; duration only) ────────────────

_conn = None
_request_count = 0
_t_start = time.time()


def _get_conn():
    global _conn
    if _conn is None:
        _conn = http.client.HTTPConnection(OSRM_HOST, OSRM_PORT, timeout=15)
    return _conn


def osrm_duration(lat1, lon1, lat2, lon2, retries=3):
    """Return total OSRM driving duration (s) from (lat1,lon1)→(lat2,lon2), or None."""
    global _conn, _request_count
    path = (f"/route/v1/driving/{lon1},{lat1};{lon2},{lat2}"
            f"?overview=false&annotations=false")
    for attempt in range(retries):
        try:
            conn = _get_conn()
            conn.request("GET", path)
            r = conn.getresponse()
            body = r.read()
            _request_count += 1
            if _request_count % 500 == 0:
                el = time.time() - _t_start
                print(f"  {_request_count} queries in {el:.0f}s ({_request_count/el:.0f} q/s)")
            data = json.loads(body)
            if data.get("code") != "Ok":
                return None
            return data["routes"][0]["duration"]
        except (http.client.HTTPException, ConnectionError, json.JSONDecodeError, KeyError):
            _conn = None
            if attempt < retries - 1:
                time.sleep(0.1)
            else:
                return None


def _check_osrm():
    print(f"\nChecking OSRM at {OSRM_HOST}:{OSRM_PORT} …")
    if osrm_duration(54.5933, -5.6960, 54.5933, -5.6960) is None:
        print(f"ERROR: Cannot reach OSRM at {OSRM_HOST}:{OSRM_PORT}.")
        print("Start OSRM with:")
        print("  docker run -t -i -p 5000:5000 -v $(pwd):/data osrm/osrm-backend \\")
        print("    osrm-routed --algorithm mld /data/northern-ireland.osrm")
        sys.exit(1)
    print("  OSRM reachable")


# ── Polygon loading + uniform point sampling ───────────────────────────────────

def load_polygons():
    """Return {(level, code): shapely geometry in WGS84} for all three levels."""
    polys = {}
    for level, path in (("SDZ", SDZ_FILE), ("DZ", DZ_FILE), ("DEA", DEA_FILE)):
        gdf = gpd.read_file(path).to_crs("EPSG:4326")
        col = CODE_COL[level]
        if col not in gdf.columns:
            print(f"ERROR: {path} has no '{col}' column (have {list(gdf.columns)})")
            sys.exit(1)
        for code, geom in zip(gdf[col], gdf.geometry):
            polys[(level, str(code))] = geom
        print(f"  {level}: {len(gdf)} polygons from {path}")
    return polys


def sample_points(geom, n, rng):
    """Rejection-sample n uniform points (lon, lat) inside a (Multi)Polygon."""
    minx, miny, maxx, maxy = geom.bounds
    pg = prep(geom)
    pts = []
    # batch to keep rejection cheap even for low acceptance rates
    while len(pts) < n:
        need = n - len(pts)
        xs = rng.uniform(minx, maxx, size=need * 4)
        ys = rng.uniform(miny, maxy, size=need * 4)
        for x, y in zip(xs, ys):
            if pg.contains(Point(x, y)):
                pts.append((x, y))   # (lon, lat)
                if len(pts) == n:
                    break
    return pts


# ── Main ────────────────────────────────────────────────────────────────────────

def main():
    M = M_DEFAULT
    if "--m" in sys.argv:
        M = int(sys.argv[sys.argv.index("--m") + 1])

    print("Loading census zones …")
    with open(CENSUS_ZONES_FILE) as f:
        external_nodes = json.load(f)["external_nodes"]
    print(f"  {len(external_nodes)} external nodes")

    print("Loading zone polygons …")
    polys = load_polygons()

    _check_osrm()

    print(f"\nSampling {M} intra-zonal point-pairs per zone "
          f"(~{len(external_nodes) * M:,} OSRM routes) …")
    rng = np.random.default_rng(SEED)

    intra_times = {}
    missing_poly = []
    degenerate = []
    for node in external_nodes:
        zid, level = node["id"], node["level"]
        geom = polys.get((level, str(zid)))
        if geom is None:
            missing_poly.append((level, zid))
            continue
        origins = sample_points(geom, M, rng)
        dests   = sample_points(geom, M, rng)
        times = []
        for (olon, olat), (dlon, dlat) in zip(origins, dests):
            dur = osrm_duration(olat, olon, dlat, dlon)
            if dur is not None and dur > 0:
                times.append(round(float(dur), 2))
        if not times:
            degenerate.append((level, zid))
            continue
        intra_times[zid] = times
        if len(times) < M:
            degenerate.append((level, zid, len(times)))

    # ── Report (loud on any zone we could not fully sample) ──────────────────────
    print(f"\nSampled {len(intra_times)}/{len(external_nodes)} external zones")
    if missing_poly:
        print(f"  WARNING: {len(missing_poly)} zones have NO polygon (no self-term applied):")
        for level, zid in missing_poly:
            print(f"    {level} {zid}")
    if degenerate:
        print(f"  WARNING: {len(degenerate)} zones had < {M} successful routes:")
        for d in degenerate:
            print(f"    {d}")
    allt = np.array([t for ts in intra_times.values() for t in ts])
    if len(allt):
        print(f"  intra-zonal times (s): median {np.median(allt):.0f}  "
              f"p10 {np.percentile(allt,10):.0f}  p90 {np.percentile(allt,90):.0f}")

    out = {
        "_meta": {
            "M": M, "seed": SEED, "n_zones": len(intra_times),
            "sampling": "uniform-in-polygon random pairs",
            "note": "intra-zonal OSRM durations (s); denominator-only self-term. "
                    "Same OSRM profile as build_external_links.py.",
        },
        **intra_times,
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(out, f)
    print(f"\nWrote {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
