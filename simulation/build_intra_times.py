#!/usr/bin/env python3
"""Mass-weighted, per-component intra-zonal self-term sampler.

Each external zone is a single centroid node, so its production-constrained denominator
`D^c_i = Σ_k a^c_k·f_c(d_ik)` loses the `k=i` diagonal (its intra-zonal trips), which
over-allocates the zone's fixed trip budget to the rest of the network (see CLAUDE.md
"external intra-zonal self-term").  This script measures that diagonal as
`D_self^c_i = a^c_i · S^c_i`, where `S^c_i` is the **producer×attractor mass-weighted mean
kernel over intra-zonal trips**:

    S^c_i = E_{o ∝ producer^c, b ∝ attractor^c, both in zone i}[ f_c(d_ob) ]

We store the *geometry* (a weighted histogram of the sampled intra-zonal times) per zone per
component; the model applies the tuned kernel `f_c` at eval time (`model.load_self_terms` /
`constrained_od_flows`).  Sampling ∝ mass, road-snapped, with real POI destinations for
retail/school captures **clustering**: a sparse rural zone whose people and jobs both sit in
the same villages reads short (strong self-suppression), while a genuinely spread zone reads
long — the opposite of the old uniform-in-polygon single-average, which sampled empty fields.
Because the `p_o·a_b·f(d_ob)` interaction is symmetric, one `S^c` serves both legs of a
component (out and return), so the previous leg-asymmetry (different diagonal per leg)
dissolves.

Method (per external zone × per component, mirrors analysis/build_n_of_t.py):
  origins ∝ producer over the zone's member small areas, each a FRESH uniform-in-polygon
  point snapped to a road on demand (live, over a thread pool — no static point cache, whose
  fixed few-points-per-area collapsed small single-area zones onto one OSRM node);
  destinations ∝ attractor — member-area road points (res/commute) or real POIs within the
  zone (parking ∝ spaces for retail, schools ∝ per-level enrolment); one OSRM `/table` per
  batch (independent draws ⇒ the S×D matrix is ∝ p⊗a, histogrammed unweighted); off-road
  endpoints (`/nearest`/`/table` snap > SNAP_TOL_M) discarded.  A zone with genuinely one road
  node collapses to all-zero times ⇒ correctly no self-term.

Reuses: build_n_of_t (osrm_table, load_area_masses, load_poi_layers, _check_osrm, SNAP_TOL_M,
NEAREST_TRIES) + the local sample_points; ingest_ni_census/ingest_roi_census +
build_census_zones' parent-map aggregation for zone→member-area membership.

Inputs:  data/census_zones.json, data/island_opportunity_table.csv, the NI/RoI boundary files,
         and the parking/school POI caches.  Needs OSRM up (localhost:5000).
Output:  data/external_intra_times.json — {"<code>": {"<component>": {"t":[s…], "w":[…]}}} + _meta.

Model-layer only: NOT in the paths-cache signature, so re-running needs no paths rebuild —
re-tune afterwards.  `--s/--d/--batches` set the (fixed-generous) sample budget per zone-component.
"""
import argparse
import http.client
import json
import os
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, "simulation")
sys.path.insert(0, "analysis")

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from shapely.prepared import prep

# OSRM endpoint — kept here because build_n_of_t imports these from this module.
OSRM_HOST = "localhost"
OSRM_PORT = 5000

CENSUS_ZONES_FILE = "data/census_zones.json"
OPP_TABLE         = "data/island_opportunity_table.csv"
OUTPUT_FILE       = "data/external_intra_times.json"

SEED       = 20260706
BIN_STEP_S = 30.0
BIN_CAP_S  = 14400.0                    # 240 min; intra-zonal times sit well below this
S_DEFAULT  = 50                         # origin points per zone-component per /table batch
D_DEFAULT  = 50                         # destination points (S+D ≤ osrm --max-table-size)
BATCHES    = 8                          # /table batches per zone-component (⇒ 20k pairs — from the
                                        # convergence probe: ~1% rel-SE on res/commute/retail and
                                        # ~2.3% on the binding case, sparse-rural school; dense/near
                                        # zones (the SDZs) converge <1% far sooner ⇒ harmless overkill)

# The six model kernel components == build_n_of_t.PURPOSES.  Per component:
#   producer  = opportunity-table column; destination = area column or POI layer/weight-column.
COMPONENTS = ("res", "commute", "retail",
              "school_primary", "school_postprimary", "school_tertiary")
COMP_PRODUCER = {
    "res": "population", "commute": "commute_producers", "retail": "population",
    "school_primary": "school_producers_primary",
    "school_postprimary": "school_producers_postprimary",
    "school_tertiary": "school_producers_tertiary",
}
# destination spec: ("area", <opp col>) or ("poi", <layer>, <weight col or None>[, <area fallback col>]).
# The optional 4th element is an opportunity-table column to fall back to (area-level, road-snapped)
# when the zone has no POIs of that layer — used for retail, whose attractor `retail_spaces` carries a
# workplace-derived fallback for no-parking zones (build_census_zones / build_opportunity_table), so the
# self-term must mirror it. School levels have no fallback: no school of a level ⇒ no intra-zonal school
# trips of that level ⇒ correctly no self-term.
COMP_DEST = {
    "res":     ("area", "population"),
    "commute": ("area", "commute_attractor"),
    "retail":  ("poi", "parking", None, "retail_spaces"),
    "school_primary":     ("poi", "school", "enrol_primary"),
    "school_postprimary": ("poi", "school", "enrol_postprimary"),
    "school_tertiary":    ("poi", "school", "enrol_tertiary"),
}


def sample_points(geom, n, rng):
    """Rejection-sample n uniform points (lon, lat) inside a (Multi)Polygon.
    (Kept in this module — build_n_of_t imports it.)"""
    minx, miny, maxx, maxy = geom.bounds
    pg = prep(geom)
    pts = []
    while len(pts) < n:
        need = n - len(pts)
        xs = rng.uniform(minx, maxx, size=need * 4)
        ys = rng.uniform(miny, maxy, size=need * 4)
        for x, y in zip(xs, ys):
            if pg.contains(Point(x, y)):
                pts.append((x, y))
                if len(pts) == n:
                    break
    return pts


# ── Concurrent b-faithful road-point sampling (live, no static point cache) ─────
# Each intra-zonal origin/destination is a FRESH uniform-in-polygon point retried
# until it snaps < SNAP_TOL_M to a road (so it preserves the drawn area's mass) —
# identical in spirit to build_n_of_t.road_point, but sampled on demand instead of
# from a fixed 3-point-per-area cache (which collapsed to one OSRM node for small
# single-area zones). The /nearest snap calls dominate the runtime and are
# independent, so they run over a thread pool. build_n_of_t's OSRM helpers share one
# HTTPConnection (not thread-safe), so this path uses its own per-thread connections
# and never touches build_n_of_t._conn.
N_SNAP_WORKERS = 16                     # concurrent /nearest workers

_tls = threading.local()


def _nearest_conn():
    c = getattr(_tls, "conn", None)
    if c is None:
        c = _tls.conn = http.client.HTTPConnection(OSRM_HOST, OSRM_PORT, timeout=60)
    return c


def _osrm_nearest_ts(lon, lat):
    """Thread-safe /nearest snap distance (m), or None — one HTTPConnection per thread."""
    for _ in range(3):
        try:
            c = _nearest_conn()
            c.request("GET", f"/nearest/v1/driving/{lon},{lat}?number=1")
            d = json.loads(c.getresponse().read())
            return d["waypoints"][0]["distance"] if d.get("code") == "Ok" else None
        except (http.client.HTTPException, ConnectionError, OSError, ValueError):
            _tls.conn = None
    return None


def sample_road_points(geoms, seeds, snap_tol, tries, pool):
    """One road-proximate point per geom (b-faithful: a fresh uniform-in-polygon point retried
    until it snaps < snap_tol, else the best-snapping candidate). Retries proceed in rounds:
    the GEOS point sampling runs on THIS thread only — shapely/GEOS is NOT thread-safe — while
    each round's independent /nearest snap calls fan out over `pool`. Deterministic: each point
    has its own rng seeded from `seeds` (assigned in caller order), independent of scheduling.
    A point is None only if it never reached OSRM in any round."""
    n = len(geoms)
    rngs = [np.random.default_rng(s) for s in seeds]
    best = [None] * n
    best_d = [np.inf] * n
    pending = list(range(n))
    for _ in range(tries):
        if not pending:
            break
        cand = [sample_points(geoms[i], 1, rngs[i])[0] for i in pending]   # GEOS — single-thread
        dists = list(pool.map(lambda ll: _osrm_nearest_ts(*ll), cand))     # /nearest — concurrent
        still = []
        for i, lonlat, d in zip(pending, cand, dists):
            if d is not None and d < best_d[i]:
                best_d[i], best[i] = d, lonlat
            if d is None or d >= snap_tol:
                still.append(i)
        pending = still
    return best


# ── Zone → member small-area membership (reconstructed like build_census_zones) ──

def build_membership():
    """Return {external_zone_id: [member small-area codes]} for every external node,
    reconstructing build_census_zones' parent-map aggregation (handles the NI SDZ→DEA
    spatial-join quirk and RoI ED/LEA-by-dissolve via the ingest loaders)."""
    from ingest_ni_census import load_ni_census
    from ingest_roi_census import load_roi_census
    dz_gdf, sdz_gdf, _ = load_ni_census()
    sa_gdf, ed_gdf, _ = load_roi_census()
    dz  = pd.concat([dz_gdf,  sa_gdf], ignore_index=True)    # small areas (DZ + SA)
    sdz = pd.concat([sdz_gdf, ed_gdf], ignore_index=True)    # intermediate (SDZ + ED)
    dz_to_sdz  = dz.set_index("area_code")["parent_code"].to_dict()   # small area → intermediate
    sdz_to_dea = sdz.set_index("area_code")["parent_code"].to_dict()  # intermediate → outer

    sa_by_int = defaultdict(list)
    for sa, p in dz_to_sdz.items():
        sa_by_int[p].append(sa)
    int_by_outer = defaultdict(list)
    for s, d in sdz_to_dea.items():
        int_by_outer[d].append(s)

    with open(CENSUS_ZONES_FILE) as f:
        nodes = json.load(f)["external_nodes"]
    membership = {}
    for node in nodes:
        zid, level = node["id"], node["level"]
        if level in ("SDZ", "ED"):                          # intermediate node → its small areas
            members = list(sa_by_int.get(zid, []))
        elif level in ("DEA", "LEA"):                       # outer node → grandchild small areas
            members = [sa for it in int_by_outer.get(zid, []) for sa in sa_by_int.get(it, [])]
        else:                                               # orphan DZ/SA node → itself
            members = [zid]
        membership[zid] = members
    return membership, nodes


# ── POI → external-zone assignment (point-in-zone) ─────────────────────────────

def assign_pois(zone_polys, pois):
    """{zone_id: {'parking': (coords, w), 'school': (coords, {col: w})}} via point-in-zone
    containment of the global POI clouds (build_n_of_t.load_poi_layers)."""
    zids = list(zone_polys)
    zgdf = gpd.GeoDataFrame({"zid": zids}, geometry=[zone_polys[z] for z in zids],
                            crs="EPSG:4326")
    out = {z: {} for z in zids}
    for layer, (coords, w) in pois.items():
        pts = gpd.GeoDataFrame(geometry=[Point(lo, la) for lo, la in coords], crs="EPSG:4326")
        j = gpd.sjoin(pts, zgdf, predicate="within", how="inner")
        for z, grp in j.groupby("zid"):
            idx = grp.index.values
            if layer == "parking":
                out[z]["parking"] = (coords[idx], w[idx])
            else:
                out[z]["school"] = (coords[idx], {c: w[c][idx] for c in w})
    return out


# ── Per zone × component sampling ──────────────────────────────────────────────

def sample_zone_component(comp, members, area_geom, area_mass, zpois,
                          S, D, batches, edges, osrm_table, rng, snap_tol, tries, pool):
    """Weighted histogram of intra-zonal times for one zone × component, or None if the
    zone has no producer/attractor for this component, or all intra-zonal pairs collapse
    onto a single road node (a zone with genuinely one road node ⇒ correctly no self-term)."""
    prod_col = COMP_PRODUCER[comp]
    pw = np.array([area_mass.get(a, {}).get(prod_col, 0.0) for a in members], dtype=float)
    if pw.sum() <= 0:
        return None
    pw = pw / pw.sum()

    def _area_weights(col):
        aw = np.array([area_mass.get(a, {}).get(col, 0.0) for a in members], dtype=float)
        return (aw / aw.sum()) if aw.sum() > 0 else None

    dspec = COMP_DEST[comp]
    dest_mode = dspec[0]                                     # "area" or "poi"
    if dest_mode == "area":
        aw = _area_weights(dspec[1])
        if aw is None:
            return None
    else:                                                   # poi, with an optional area fallback (dspec[3])
        layer, wcol = dspec[1], dspec[2]
        entry = zpois.get(layer)
        pw_raw = None
        if entry is not None:
            pcoords, pw_raw = entry
            pw_raw = pw_raw if wcol is None else pw_raw[wcol]
        if pw_raw is not None and (pw_raw > 0).any():
            m = pw_raw > 0
            poi_coords, poi_w = pcoords[m], pw_raw[m] / pw_raw[m].sum()
        elif len(dspec) > 3:                                # no POIs in zone → area-level fallback column
            aw = _area_weights(dspec[3])                    # retail: workplace-derived retail_spaces
            if aw is None:
                return None
            dest_mode = "area"
        else:
            return None                                     # e.g. no school of this level ⇒ no self-term

    def _geoms(idx):                                    # geometries of the drawn member areas
        return [area_geom[members[a]] for a in idx]

    counts = np.zeros(len(edges) - 1, dtype=np.float64)
    for _ in range(batches):
        oi = rng.choice(len(members), size=S, p=pw)
        src = [p for p in sample_road_points(_geoms(oi), rng.integers(0, 2**63 - 1, size=S),
                                             snap_tol, tries, pool) if p is not None]
        if not src:
            continue
        if dest_mode == "area":
            di = rng.choice(len(members), size=D, p=aw)
            dst = [p for p in sample_road_points(_geoms(di), rng.integers(0, 2**63 - 1, size=D),
                                                 snap_tol, tries, pool) if p is not None]
        else:
            di = rng.choice(len(poi_coords), size=D, p=poi_w)
            dst = [tuple(poi_coords[k]) for k in di]
        if not dst:
            continue
        res = osrm_table(src, dst)
        if res is None:
            continue
        dur, ssnap, dsnap = res
        valid = (np.isfinite(dur) & (dur > 0)          # drop degenerate same-node self-pairs
                 & (ssnap < snap_tol)[:, None]
                 & (dsnap < snap_tol)[None, :])
        counts += np.histogram(dur[valid], bins=edges)[0]

    if counts.sum() <= 0:
        return None
    centers = 0.5 * (edges[:-1] + np.minimum(edges[1:], BIN_CAP_S + BIN_STEP_S))
    nz = counts > 0
    w = counts[nz] / counts.sum()
    return {"t": [round(float(t), 1) for t in centers[nz]],
            "w": [round(float(x), 6) for x in w]}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s", type=int, default=S_DEFAULT, help="origin points per /table batch")
    ap.add_argument("--d", type=int, default=D_DEFAULT, help="destination points per /table batch")
    ap.add_argument("--batches", type=int, default=BATCHES, help="/table batches per zone-component")
    ap.add_argument("--component", choices=list(COMPONENTS),
                    help="sample ONLY this component and MERGE it into the existing "
                         "external_intra_times.json (updates that component per zone, keeps the "
                         "others). Default: all 6, overwrite.")
    args = ap.parse_args()
    comps_to_run = [args.component] if args.component else list(COMPONENTS)

    from build_n_of_t import (osrm_table, _check_osrm, load_area_masses,
                              load_poi_layers, SNAP_TOL_M, NEAREST_TRIES)
    _check_osrm()

    print("Reconstructing zone → member-area membership …")
    membership, nodes = build_membership()
    print(f"  {len(membership)} external zones")

    print("Loading area masses + geometries + POIs …")
    df, geoms = load_area_masses()                          # df ∝ opportunity table, geoms WGS84
    code_to_row = {c: i for i, c in enumerate(df["area_code"].values)}
    area_mass = df.set_index("area_code").to_dict("index")
    pois = load_poi_layers()

    rng = np.random.default_rng(SEED)
    area_geom = {c: geoms[i] for c, i in code_to_row.items()}   # area code → WGS84 polygon

    print("Assigning POIs to zones (point-in-zone) …")
    zone_polys = {}
    for zid, members in membership.items():
        gs = [geoms[code_to_row[c]] for c in members if c in code_to_row]
        if gs:
            zone_polys[zid] = gs[0] if len(gs) == 1 else gpd.GeoSeries(gs).unary_union
    zpois = assign_pois(zone_polys, pois)

    edges = np.append(np.arange(0.0, BIN_CAP_S + BIN_STEP_S, BIN_STEP_S), np.inf)
    print(f"\nSampling {args.s}×{args.d}×{args.batches} intra-zonal pairs per zone × "
          f"{len(comps_to_run)} component(s): {', '.join(comps_to_run)} "
          f"(live road-point sampling, {N_SNAP_WORKERS} snap workers) …")
    t0 = time.time()
    out_zones = {}
    empty_zone, skipped = [], defaultdict(int)
    pool = ThreadPoolExecutor(max_workers=N_SNAP_WORKERS)
    try:
        for n, zid in enumerate(membership):
            members = [c for c in membership[zid] if c in code_to_row]
            if not members:
                empty_zone.append(zid)
                continue
            byc = {}
            for comp in comps_to_run:
                h = sample_zone_component(comp, members, area_geom, area_mass,
                                          zpois.get(zid, {}), args.s, args.d, args.batches,
                                          edges, osrm_table, rng, SNAP_TOL_M, NEAREST_TRIES, pool)
                if h is not None:
                    byc[comp] = h
                else:
                    skipped[comp] += 1
            if byc:
                out_zones[zid] = byc
            if (n + 1) % 25 == 0:
                print(f"  {n+1}/{len(membership)} zones ({time.time()-t0:.0f}s)", flush=True)
    finally:
        pool.shutdown(wait=True)

    print(f"\nSampled {len(out_zones)}/{len(membership)} zones in {time.time()-t0:.0f}s")
    if empty_zone:
        print(f"  {len(empty_zone)} zones had no member areas in the opportunity table (no self-term)")
    for comp, k in skipped.items():
        if k:
            print(f"  component {comp}: {k} zones with no self-term "
                  f"(no producer/attractor, or intra-zonal pairs collapse to a single road node)")

    per_comp = {c: {"s": args.s, "d": args.d, "batches": args.batches} for c in comps_to_run}
    if args.component and os.path.exists(OUTPUT_FILE):
        # MERGE: replace only this component per zone, keep the others from the existing file.
        with open(OUTPUT_FILE) as f:
            prev = json.load(f)
        prev_meta = prev.pop("_meta", {})
        c = args.component
        for zid in list(prev):
            prev[zid].pop(c, None)                          # drop the stale target-component entry
        for zid, byc in out_zones.items():
            prev.setdefault(zid, {}).update(byc)            # add the freshly-sampled one
        zones = {z: d for z, d in prev.items() if d}        # drop zones left with no components
        per_comp = {**prev_meta.get("per_component", {}), **per_comp}  # keep other comps' provenance
        print(f"  MERGED component '{c}' into existing {OUTPUT_FILE}")
    else:
        if args.component:
            print(f"  NOTE: {OUTPUT_FILE} absent — writing '{args.component}' only "
                  f"(other components missing until a full run)")
        zones = out_zones

    components = sorted({k for d in zones.values() for k in d})
    out = {
        "_meta": {
            "seed": SEED, "components": components, "per_component": per_comp,
            "bin_step_s": BIN_STEP_S, "bin_cap_s": BIN_CAP_S, "n_zones": len(zones),
            "sampling": "mass-weighted (producer×attractor) intra-zonal, road-snapped, real POIs",
            "note": "per-zone per-component weighted time histograms; denominator-only self-term. "
                    "Model applies the tuned kernel f_c to the bin centres (see model.load_self_terms).",
        },
        **zones,
    }
    with open(OUTPUT_FILE, "w") as f:
        json.dump(out, f)
    print(f"Wrote {OUTPUT_FILE} ({len(zones)} zones, components: {', '.join(components)})")


if __name__ == "__main__":
    main()
