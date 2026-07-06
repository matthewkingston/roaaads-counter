"""National Ireland opportunity-density-in-cost `n_Ire(t)` — point-cloud Monte-Carlo sampler.

Per purpose, estimates `n(t) = Σ_{i,j} P_i·A_j·δ(c_ij − t)` (P=producer/origin mass, A=attractor/
destination mass, c=OSRM car time on the deployed car_roaaads.lua) by drawing origin points ∝ P and
destination points ∝ A, routing them on the local island OSRM, and histogramming the times.

**This is the NAIVE outer-product sampler (v1) + a one-purpose PILOT** — the agreed first checkpoint.
Because origin/dest draws are independent, the full B×B cross-product of an origin-batch × dest-batch is
already distributed ∝ P_i·A_j, so one OSRM `/table` call yields B² correctly-distributed pairs
(accumulated unweighted). The short-range **head** is the known weak spot; the pilot measures its
raggedness so we can decide whether distance-band stratification (deferred) is needed before scaling up.

Six purposes (outbound leg only, unconstrained N(t) — no 1/D_i; both v1 simplifications):
  res              P=population           A=population            (area→area)
  commute          P=commute_producers    A=commute_attractor     (area→area)
  retail           P=population           A=retail parking POIs ∝ parking_spaces
  school_primary   P=school_producers_primary      A=school POIs ∝ enrol_primary
  school_postprimary  …postprimary        …enrol_postprimary
  school_tertiary  …tertiary              …enrol_tertiary

Destinations for area→area purposes are road-proximate points sampled inside the destination area;
retail/school destinations are the real POIs directly (finer than an area total). Off-road points are
handled by an OSRM `/nearest` pre-filter (resample until snap < tolerance, preserving the area's mass);
POIs are QA'd via the `/table` snap distance. Deterministic (fixed SEED). Caveat: car_roaaads.lua was
calibrated on Newtownards Google data — a national approximation.

Run:
  python3 analysis/build_n_of_t.py --pilot --purpose res --pairs 1000000   # pilot (needs OSRM up)
  python3 analysis/build_n_of_t.py --pairs 10000000                        # all six, scale-up
Needs OSRM up (island extract, localhost:5000). `--batch B` needs `osrm-routed --max-table-size >= 2B`
(default 100 ⇒ B<=50, no restart needed for the pilot).
"""
import argparse
import glob
import hashlib
import http.client
import json
import os
import sys
import time
from datetime import date

import numpy as np
import pyproj
from scipy.spatial import cKDTree

sys.path.insert(0, "simulation")
import geopandas as gpd
from build_intra_times import sample_points, OSRM_HOST, OSRM_PORT
from school_attractor import add_level_enrolments, LEVEL_ENROL_COLS
from parking_demand import parking_spaces
from demographics_config import PROJECTED_CRS, PARKING_ISLAND_CACHE, SCHOOL_ISLAND_CACHE, OSRM_LUA

OPP_TABLE   = "data/island_opportunity_table.csv"
OUTPUT_FILE = "data/national_n_of_t.json"
# OSRM_LUA (single-source, from demographics_config) — used for a profile hash (best-effort)

SEED         = 20260703
SNAP_TOL_M   = 250.0        # reject a trip end snapping further than this from a road
NEAREST_TRIES = 12          # resample attempts to land an area point near a road
BIN_STEP_S   = 30.0         # histogram bin width (seconds) — fine in the head
BIN_CAP_S    = 14400.0      # last explicit edge (240 min); everything beyond → overflow bin
BATCH        = 45           # sources = dests = B; needs OSRM --max-table-size >= 2B

# stratification (--stratified): near-field distance bands + far tail + cached area points.
BANDS_KM    = (0.0, 3.0, 10.0, 25.0)   # near-band edges (km); far band = [last, inf)
K_CACHE     = 3                         # cached road points per area (purpose-independent, reused)
NEAR_BUDGET = 150_000                   # sampled pairs per near band
FAR_BUDGET  = 1_000_000
DESTS_PER_CALL = 49                     # D in the near 1×D /table calls (1+D <= max-table-size)
POINT_CACHE = "data/_area_road_points.json"    # gitignored; resumable

# purpose → (producer column, destination spec). dest "area:<col>" or "poi:<layer>[:<col>]".
PURPOSES = {
    "res":                ("population",                 "area:population"),
    "commute":            ("commute_producers",          "area:commute_attractor"),
    "retail":             ("population",                 "poi:parking"),
    "school_primary":     ("school_producers_primary",    "poi:school:enrol_primary"),
    "school_postprimary": ("school_producers_postprimary", "poi:school:enrol_postprimary"),
    "school_tertiary":    ("school_producers_tertiary",   "poi:school:enrol_tertiary"),
}


# ── OSRM helpers (/table + /nearest) ──────────────────────────────────────────
_conn = None


def _get_conn():
    global _conn
    if _conn is None:
        _conn = http.client.HTTPConnection(OSRM_HOST, OSRM_PORT, timeout=60)
    return _conn


def _osrm_get(path, retries=3):
    global _conn
    for attempt in range(retries):
        try:
            c = _get_conn()
            c.request("GET", path)
            data = json.loads(c.getresponse().read())
            if data.get("code") != "Ok":
                return None
            return data
        except (http.client.HTTPException, ConnectionError, json.JSONDecodeError, OSError):
            _conn = None
            if attempt < retries - 1:
                time.sleep(0.1)
    return None


def osrm_nearest(lon, lat):
    """Snap distance (m) of (lon,lat) to the nearest road, or None."""
    d = _osrm_get(f"/nearest/v1/driving/{lon},{lat}?number=1")
    if d is None:
        return None
    return d["waypoints"][0]["distance"]


def osrm_table(src, dst):
    """src, dst = lists of (lon,lat). Returns (durations BxD, src_snap[], dst_snap[]) or None."""
    coords = ";".join(f"{lo},{la}" for lo, la in (*src, *dst))
    ns = len(src)
    srcs = ";".join(str(i) for i in range(ns))
    dsts = ";".join(str(i) for i in range(ns, ns + len(dst)))
    d = _osrm_get(f"/table/v1/driving/{coords}?sources={srcs}&destinations={dsts}"
                  f"&annotations=duration")
    if d is None:
        return None
    dur = np.array(d["durations"], dtype=float)               # None → nan
    ssnap = np.array([w["distance"] for w in d["sources"]], dtype=float)
    dsnap = np.array([w["distance"] for w in d["destinations"]], dtype=float)
    return dur, ssnap, dsnap


def _check_osrm():
    if osrm_nearest(-5.696, 54.593) is None:
        sys.exit(f"ERROR: OSRM not reachable at {OSRM_HOST}:{OSRM_PORT}. Start the island instance "
                 f"(osrm-routed --algorithm mld ...) — see build_intra_times._check_osrm().")


# ── Data loading ──────────────────────────────────────────────────────────────
DZ_FILE = "simulation/dz2021/DZ2021.geojson"
SA_GLOB = "data/ireland_data/Small_Area_National_Statistical_Boundaries_2022_*.geojson"


def _load_dz_sa_polygons():
    """{area_code: geometry WGS84} for NI DZ + RoI SA only — lighter than
    build_intra_times.load_polygons (which also loads SDZ/DEA/ED/LEA + a dissolve)."""
    polys = {}
    dz = gpd.read_file(DZ_FILE)[["DZ2021_cd", "geometry"]].to_crs("EPSG:4326")
    for code, g in zip(dz["DZ2021_cd"].astype(str), dz.geometry):
        polys[code] = g
    sa_files = glob.glob(SA_GLOB)
    if not sa_files:
        sys.exit(f"ERROR: RoI SA boundary not found ({SA_GLOB})")
    sa = gpd.read_file(sa_files[0])[["SA_PUB2022", "geometry"]].to_crs("EPSG:4326")
    for code, g in zip(sa["SA_PUB2022"].astype(str), sa.geometry):
        polys[code] = g
    return polys


def load_area_masses():
    """Opportunity table joined to DZ+SA polygons (WGS84). Returns (df, geoms aligned to df rows)."""
    if not os.path.exists(OPP_TABLE):
        sys.exit(f"ERROR: {OPP_TABLE} not found — run simulation/build_opportunity_table.py first.")
    df = __import__("pandas").read_csv(OPP_TABLE, dtype={"area_code": str})
    polys = _load_dz_sa_polygons()
    geoms, keep = [], []
    for i, code in enumerate(df["area_code"].values):
        g = polys.get(code)
        if g is not None:
            geoms.append(g)
            keep.append(i)
    df = df.iloc[keep].reset_index(drop=True)
    print(f"Area masses: {len(df):,}/{len(polys):,} areas joined to DZ+SA polygons")
    return df, geoms


def load_poi_layers():
    """{'parking': (coords Nx2 lon/lat, weights), 'school': (coords, {enrol_col: weights})}."""
    park = gpd.read_file(PARKING_ISLAND_CACHE).to_crs(PROJECTED_CRS)
    park = park[park.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    tags = [c for c in park.columns if c != "geometry"]
    park["w"] = [parking_spaces({c: r[c] for c in tags}, g.area)
                 for r, g in zip(park.to_dict("records"), park.geometry)]
    park = park[park["w"] > 0].copy()
    pc = park.geometry.centroid.to_crs("EPSG:4326")
    parking = (np.column_stack([pc.x.values, pc.y.values]), park["w"].to_numpy(float))

    sch = add_level_enrolments(gpd.read_file(SCHOOL_ISLAND_CACHE).to_crs(PROJECTED_CRS))
    sc = sch.geometry.centroid.to_crs("EPSG:4326")
    scoords = np.column_stack([sc.x.values, sc.y.values])
    school = (scoords, {c: sch[c].to_numpy(float) for c in LEVEL_ENROL_COLS})
    print(f"POIs: {len(parking[1]):,} parking lots, {len(scoords):,} schools")
    return {"parking": parking, "school": school}


# ── Sampling ──────────────────────────────────────────────────────────────────
def road_point(geom, rng):
    """Uniform-in-polygon point (lon,lat) resampled until it snaps < SNAP_TOL_M to a road.
    Preserves the area's mass (always returns an in-area, road-proximate point)."""
    best, best_d = None, np.inf
    for _ in range(NEAREST_TRIES):
        lon, lat = sample_points(geom, 1, rng)[0]
        d = osrm_nearest(lon, lat)
        if d is None:
            continue
        if d < best_d:
            best, best_d = (lon, lat), d
        if d < SNAP_TOL_M:
            return (lon, lat), d
    return best, best_d                                        # best effort (flagged via snap dist)


def _draw(idx_rng, weights, n):
    return idx_rng.choice(len(weights), size=n, p=weights)


def run_purpose(name, df, geoms, pois, n_target, batch, edges, rng):
    prod_col, dest_spec = PURPOSES[name]
    pw = df[prod_col].to_numpy(float)
    if pw.sum() <= 0:
        sys.exit(f"ERROR: producer column {prod_col!r} sums to 0 for {name}.")
    pw = pw / pw.sum()

    dest_kind = dest_spec.split(":")[0]
    if dest_kind == "area":
        aw = df[dest_spec.split(":")[1]].to_numpy(float)
        aw = aw / aw.sum()
    else:                                                      # poi[:col]
        layer = dest_spec.split(":")[1]
        coords, w = pois[layer]
        if layer == "school":
            w = w[dest_spec.split(":")[2]]
        pmask = w > 0
        poi_coords, poi_w = coords[pmask], (w[pmask] / w[pmask].sum())

    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    n_pairs = n_kept = n_disc = 0
    t0 = time.time()
    while n_pairs < n_target:
        # origins: draw B areas ∝ producer, road-proximate point each
        oi = _draw(rng, pw, batch)
        src, ssnap = [], []
        for a in oi:
            (lo, la), d = road_point(geoms[a], rng)
            src.append((lo, la)); ssnap.append(d)
        # destinations
        if dest_kind == "area":
            di = _draw(rng, aw, batch)
            dst, dsnap = [], []
            for a in di:
                (lo, la), d = road_point(geoms[a], rng)
                dst.append((lo, la)); dsnap.append(d)
            dsnap = np.array(dsnap)
        else:
            di = _draw(rng, poi_w, batch)
            dst = [tuple(poi_coords[k]) for k in di]
            dsnap = None                                       # POI snap comes from /table

        res = osrm_table(src, dst)
        if res is None:
            continue
        dur, tsrc, tdst = res
        ssnap = np.array(ssnap)
        dsnap = tdst if dsnap is None else dsnap
        ok_s = ssnap < SNAP_TOL_M
        ok_d = dsnap < SNAP_TOL_M
        valid = np.isfinite(dur) & ok_s[:, None] & ok_d[None, :]
        n_pairs += batch * batch
        n_disc += int((~valid).sum())
        good = dur[valid]
        n_kept += good.size
        counts += np.histogram(good, bins=edges)[0]
        if n_pairs % (batch * batch * 40) == 0:
            el = time.time() - t0
            print(f"  [{name}] {n_pairs:,} pairs, {n_kept:,} kept "
                  f"({100*n_disc/max(n_pairs,1):.1f}% discarded), {el:.0f}s", flush=True)
    return {
        "producer": prod_col, "destination": dest_spec,
        "n_pairs": int(n_pairs), "n_kept": int(n_kept),
        "discard_rate": round(n_disc / max(n_pairs, 1), 4),
        "counts": counts.tolist(),
    }


# ── Stratified sampler (near-field distance bands + far tail) ──────────────────
# Partition pair-space by centroid haversine into a few near bands + a far tail. Sample each band
# with its own budget (near bands generous ⇒ dense head), reweight by the band's exact P·A mass, and
# reconstruct n(t) = Σ_b M_b·ŝ_b(t). Unbiased (within-band draws ∝ P_i·A_j), and the banding only
# allocates samples — final heights come from routed times × exact M_b, so the haversine proxy can't
# bias n(t). Cached road points per area make the ring sampling cheap (no per-pair /nearest).
_TO_ITM = pyproj.Transformer.from_crs("EPSG:4326", PROJECTED_CRS, always_xy=True)


def _itm(lon, lat):
    x, y = _TO_ITM.transform(np.asarray(lon, float), np.asarray(lat, float))
    return np.column_stack([x, y])


def build_point_cache(df, geoms, rng, k=K_CACHE):
    """k road-proximate points (lon,lat) per area — purpose-independent, cached to disk (resumable)."""
    codes = df["area_code"].tolist()
    if os.path.exists(POINT_CACHE):
        obj = json.load(open(POINT_CACHE))
        if obj.get("area_codes") == codes and obj.get("k") == k:
            print(f"Point cache: loaded {len(codes):,} areas × {k} from {POINT_CACHE}")
            return [[tuple(p) for p in a] for a in obj["points"]]
    print(f"Point cache: sampling {k} road points × {len(geoms):,} areas (one-off)…", flush=True)
    t0 = time.time(); pts = []
    for i, g in enumerate(geoms):
        pts.append([road_point(g, rng)[0] for _ in range(k)])
        if (i + 1) % 4000 == 0:
            print(f"  cached {i+1:,}/{len(geoms):,} ({time.time()-t0:.0f}s)", flush=True)
    json.dump({"area_codes": codes, "k": k, "points": pts}, open(POINT_CACHE, "w"))
    print(f"  saved {POINT_CACHE} ({time.time()-t0:.0f}s)")
    return pts


def _geo(name, df, pois):
    """Origin (area) + destination geometry/masses for a purpose (dest = 'area' or 'poi')."""
    prod_col, dest_spec = PURPOSES[name]
    ocent = _itm(df["centroid_lon"].values, df["centroid_lat"].values)
    prod = df[prod_col].to_numpy(float)
    kind = dest_spec.split(":")[0]
    if kind == "area":
        dw = df[dest_spec.split(":")[1]].to_numpy(float)
        dloc = ocent                                # dest index space == origin (area) index space
        dll = None
    else:
        layer = dest_spec.split(":")[1]
        coords, w = pois[layer]
        if layer == "school":
            w = w[dest_spec.split(":")[2]]
        m = w > 0
        dll = coords[m]; dw = w[m]
        dloc = _itm(dll[:, 0], dll[:, 1])
    return prod_col, dest_spec, kind, ocent, prod, dloc, dw, dll


def _band_masses(ocent, prod, dloc, dw, edges_m):
    """Per-origin in-band dest mass S[i,b] and band masses M_b=Σ_i P_i·S[i,b] (+ far tail)."""
    tree = cKDTree(dloc)
    nb = len(edges_m) - 1
    S = np.zeros((len(ocent), nb))
    nbrs = tree.query_ball_point(ocent, edges_m[-1])
    for i, idx in enumerate(nbrs):
        if not idx:
            continue
        idx = np.asarray(idx)
        d = np.hypot(dloc[idx, 0] - ocent[i, 0], dloc[idx, 1] - ocent[i, 1])
        b = np.clip(np.searchsorted(edges_m, d, side="right") - 1, 0, nb - 1)
        for bb in range(nb):
            S[i, bb] = dw[idx[b == bb]].sum()
    Sfar = np.maximum(0.0, dw.sum() - S.sum(axis=1))
    return tree, S, (prod[:, None] * S).sum(axis=0), float((prod * Sfar).sum())


def _dpt(kind, k, cache, dll, rng):
    return (cache[k][rng.integers(len(cache[k]))] if kind == "area" else tuple(dll[k]))


def _sample_near(b, edges_m, tree, S, prod, ocent, dloc, dw, kind, cache, dll, edges, budget, D, rng):
    """Band b: draw origin ∝ P·S[:,b], D dests ∝ dw from its ring, route 1×D via /table."""
    r_lo, r_hi = edges_m[b], edges_m[b + 1]
    ow = prod * S[:, b]
    if ow.sum() <= 0:
        return np.zeros(len(edges) - 1), 0
    ow = ow / ow.sum()
    hist = np.zeros(len(edges) - 1); n = 0
    while n < budget:
        i = int(rng.choice(len(ow), p=ow))
        idx = np.asarray(tree.query_ball_point(ocent[i], r_hi))
        if idx.size:
            d = np.hypot(dloc[idx, 0] - ocent[i, 0], dloc[idx, 1] - ocent[i, 1])
            idx = idx[(d >= r_lo) & (d < r_hi)]
        if idx.size == 0:
            n += D; continue
        w = dw[idx]; w = w / w.sum()
        pick = idx[rng.choice(len(idx), size=D, p=w)]
        src = [cache[i][rng.integers(len(cache[i]))]]
        dst = [_dpt(kind, int(k), cache, dll, rng) for k in pick]
        n += D
        res = osrm_table(src, dst)
        if res is None:
            continue
        dur, ts, td = res
        ok = np.isfinite(dur[0]) & (td < SNAP_TOL_M) & (ts[0] < SNAP_TOL_M)
        hist += np.histogram(dur[0][ok], bins=edges)[0]
    return hist, n


def _sample_far(r_far_m, prod, ocent, dloc, dw, kind, cache, dll, edges, budget, B, rng):
    """Far tail: naive outer-product with a centroid-distance ≥ r_far mask (cached points)."""
    ow = prod / prod.sum(); dwn = dw / dw.sum()
    hist = np.zeros(len(edges) - 1); n = 0
    while n < budget:
        oi = rng.choice(len(ow), size=B, p=ow)
        di = rng.choice(len(dwn), size=B, p=dwn)
        src = [cache[int(a)][rng.integers(len(cache[int(a)]))] for a in oi]
        dst = [_dpt(kind, int(k), cache, dll, rng) for k in di]
        n += B * B
        res = osrm_table(src, dst)
        if res is None:
            continue
        dur, ts, td = res
        dd = np.hypot(ocent[oi][:, None, 0] - dloc[di][None, :, 0],
                      ocent[oi][:, None, 1] - dloc[di][None, :, 1])
        ok = (np.isfinite(dur) & (ts[:, None] < SNAP_TOL_M) & (td[None, :] < SNAP_TOL_M)
              & (dd >= r_far_m))
        hist += np.histogram(dur[ok], bins=edges)[0]
    return hist, n


def run_purpose_stratified(name, df, pois, cache, edges, rng, near_budget, far_budget, D, B):
    prod_col, dest_spec, kind, ocent, prod, dloc, dw, dll = _geo(name, df, pois)
    if prod.sum() <= 0 or dw.sum() <= 0:
        sys.exit(f"ERROR: {name} has zero producer or attractor mass.")
    edges_m = np.array(BANDS_KM) * 1000.0
    tree, S, Mnear, Mfar = _band_masses(ocent, prod, dloc, dw, edges_m)
    rec = np.zeros(len(edges) - 1); bands = []
    for b in range(len(edges_m) - 1):
        hist, n = _sample_near(b, edges_m, tree, S, prod, ocent, dloc, dw, kind, cache, dll,
                               edges, near_budget, D, rng)
        s = hist.sum()
        if s > 0:
            rec += Mnear[b] * hist / s
        bands.append({"band_km": [BANDS_KM[b], BANDS_KM[b + 1]], "mass": float(Mnear[b]),
                      "kept": int(s), "attempts": int(n)})
        print(f"  [{name} {BANDS_KM[b]:g}-{BANDS_KM[b+1]:g}km] mass={Mnear[b]:.3e} kept={int(s):,}",
              flush=True)
    hist_f, nf = _sample_far(edges_m[-1], prod, ocent, dloc, dw, kind, cache, dll,
                             edges, far_budget, B, rng)
    sf = hist_f.sum()
    if sf > 0:
        rec += Mfar * hist_f / sf
    print(f"  [{name} >{BANDS_KM[-1]:g}km] mass={Mfar:.3e} kept={int(sf):,}", flush=True)
    return {
        "producer": prod_col, "destination": dest_spec, "mode": "stratified",
        "total_mass": float(Mnear.sum() + Mfar), "bands": bands,
        "far": {"mass": Mfar, "kept": int(sf), "attempts": int(nf)},
        "n_of_t": rec.tolist(),
    }


# ── Output + pilot plot ───────────────────────────────────────────────────────
def _profile_hash():
    try:
        return hashlib.sha1(open(OSRM_LUA, "rb").read()).hexdigest()[:12]
    except OSError:
        return "unknown"


def save(results, edges, args):
    out = {
        "_meta": {
            "date": date.today().isoformat(), "seed": SEED,
            "snap_tol_m": SNAP_TOL_M, "bin_step_s": BIN_STEP_S, "bin_cap_s": BIN_CAP_S,
            "batch": args.batch, "osrm_profile_sha1": _profile_hash(),
            "unconstrained": True, "leg": "outbound",
            "mode": "stratified" if getattr(args, "stratified", False) else "naive",
        },
        "bin_edges_s": edges.tolist(),
        "n_of_t": results,
    }
    path = args.out or (OUTPUT_FILE if not args.pilot else
                        "data/national_n_of_t_pilot.json")
    json.dump(out, open(path, "w"), indent=1)
    print(f"\nSaved {path}")
    return path


def pilot_plot(name, values, edges):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    values = np.array(values, float)
    centers = 0.5 * (edges[:-1] + edges[1:]) / 60.0            # minutes
    dens = values / max(values.sum(), 1)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    head = centers <= 10
    ax[0].bar(centers[head], values[head], width=BIN_STEP_S/60*0.9)
    ax[0].set(title=f"{name}: head (0–10 min)", xlabel="min", ylabel="n(t) [mass/bin]")
    m = centers <= 90
    ax[1].plot(centers[m], dens[m], label="empirical n(t)")
    lin = centers[m] * (centers[m] <= centers[m].max())
    lin = lin / lin.sum() * dens[m].sum()
    ax[1].plot(centers[m], lin, "--", label="n∝t (old assumption)")
    ax[1].set(title=f"{name}: n(t) vs n∝t", xlabel="min", ylabel="norm. density"); ax[1].legend()
    fig.tight_layout()
    p = f"reports/n_of_t_pilot_{name}.png"
    os.makedirs("reports", exist_ok=True)
    fig.savefig(p, dpi=110); print(f"Saved {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true", help="single purpose + diagnostic plot")
    ap.add_argument("--purpose", default="res", choices=list(PURPOSES))
    ap.add_argument("--pairs", type=int, default=1_000_000, help="naive mode: pairs/purpose")
    ap.add_argument("--batch", type=int, default=BATCH)
    ap.add_argument("--stratified", action="store_true",
                    help="near-field distance bands + far tail (dense head)")
    ap.add_argument("--near", type=int, default=NEAR_BUDGET, help="stratified: pairs per near band")
    ap.add_argument("--far", type=int, default=FAR_BUDGET, help="stratified: far-tail pairs")
    ap.add_argument("--dests", type=int, default=DESTS_PER_CALL, help="stratified: D per near /table")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    _check_osrm()
    edges = np.append(np.arange(0.0, BIN_CAP_S + BIN_STEP_S, BIN_STEP_S), np.inf)
    df, geoms = load_area_masses()
    pois = load_poi_layers()
    rng = np.random.default_rng(SEED)
    cache = build_point_cache(df, geoms, rng) if args.stratified else None

    key = "n_of_t" if args.stratified else "counts"
    names = [args.purpose] if args.pilot else list(PURPOSES)
    results = {}
    for nm in names:
        if args.stratified:
            print(f"\n=== {nm} stratified (near {args.near:,}/band, far {args.far:,}) ===")
            results[nm] = run_purpose_stratified(nm, df, pois, cache, edges, rng,
                                                 args.near, args.far, args.dests, args.batch)
        else:
            print(f"\n=== {nm} (target {args.pairs:,} pairs, batch {args.batch}) ===")
            results[nm] = run_purpose(nm, df, geoms, pois, args.pairs, args.batch, edges, rng)
        results[nm]["_edges_ref"] = "bin_edges_s"
    save(results, edges, args)
    if args.pilot:
        pilot_plot(args.purpose, results[args.purpose][key], edges)


if __name__ == "__main__":
    main()
