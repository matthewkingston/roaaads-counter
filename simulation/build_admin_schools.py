"""Build the island-wide **admin-roll** school cache (school-age: primary / post-primary /
special), the data-anchored replacement for the OSM school-age estimate in `build_schools.py`.

Output (gitignored): `data/cache_admin_schools_island.geojson` — one Point feature per school:
  properties: jurisdiction (NI|RoI), level (primary|post_primary|special), name,
              enrolment (school-age pupils), geocode_method, geocode_score, needs_review
  geometry:   Point(lon,lat) in WGS84, or null for the unmatched NI tail (manual/Nominatim later)

This is **school-age only** — third-level stays in the OSM/curated path for now (the lumped
Stage-1 milestone; tertiary is split out and FT-curated in Stage 2). Level is tagged per school
so the eventual primary/post-primary/tertiary split needs no re-ingest.

Data sources (gitignored, on disk):
  NI  `data/ni_data/School level - {Primary…,post primary…,special…}.{xlsx,XLSX}`
      school-age primary = `Total enrolment − Total: Nursery FT − Total: Nursery PT
                            − Total: Pre-school age`; post-primary = `Total enrolment`;
      special total from the `Sex` sheet. No coords → geocoded (below).
  RoI `data/ireland_data/Data_on_Individual_Schools_{Primary_Mainstream,post_primary}.xlsx`
      enrolment = `Enrolment per Return` (primary + `Special` tab) / `Total 2025-2026`
      (post-primary); **in-file `School Latitude/Longitude`**. Pure-`Boarding` post-primary
      schools are dropped (no daily home→school car trip); `Mixed` kept (no per-pupil split).

NI geocoding (no coords in the rolls) — uniform island-wide, **prototype-validated** cascade,
zero external queries here (Nominatim tail handled separately/manually):
  match each school to a same-jurisdiction OSM school POI (`cache_osm_schools_island.geojson`)
  by name, gated by (a) `amenity` vs level compatibility and (b) **DZ containment** (the school's
  `Datazone` must contain the OSM POI). Type-preserving name normalisation. Score banding:
  exact or fuzzy ≥ 0.85 → confident; 0.70–0.85 → accepted but `needs_review`; < 0.70 or no
  compatible POI in the DZ → unmatched tail (null geometry, `needs_review`).

Run: `python3 simulation/build_admin_schools.py`  (needs the OSM school cache + DZ boundary).
"""
import json
import re
import difflib
import warnings
from collections import defaultdict

import unicodedata
import openpyxl
import geopandas as gpd
from shapely.geometry import Point, mapping
from pyproj import Transformer

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
NI_PRIMARY  = "data/ni_data/School level - Primary schools data 202526 - Revised 3 June 2026.XLSX"
NI_POST     = "data/ni_data/School level - post primary schools 202526.xlsx"
NI_SPECIAL  = "data/ni_data/School level - special schools 202526.xlsx"
ROI_PRIMARY = "data/ireland_data/Data_on_Individual_Schools_Primary_Mainstream.xlsx"
ROI_POST    = "data/ireland_data/Data_on_Individual_Schools_post_primary.xlsx"
OSM_CACHE   = "data/cache_osm_schools_island.geojson"
DZ_BOUNDARY = "simulation/dz2021/DZ2021.geojson"
OUT         = "data/cache_admin_schools_island.geojson"

NI_LON0, NI_LON1, NI_LAT0, NI_LAT1 = -8.3, -5.3, 54.0, 55.4
FUZZY_CONFIDENT = 0.85   # exact or ≥ this → no review flag
FUZZY_ACCEPT    = 0.70   # below this (or no compatible POI in DZ) → unmatched tail
# admin level → compatible OSM amenity tags (school-age never matches college/university)
COMPAT = {"primary": {"school", "kindergarten"},
          "post_primary": {"school"},
          "special": {"school"}}

# --- "same school?" gate (beyond fuzzy score) -------------------------------------
# A high fuzzy score can still be a *different* school sharing tokens (St Matthew's vs
# St Malachy's; a main school vs its prep dept). So a fuzzy candidate is accepted only if its
# DISTINCTIVE name agrees: type words carry no identity and are stripped; the remaining "core"
# tokens must be equal; and "role" modifiers (which denote a separate establishment) must match.
# Anything else → demoted to the tail (better unmatched than wrongly placed; some schools simply
# aren't in OSM).
_NAME_TYPE = set("school schools primary post secondary college high grammar academy integrated "
                 "controlled maintained community national ns ps gaelscoil bunscoil naiscoil "
                 "scoil colaiste coliste pobalscoil nursery unit naomh special".split())
_NAME_ROLE = set("prep preparatory infant infants junior senior assessment resource centre "
                 "annexe model".split())
def _core_roles(nrm):
    toks = set(nrm.split())
    return (toks - _NAME_TYPE - _NAME_ROLE), (toks & _NAME_ROLE)
def _stem(t):
    return t[:-1] if len(t) > 4 and t.endswith("s") else t   # drop possessive/plural 's
def _sig(core):
    return "".join(sorted(_stem(t) for t in core))           # order/spacing-independent signature
def _same_school(a_norm, b_norm):
    """Same establishment? Core distinctive tokens must agree (modulo possessive-s, spacing,
    word order — accents already removed in _norm) and role modifiers must match exactly.
    Rejects different schools sharing tokens (St Matthew's vs St Malachy's) and disambiguated
    siblings (St Patrick's vs St Patrick's (Glen))."""
    ca, ra = _core_roles(a_norm); cb, rb = _core_roles(b_norm)
    if not ca or not cb or ra != rb:
        return False
    sa, sb = _sig(ca), _sig(cb)
    return sa == sb or difflib.SequenceMatcher(None, sa, sb).ratio() >= 0.90

# ── Name / DZ normalisation + similarity (type-PRESERVING — keep grammar/college/high/…) ──
_GENERIC = set("the of a and national co".split())
def _deaccent(s):
    return unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode()
def _norm(s):
    s = _deaccent(s).split(",")[0].lower().replace("&", " and ").replace("'", "")
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return " ".join("st" if t == "saint" else t for t in s.split() if t not in _GENERIC)
def _dznorm(s):
    return re.sub(r"\s+", "", str(s).lower())
def _sim(a, b):
    ta, tb = set(a.split()), set(b.split())
    jac = len(ta & tb) / len(ta | tb) if (ta | tb) else 0.0
    return max(jac, difflib.SequenceMatcher(None, a, b).ratio())


# ── Generic xlsx sheet reader (auto-detect header row by a sentinel column) ─────
def _load_sheet(path, sheet, sentinel):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    rows = list(ws.iter_rows(values_only=True))
    hi = None
    for i, r in enumerate(rows[:15]):
        if r and any(str(c).strip() == sentinel for c in r if c is not None):
            hi = i
            break
    if hi is None:
        raise ValueError(f"{path}::{sheet}: no header row containing {sentinel!r}")
    hdr = [str(c).strip() if c is not None else "" for c in rows[hi]]
    out = []
    for r in rows[hi + 1:]:
        if not r or all(c is None for c in r):
            continue
        d = dict(zip(hdr, r))
        if not d.get(sentinel):
            continue
        out.append(d)
    return out

def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


# ── NI ingest (school-age enrolment + name + datazone, joined on DENI ref) ──────
def _ni_join(enrol_path, enrol_sheet, enrol_fn, ref_path=None):
    """Return list of dict(name, datazone, enrolment) for one NI file."""
    ref_path = ref_path or enrol_path
    enrol = {d["DENI ref"]: enrol_fn(d) for d in _load_sheet(enrol_path, enrol_sheet, "DENI ref")}
    ref   = {d["DENI ref"]: d for d in _load_sheet(ref_path, "Reference Data", "DENI ref")}
    out = []
    for ref_id, rd in ref.items():
        out.append({"name": rd["School name"], "datazone": rd.get("Datazone"),
                    "enrolment": enrol.get(ref_id, 0.0)})
    return out

def load_ni_schools():
    """{level: [ {name, datazone, enrolment} ]} for NI primary / post_primary / special."""
    primary = _ni_join(
        NI_PRIMARY, "Enrolments",
        lambda d: _num(d.get("Total enrolment")) - _num(d.get("Total: Nursery FT"))
                  - _num(d.get("Total: Nursery PT")) - _num(d.get("Total: Pre-school age")))
    post = _ni_join(NI_POST, "Enrolments", lambda d: _num(d.get("Total enrolment")))
    special = _ni_join(NI_SPECIAL, "Sex", lambda d: _num(d.get("Total enrolment")))
    return {"primary": primary, "post_primary": post, "special": special}


# ── RoI ingest (enrolment + in-file coords; drop pure boarders) ─────────────────
# Coordinate-column quirks in the RoI post-primary file (the primary file is clean). Two
# malformed patterns, both recovered exactly and validated; anything else → null (flagged):
#   (A) ~108 rows: `School Latitude` = ITM easting, `School Longitude` = WGS84 latitude
#       (WGS longitude absent) → recover longitude by inverting the easting at that latitude.
#   (B)  ~11 rows (Gaeltacht schools): both columns are a full ITM pair — `School Latitude`
#       = Northing, `School Longitude` = Easting → inverse-transform the pair to WGS84.
_ITM_FWD = Transformer.from_crs("EPSG:4326", "EPSG:2157", always_xy=True)  # lon,lat → E,N
_ITM_INV = Transformer.from_crs("EPSG:2157", "EPSG:4326", always_xy=True)  # E,N → lon,lat

def _recover_lon_from_itm_easting(easting, lat):
    lo, hi = -11.0, -5.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if _ITM_FWD.transform(mid, lat)[0] < easting:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def _on_island(lon, lat):
    return lon is not None and -11 < lon < -5.3 and 51.2 < lat < 55.5

def _roi_coords(lat_raw, lon_raw):
    """Resolve one RoI row's coordinate columns → (lon, lat) WGS84, or (None, None)."""
    try:
        a, b = float(lat_raw), float(lon_raw)
    except (TypeError, ValueError):
        return (None, None)
    if 51 < a < 56 and -11 < b < -5:                 # normal: a=lat, b=lon
        return (b, a)
    if a > 400_000 and 51 < b < 56:                  # (A) a=ITM easting, b=lat
        lon = _recover_lon_from_itm_easting(a, b)
        return (lon, b) if _on_island(lon, b) else (None, None)
    if a > 400_000 and b > 400_000:                  # (B) full ITM: a=Northing, b=Easting
        lon, lat = _ITM_INV.transform(b, a)
        return (lon, lat) if _on_island(lon, lat) else (None, None)
    return (None, None)                              # uninterpretable → flag, don't guess

def _roi_rows(path, sheet, name_col, enrol_col):
    out = []
    for d in _load_sheet(path, sheet, "Roll Number"):
        if str(d.get("Pupil Attendance Type", "")).strip().lower() == "boarding":
            continue  # boarders generate no daily home→school car trip
        lon, lat = _roi_coords(d.get("School Latitude"), d.get("School Longitude"))
        out.append({"name": d.get(name_col), "lat": lat, "lon": lon,
                    "enrolment": _num(d.get(enrol_col))})
    return out

def load_roi_schools():
    """{level: [ {name, lat, lon, enrolment} ]} for RoI primary / post_primary / special."""
    return {
        "primary":      _roi_rows(ROI_PRIMARY, "Mainstream", "Official Name", "Enrolment per Return"),
        "special":      _roi_rows(ROI_PRIMARY, "Special",    "Official Name", "Enrolment per Return"),
        "post_primary": _roi_rows(ROI_POST,    "School Lists", "Official School Name", "Total 2025-2026"),
    }


# ── NI geocoder: OSM cache → DZ-contained, amenity-gated, name-matched ──────────
def build_ni_geocoder():
    """Return by_dz: {dz_norm: [(norm_name, amenity, lon, lat)]} for NI OSM school POIs."""
    osm_raw = json.load(open(OSM_CACHE))
    pts, nn, raw, ams = [], [], [], []
    for f in osm_raw["features"]:
        g = f["geometry"]; nm = f["properties"].get("name"); am = f["properties"].get("amenity")
        if g["type"] != "Point" or not nm:
            continue
        lon, lat = g["coordinates"]
        if NI_LON0 < lon < NI_LON1 and NI_LAT0 < lat < NI_LAT1:
            pts.append(Point(lon, lat)); nn.append(_norm(nm)); raw.append(nm); ams.append(am)
    osm = gpd.GeoDataFrame({"nn": nn, "raw": raw, "am": ams}, geometry=pts, crs="EPSG:4326")
    dz = gpd.read_file(DZ_BOUNDARY)[["DZ2021_nm", "geometry"]].to_crs("EPSG:4326")
    dz["dzn"] = dz["DZ2021_nm"].map(_dznorm)
    osm = gpd.sjoin(osm, dz[["dzn", "geometry"]], how="left", predicate="within")
    by_dz = defaultdict(list)
    for _, r in osm.iterrows():
        by_dz[r["dzn"]].append((r["nn"], r["raw"], r["am"], r.geometry.x, r.geometry.y))
    return by_dz

def geocode_ni(name, datazone, level, by_dz):
    """→ (lon, lat, method, score, matched_osm_name).

    Tail (None coords) when: no compatible POI in the DZ; best fuzzy < FUZZY_ACCEPT; or the
    best fuzzy candidate fails the `_same_school` identity check (different school / prep dept /
    etc.) — demoted on doubt rather than placed wrongly.
    """
    nn = _norm(name); dzt = _dznorm(datazone)
    cands = [c for c in by_dz.get(dzt, []) if c[2] in COMPAT[level]]
    if not cands:
        return (None, None, "unmatched", 0.0, None)
    exact = [c for c in cands if c[0] == nn]
    if exact:
        c = exact[0]
        return (c[3], c[4], "exact_dz_type", 1.0, c[1])
    best = max(cands, key=lambda c: _sim(nn, c[0])); score = round(_sim(nn, best[0]), 3)
    if score >= FUZZY_ACCEPT and _same_school(nn, best[0]):
        return (best[3], best[4], "fuzzy_dz_type", score, best[1])
    reason = "doubtful_match" if score >= FUZZY_ACCEPT else "unmatched"
    return (None, None, reason, score, best[1])


# ── Assemble + write ───────────────────────────────────────────────────────────
def main():
    print("Building NI geocoder from OSM cache + DZ boundaries …")
    by_dz = build_ni_geocoder()

    features = []
    stats = defaultdict(lambda: defaultdict(float))  # (juris,level) -> {n, enrol, unmatched, review}

    print("Ingesting NI rolls + geocoding …")
    for level, schools in load_ni_schools().items():
        for s in schools:
            lon, lat, method, score, matched = geocode_ni(s["name"], s["datazone"], level, by_dz)
            has = lon is not None
            geom = mapping(Point(lon, lat)) if has else None
            features.append({"type": "Feature", "geometry": geom, "properties": {
                "jurisdiction": "NI", "level": level, "name": s["name"],
                "enrolment": round(s["enrolment"], 1), "geocode_method": method,
                "geocode_score": score, "matched_osm_name": matched,
                "needs_review": method == "fuzzy_dz_type"}})
            k = ("NI", level)
            stats[k]["n"] += 1; stats[k]["enrol"] += s["enrolment"]
            stats[k]["fuzzy"] += method == "fuzzy_dz_type"
            stats[k]["tail"] += not has

    print("Ingesting RoI rolls (in-file coords) …")
    for level, schools in load_roi_schools().items():
        for s in schools:
            has = s["lon"] is not None
            geom = mapping(Point(s["lon"], s["lat"])) if has else None
            features.append({"type": "Feature", "geometry": geom, "properties": {
                "jurisdiction": "RoI", "level": level, "name": s["name"],
                "enrolment": round(s["enrolment"], 1),
                "geocode_method": "roi_infile" if has else "unmatched",
                "geocode_score": 1.0 if has else 0.0, "matched_osm_name": None,
                "needs_review": not has}})
            k = ("RoI", level)
            stats[k]["n"] += 1; stats[k]["enrol"] += s["enrolment"]
            stats[k]["tail"] += not has

    with open(OUT, "w") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)

    # ── Summary ──
    print(f"\nWrote {OUT}  ({len(features)} schools)")
    print(f"{'juris/level':22} {'n':>5} {'enrol':>9} {'fuzzy':>6} {'tail':>5}")
    tot_n = tot_e = tot_t = 0
    for k in sorted(stats):
        st = stats[k]
        print(f"{k[0]+'/'+k[1]:22} {int(st['n']):5} {int(st['enrol']):9,} "
              f"{int(st['fuzzy']):6} {int(st['tail']):5}")
        tot_n += st["n"]; tot_e += st["enrol"]; tot_t += st["tail"]
    print(f"{'TOTAL':22} {int(tot_n):5} {int(tot_e):9,} {'':6} {int(tot_t):5}")
    print(f"\nGeocoded: {int(tot_n - tot_t)}/{int(tot_n)} ({100*(tot_n-tot_t)/tot_n:.1f}%); "
          f"{int(tot_t)} tail (null geometry — manual/Nominatim). "
          f"'fuzzy' = accepted same-school variants (audit via matched_osm_name).")


if __name__ == "__main__":
    main()
