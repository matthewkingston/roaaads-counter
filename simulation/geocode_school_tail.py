"""Geocode the null-geometry tail of `cache_admin_schools_island.geojson` via Nominatim.

The offline cascade in `build_admin_schools.py` geocodes ~98% of schools from on-disk data
(OSM cache + DZ containment). This fills the remainder — the `doubtful_match` + `unmatched`
features (no trusted coordinate) — by querying Nominatim on each school's address, trying a
cascade of progressively coarser queries (name+town → street+town+postcode → postcode/eircode
centroid) and taking the first hit.

**External queries.** Polite single-threaded use: ≤1 request/second, identifying User-Agent,
results cached to `data/cache_nominatim_schools.json` (resumable — re-runs skip cached queries
and already-placed schools). NI results are validated against the school's stated DZ
(point-in-polygon); a hit outside its DZ, or a coarse postcode-centroid hit, is kept but flagged
`needs_review`. Recovered coords are written back with `geocode_method` =
`nominatim` (name/street hit) or `nominatim_postcode` (centroid fallback).

Run: `python3 simulation/geocode_school_tail.py`  (needs network + the admin cache + DZ boundary).
"""
import json
import os
import time
import urllib.request
import urllib.parse
import warnings

import geopandas as gpd
from shapely.geometry import Point

import build_admin_schools as B

warnings.filterwarnings("ignore")

CACHE    = B.OUT
NOMCACHE = "data/cache_nominatim_schools.json"
UA = "roaaads-traffic-model/1.0 (all-Ireland gravity model; school geocode one-off)"
SLEEP = 1.1   # Nominatim usage policy: ≤ 1 req/s


def _ni_addr():
    out = {}
    for path in (B.NI_PRIMARY, B.NI_POST, B.NI_SPECIAL):
        for d in B._load_sheet(path, "Reference Data", "DENI ref"):
            nm = d.get("School name")
            if nm:
                out[nm] = {"line": str(d.get("Address Line", "") or "").strip(),
                           "town": str(d.get("Town", "") or "").strip(),
                           "postcode": str(d.get("Postcode", "") or "").strip(),
                           "dz": d.get("Datazone")}
    return out

def _roi_addr():
    out = {}
    for path, sheet, namec in [(B.ROI_PRIMARY, "Special", "Official Name"),
                               (B.ROI_PRIMARY, "Mainstream", "Official Name"),
                               (B.ROI_POST, "School Lists", "Official School Name")]:
        for d in B._load_sheet(path, sheet, "Roll Number"):
            nm = d.get(namec)
            if not nm:
                continue
            line = ", ".join(str(d.get(c, "") or "").strip() for c in
                             ("Address (Line 1)", "Address (Line 2)", "Address 1", "Address 2")
                             if d.get(c))
            out[nm] = {"line": line,
                       "county": str(d.get("County Description", d.get("County", "")) or "").strip(),
                       "eircode": str(d.get("Eircode", "") or "").strip()}
    return out


def _nominatim(query, cache):
    if query in cache:
        return cache[query]
    qs = urllib.parse.urlencode({"q": query, "format": "json", "limit": 1,
                                 "countrycodes": "gb,ie"})
    req = urllib.request.Request("https://nominatim.openstreetmap.org/search?" + qs,
                                 headers={"User-Agent": UA})
    try:
        res = json.load(urllib.request.urlopen(req, timeout=20))
        hit = {"lat": float(res[0]["lat"]), "lon": float(res[0]["lon"]),
               "display": res[0]["display_name"]} if res else None
    except Exception as e:
        hit = {"error": f"{type(e).__name__}: {e}"}
    cache[query] = hit
    time.sleep(SLEEP)
    return hit


def _queries(nm, juris, a):
    """Cascade of query strings, precise → coarse. Last element is the centroid fallback;
    its index is needed so we can flag it `needs_review`."""
    if juris == "NI":
        q = []
        if a.get("town"):
            q.append(f"{nm}, {a['town']}, Northern Ireland")
        if a.get("line") and a.get("town"):
            q.append(f"{a['line']}, {a['town']}, {a.get('postcode','')}, Northern Ireland")
        if a.get("postcode"):
            q.append(a["postcode"])           # postcode centroid (coarse fallback)
        return q, (len(q) - 1 if a.get("postcode") else -2)
    q = []
    if a.get("county"):
        q.append(f"{nm}, {a['county']}, Ireland")
    if a.get("line") and a.get("county"):
        q.append(f"{a['line']}, {a['county']}, Ireland")
    if a.get("eircode"):
        q.append(f"{a['eircode']}, Ireland")
    return q, (len(q) - 1 if a.get("eircode") else -2)


def main():
    fc = json.load(open(CACHE))
    feats = fc["features"]
    nomcache = json.load(open(NOMCACHE)) if os.path.exists(NOMCACHE) else {}
    ni_addr, roi_addr = _ni_addr(), _roi_addr()
    dz = gpd.read_file(B.DZ_BOUNDARY)[["DZ2021_nm", "geometry"]].to_crs("EPSG:2157")
    dz["dzn"] = dz["DZ2021_nm"].map(B._dznorm)
    dz_geom = dz.set_index("dzn")["geometry"].to_dict()  # metric (ITM) for distance
    FAR_M = 3000   # NI hit > this from its stated DZ ⇒ likely a same-name wrong match → reject

    tail = [f for f in feats if f["geometry"] is None]
    print(f"Tail to geocode: {len(tail)}")
    placed = coarse = dz_ok = dz_bad = far = fail = 0
    for f in tail:
        p = f["properties"]; nm = p["name"]; juris = p["jurisdiction"]
        a = (ni_addr if juris == "NI" else roi_addr).get(nm, {})
        queries, fb_idx = _queries(nm, juris, a) if a else ([nm], -2)
        hit = None; used = -1
        for i, q in enumerate(queries):
            h = _nominatim(q, nomcache)
            if h and "lat" in h and -11 < h["lon"] < -5.3 and 51.2 < h["lat"] < 55.5:
                hit = h; used = i; break
        if not hit:
            fail += 1
            p["geocode_method"] = "nominatim_fail"; p["needs_review"] = True
            continue
        lon, lat = hit["lon"], hit["lat"]
        is_centroid = used == fb_idx
        # NI: validate against the stated DZ; reject a hit too far from it (same-name wrong match).
        dist = None
        if juris == "NI" and a.get("dz"):
            poly = dz_geom.get(B._dznorm(a["dz"]))
            if poly is not None:
                E, N = B._ITM_FWD.transform(lon, lat)
                dist = poly.distance(Point(E, N))
        if dist is not None and dist > FAR_M:
            far += 1
            p["geocode_method"] = "nominatim_far"; p["needs_review"] = True
            p["matched_osm_name"] = hit["display"][:80]
            continue                              # leave geometry None → manual
        f["geometry"] = {"type": "Point", "coordinates": [lon, lat]}
        p["geocode_method"] = "nominatim_postcode" if is_centroid else "nominatim"
        p["geocode_score"] = None; p["matched_osm_name"] = hit["display"][:80]
        placed += 1; coarse += is_centroid
        if dist == 0:
            dz_ok += 1; p["needs_review"] = is_centroid   # in-DZ precise → trusted
        elif dist is not None:
            dz_bad += 1; p["needs_review"] = True          # near but outside DZ → flag
        else:
            p["needs_review"] = True                        # RoI / no DZ → flag

    json.dump(nomcache, open(NOMCACHE, "w"))
    json.dump(fc, open(CACHE, "w"))
    geocoded = sum(1 for f in feats if f["geometry"])
    print(f"\nNominatim: {placed} placed ({coarse} via postcode/eircode centroid), "
          f"{far} rejected >3km from DZ, {fail} failed")
    print(f"  NI DZ-validated: {dz_ok} | NI near but outside DZ (flagged): {dz_bad}")
    print(f"Cache now geocoded: {geocoded}/{len(feats)} ({100*geocoded/len(feats):.1f}%); "
          f"{len(feats)-geocoded} still missing.")


if __name__ == "__main__":
    main()
