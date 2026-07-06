"""Build the island-wide OSM school cache → data/cache_osm_schools_island.geojson.

Single school source shared by build_demographics.py (internal core, clipped to the
core polygon) and build_census_zones.py (external zones, sjoin), mirroring build_parking.py.
Extracts amenity=school/college/university/kindergarten ways+nodes from the all-island pbf
via osmctools (reuses ni.o5m; RAM-light), then applies school_demand.assign_enrolments
GLOBALLY — light dedup for schools/kindergartens, curated-institution split for
third-level — so each institution's SOURCED total is split consistently across all its
POIs. Writes one point per KEPT feature with its enrolment, class and name.

    python3 simulation/build_schools.py

Output (gitignored, data/cache_*): point geometry + enrolment + amenity + name. **Needs
Docker + the pbf / ni.o5m.** One-off; re-run only when OSM schools change materially.
"""
import os
import subprocess
import xml.etree.ElementTree as ET

import geopandas as gpd
import pyproj
from shapely.geometry import Point
from shapely.prepared import prep

from school_demand import assign_enrolments

from demographics_config import OSRM_DIR as OSRM_ROOT   # single-source OSRM path
PBF_NAME  = "ireland-and-northern-ireland-latest.osm.pbf"
WORK_DIR  = os.path.join(OSRM_ROOT, "edge_index")
OSMCTOOLS_IMAGE = "osmctools-roaaads"
OUT_GEOJSON = "data/cache_osm_schools_island.geojson"
DEA_BOUNDARY = "simulation/dea2021/DEA2021.geojson"   # NI extent for jurisdiction tagging
PROJECTED_CRS = "EPSG:2157"
AMENITIES = {"school", "college", "university", "kindergarten"}


def _docker(cmd):
    uidgid = f"{os.getuid()}:{os.getgid()}"
    subprocess.run(["docker", "run", "--rm", "--user", uidgid,
                    "-v", f"{OSRM_ROOT}:/data", "-v", f"{WORK_DIR}:/out",
                    OSMCTOOLS_IMAGE, "sh", "-c", cmd], check=True)


def _extract():
    os.makedirs(WORK_DIR, exist_ok=True)
    o5m = os.path.join(WORK_DIR, "ni.o5m")
    osm = os.path.join(WORK_DIR, "schools.osm")
    if not os.path.exists(o5m):
        print("Converting pbf → o5m (osmconvert, streaming) …")
        _docker(f"osmconvert /data/{PBF_NAME} -t=/out/_osmconvert_tmp -o=/out/ni.o5m")
    else:
        print(f"Reusing existing {o5m}")
    print("Filtering amenity=school/college/university/kindergarten (osmfilter) …")
    _docker('osmfilter /out/ni.o5m -t=/out/_osmfilter_tmp '
            '--keep="amenity=school =college =university =kindergarten" -o=/out/schools.osm')
    print(f"  Wrote {osm}  ({os.path.getsize(osm) / 1e6:.0f} MB)")
    return osm


def _parse(osm):
    """Return list of feature dicts (amenity, name, lon, lat) for school POIs
    (nodes + ways→centroid)."""
    nodes = {}
    feats = []
    for _ev, el in ET.iterparse(osm, events=("end",)):
        if el.tag == "node":
            nodes[el.get("id")] = (float(el.get("lon")), float(el.get("lat")))
            tags = {t.get("k"): t.get("v") for t in el.findall("tag")}
            if tags.get("amenity") in AMENITIES:
                feats.append(dict(amenity=tags["amenity"], name=tags.get("name"),
                                  school=tags.get("school"),
                                  lon=float(el.get("lon")), lat=float(el.get("lat"))))
            el.clear()
        elif el.tag == "way":
            tags = {t.get("k"): t.get("v") for t in el.findall("tag")}
            if tags.get("amenity") in AMENITIES:
                refs = [r.get("ref") for r in el.findall("nd") if r.get("ref") in nodes]
                if refs:
                    lon = sum(nodes[r][0] for r in refs) / len(refs)
                    lat = sum(nodes[r][1] for r in refs) / len(refs)
                    feats.append(dict(amenity=tags["amenity"], name=tags.get("name"),
                                      school=tags.get("school"), lon=lon, lat=lat))
            el.clear()
    return feats


def main():
    osm = _extract()
    print("Parsing school features …")
    feats = _parse(osm)
    print(f"  {len(feats)} raw school POIs")

    to_itm = pyproj.Transformer.from_crs("EPSG:4326", PROJECTED_CRS, always_xy=True).transform
    for f in feats:
        f["x"], f["y"] = to_itm(f["lon"], f["lat"])

    # Tag jurisdiction (NI vs RoI) for jurisdiction-aware primary/secondary enrolment.
    print("Tagging jurisdiction from NI boundary …")
    ni = gpd.read_file(DEA_BOUNDARY).to_crs(PROJECTED_CRS)
    ni_prep = prep(ni.geometry.unary_union.simplify(200))
    del ni
    for f in feats:
        f["juris"] = "NI" if ni_prep.contains(Point(f["x"], f["y"])) else "RoI"

    kept = assign_enrolments(feats)
    print(f"  {len(kept)} features kept after dedup/curation; "
          f"total enrolment {sum(e for _f, e in kept):,.0f} pupils")

    gdf = gpd.GeoDataFrame(
        [{"amenity": f["amenity"], "name": f["name"], "enrolment": round(e, 1)}
         for f, e in kept],
        geometry=[Point(f["lon"], f["lat"]) for f, _e in kept],
        crs="EPSG:4326",
    )
    gdf.to_file(OUT_GEOJSON, driver="GeoJSON")
    print(f"Saved {len(gdf)} school POIs → {OUT_GEOJSON}")
    print("Next: re-run build_census_zones.py (external) + build_demographics.py (internal).")


if __name__ == "__main__":
    main()
