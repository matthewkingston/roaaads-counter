"""Build the island-wide OSM parking cache → data/cache_osm_parking_island.geojson.

Single parking source shared by build_demographics.py (internal core nodes, filtered
to the core polygon) and build_census_zones.py (external census zones, sjoin to each
zone polygon), so internal and external retail demand use ONE estimator with identical
tag handling (see simulation/parking_demand.py). Replaces the old per-CENTRE Overpass
download (data/cache_osm_parking.geojson) which stored only amenity/landuse/access/name
and so lacked the capacity/parking/building tags the spaces estimator needs.

Streaming, RAM-light (osmctools, ~0.5 GB — NOT osmium, whose referenced-node id-set
needs several GB; same rationale as build_network.py / build_edge_index.py). Reuses the
cached whole-island ni.o5m written by build_edge_index.py if present (skips the slow
pbf→o5m conversion). One-off: re-run only when OSM parking changes materially.

    python3 simulation/build_parking.py

Output (gitignored, matches data/cache_*): every parking polygon on the island with the
tags the estimator reads (access, parking, building, building:levels, parking:levels,
capacity, fee, amenity, landuse, name) in WGS84. **Needs Docker + the pbf / ni.o5m.**
"""
import json, os, subprocess, sys
import xml.etree.ElementTree as ET

import geopandas as gpd
from shapely.geometry import Polygon

# Reuse build_edge_index.py's o5m scratch + image (one OSM snapshot, one cache).
OSRM_ROOT = "/home/matthew/Documents/CodingFun/osrm"
PBF_NAME  = "ireland-and-northern-ireland-latest.osm.pbf"
WORK_DIR  = os.path.join(OSRM_ROOT, "edge_index")     # ni.o5m lives here
OSMCTOOLS_IMAGE = "osmctools-roaaads"                 # built by build_network.py
OUT_GEOJSON = "data/cache_osm_parking_island.geojson"

# Tags the estimator (parking_demand.parking_spaces) reads + a few for provenance.
KEEP_TAGS = ["access", "parking", "building", "building:levels", "parking:levels",
             "capacity", "fee", "amenity", "landuse", "name"]


def _docker(cmd):
    uidgid = f"{os.getuid()}:{os.getgid()}"
    subprocess.run(["docker", "run", "--rm", "--user", uidgid,
                    "-v", f"{OSRM_ROOT}:/data", "-v", f"{WORK_DIR}:/out",
                    OSMCTOOLS_IMAGE, "sh", "-c", cmd], check=True)


def _osmctools_parking_extract():
    """ni.o5m → parking.osm (amenity=parking OR landuse=parking ways + dependent
    nodes, with tags + coords). Returns the host path of the .osm."""
    os.makedirs(WORK_DIR, exist_ok=True)
    o5m = os.path.join(WORK_DIR, "ni.o5m")
    osm = os.path.join(WORK_DIR, "parking.osm")
    if not os.path.exists(o5m):
        print("Converting pbf → o5m (osmconvert, streaming) …")
        _docker(f"osmconvert /data/{PBF_NAME} -t=/out/_osmconvert_tmp -o=/out/ni.o5m")
    else:
        print(f"Reusing existing {o5m}")
    print("Filtering amenity=parking / landuse=parking ways + dependent nodes (osmfilter) …")
    _docker('osmfilter /out/ni.o5m -t=/out/_osmfilter_tmp '
            '--keep="amenity=parking landuse=parking" -o=/out/parking.osm')
    print(f"  Wrote {osm}  ({os.path.getsize(osm) / 1e6:.0f} MB)")
    return osm


def _iter_clear(osm_path):
    """Stream top-level <node>/<way> elements with bounded RAM (OSM XML lists all
    nodes before any way, so dropping each parsed top-level child keeps the root
    child-list O(1))."""
    context = ET.iterparse(osm_path, events=("start", "end"))
    _event, root = next(context)
    for event, elem in context:
        if event != "end":
            continue
        if elem.tag in ("node", "way"):
            yield elem
            root.clear()


def main():
    osm = _osmctools_parking_extract()

    print("Assembling parking polygons …")
    nodes = {}        # id → (lon, lat)
    feats = []        # GeoJSON-ish dicts for geopandas
    n_ways = n_poly = 0
    for elem in _iter_clear(osm):
        if elem.tag == "node":
            nodes[elem.get("id")] = (float(elem.get("lon")), float(elem.get("lat")))
            continue
        # way
        n_ways += 1
        refs = [nd.get("ref") for nd in elem.findall("nd")]
        if len(refs) < 4 or refs[0] != refs[-1]:
            continue                      # not a closed ring → not an area
        coords = [nodes[r] for r in refs if r in nodes]
        if len(coords) < 4:
            continue                      # missing coords (shouldn't happen post-osmfilter)
        tags = {t.get("k"): t.get("v") for t in elem.findall("tag")}
        props = {k: tags.get(k) for k in KEEP_TAGS}
        feats.append({"geometry": Polygon(coords), "props": props})
        n_poly += 1

    print(f"  {n_ways} parking ways scanned → {n_poly} closed polygons")
    gdf = gpd.GeoDataFrame(
        [f["props"] for f in feats],
        geometry=[f["geometry"] for f in feats],
        crs="EPSG:4326",
    )
    gdf.to_file(OUT_GEOJSON, driver="GeoJSON")
    print(f"Saved {len(gdf)} parking polygons → {OUT_GEOJSON}")
    print("Next: re-run build_census_zones.py (external retail) and build_demographics.py (internal retail).")


if __name__ == "__main__":
    main()
