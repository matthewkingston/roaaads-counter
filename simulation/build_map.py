"""
Build the interactive Newtownards model map (newtownards_map.html).

Reads the artifacts written by build_demographics.py (node_weights.json,
newtownards_demographics.geojson), the road graphs, the cached OSM POI/parking
layers, and — if present — newtownards_flows.json, and renders a folium map with
road, demographic, demand, and flow layers.

This used to be the `--map-only` path of build_demographics.py; it is now a
separate step. Run it after build_demographics.py (and after build_assignment.py
when you want the flow layers refreshed):

  python3 simulation/build_map.py
"""

import json, os, argparse, math
import geopandas as gpd
import pandas as pd
import osmnx as ox
import folium
from folium.plugins import FeatureGroupSubGroup
import pyproj
import odf.opendocument, odf.table, odf.text

from demographics_config import (
    CENTRE, OUT_DIR, GRAPH_PATH, POI_CACHE, PARKING_CACHE,
    CENSUS_ZONES_FILE, EXCLUDE_AMENITY, SCHOOL_ISLAND_CACHE,
    HIGHWAY_STYLE, ROAD_TYPE_LABELS, PROJECTED_CRS,
)
from model import walking_session_residual, EXCLUDE_LINKS
import branca.colormap as bcm

_SCHOOL_TAGS = {"school", "college", "university", "kindergarten"}

argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
).parse_args()

# Load core polygon / external zone metadata for the core-DZ map layer.
_census_zones = None
if os.path.exists(CENSUS_ZONES_FILE):
    with open(CENSUS_ZONES_FILE) as _f:
        _census_zones = json.load(_f)

_ext_node_trips = {}  # populated from flows JSON when available

# POI/parking GeoDataFrames are reloaded from cache below (see map layers).
_park_utm = None
pois_utm  = None

# ── Load saved outputs ──────────────────────────────────────────────────────────
print("Loading saved outputs …")
transformer_to_utm = pyproj.Transformer.from_crs("EPSG:4326", PROJECTED_CRS, always_xy=True)
transformer_to_wgs = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)
centre_utm_x, centre_utm_y = transformer_to_utm.transform(CENTRE[1], CENTRE[0])
print("Loading graph …")
G_cons = ox.load_graphml(GRAPH_PATH)
node_ids = list(G_cons.nodes())
node_coords_utm = [(G_cons.nodes[n]["x"], G_cons.nodes[n]["y"]) for n in node_ids]
print("Loading node weights …")
with open(f"{OUT_DIR}/node_weights.json") as _f:
    _w = json.load(_f)
_pnid = lambda k: (int(k) if k.lstrip("-").isdigit() else k)
node_population      = {_pnid(k): v for k, v in _w["node_population"].items()}
node_workplace       = {_pnid(k): v for k, v in _w["node_workplace"].items()}
node_retail_spaces   = {_pnid(k): v for k, v in _w.get("node_retail_spaces", {}).items()}
node_school_demand   = {_pnid(k): v for k, v in _w.get("node_school_demand", {}).items()}
_boundary_ids        = set(int(x) for x in _w.get("boundary_node_ids", []))
_boundary_ids_cons   = set(int(x) for x in _w.get("boundary_node_ids_cons", _w.get("boundary_node_ids", [])))
print("Loading DZ boundaries …")
dz_final = gpd.read_file(f"{OUT_DIR}/newtownards_demographics.geojson")
print(f"  {len(dz_final)} Data Zones · {len(node_ids)} nodes")

# ── Build map ───────────────────────────────────────────────────────────────────

print("Building map …")
G_raw = ox.load_graphml(f"{OUT_DIR}/newtownards_network.graphml")

m = folium.Map(location=list(CENTRE), zoom_start=14, tiles="CartoDB positron")

# ── Roads (master "Roads" group + per-type subgroups) ─────────────────────────
# A single parent FeatureGroup toggles every road type at once; each road type is
# a FeatureGroupSubGroup so it can also be toggled individually. (Folium has no
# true collapsible tree in the layer control, so the subgroups appear as their own
# entries; toggling the parent shows/hides them collectively.)
from collections import defaultdict
by_type = defaultdict(list)
for u, v, data in G_raw.edges(data=True):
    htype = data.get("highway", "unclassified")
    if isinstance(htype, list): htype = htype[0]
    by_type[htype].append((u, v, data))

type_order = [
    "living_street", "unclassified", "residential",
    "tertiary_link", "tertiary", "secondary",
    "primary_link", "primary", "trunk_link", "trunk", "motorway",
]
all_types = type_order + [t for t in by_type if t not in type_order]

roads_parent = folium.FeatureGroup(name="Roads", show=False)
roads_parent.add_to(m)
for htype in all_types:
    edges = by_type.get(htype)
    if not edges:
        continue
    style = HIGHWAY_STYLE.get(htype, {"color": "#aaaaaa", "weight": 1})
    label = ROAD_TYPE_LABELS.get(htype, f"Roads · {htype}")
    sub = FeatureGroupSubGroup(roads_parent, label, show=True)
    for u, v, data in edges:
        geom = data.get("geometry")
        if geom and hasattr(geom, "coords"):
            coords = [(lat, lon) for lon, lat in geom.coords]
        else:
            ud, vd = G_raw.nodes[u], G_raw.nodes[v]
            coords = [(ud["y"], ud["x"]), (vd["y"], vd["x"])]
        name = data.get("name", "")
        if isinstance(name, list):
            name = name[0]
        length = float(data.get("length", 0))
        folium.PolyLine(
            coords, color=style["color"], weight=style["weight"],
            opacity=0.7,
            tooltip=f"{name or '(unnamed)'} [{htype}] · {length:.0f}m",
        ).add_to(sub)
    sub.add_to(m)

# ── Node classification (boundary vs interior) ────────────────────────────────
BOUNDARY_NODE_IDS = _boundary_ids_cons  # consolidated IDs for map display

boundary_nodes_map = {}   # node_id → (wgs_lat, wgs_lon, dist, degree)
interior_nodes_map = {}

for node_id, (nx_utm, ny_utm) in zip(node_ids, node_coords_utm):
    dist = math.sqrt((nx_utm - centre_utm_x)**2 + (ny_utm - centre_utm_y)**2)
    deg  = G_cons.degree(node_id)
    nlon, nlat = transformer_to_wgs.transform(nx_utm, ny_utm)
    if node_id in BOUNDARY_NODE_IDS:
        boundary_nodes_map[node_id] = (nlat, nlon, dist, deg)
    else:
        interior_nodes_map[node_id] = (nlat, nlon, dist, deg)

all_nodes_map = {**boundary_nodes_map, **interior_nodes_map}
print(f"  Node classification: {len(boundary_nodes_map)} boundary, {len(interior_nodes_map)} interior")

# ── Nodes — single layer; boundary nodes drawn as orange squares ──────────────
nodes_fg = folium.FeatureGroup(name="Nodes", show=True)
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    is_boundary = node_id in BOUNDARY_NODE_IDS
    node_pop  = node_population.get(node_id, 0)
    node_wp   = node_workplace.get(node_id, 0)
    node_ret  = node_retail_spaces.get(node_id, 0)
    node_sch  = node_school_demand.get(node_id, 0)
    tooltip = (
        f"<b>Node {node_id}</b>{' [boundary]' if is_boundary else ''}<br>"
        f"population: {node_pop:,.0f}<br>"
        f"workplace jobs: {node_wp:,.1f}<br>"
        f"retail parking spaces: {node_ret:,.1f}<br>"
        f"school capacity: {node_sch:,.1f}"
    )
    if is_boundary:
        folium.RegularPolygonMarker(
            location=[nlat, nlon],
            number_of_sides=4, radius=6, rotation=45,
            color="#e05c00", fill=True, fill_color="#ff7c20", fill_opacity=0.9, weight=1.5,
            tooltip=tooltip,
        ).add_to(nodes_fg)
    else:
        folium.CircleMarker(
            location=[nlat, nlon], radius=3,
            color="#1a73e8", fill=True, fill_color="#1a73e8", fill_opacity=0.8, weight=1,
            tooltip=tooltip,
        ).add_to(nodes_fg)
nodes_fg.add_to(m)

# ── Metric node layers — size scaled by the metric, OFF by default ────────────
# Ranges are set from internal nodes only so boundary giants (e.g. Belfast) don't
# collapse the visible size range.
def _scaled_radius(value, vmax):
    return max(3, 14 * (value / vmax) ** 0.5)

# Population nodes
max_pop_internal = max((node_population.get(n, 0) for n in all_nodes_map), default=1) or 1
pop_fg = folium.FeatureGroup(name="Population nodes", show=False)
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    pop = node_population.get(node_id, 0)
    if pop <= 0:
        continue
    folium.CircleMarker(
        location=[nlat, nlon], radius=_scaled_radius(pop, max_pop_internal),
        color="#1a4fa0", fill=True, fill_color="#5b9bf3", fill_opacity=0.65, weight=1,
        tooltip=f"<b>Node {node_id}</b><br>population: {pop:,.1f}",
    ).add_to(pop_fg)
pop_fg.add_to(m)

# Workplace nodes (commute attractor)
max_wp_internal = max((node_workplace.get(n, 0) for n in interior_nodes_map), default=1) or 1
wp_fg = folium.FeatureGroup(name="Workplace nodes", show=False)
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    wp = node_workplace.get(node_id, 0)
    if wp <= 0:
        continue
    folium.CircleMarker(
        location=[nlat, nlon], radius=_scaled_radius(wp, max_wp_internal),
        color="#7b2d8b", fill=True, fill_color="#b05ec0", fill_opacity=0.65, weight=1,
        tooltip=f"<b>Node {node_id}</b><br>workplace jobs: {wp:,.1f}",
    ).add_to(wp_fg)
wp_fg.add_to(m)

# Retail nodes (retail attractor — parking spaces)
max_ret_internal = max((node_retail_spaces.get(n, 0) for n in interior_nodes_map), default=1) or 1
ret_fg = folium.FeatureGroup(name="Retail nodes", show=False)
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    ret = node_retail_spaces.get(node_id, 0)
    if ret <= 0:
        continue
    folium.CircleMarker(
        location=[nlat, nlon], radius=_scaled_radius(ret, max_ret_internal),
        color="#b06a00", fill=True, fill_color="#e8a020", fill_opacity=0.65, weight=1,
        tooltip=f"<b>Node {node_id}</b><br>retail parking spaces: {ret:,.1f}",
    ).add_to(ret_fg)
ret_fg.add_to(m)

# School nodes
max_sch_internal = max((node_school_demand.get(n, 0) for n in all_nodes_map), default=1) or 1
sch_fg = folium.FeatureGroup(name="School nodes", show=False)
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    sch = node_school_demand.get(node_id, 0)
    if sch <= 0:
        continue
    folium.CircleMarker(
        location=[nlat, nlon], radius=_scaled_radius(sch, max_sch_internal),
        color="#1a7a3c", fill=True, fill_color="#2ecc71", fill_opacity=0.65, weight=1,
        tooltip=f"<b>Node {node_id}</b><br>school capacity: {sch:,.1f}",
    ).add_to(sch_fg)
sch_fg.add_to(m)

# ── Core Data Zones — pop choropleth, OFF by default ──────────────────────────
dz_plot = dz_final.to_crs("EPSG:4326")
pop_vals = dz_plot["pop_estimated"].dropna().astype(float)
pop_min, pop_max = pop_vals.min(), pop_vals.max()

def pop_color(pop):
    if pd.isna(pop): return "#eeeeee"
    t = (float(pop) - pop_min) / max(pop_max - pop_min, 1)
    r = int(255 * (1 - t * 0.8))
    g = int(255 * (1 - t * 0.8))
    b = int(255 * (1 - t * 0.3))
    return f"#{r:02x}{g:02x}{b:02x}"

dz_fg = folium.FeatureGroup(name="Core Data Zones", show=False)
for _, row in dz_plot.iterrows():
    geojson_str = json.dumps(row.geometry.__geo_interface__)
    pop_full = int(row["population"]) if pd.notna(row["population"]) else None
    folium.GeoJson(
        geojson_str,
        style_function=lambda f, pop=row["pop_estimated"]: {
            "fillColor":   pop_color(pop),
            "color":       "#444444",
            "weight":      1.5,
            "fillOpacity": 0.45,
        },
        tooltip=folium.Tooltip(
            f"<b>{row['DZ2021_nm']}</b><br>"
            f"{row['DZ2021_cd']}<br>"
            f"population: {pop_full if pop_full is not None else 'N/A'}"
        ),
    ).add_to(dz_fg)
dz_fg.add_to(m)

# ── Car parks — realism check layer (red = private, blue = public/untagged) ────
_park_wgs = None
if _park_utm is not None:  # full run: set in else branch above
    _park_wgs = _park_utm.to_crs("EPSG:4326")
elif os.path.exists(PARKING_CACHE):  # --map-only: reload and reprocess from cache
    _p = gpd.read_file(PARKING_CACHE).to_crs(PROJECTED_CRS)
    _p = _p[_p.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    _p["area_m2"] = _p.geometry.area
    _p["is_private"] = _p["access"].isin(["private"]) if "access" in _p.columns else False
    _park_wgs = _p.to_crs("EPSG:4326")

if _park_wgs is not None:
    park_fg = folium.FeatureGroup(name="Car parks", show=False)
    for _, _row in _park_wgs.iterrows():
        _access_col = "access" in _park_wgs.columns
        _access_val = _row["access"] if _access_col else None
        _access_label = str(_access_val) if (_access_col and pd.notna(_access_val)) else "public/untagged"
        _area = _row["area_m2"]
        _name_val = _row["name"] if "name" in _park_wgs.columns and pd.notna(_row.get("name")) else ""
        _is_priv = bool(_row.get("is_private", False))
        _color = "#cc4444" if _is_priv else "#4488cc"
        _tip = (
            f"<b>{_name_val or 'Car park'}</b><br>"
            f"access: {_access_label}<br>"
            f"area: {_area:,.0f} m²"
        )
        folium.GeoJson(
            _row.geometry.__geo_interface__,
            style_function=lambda f, c=_color: {
                "fillColor": c, "color": c, "weight": 1.5, "fillOpacity": 0.4,
            },
            tooltip=folium.Tooltip(_tip),
        ).add_to(park_fg)
    park_fg.add_to(m)
    print(f"Added car parks layer ({len(_park_wgs)} polygons)")

# ── Businesses — amenity/shop/office POIs used for workplace allocation ────────
#     (school POIs excluded — shown in their own layer below)
_poi_wgs = None
if pois_utm is not None:  # full run: set in else branch above
    _poi_wgs = pois_utm.to_crs("EPSG:4326")
elif os.path.exists(POI_CACHE):  # --map-only: reload and re-filter from cache
    _pp = gpd.read_file(POI_CACHE)
    if "amenity" in _pp.columns:
        _pp = _pp[_pp["amenity"].isna() | ~_pp["amenity"].isin(EXCLUDE_AMENITY)]
    _pp = _pp[_pp.geometry.notna()].copy()
    _pp["geometry"] = _pp.geometry.centroid
    _poi_wgs = _pp

if _poi_wgs is not None:
    # Separate school POIs before building layers
    _school_mask = (_poi_wgs["amenity"].isin(_SCHOOL_TAGS)
                    if "amenity" in _poi_wgs.columns
                    else pd.Series(False, index=_poi_wgs.index))
    _school_pois_wgs = _poi_wgs[_school_mask].copy()
    _biz_pois_wgs    = _poi_wgs[~_school_mask].copy()

    _POI_COLOURS = {"amenity": "#e67e22", "shop": "#27ae60", "office": "#2980b9"}
    poi_fg = folium.FeatureGroup(name="Businesses", show=False)
    for _, _row in _biz_pois_wgs.iterrows():
        _kind, _val, _color = None, None, "#888888"
        for _col in ("amenity", "shop", "office"):
            if _col in _biz_pois_wgs.columns and pd.notna(_row.get(_col)):
                _kind, _val, _color = _col, _row[_col], _POI_COLOURS[_col]
                break
        _name = ""
        if "name" in _biz_pois_wgs.columns and pd.notna(_row.get("name")):
            _name = str(_row["name"])
        _type_str = f"{_kind}: {_val}" if _kind else "unknown"
        _tip = f"<b>{_name}</b><br>{_type_str}" if _name else f"<b>{_type_str}</b>"
        folium.CircleMarker(
            location=[_row.geometry.y, _row.geometry.x],
            radius=5,
            color=_color, fill=True, fill_color=_color, fill_opacity=0.85, weight=1,
            tooltip=folium.Tooltip(_tip),
        ).add_to(poi_fg)
    poi_fg.add_to(m)
    print(f"Added businesses layer ({len(_biz_pois_wgs)} non-school POIs: orange=amenity, green=shop, blue=office)")

    # ── Schools — markers from the island school cache (per-POI enrolment), OFF by default ──
    # Uses the same enrolment build_schools.py / school_demand.py compute for the model,
    # clipped to the core polygon, so the map matches node_school_demand.
    if os.path.exists(SCHOOL_ISLAND_CACHE) and _census_zones is not None:
        from shapely.geometry import Polygon as _Poly
        _core_wgs = _Poly(_census_zones["core_polygon"])
        _sch = gpd.read_file(SCHOOL_ISLAND_CACHE)
        _sch = _sch[_sch.geometry.within(_core_wgs)].copy()
        if len(_sch) > 0:
            school_poi_fg = folium.FeatureGroup(name="Schools", show=False)
            for _, _row in _sch.iterrows():
                _amenity = _row.get("amenity")
                _name = str(_row["name"]) if pd.notna(_row.get("name")) else ""
                _enroll = float(_row.get("enrolment") or 0.0)
                _label = _name or str(_amenity or "school")
                _tip = (f"<b>{_label}</b><br>type: {_amenity}<br>"
                        f"enrolment: {_enroll:,.0f} pupils")
                folium.CircleMarker(
                    location=[_row.geometry.y, _row.geometry.x],
                    radius=7,
                    color="#1a7a3c", fill=True, fill_color="#2ecc71", fill_opacity=0.85, weight=1.5,
                    tooltip=folium.Tooltip(_tip),
                ).add_to(school_poi_fg)
            school_poi_fg.add_to(m)
            print(f"Added schools layer ({len(_sch)} schools from island cache)")

# ── Optional flow layers (loaded from newtownards_flows.json if it exists) ────
_flows_path = f"{OUT_DIR}/newtownards_flows.json"
if os.path.exists(_flows_path):
    import pyproj as _pyproj
    _tr_flow = _pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)

    with open(_flows_path) as _f:
        _flows_data = json.load(_f)

    def _parse_flows(key):
        return {tuple(int(x) for x in k.split(",")): v
                for k, v in _flows_data.get(key, {}).items()}

    _link_flow         = _parse_flows("flows")
    _link_flow_res     = _parse_flows("flows_res")
    _link_flow_commute = _parse_flows("flows_commute")
    _link_flow_retail  = _parse_flows("flows_retail")
    _link_flow_school  = _parse_flows("flows_school")
    _has_components    = bool(_link_flow_res)
    _has_school_layer  = bool(_link_flow_school)
    _ext_node_trips    = _flows_data.get("ext_node_trips", {})

    # ── colour helpers ────────────────────────────────────────────────────────
    def _log_scale(flow_dict):
        """Return (lmin, lmax) from P10/P90 of a flow dict."""
        vals = sorted(v for v in flow_dict.values() if v > 0)
        if not vals:
            return 0.0, 1.0
        p10 = vals[max(0, int(len(vals) * 0.10))]
        p90 = vals[min(len(vals) - 1, int(len(vals) * 0.90))]
        return math.log10(max(p10, 1)), math.log10(max(p90, 1))

    def _t(flow, lmin, lmax):
        if flow <= 0:
            return 0.0
        return max(0.0, min(1.0, (math.log10(max(flow, 1)) - lmin) / max(lmax - lmin, 1e-6)))

    def _weight(t):
        return 1 + 7 * t

    def _color_combined(t):
        """Blue → yellow → red."""
        if t < 0.33:
            r, g, b = 0, int(180 * (t / 0.33)), int(200 * (1 - t / 0.33))
        elif t < 0.66:
            s = (t - 0.33) / 0.33
            r, g, b = int(220 * s), 180, 0
        else:
            s = (t - 0.66) / 0.34
            r, g, b = 220 + int(35 * s), int(180 * (1 - s)), 0
        return f"#{r:02x}{g:02x}{b:02x}"

    def _color_res(t):
        """Light teal → dark forest green."""
        r = int(30  + 80  * (1 - t))
        g = int(140 + 80  * (1 - t * 0.5))
        b = int(80  + 80  * (1 - t))
        return f"#{r:02x}{g:02x}{b:02x}"

    def _color_commute(t):
        """Amber → dark orange-red (commute trips)."""
        r = min(255, int(200 + 55 * t))
        g = int(160 * (1 - t * 0.85))
        b = int(20  * (1 - t))
        return f"#{r:02x}{g:02x}{b:02x}"

    def _color_retail(t):
        """Light cyan → deep teal-blue (retail trips)."""
        r = int(80  * (1 - t))
        g = int(150 + 50 * (1 - t))
        b = min(255, int(160 + 95 * t))
        return f"#{r:02x}{g:02x}{b:02x}"

    def _color_school(t):
        """Light violet → deep purple (school trips)."""
        r = int(160 + 60 * (1 - t))
        g = int(80  * (1 - t * 0.6))
        b = min(255, int(180 + 75 * t))
        return f"#{r:02x}{g:02x}{b:02x}"

    # ── helper to build one FeatureGroup ─────────────────────────────────────
    def _add_flow_fg(flow_dict, name, color_fn, tooltip_fn, show):
        lmin, lmax = _log_scale(flow_dict)
        fg = folium.FeatureGroup(name=name, show=show)
        for u, v, data in G_cons.edges(data=True):
            flow = flow_dict.get((u, v), 0) + flow_dict.get((v, u), 0)
            geom = data.get("geometry")
            if geom and hasattr(geom, "coords"):
                coords = [_tr_flow.transform(x, y)[::-1] for x, y in geom.coords]
            else:
                ud, vd = G_cons.nodes[u], G_cons.nodes[v]
                lon_u, lat_u = _tr_flow.transform(float(ud["x"]), float(ud["y"]))
                lon_v, lat_v = _tr_flow.transform(float(vd["x"]), float(vd["y"]))
                coords = [(lat_u, lon_u), (lat_v, lon_v)]
            ti = _t(flow, lmin, lmax)
            folium.PolyLine(
                coords,
                color=color_fn(ti) if flow > 0 else "#cccccc",
                weight=_weight(ti),
                opacity=0.85,
                tooltip=tooltip_fn(data.get("name", ""), flow,
                                   float(data.get("length", 0))),
            ).add_to(fg)
        fg.add_to(m)

    # ── Combined layer (always shown) ─────────────────────────────────────────
    # Tooltip built inline because it needs per-edge (u,v) to look up res/biz components.
    lmin_c, lmax_c = _log_scale(_link_flow)
    fg_combined = folium.FeatureGroup(name="Flows", show=True)
    for _u, _v, _data in G_cons.edges(data=True):
        _flow = _link_flow.get((_u, _v), 0) + _link_flow.get((_v, _u), 0)
        _geom = _data.get("geometry")
        if _geom and hasattr(_geom, "coords"):
            _coords = [_tr_flow.transform(x, y)[::-1] for x, y in _geom.coords]
        else:
            _ud, _vd = G_cons.nodes[_u], G_cons.nodes[_v]
            _lon_u, _lat_u = _tr_flow.transform(float(_ud["x"]), float(_ud["y"]))
            _lon_v, _lat_v = _tr_flow.transform(float(_vd["x"]), float(_vd["y"]))
            _coords = [(_lat_u, _lon_u), (_lat_v, _lon_v)]
        _ti = _t(_flow, lmin_c, lmax_c)
        _tip = f"{_data.get('name', '') or 'link'}<br>est. AADT: {_flow:,.0f}"
        if _has_components:
            _r = _link_flow_res.get((_u, _v), 0) + _link_flow_res.get((_v, _u), 0)
            _c = _link_flow_commute.get((_u, _v), 0) + _link_flow_commute.get((_v, _u), 0)
            _rt = _link_flow_retail.get((_u, _v), 0) + _link_flow_retail.get((_v, _u), 0)
            _s = (_link_flow_school.get((_u, _v), 0) + _link_flow_school.get((_v, _u), 0)
                  if _has_school_layer else 0)
            _tot = _r + _c + _rt + _s
            if _tot > 0:
                _tip += f"<br>&nbsp;&nbsp;residential: {_r:,.0f} ({100*_r/_tot:.0f}%)"
                _tip += f"<br>&nbsp;&nbsp;commute: {_c:,.0f} ({100*_c/_tot:.0f}%)"
                _tip += f"<br>&nbsp;&nbsp;retail: {_rt:,.0f} ({100*_rt/_tot:.0f}%)"
                if _has_school_layer:
                    _tip += f"<br>&nbsp;&nbsp;school: {_s:,.0f} ({100*_s/_tot:.0f}%)"
        _tip += f"<br>length: {float(_data.get('length', 0)):.0f}m"
        folium.PolyLine(
            _coords,
            color=_color_combined(_ti) if _flow > 0 else "#cccccc",
            weight=_weight(_ti), opacity=0.85, tooltip=_tip,
        ).add_to(fg_combined)
    fg_combined.add_to(m)

    # ── Component layers (off by default, only added when data available) ─────
    if _has_components:
        _add_flow_fg(
            _link_flow_res, "Flows - residential", _color_res,
            lambda name, flow, length: (
                f"{name or 'link'}<br>residential AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        _add_flow_fg(
            _link_flow_commute, "Flows - commute", _color_commute,
            lambda name, flow, length: (
                f"{name or 'link'}<br>commute AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        _add_flow_fg(
            _link_flow_retail, "Flows - retail", _color_retail,
            lambda name, flow, length: (
                f"{name or 'link'}<br>retail AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        if _has_school_layer:
            _add_flow_fg(
                _link_flow_school, "Flows - school", _color_school,
                lambda name, flow, length: (
                    f"{name or 'link'}<br>school AADT: {flow:,.0f}<br>length: {length:.0f}m"
                ),
                show=False,
            )
        _layer_str = "combined + res + commute + retail" + (" + school" if _has_school_layer else "")
        print(f"Added flow layers from {_flows_path} ({_layer_str}, {len(_link_flow)} links)")
    else:
        print(f"Added flow layer from {_flows_path} ({len(_link_flow)} links)")

    # ── Residuals layer (observed vs modelled; tuner-faithful per-session z) ──
    # Colours each observed directed link by the signed RMS of the exact per-session
    # Poisson z's the calibration minimises (model.walking_session_residual), so the
    # map and the objective agree by construction. Two directions of a street render
    # as side-by-side offset lines.
    _TUNED_PARAMS = f"{OUT_DIR}/tuned_params.json"
    _LINK_AADT    = "data/link_aadt.json"
    if not _has_components:
        pass  # res/biz component flows required
    elif not (os.path.exists(_TUNED_PARAMS) and os.path.exists(_LINK_AADT)):
        print("Residuals layer skipped: tuned_params.json or link_aadt.json not found")
    else:
        with open(_TUNED_PARAMS) as _f:
            _tp = json.load(_f)

        def _parse_sf(key):
            return {tuple(int(x) for x in k.split(",")): v
                    for k, v in _tp.get(key, {}).items()}

        _sf_res     = _parse_sf("slot_fracs_res")
        _sf_commute = _parse_sf("slot_fracs_commute")
        _sf_retail  = _parse_sf("slot_fracs_retail")
        _sf_sch     = _parse_sf("slot_fracs_school")

        if not (_sf_res and _sf_commute and _sf_retail):
            print("Residuals layer skipped: slot_fracs_res/commute/retail missing from tuned_params.json")
        else:
            with open(_LINK_AADT) as _f:
                _link_aadt = json.load(_f)["links"]

            _ZMAX = 4.0

            def _resid_color(z):
                # diverging: blue = model over-predicts (z>0), red = under-predicts (z<0)
                t = max(-1.0, min(1.0, z / _ZMAX))
                if t >= 0:
                    r, g, b = 220, int(220 * (1 - t)), int(220 * (1 - t))
                else:
                    s = -t
                    r, g, b = int(220 * (1 - s)), int(220 * (1 - s)), 220
                return f"#{r:02x}{g:02x}{b:02x}"

            def _resid_coords(u, v, data):
                """Edge polyline (lat,lon) offset ~6 m to the right of u→v travel,
                so (u,v) and (v,u) sit side by side."""
                OFFSET_M = 6.0
                geom = data.get("geometry") if data else None
                if geom is not None and hasattr(geom, "coords"):
                    pts = [(x, y) for x, y in geom.coords]
                else:
                    ud, vd = G_cons.nodes[u], G_cons.nodes[v]
                    pts = [(float(ud["x"]), float(ud["y"])), (float(vd["x"]), float(vd["y"]))]
                out = []
                n = len(pts)
                for i, (x, y) in enumerate(pts):
                    if i < n - 1:
                        dx, dy = pts[i + 1][0] - x, pts[i + 1][1] - y
                    else:
                        dx, dy = x - pts[i - 1][0], y - pts[i - 1][1]
                    L = math.hypot(dx, dy) or 1.0
                    nx, ny = dy / L, -dx / L          # right-hand normal of travel
                    ox, oy = x + OFFSET_M * nx, y + OFFSET_M * ny
                    lon, lat = _tr_flow.transform(ox, oy)
                    out.append((lat, lon))
                return out

            resid_fg = folium.FeatureGroup(name="Residuals (obs vs model)", show=False)
            _n_drawn = _n_excl = _n_skip = 0
            for _key, _entry in _link_aadt.items():
                _u, _v = (int(x) for x in _key.split(","))
                if G_cons.has_edge(_u, _v):
                    _data = G_cons.get_edge_data(_u, _v)[0]
                elif G_cons.has_edge(_v, _u):
                    _data = G_cons.get_edge_data(_v, _u)[0]
                elif _u in G_cons.nodes and _v in G_cons.nodes:
                    _data = None
                else:
                    _n_skip += 1
                    continue
                _lbl = (_data or {}).get("name", "")
                if isinstance(_lbl, list):
                    _lbl = _lbl[0]
                _coords = _resid_coords(_u, _v, _data)
                if (_u, _v) in EXCLUDE_LINKS:
                    folium.PolyLine(
                        _coords, color="#999999", weight=3, opacity=0.7,
                        tooltip=f"{_lbl or 'link'} ({_u}→{_v})<br>excluded from calibration",
                    ).add_to(resid_fg)
                    _n_excl += 1
                    continue
                _m_r  = _link_flow_res.get((_u, _v), 0.0)
                _m_c  = _link_flow_commute.get((_u, _v), 0.0)
                _m_rt = _link_flow_retail.get((_u, _v), 0.0)
                _m_s  = _link_flow_school.get((_u, _v), 0.0) if _has_school_layer else 0.0
                _comps = [(_m_r, _sf_res), (_m_c, _sf_commute), (_m_rt, _sf_retail)]
                if _has_school_layer:
                    _comps.append((_m_s, _sf_sch))
                _zs = []
                for _sess in _entry.get("observations", []):
                    _r = walking_session_residual(_comps, _sess)
                    if _r is not None:
                        _zs.append(_r["z"])
                if not _zs:
                    _n_skip += 1
                    continue
                _mean_z = sum(_zs) / len(_zs)
                _rms    = math.sqrt(sum(z * z for z in _zs) / len(_zs))
                _signed = math.copysign(_rms, _mean_z) if _mean_z != 0 else 0.0
                _model_aadt = _m_r + _m_c + _m_rt + _m_s
                _obs_aadt   = _entry.get("aadt", 0)
                _obs_sig    = _entry.get("aadt_uncertainty", 0)
                _nobs       = _entry.get("n_observations", len(_zs))
                _tip = (
                    f"<b>{_lbl or 'link'}</b> ({_u}→{_v})<br>"
                    f"model AADT: {_model_aadt:,.0f}<br>"
                    f"observed AADT: {_obs_aadt:,.0f} ± {_obs_sig:,.0f}<br>"
                    f"signed RMS z: {_signed:+.2f}<br>"
                    f"observations: {_nobs}"
                )
                folium.PolyLine(
                    _coords, color=_resid_color(_signed), weight=4, opacity=0.9,
                    tooltip=_tip,
                ).add_to(resid_fg)
                _n_drawn += 1
            resid_fg.add_to(m)

            _cmap = bcm.LinearColormap(
                ["#1414dc", "#ffffff", "#dc1414"], vmin=-_ZMAX, vmax=_ZMAX,
                caption="Residual signed RMS z (blue = model over-predicts, red = under-predicts)",
            )
            _cmap.add_to(m)
            print(f"Added residuals layer ({_n_drawn} links, {_n_excl} excluded, "
                  f"{_n_skip} skipped)")
else:
    print(f"No flow data found at {_flows_path} — skipping flow layer")

# ── Count sites (from ODS, Irish Grid → WGS84) ───────────────────────────────
_ODS_PATH = os.path.join(os.path.dirname(__file__), "..", "data",
                         "2023-northern-ireland-traffic-count-data-in-ods-format.ods")
# Sites used in calibration (shown with a distinct star marker)
_CALIBRATION_SITES = {"507", "508", "444"}

def _cell_text(cell):
    return " ".join(str(p) for p in cell.getElementsByType(odf.text.P)).strip()

if os.path.exists(_ODS_PATH):
    _ods_doc  = odf.opendocument.load(_ODS_PATH)
    _ods_sheets = _ods_doc.spreadsheet.getElementsByType(odf.table.Table)
    _ig_to_wgs = pyproj.Transformer.from_crs("EPSG:29902", "EPSG:4326", always_xy=True)
    _ods_sites = []
    for _sheet in _ods_sheets:
        _sname = _sheet.getAttribute("name")
        if _sname == "Mastersheet":
            continue
        _rows = _sheet.getElementsByType(odf.table.TableRow)
        _sid, _grid, _desc = None, None, None
        for _row in _rows[:8]:
            _cells = _row.getElementsByType(odf.table.TableCell)
            _vals = [_cell_text(_c) for _c in _cells]
            if _vals and _vals[0] == "Site ID":
                _sid = _vals[1] if len(_vals) > 1 else None
            elif _vals and _vals[0] == "Grid":
                _grid = _vals[1] if len(_vals) > 1 else None
            elif _vals and _vals[0] == "Description":
                _desc = _vals[1] if len(_vals) > 1 else None
        if _grid and len(_grid) == 12 and _grid.isdigit():
            _e, _n = int(_grid[:6]), int(_grid[6:])
            _lon, _lat = _ig_to_wgs.transform(_e, _n)
            _ods_sites.append({
                "id": _sname, "site_id": _sid or _sname, "desc": _desc or "",
                "lat": _lat, "lon": _lon, "e": _e, "n": _n,
            })

    count_fg = folium.FeatureGroup(name="Count sites", show=True)
    for _s in _ods_sites:
        _is_cal = _s["id"] in _CALIBRATION_SITES
        _tip = (
            f"<b>Site {_s['id']}</b>"
            + (" ★ calibration" if _is_cal else "")
            + f"<br>{_s['desc']}"
            f"<br>Irish Grid E {_s['e']:,} N {_s['n']:,}"
            f"<br>{_s['lat']:.5f}°N {abs(_s['lon']):.5f}°W"
        )
        if _is_cal:
            folium.RegularPolygonMarker(
                location=[_s["lat"], _s["lon"]],
                number_of_sides=5, radius=10, rotation=54,
                color="#8b0000", fill=True, fill_color="#e00000",
                fill_opacity=0.9, weight=2,
                tooltip=folium.Tooltip(_tip),
            ).add_to(count_fg)
        else:
            folium.CircleMarker(
                location=[_s["lat"], _s["lon"]], radius=7,
                color="#5a0080", fill=True, fill_color="#a040c0",
                fill_opacity=0.75, weight=1.5,
                tooltip=folium.Tooltip(_tip),
            ).add_to(count_fg)
    count_fg.add_to(m)
    print(f"Added count sites layer ({len(_ods_sites)} sites, "
          f"{sum(1 for s in _ods_sites if s['id'] in _CALIBRATION_SITES)} calibration)")
else:
    print(f"ODS file not found at {_ODS_PATH} — skipping count sites layer")

# ── External nodes ────────────────────────────────────────────────────────────
if _census_zones is not None:
    _ext_zones = _census_zones["external_nodes"]
    _max_ext_tot = max(
        (_ext_node_trips.get(_ez["id"], {}).get("trips_through", 0)
         + _ext_node_trips.get(_ez["id"], {}).get("trips_internal", 0))
        for _ez in _ext_zones
    ) if _ext_node_trips else 1.0
    ext_fg = folium.FeatureGroup(name="External nodes", show=False)
    for _ez in _ext_zones:
        _nid    = _ez["id"]
        _ez_pop = node_population.get(_nid, 0)
        _ez_wp  = node_workplace.get(_nid, 0)
        _ez_ret = node_retail_spaces.get(_nid, 0)
        _ez_sch = node_school_demand.get(_nid, 0)
        _trips  = _ext_node_trips.get(_nid, {})
        _t_int  = _trips.get("trips_internal", 0)
        _t_thr  = _trips.get("trips_through",  0)
        _t_tot  = _t_int + _t_thr
        _radius = max(4, 12 * (_t_tot / max(_max_ext_tot, 1)) ** 0.4) if _t_tot > 0 else 4
        _tip = (
            f"<b>{_nid}</b> [{_ez['level']}]<br>"
            f"population: {_ez_pop:,.0f}<br>"
            f"workplace jobs: {_ez_wp:,.1f}<br>"
            f"retail parking spaces: {_ez_ret:,.1f}<br>"
            f"school demand: {_ez_sch:,.1f}"
        )
        if _trips:
            _tip += (
                f"<br>trips sent: {_t_tot:,.0f}"
                f"<br>&nbsp;&nbsp;to internal: {_t_int:,.0f}"
                f"<br>&nbsp;&nbsp;through-trips: {_t_thr:,.0f}"
            )
        folium.CircleMarker(
            location=[_ez["centroid_lat"], _ez["centroid_lon"]],
            radius=_radius,
            color="#2c6e49", fill=True, fill_color="#40916c", fill_opacity=0.7, weight=1.5,
            tooltip=folium.Tooltip(_tip),
        ).add_to(ext_fg)
    ext_fg.add_to(m)
    print(f"Added external nodes layer ({len(_ext_zones)} nodes"
          + (f", {len(_ext_node_trips)} with trip data" if _ext_node_trips else ", no trip data")
          + ")")

folium.LayerControl(collapsed=False).add_to(m)

out_path = f"{OUT_DIR}/newtownards_map.html"
m.save(out_path)
total_node_pop = sum(node_population.values())
print(f"\nSaved: {out_path}")
print(f"  {len(dz_plot)} Data Zones · est. pop range {pop_min:.0f}–{pop_max:.0f}")
print(f"  {len(node_population)} nodes with population assigned · total {total_node_pop:.0f}")
