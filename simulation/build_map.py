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
import pyproj
import odf.opendocument, odf.table, odf.text

from demographics_config import (
    CENTRE, OUT_DIR, GRAPH_PATH, DZ_BOUNDARY_FILE, POI_CACHE, PARKING_CACHE,
    CENSUS_ZONES_FILE, EXCLUDE_AMENITY, SCHOOL_ENROLL_FALLBACK,
    HIGHWAY_STYLE, ROAD_TYPE_LABELS, PROJECTED_CRS,
)

_SCHOOL_TAGS = set(SCHOOL_ENROLL_FALLBACK)

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
node_business_demand = {_pnid(k): v for k, v in _w["node_business_demand"].items()}
node_parking_equiv   = {_pnid(k): v for k, v in _w.get("node_parking_equiv", {}).items()}
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

# 5a. Road network edges (from raw geographic graph)
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

for htype in all_types:
    edges = by_type.get(htype)
    if not edges:
        continue
    style = HIGHWAY_STYLE.get(htype, {"color": "#aaaaaa", "weight": 1})
    label = ROAD_TYPE_LABELS.get(htype, f"Roads · {htype}")
    fg = folium.FeatureGroup(name=label, show=False)
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
        ).add_to(fg)
    fg.add_to(m)

# 5b. Classify nodes as boundary (in/out flow) vs interior.
#
# Set established by auto-detection followed by manual inspection.
# To update: add/remove node IDs from BOUNDARY_NODE_IDS.

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

print(f"  Node classification: {len(boundary_nodes_map)} boundary, {len(interior_nodes_map)} interior")

# Boundary nodes — ON by default, orange diamonds
boundary_fg = folium.FeatureGroup(name=f"In/Out flow nodes — auto-detected ({len(boundary_nodes_map)})", show=True)
for node_id, (nlat, nlon, dist, deg) in boundary_nodes_map.items():
    node_pop = node_population.get(node_id, 0)
    node_biz = node_business_demand.get(node_id, 0)
    node_sch = node_school_demand.get(node_id, 0)
    _sch_str = f"<br>school pupils: {node_sch:.0f}" if node_sch > 0 else ""
    folium.RegularPolygonMarker(
        location=[nlat, nlon],
        number_of_sides=4, radius=6, rotation=45,
        color="#e05c00", fill=True, fill_color="#ff7c20", fill_opacity=0.9, weight=1.5,
        tooltip=(
            f"<b>Node {node_id}</b> [boundary]<br>"
            f"degree={deg} · {dist:.0f}m from centre<br>"
            f"pop: {node_pop:,.0f}<br>"
            f"workplace: {node_biz:,.0f}"
            f"{_sch_str}"
        ),
    ).add_to(boundary_fg)
boundary_fg.add_to(m)

# Interior nodes — ON by default, small blue dots
interior_fg = folium.FeatureGroup(name=f"Interior nodes ({len(interior_nodes_map)})", show=True)
for node_id, (nlat, nlon, dist, deg) in interior_nodes_map.items():
    node_pop = node_population.get(node_id, 0)
    node_biz = node_business_demand.get(node_id, 0)
    node_sch = node_school_demand.get(node_id, 0)
    _sch_str = f"<br>school pupils: {node_sch:.1f}" if node_sch > 0 else ""
    folium.CircleMarker(
        location=[nlat, nlon], radius=3,
        color="#1a73e8", fill=True, fill_color="#1a73e8", fill_opacity=0.8, weight=1,
        tooltip=(
            f"Node {node_id} · degree={deg} · {dist:.0f}m from centre<br>"
            f"est. pop: {node_pop:.1f}<br>"
            f"workplace pop: {node_biz:.1f}"
            f"{_sch_str}"
        ),
    ).add_to(interior_fg)
interior_fg.add_to(m)

# 5d-pre. Core area DZ boundaries — actual DZs in the core polygon from census_zones.json,
# not the old 3km-clipped set.  Useful for seeing exactly which DZs were pulled into core.
if _census_zones is not None:
    from shapely.geometry import Polygon as _CorePoly
    _core_poly_wgs_map = _CorePoly(_census_zones["core_polygon"])
    _dz_all_wgs = gpd.read_file(DZ_BOUNDARY_FILE).to_crs("EPSG:4326")
    _dz_core_map = _dz_all_wgs[_dz_all_wgs.geometry.centroid.within(_core_poly_wgs_map)].copy()
    core_dz_fg = folium.FeatureGroup(
        name=f"Core area DZs — from census_zones.json ({len(_dz_core_map)})", show=False)
    for _, _row in _dz_core_map.iterrows():
        folium.GeoJson(
            _row.geometry.__geo_interface__,
            style_function=lambda f: {
                "fillColor": "#ff6600", "color": "#cc4400",
                "weight": 2, "fillOpacity": 0.12,
            },
            tooltip=folium.Tooltip(
                f"<b>{_row.get('DZ2021_nm','')}</b><br>{_row.get('DZ2021_cd','')}"),
        ).add_to(core_dz_fg)
    core_dz_fg.add_to(m)

# 5d. DZ choropleth (estimated population in clipped area)
# Convert clipped DZs back to WGS84 for Leaflet rendering
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

dz_fg = folium.FeatureGroup(name=f"Data Zones — pop estimate ({len(dz_plot)})", show=False)
for _, row in dz_plot.iterrows():
    geojson_str = json.dumps(row.geometry.__geo_interface__)
    pop_est = int(row["pop_estimated"]) if pd.notna(row["pop_estimated"]) else None
    pop_full = int(row["population"]) if pd.notna(row["population"]) else None
    area_pct = row["area_pct"]
    clipped = area_pct < 0.99
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
            f"Code: {row['DZ2021_cd']}<br>"
            f"Est. pop (clipped): {pop_est if pop_est else 'N/A'}<br>"
            + (f"Census pop (full DZ): {pop_full} · {area_pct*100:.0f}% of DZ in study area<br>" if clipped else
               f"Census pop: {pop_full}<br>")
            + f"Area in study: {row['area_clipped_m2']/10000:.1f} ha"
        ),
    ).add_to(dz_fg)
dz_fg.add_to(m)

# 5e. Business demand nodes — proportional purple circles, OFF by default
# Scale based on internal nodes only so boundary giants (Belfast) don't collapse the range.
max_biz_internal = max((node_business_demand.get(n, 0) for n in interior_nodes_map), default=1)
biz_fg = folium.FeatureGroup(name="Business demand nodes", show=False)
all_nodes_map = {**boundary_nodes_map, **interior_nodes_map}
for node_id, (nlat, nlon, dist, deg) in all_nodes_map.items():
    biz = node_business_demand.get(node_id, 0)
    if biz <= 0:
        continue
    radius = max(3, 14 * (biz / max_biz_internal) ** 0.5)
    node_pop  = node_population.get(node_id, 0)
    park_eq   = node_parking_equiv.get(node_id, 0)
    wp_demand = biz - park_eq
    folium.CircleMarker(
        location=[nlat, nlon], radius=radius,
        color="#7b2d8b", fill=True, fill_color="#b05ec0", fill_opacity=0.65, weight=1,
        tooltip=(
            f"<b>Node {node_id}</b><br>"
            f"workplace population: {wp_demand:.1f}<br>"
            f"parking eq.: {park_eq:.1f}<br>"
            f"business demand: {biz:.1f}<br>"
            f"est. pop: {node_pop:.1f}"
        ),
    ).add_to(biz_fg)
biz_fg.add_to(m)

# 5f. Parking polygons — realism check layer (red = private, blue = public/untagged)
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
    park_fg = folium.FeatureGroup(name=f"Car parks — OSM ({len(_park_wgs)})", show=False)
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
    print(f"Added parking layer ({len(_park_wgs)} polygons)")

# 5g. POI layer — amenity/shop/office features used for workplace allocation
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
    poi_fg = folium.FeatureGroup(name=f"POIs — workplace allocation ({len(_biz_pois_wgs)})", show=False)
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
    print(f"Added POI layer ({len(_biz_pois_wgs)} non-school POIs: orange=amenity, green=shop, blue=office)")

    # 5h. School POI layer — schools with enrollment, green circles, OFF by default
    if len(_school_pois_wgs) > 0:
        school_poi_fg = folium.FeatureGroup(name=f"Schools — OSM ({len(_school_pois_wgs)})", show=False)
        for _, _row in _school_pois_wgs.iterrows():
            _amenity = _row.get("amenity") if "amenity" in _school_pois_wgs.columns else None
            _name = ""
            if "name" in _school_pois_wgs.columns and pd.notna(_row.get("name")):
                _name = str(_row["name"])
            # Resolve enrollment: OSM capacity if present and numeric, else fallback
            _cap_raw = _row.get("capacity") if "capacity" in _school_pois_wgs.columns else None
            if pd.notna(_cap_raw):
                try:
                    _enroll = int(float(_cap_raw))
                    _enroll_src = "OSM capacity"
                except (ValueError, TypeError):
                    _enroll = SCHOOL_ENROLL_FALLBACK.get(_amenity, 300)
                    _enroll_src = "fallback"
            else:
                _enroll = SCHOOL_ENROLL_FALLBACK.get(_amenity, 300)
                _enroll_src = "fallback"
            _label = _name or str(_amenity or "school")
            _tip = (
                f"<b>{_label}</b><br>"
                f"type: {_amenity}<br>"
                f"enrollment: {_enroll} pupils ({_enroll_src})"
            )
            folium.CircleMarker(
                location=[_row.geometry.y, _row.geometry.x],
                radius=7,
                color="#1a7a3c", fill=True, fill_color="#2ecc71", fill_opacity=0.85, weight=1.5,
                tooltip=folium.Tooltip(_tip),
            ).add_to(school_poi_fg)
        school_poi_fg.add_to(m)
        print(f"Added school POI layer ({len(_school_pois_wgs)} schools)")

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

    _link_flow        = _parse_flows("flows")
    _link_flow_res    = _parse_flows("flows_res")
    _link_flow_biz    = _parse_flows("flows_biz")
    _link_flow_school = _parse_flows("flows_school")
    _has_components   = bool(_link_flow_res)
    _has_school_layer = bool(_link_flow_school)
    _ext_node_trips   = _flows_data.get("ext_node_trips", {})

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

    def _color_biz(t):
        """Amber → dark orange-red."""
        r = min(255, int(200 + 55 * t))
        g = int(160 * (1 - t * 0.85))
        b = int(20  * (1 - t))
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
    fg_combined = folium.FeatureGroup(name="Road flows — est. AADT", show=True)
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
            _b = _link_flow_biz.get((_u, _v), 0) + _link_flow_biz.get((_v, _u), 0)
            _s = (_link_flow_school.get((_u, _v), 0) + _link_flow_school.get((_v, _u), 0)
                  if _has_school_layer else 0)
            _tot = _r + _b + _s
            if _tot > 0:
                _tip += f"<br>&nbsp;&nbsp;residential: {_r:,.0f} ({100*_r/_tot:.0f}%)"
                _tip += f"<br>&nbsp;&nbsp;business: {_b:,.0f} ({100*_b/_tot:.0f}%)"
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
            _link_flow_res, "Road flows — residential (pop→pop)", _color_res,
            lambda name, flow, length: (
                f"{name or 'link'}<br>residential AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        _add_flow_fg(
            _link_flow_biz, "Road flows — business-adjacent (pb+bb)", _color_biz,
            lambda name, flow, length: (
                f"{name or 'link'}<br>business AADT: {flow:,.0f}<br>length: {length:.0f}m"
            ),
            show=False,
        )
        if _has_school_layer:
            _add_flow_fg(
                _link_flow_school, "Road flows — school (pop→school)", _color_school,
                lambda name, flow, length: (
                    f"{name or 'link'}<br>school AADT: {flow:,.0f}<br>length: {length:.0f}m"
                ),
                show=False,
            )
        _layer_str = "combined + res + biz" + (" + school" if _has_school_layer else "")
        print(f"Added flow layers from {_flows_path} ({_layer_str}, {len(_link_flow)} links)")
    else:
        print(f"Added flow layer from {_flows_path} ({len(_link_flow)} links)")
else:
    print(f"No flow data found at {_flows_path} — skipping flow layer")

# ── NI traffic count sites (from ODS, Irish Grid → WGS84) ────────────────────
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

    count_fg = folium.FeatureGroup(
        name=f"NI traffic count sites — 2023 ODS ({len(_ods_sites)})", show=True)
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
    print(f"Added NI count sites layer ({len(_ods_sites)} sites, "
          f"{sum(1 for s in _ods_sites if s['id'] in _CALIBRATION_SITES)} calibration)")
else:
    print(f"ODS file not found at {_ODS_PATH} — skipping count sites layer")

# ── External zone markers ─────────────────────────────────────────────────────
if _census_zones is not None:
    _ext_zones = _census_zones["external_nodes"]
    _max_ext_tot = max(
        (_ext_node_trips.get(_ez["id"], {}).get("trips_through", 0)
         + _ext_node_trips.get(_ez["id"], {}).get("trips_internal", 0))
        for _ez in _ext_zones
    ) if _ext_node_trips else 1.0
    ext_fg = folium.FeatureGroup(name=f"External zones ({len(_ext_zones)})", show=False)
    for _ez in _ext_zones:
        _nid    = _ez["id"]
        _ez_pop = node_population.get(_nid, 0)
        _ez_biz = node_business_demand.get(_nid, 0)
        _ez_sch = node_school_demand.get(_nid, 0)
        _trips  = _ext_node_trips.get(_nid, {})
        _t_int  = _trips.get("trips_internal", 0)
        _t_thr  = _trips.get("trips_through",  0)
        _t_tot  = _t_int + _t_thr
        _radius = max(4, 12 * (_t_tot / max(_max_ext_tot, 1)) ** 0.4) if _t_tot > 0 else 4
        _tip = (
            f"<b>{_nid}</b> [{_ez['level']}]<br>"
            f"population: {_ez_pop:,.0f}<br>"
            f"business demand: {_ez_biz:,.1f}<br>"
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
    print(f"Added external zones layer ({len(_ext_zones)} nodes"
          + (f", {len(_ext_node_trips)} with trip data" if _ext_node_trips else ", no trip data")
          + ")")

folium.LayerControl(collapsed=False).add_to(m)

out_path = f"{OUT_DIR}/newtownards_map.html"
m.save(out_path)
total_node_pop = sum(node_population.values())
print(f"\nSaved: {out_path}")
print(f"  {len(dz_plot)} Data Zones · est. pop range {pop_min:.0f}–{pop_max:.0f}")
print(f"  {len(node_population)} nodes with population assigned · total {total_node_pop:.0f}")
