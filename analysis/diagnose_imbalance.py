#!/usr/bin/env python3
"""Generation vs consumption imbalance per node × component (fit-diagnosis tool).

For each gravity component and each node i this compares two projections of the
tuned OD matrix T^c (per-pair flows from model.constrained_od_flows, scaled by the
tuned K_c and the AADT weight W_c → daily veh):

    gen_i = Σ_j T^c_ij   (row sum) — trips GENERATED / exported at i
    con_i = Σ_j T^c_ji   (col sum) — trips CONSUMED  / imported at i
    imbalance_i = con_i − gen_i   (signed daily veh);   ratio = con_i / gen_i

The model is generation- (production-) constrained, so gen_i is a FIXED input:
    • internal node  →  gen_i = producer weight × K_c   (its whole production)
    • external node  →  gen_i = the EXPORTED production only (the intra-zonal
      self-term keeps the rest at home and never materialises it as flow).
con_i is EMERGENT.  Σ_i gen = Σ_i con per component (every materialised trip has one
origin and one destination), so the imbalance is a pure spatial REDISTRIBUTION:
    imbalance < 0  → net SOURCE  (more leaves than returns — "leave, fewer come back")
    imbalance > 0  → net SINK    (more arrives than leaves)
The intra-zonal self-term is denominator-only (never a materialised trip), so it does
NOT appear in gen or con — external self-retention can't masquerade as an imbalance.

Each two-leg component is a symmetric producer↔attractor round-trip, so under uniform
accessibility gen = con at every node — the imbalance is therefore PURELY the
accessibility-driven balancing distortion (accessible nodes over-served, remote nodes
under-served).  The headline per component is the

    relocatable fraction = ½ Σ_i |imbalance_i| / Σ_i gen_i

— the share of that component's flow a doubly (production+consumption)-constrained model
would move.  Small ⇒ the single constraint is already self-consistent (consumption
constraint buys little); large ⇒ consumption constraint would matter.  Reported overall
and split internal / external (external carries extra imbalance from world-edge truncation,
a different effect from accessibility — read it separately).

Read-only: no re-tune, no OSRM, no paths rebuild.  Reads the CURRENT tuned_params.json,
so run it after a re-tune for the definitive picture (a stale param file still shows the
imbalance STRUCTURE, just not the final magnitudes).

Outputs (in reports/): per-component source/sink tables (stdout), imbalance.csv,
imbalance_scatter.png (gen-vs-con scatter), and imbalance_map.html — an interactive
folium map with one marker per node (colour = con/gen, blue source ↔ red sink; size =
generation) and a toggleable layer per component.  Open the HTML locally (it needs CDN
map tiles, so it won't render inside a sandboxed artifact).

Usage:
  python3 analysis/diagnose_imbalance.py                 # tables + CSV + scatter + map (all 6 components)
  python3 analysis/diagnose_imbalance.py --top 25        # longer source/sink tables
  python3 analysis/diagnose_imbalance.py --component commute
  python3 analysis/diagnose_imbalance.py --no-plot --no-map   # tables + CSV only
"""
import os, sys, json, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "simulation"))
from model import (PATHS_CACHE, WEIGHTS_FILE, TUNED_PARAMS, SCHOOL_LEVELS,
                   constrained_od_flows, load_self_terms, aadt_weights,
                   load_generation_rates, compute_generation_scales,
                   assert_paths_cache_fresh, willingness_keys, willingness_from_flat)

CENSUS_ZONES = "data/census_zones.json"
REPORTS_DIR  = "reports"

_pnid = lambda k: (int(k) if str(k).lstrip("-").isdigit() else k)


def load_model_flows():
    """Replicate build_assignment's forward pass and return the per-pair daily-veh flow
    (K_c·W_c·t_c) for each of the six components, plus node ids and od src/dst."""
    with open(WEIGHTS_FILE) as f:
        weights = json.load(f)
    with open(TUNED_PARAMS) as f:
        tp = json.load(f)
    if not all(k in tp for k in willingness_keys()):
        raise SystemExit("tuned_params.json lacks the double-exp willingness kernels — "
                         "run reset_gravity_params.py then tune_assignment.py")
    if not all(k in tp for k in ("K_res", "K_commute", "K_retail")):
        raise SystemExit("tuned_params.json lacks the multi-component K's — re-tune first")
    willingness = willingness_from_flat(tp)
    K_res, K_commute, K_retail = tp["K_res"], tp["K_commute"], tp["K_retail"]
    K_school = {lvl: tp.get(f"K_{lvl}", 0.0) for lvl in SCHOOL_LEVELS}

    _sf = lambda key: {tuple(int(x) for x in k.split(",")): v for k, v in tp.get(key, {}).items()}
    W_res, W_commute, W_retail, W_school = aadt_weights(
        _sf("slot_fracs_res"), _sf("slot_fracs_commute"), _sf("slot_fracs_retail"),
        {lvl: _sf(f"slot_fracs_school_{lvl}") for lvl in SCHOOL_LEVELS})

    cache = np.load(PATHS_CACHE, allow_pickle=True)
    assert_paths_cache_fresh(cache)
    node_ids = cache["node_ids"]
    od_src, od_dst = cache["od_src"], cache["od_dst"]
    od_dist = cache["od_dist"].astype(np.float64)
    N_nodes = len(node_ids)

    nw = lambda layer, sub=None: (
        {_pnid(k): v for k, v in weights.get(layer, {}).items()} if sub is None
        else {_pnid(k): v for k, v in weights.get(layer, {}).items()})
    arr = lambda d: np.array([d.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
    w_pop        = arr(nw("node_population"))
    w_commute_at = arr(nw("node_commute_attractor"))
    w_retail     = arr(nw("node_retail_spaces"))
    w_commute_pr = arr(nw("node_commute_producers"))
    w_sch  = {lvl: arr(nw(f"node_school_demand_{lvl}"))    for lvl in SCHOOL_LEVELS}
    w_schp = {lvl: arr(nw(f"node_school_producers_{lvl}")) for lvl in SCHOOL_LEVELS}

    gen_rates = load_generation_rates()
    gen_scale = (compute_generation_scales(weights, gen_rates) if gen_rates is not None else None)
    self_src, self_dist, self_w = load_self_terms(list(node_ids))
    active = [lvl for lvl in SCHOOL_LEVELS if K_school.get(lvl, 0.0) > 0 and w_sch[lvl].sum() > 0]

    t_res, t_com, t_ret, t_sch = constrained_od_flows(
        od_src, od_dst, od_dist, N_nodes,
        w_pop, w_commute_at, w_retail, willingness,
        with_school=len(active) > 0,
        w_school_levels=w_sch, w_school_prod_levels=w_schp,
        self_src=self_src, self_dist=self_dist, self_w=self_w,
        w_commute_prod=w_commute_pr, gen_scale=gen_scale)

    comps = {
        "res":     t_res * (K_res     * W_res),
        "commute": t_com * (K_commute * W_commute),
        "retail":  t_ret * (K_retail  * W_retail),
    }
    for lvl in SCHOOL_LEVELS:
        comps[f"school_{lvl}"] = t_sch[lvl] * (K_school.get(lvl, 0.0) * W_school.get(lvl, 0.0))
    return comps, node_ids, od_src, od_dst


def load_ext_meta():
    """{census_code: (level, lat, lon, population)} for external nodes (labelling only)."""
    if not os.path.exists(CENSUS_ZONES):
        return {}
    with open(CENSUS_ZONES) as f:
        zones = json.load(f)
    return {e["id"]: (e.get("level", "?"), e.get("centroid_lat"), e.get("centroid_lon"),
                      e.get("population")) for e in zones.get("external_nodes", [])}


def node_label(nid, is_ext, ext_meta):
    if is_ext:
        lvl = ext_meta.get(nid, ("?",))[0]
        return f"{nid} [{lvl}]"
    return str(nid)


def report_component(name, flow, node_ids, od_src, od_dst, is_ext, ext_meta, top):
    N = len(node_ids)
    gen = np.bincount(od_src, weights=flow, minlength=N)
    con = np.bincount(od_dst, weights=flow, minlength=N)
    imb = con - gen
    tot = gen.sum()

    print(f"\n{'='*78}\n{name.upper()}")
    if tot <= 0:
        print("  (no flow — component inactive)")
        return None
    # conservation sanity: Σgen must equal Σcon (identical materialised-pair set)
    cons_err = abs(gen.sum() - con.sum()) / tot
    print(f"  Σgen = {gen.sum():,.0f}   Σcon = {con.sum():,.0f}   (conservation resid {cons_err:.2e})")

    def reloc(mask):
        g = gen[mask].sum()
        return (0.5 * np.abs(imb[mask]).sum() / g) if g > 0 else float("nan")
    all_m = (gen > 0) | (con > 0)
    int_m = all_m & ~is_ext
    ext_m = all_m & is_ext
    print(f"  relocatable fraction  overall {reloc(all_m):6.1%}   "
          f"internal {reloc(int_m):6.1%}   external {reloc(ext_m):6.1%}")
    print(f"    (= ½Σ|con−gen|/Σgen — share of flow a consumption-constraint would move)")

    order = np.argsort(imb)                              # most negative (sources) first
    def show(idxs, header):
        print(f"  {header}")
        print(f"    {'node':<22} {'type':<4} {'gen':>10} {'con':>10} {'imbalance':>11} {'con/gen':>8}")
        for i in idxs:
            if gen[i] <= 0 and con[i] <= 0:
                continue
            r = (con[i] / gen[i]) if gen[i] > 0 else float("inf")
            typ = "ext" if is_ext[i] else "int"
            print(f"    {node_label(node_ids[i], is_ext[i], ext_meta):<22} {typ:<4} "
                  f"{gen[i]:>10,.0f} {con[i]:>10,.0f} {imb[i]:>11,.0f} {r:>8.2f}")
    show(order[:top],       f"top {top} net SOURCES (con < gen — trips leave, fewer return):")
    show(order[::-1][:top], f"top {top} net SINKS   (con > gen — more arrive than leave):")
    return gen, con, imb


def make_plot(results, node_ids, is_ext):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"\n[plot] matplotlib unavailable ({e}) — skipping scatter")
        return
    names = [n for n in results if results[n] is not None]
    ncol = 3
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.5 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, name in zip(axes, names):
        gen, con, imb = results[name]
        m = (gen > 0) & (con > 0)
        for mask, col, lab in ((m & ~is_ext, "#2166ac", "internal"),
                               (m & is_ext,  "#b2182b", "external")):
            ax.scatter(gen[mask], con[mask], s=8, alpha=0.5, c=col, label=lab, edgecolors="none")
        lo = min(gen[m].min(), con[m].min()) if m.any() else 1
        hi = max(gen[m].max(), con[m].max()) if m.any() else 10
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(name); ax.set_xlabel("generation (daily veh)"); ax.set_ylabel("consumption")
        ax.legend(fontsize=7, loc="upper left")
    for ax in axes[len(names):]:
        ax.axis("off")
    fig.suptitle("Generation vs consumption per node — off-diagonal = imbalance", y=1.0)
    fig.tight_layout()
    os.makedirs(REPORTS_DIR, exist_ok=True)
    out = os.path.join(REPORTS_DIR, "imbalance_scatter.png")
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"\n[plot] wrote {out}")


def write_csv(results, node_ids, is_ext, ext_meta):
    os.makedirs(REPORTS_DIR, exist_ok=True)
    out = os.path.join(REPORTS_DIR, "imbalance.csv")
    with open(out, "w") as f:
        f.write("component,node_id,type,gen,con,imbalance,ratio\n")
        for name, res in results.items():
            if res is None:
                continue
            gen, con, imb = res
            for i in range(len(node_ids)):
                if gen[i] <= 0 and con[i] <= 0:
                    continue
                r = (con[i] / gen[i]) if gen[i] > 0 else ""
                typ = "ext" if is_ext[i] else "int"
                f.write(f"{name},{node_ids[i]},{typ},{gen[i]:.3f},{con[i]:.3f},"
                        f"{imb[i]:.3f},{r if r == '' else f'{r:.4f}'}\n")
    print(f"[csv]  wrote {out}")


def _internal_coords():
    """{int node_id: (lat, lon)} for internal nodes, from the reduced routing graph (ITM→WGS84)."""
    import osmnx as ox, pyproj
    from model import ROUTING_GRAPH
    from demographics_config import PROJECTED_CRS
    G = ox.load_graphml(ROUTING_GRAPH)
    to_wgs = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)
    coords = {}
    for n, dat in G.nodes(data=True):
        try:
            x, y = float(dat["x"]), float(dat["y"])
        except (KeyError, ValueError, TypeError):
            continue
        lon, lat = to_wgs.transform(x, y)
        coords[_pnid(n)] = (lat, lon)
    return coords


def make_map(results, node_ids, is_ext, ext_meta):
    """Folium HTML: one CircleMarker per node, colour = con/gen (blue source ↔ red sink),
    size = generation, per-component toggleable layers.  Open locally (needs CDN tiles)."""
    try:
        import folium
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
    except Exception as e:
        print(f"\n[map] folium/matplotlib unavailable ({e}) — skipping map")
        return
    int_coords = _internal_coords()
    latlon = {}
    for i, nid in enumerate(node_ids):
        if is_ext[i]:
            meta = ext_meta.get(nid)
            if meta and meta[1] is not None and meta[2] is not None:
                latlon[i] = (meta[1], meta[2])
        else:
            c = int_coords.get(_pnid(nid))
            if c is not None:
                latlon[i] = c
    if not latlon:
        print("\n[map] no node coordinates found — skipping map")
        return

    cmap = cm.get_cmap("coolwarm")

    def color_for(gen, con):
        r = 4.0 if gen <= 0 else (0.25 if con <= 0 else con / gen)   # con/gen; clamp to [1/4,4]
        return mcolors.to_hex(cmap((np.clip(np.log2(r), -2, 2) + 2) / 4))

    def radius_for(gen):
        return 2.0 + 2.3 * np.log10(max(gen, 0.0) + 1.0)

    lats = [ll[0] for ll in latlon.values()]
    lons = [ll[1] for ll in latlon.values()]
    m = folium.Map(location=[float(np.mean(lats)), float(np.mean(lons))],
                   zoom_start=8, tiles="CartoDB positron")
    names = [n for n in results if results[n] is not None]
    for name in names:
        gen, con, imb = results[name]
        fg = folium.FeatureGroup(name=name, show=(name == "commute"))
        for i, (lat, lon) in latlon.items():
            g, c, d = gen[i], con[i], imb[i]
            if g <= 0 and c <= 0:
                continue
            ratio = (c / g) if g > 0 else float("inf")
            typ = "ext" if is_ext[i] else "int"
            lvl = ext_meta.get(node_ids[i], ("core",))[0] if is_ext[i] else "core"
            popup = (f"<b>{node_ids[i]}</b> ({typ}, {lvl})<br>"
                     f"gen {g:,.0f} &rarr; con {c:,.0f}<br>"
                     f"imbalance {d:+,.0f} veh/day<br>con/gen {ratio:.2f}")
            folium.CircleMarker(
                location=[lat, lon], radius=radius_for(g),
                weight=0, fill=True, fill_color=color_for(g, c), fill_opacity=0.75,
                tooltip=f"{name}: con/gen {ratio:.2f} ({d:+,.0f})",
                popup=folium.Popup(popup, max_width=260)).add_to(fg)
        fg.add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)
    m.fit_bounds([[min(lats), min(lons)], [max(lats), max(lons)]])

    legend = ('<div style="position:fixed;bottom:24px;left:24px;z-index:9999;background:white;'
              'padding:10px 12px;border:1px solid #999;border-radius:5px;font:12px sans-serif;'
              'max-width:235px"><b>Generation–consumption imbalance</b><br>'
              'colour = con/gen &nbsp;'
              '<span style="color:#3b4cc0">&#9679;</span>&nbsp;source (con&lt;gen) &nbsp;'
              '<span style="color:#dddddd">&#9679;</span>&nbsp;balanced &nbsp;'
              '<span style="color:#b40426">&#9679;</span>&nbsp;sink (con&gt;gen)<br>'
              'size = generation (daily veh)<br>'
              '<i>toggle components top-right; default = commute</i></div>')
    m.get_root().html.add_child(folium.Element(legend))

    os.makedirs(REPORTS_DIR, exist_ok=True)
    out = os.path.join(REPORTS_DIR, "imbalance_map.html")
    m.save(out)
    print(f"\n[map]  wrote {out}  ({len(latlon)} located nodes, {len(names)} component layers)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=15, help="rows per source/sink table (default 15)")
    ap.add_argument("--component", help="restrict to one component (res/commute/retail/school_<lvl>)")
    ap.add_argument("--no-plot", action="store_true", help="skip the gen-vs-con scatter PNG")
    ap.add_argument("--no-map", action="store_true", help="skip the interactive folium map")
    ap.add_argument("--no-csv", action="store_true")
    args = ap.parse_args()

    comps, node_ids, od_src, od_dst = load_model_flows()
    ext_meta = load_ext_meta()
    is_ext = np.array([not str(nid).lstrip("-").isdigit() for nid in node_ids])

    print(f"Loaded {len(node_ids):,} nodes  ({is_ext.sum():,} external / {(~is_ext).sum():,} internal), "
          f"{len(od_src):,} OD pairs.")
    print("gen = exported (production-constrained input); con = imported (emergent). "
          "imbalance = con − gen.")

    names = [args.component] if args.component else list(comps)
    results = {}
    for name in names:
        if name not in comps:
            raise SystemExit(f"unknown component '{name}' — choose from {list(comps)}")
        results[name] = report_component(name, comps[name], node_ids, od_src, od_dst,
                                         is_ext, ext_meta, args.top)

    if not args.no_csv:
        write_csv(results, node_ids, is_ext, ext_meta)
    if not args.no_plot:
        make_plot(results, node_ids, is_ext)
    if not args.no_map:
        make_map(results, node_ids, is_ext, ext_meta)


if __name__ == "__main__":
    main()
