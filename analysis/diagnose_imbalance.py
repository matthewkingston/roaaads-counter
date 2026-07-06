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

Intra-zonal self-flow (default ON): an external zone collapsed to a centroid keeps a
fraction s_i = a_i·E[f_intra]/D_i of its production intra-zonally via the denominator-only
self-term, which is never materialised as a trip — so its exported gen/con show only the
inter-zonal slice and self-contained coarse zones read with wild ratios.  By default the
implied self-flow (K·W·scale·p_i·s_i) is materialised as an i→i trip added to both that
zone's gen and con, so ratios/relocatable are measured against FULL production.  This does
NOT change imbalance (con−gen adds equally to both) — it only rescales ratios toward 1 for
self-contained zones; genuine inter-zonal imbalances (e.g. Belfast) are preserved, and
internal nodes (s_i=0) are untouched.  --exported-only restores the old inter-zonal-slice view.

--sides splits each two-leg component into its producer-role and attractor-role side
(e.g. commute → commute·worker {residents sent out vs returning, ref = commute_producers}
and commute·job {workers arriving vs leaving, ref = commute_attractor}), because the
combined gen/con SUMS both roles at a node and can hide a cancelling imbalance
(worker-side + job-side imbalance = the combined con−gen).  res is single-leg (pop↔pop) so
it stays one signal → 11 side-signals.  Each side is isolated by zeroing every other leg's
gen_scale (denominators are scale-free), so it costs ~11 model passes and writes
reports/imbalance_sides_{csv,scatter.png,map.html}.

Usage:
  python3 analysis/diagnose_imbalance.py                 # 6-component view: tables + CSV + scatter + map
  python3 analysis/diagnose_imbalance.py --sides         # 11 producer/attractor side-signals
  python3 analysis/diagnose_imbalance.py --exported-only # inter-zonal slice only (old external readings)
  python3 analysis/diagnose_imbalance.py --top 25        # longer source/sink tables
  python3 analysis/diagnose_imbalance.py --component commute          # one component
  python3 analysis/diagnose_imbalance.py --sides --component commute  # commute·worker + commute·job
  python3 analysis/diagnose_imbalance.py --sides --doubly             # SANITY CHECK: apply the doubly_constrained set (ratios → ~1)
  python3 analysis/diagnose_imbalance.py --no-plot --no-map   # tables + CSV only
"""
import os, sys, json, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "simulation"))
from model import (PATHS_CACHE, WEIGHTS_FILE, TUNED_PARAMS, SCHOOL_LEVELS,
                   constrained_od_flows, load_self_terms, aadt_weights,
                   load_generation_rates, compute_generation_scales,
                   assert_paths_cache_fresh, willingness_keys, willingness_from_flat,
                   _modesub_kernel)

CENSUS_ZONES = "data/census_zones.json"
REPORTS_DIR  = "reports"

_pnid = lambda k: (int(k) if str(k).lstrip("-").isdigit() else k)


def load_context():
    """Load weights, tuned params, paths cache and the generation scales once; return a
    context dict reused by both the combined (6-component) and per-leg (--sides) views."""
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
    doubly_set = set(tp.get("doubly_constrained") or [])     # honoured only under --doubly
    furness_sweeps = tp.get("furness_max_sweeps")
    K = {"res": tp["K_res"], "commute": tp["K_commute"], "retail": tp["K_retail"]}
    for lvl in SCHOOL_LEVELS:
        K[f"school_{lvl}"] = tp.get(f"K_{lvl}", 0.0)

    _sf = lambda key: {tuple(int(x) for x in k.split(",")): v for k, v in tp.get(key, {}).items()}
    W_res, W_commute, W_retail, W_school = aadt_weights(
        _sf("slot_fracs_res"), _sf("slot_fracs_commute"), _sf("slot_fracs_retail"),
        {lvl: _sf(f"slot_fracs_school_{lvl}") for lvl in SCHOOL_LEVELS})
    W = {"res": W_res, "commute": W_commute, "retail": W_retail}
    for lvl in SCHOOL_LEVELS:
        W[f"school_{lvl}"] = W_school.get(lvl, 0.0)

    cache = np.load(PATHS_CACHE, allow_pickle=True)
    assert_paths_cache_fresh(cache)
    node_ids = cache["node_ids"]
    od_src, od_dst = cache["od_src"], cache["od_dst"]
    od_dist = cache["od_dist"].astype(np.float64)
    N = len(node_ids)

    nw = lambda layer: {_pnid(k): v for k, v in weights.get(layer, {}).items()}
    arr = lambda d: np.array([d.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
    w_sch  = {lvl: arr(nw(f"node_school_demand_{lvl}"))    for lvl in SCHOOL_LEVELS}
    gen_rates = load_generation_rates()
    gen_scale = (compute_generation_scales(weights, gen_rates) if gen_rates is not None else None)
    self_terms = load_self_terms(list(node_ids))            # {component: (src, dist, w)} or None
    active = [lvl for lvl in SCHOOL_LEVELS if K[f"school_{lvl}"] > 0 and w_sch[lvl].sum() > 0]

    return dict(node_ids=node_ids, od_src=od_src, od_dst=od_dst, od_dist=od_dist, N=N,
                w_pop=arr(nw("node_population")), w_commute_at=arr(nw("node_commute_attractor")),
                w_retail=arr(nw("node_retail_spaces")), w_commute_pr=arr(nw("node_commute_producers")),
                w_sch=w_sch, w_schp={lvl: arr(nw(f"node_school_producers_{lvl}")) for lvl in SCHOOL_LEVELS},
                willingness=willingness, self_terms=self_terms,
                active=active, gen_scale=gen_scale, K=K, W=W,
                doubly=False, doubly_set=doubly_set, furness_sweeps=furness_sweeps)


def _run(ctx, gen_scale, with_school=None):
    """One constrained_od_flows pass with a given (possibly leg-masked) gen_scale."""
    if with_school is None:
        with_school = len(ctx["active"]) > 0
    extra = (dict(doubly_constrained=ctx["doubly_set"], furness_max_sweeps=ctx.get("furness_sweeps"))
             if ctx.get("doubly") else {})
    return constrained_od_flows(
        ctx["od_src"], ctx["od_dst"], ctx["od_dist"], ctx["N"],
        ctx["w_pop"], ctx["w_commute_at"], ctx["w_retail"], ctx["willingness"],
        with_school=with_school,
        w_school_levels=ctx["w_sch"], w_school_prod_levels=ctx["w_schp"],
        self_terms=ctx["self_terms"],
        w_commute_prod=ctx["w_commute_pr"], gen_scale=gen_scale, **extra)


def _select(comp, t_res, t_com, t_ret, t_sch):
    if comp == "res":     return t_res
    if comp == "commute": return t_com
    if comp == "retail":  return t_ret
    return t_sch[comp.split("_", 1)[1]]                       # school_<lvl>


ALL_COMPONENTS = ["res", "commute", "retail"] + [f"school_{lvl}" for lvl in SCHOOL_LEVELS]


def combined_flows(ctx):
    """Per-component daily-veh per-pair flow (K_c·W_c·t_c) — the 6-component view."""
    t = _run(ctx, ctx["gen_scale"])
    return {c: _select(c, *t) * (ctx["K"][c] * ctx["W"][c]) for c in ALL_COMPONENTS}


# gen_scale leg key → the component it belongs to
LEG_COMPONENT = {"res": "res", "com_out": "commute", "com_ret": "commute",
                 "ret_out": "retail", "ret_ret": "retail"}
for _lvl in SCHOOL_LEVELS:
    LEG_COMPONENT[f"sch_{_lvl}_out"] = f"school_{_lvl}"
    LEG_COMPONENT[f"sch_{_lvl}_ret"] = f"school_{_lvl}"

# component → its generation legs (inverse of LEG_COMPONENT)
COMPONENT_LEGS = {c: [k for k, v in LEG_COMPONENT.items() if v == c] for c in ALL_COMPONENTS}


def _leg_prod_attr(ctx):
    """{leg_key: (producer_array, attractor_array)} for every active generation leg."""
    pa = {"res": (ctx["w_pop"], ctx["w_pop"]),
          "com_out": (ctx["w_commute_pr"], ctx["w_commute_at"]),
          "com_ret": (ctx["w_commute_at"], ctx["w_commute_pr"]),
          "ret_out": (ctx["w_pop"], ctx["w_retail"]),
          "ret_ret": (ctx["w_retail"], ctx["w_pop"])}
    for lvl in SCHOOL_LEVELS:
        pa[f"sch_{lvl}_out"] = (ctx["w_schp"][lvl], ctx["w_sch"][lvl])
        pa[f"sch_{lvl}_ret"] = (ctx["w_sch"][lvl], ctx["w_schp"][lvl])
    return pa


def self_flows(ctx):
    """Per-leg implied intra-zonal self-flow — the withheld production the denominator-only
    self-term (a_i·E[f_intra]) accounts for but never materialises as a trip.  For a leg with
    producer p, attractor a and kernel f:  self_i = K·W·gen_scale·p_i·D_self_i/D_total_i,
    D_self_i = a_i·E[f(sampled intra times)],  D_total_i = Σ_j a_j f(d_ij) + D_self_i.
    Nonzero only for external zones (D_self=0 elsewhere).  Materialising it as an i→i trip
    (added to both the leg's gen and con) makes external zones read against their FULL
    production; imbalance (con−gen) is unchanged — only ratios/relocatable move toward 1.
    Returns {leg_key: self_flow_array} (empty if generation pinning / self-term absent)."""
    gs = ctx["gen_scale"]
    src, dst, dist, Nn = ctx["od_src"], ctx["od_dst"], ctx["od_dist"], ctx["N"]
    self_terms = ctx["self_terms"]
    if not gs or not self_terms:
        return {}
    pa = _leg_prod_attr(ctx)
    kern = {}
    out = {}
    for key in [k for k in gs if gs[k] != 0.0 and k in LEG_COMPONENT]:
        comp = LEG_COMPONENT[key]
        p, a = pa[key]
        if comp not in kern:
            kern[comp] = _modesub_kernel(dist, ctx["willingness"][comp], comp)
        f = kern[comp]
        D_inter = np.bincount(src, weights=a[dst] * f, minlength=Nn)
        st = self_terms.get(comp)                                  # component's mass-weighted intra-zonal histogram
        if st is not None:
            s_src, s_dist, s_w = st
            F_self = _modesub_kernel(s_dist, ctx["willingness"][comp], comp)
            D_self = np.bincount(s_src, weights=a[s_src] * F_self * s_w, minlength=Nn)
        else:
            D_self = np.zeros(Nn)
        D_tot = D_inter + D_self
        scale = gs[key] * ctx["K"][comp] * ctx["W"][comp]
        out[key] = np.where(D_tot > 0, scale * p * D_self / D_tot, 0.0)
    return out


def leg_rowcol(ctx):
    """{leg_key: (rowsum, colsum)} of the K·W-scaled flow for each generation leg, isolated
    by zeroing every OTHER leg's gen_scale (denominators are scale-free, so one masked call =
    one exact leg).  Row sum = the leg's constrained producer generation; col sum = its
    emergent attractor consumption."""
    gs_full = ctx["gen_scale"]
    if not gs_full:
        raise SystemExit("--sides needs generation pinning (run analysis/derive_generation_rates.py)")
    keys = [k for k in gs_full if gs_full[k] != 0.0 and k in LEG_COMPONENT]
    legs = {}
    for n, leg_key in enumerate(keys, 1):
        comp = LEG_COMPONENT[leg_key]
        print(f"  [leg {n}/{len(keys)}] isolating {leg_key} ({comp}) …", flush=True)
        mask = {k: 0.0 for k in gs_full}
        mask[leg_key] = gs_full[leg_key]
        t = _run(ctx, mask, with_school=comp.startswith("school_"))
        flow = _select(comp, *t) * (ctx["K"][comp] * ctx["W"][comp])
        legs[leg_key] = (np.bincount(ctx["od_src"], weights=flow, minlength=ctx["N"]),
                         np.bincount(ctx["od_dst"], weights=flow, minlength=ctx["N"]))
    return legs


def build_sides(legs, sf=None):
    """{side_label: (out, in)} — out = constrained generation of that role (row sum of its
    out-leg), in = emergent return consumption (col sum of the paired return-leg).  Each
    two-leg component splits into a producer-role and an attractor-role side; res is
    single-leg (pop↔pop) so it stays one side.  sf (self_flows dict) adds each leg's implied
    intra-zonal self-flow to that role's out and in (readability materialisation)."""
    sf = sf or {}
    sides = {}
    def add(label, out_leg, in_leg):
        if out_leg in legs and in_leg in legs:
            out = legs[out_leg][0] + sf.get(out_leg, 0.0)
            inn = legs[in_leg][1] + sf.get(in_leg, 0.0)
            sides[label] = (out, inn)
    add("res", "res", "res")
    add("commute·worker", "com_out", "com_ret")   # residents: sent out (AM) vs returning (PM), ref = commute_producers
    add("commute·job",    "com_ret", "com_out")   # jobs: workers leaving (PM) vs arriving (AM), ref = commute_attractor
    add("retail·population", "ret_out", "ret_ret")
    add("retail·retail",     "ret_ret", "ret_out")
    for lvl in SCHOOL_LEVELS:
        add(f"school_{lvl}·student", f"sch_{lvl}_out", f"sch_{lvl}_ret")
        add(f"school_{lvl}·school",  f"sch_{lvl}_ret", f"sch_{lvl}_out")
    return sides


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


def report_component(name, flow, node_ids, od_src, od_dst, is_ext, ext_meta, top, self_add=None):
    N = len(node_ids)
    gen = np.bincount(od_src, weights=flow, minlength=N)
    con = np.bincount(od_dst, weights=flow, minlength=N)
    if self_add is not None:                       # materialise intra-zonal self-flow (both sides)
        gen = gen + self_add
        con = con + self_add
    return report_from_gencon(name, gen, con, node_ids, is_ext, ext_meta, top)


def report_from_gencon(name, gen, con, node_ids, is_ext, ext_meta, top, balanced=True):
    """Report imbalance for one series. gen = generated/exported (or a side's constrained
    OUT); con = consumed/imported (or a side's emergent IN). imbalance = con − gen.
    balanced=True (combined component / res): gen and con are the row/col of one flow, so
    Σgen ≡ Σcon (a conservation check).  balanced=False (a side): out and in are DIFFERENT
    legs, so Σout ≠ Σin in general — the small gap is the two legs' exported-total
    difference (self-term + finite-world truncation), not a violation."""
    imb = con - gen
    tot = gen.sum()

    print(f"\n{'='*78}\n{name.upper()}")
    if tot <= 0:
        print("  (no flow — inactive)")
        return None
    if balanced:
        cons_err = abs(gen.sum() - con.sum()) / tot
        print(f"  Σgen = {gen.sum():,.0f}   Σcon = {con.sum():,.0f}   (conservation resid {cons_err:.2e})")
    else:
        print(f"  Σout = {gen.sum():,.0f}   Σin = {con.sum():,.0f}   "
              f"(Σin/Σout = {con.sum()/tot:.3f}; the two legs' exported-total gap — self-term/truncation)")

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


def make_plot(results, node_ids, is_ext, out_name="imbalance_scatter.png"):
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
    fig.suptitle("out (generation) vs in (consumption) per node — off-diagonal = imbalance", y=1.0)
    fig.tight_layout()
    os.makedirs(REPORTS_DIR, exist_ok=True)
    out = os.path.join(REPORTS_DIR, out_name)
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"\n[plot] wrote {out}")


def write_csv(results, node_ids, is_ext, ext_meta, out_name="imbalance.csv"):
    os.makedirs(REPORTS_DIR, exist_ok=True)
    out = os.path.join(REPORTS_DIR, out_name)
    with open(out, "w") as f:
        f.write("series,node_id,type,gen,con,imbalance,ratio\n")
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


def make_map(results, node_ids, is_ext, ext_meta,
             out_name="imbalance_map.html", default_show="commute"):
    """Folium HTML: one CircleMarker per node, colour = in/out (blue source ↔ red sink),
    size = generation, per-series toggleable layers.  Open locally (needs CDN tiles)."""
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
        fg = folium.FeatureGroup(name=name, show=(name == default_show))
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
              f'<i>toggle layers top-right; default = {default_show}</i></div>')
    m.get_root().html.add_child(folium.Element(legend))

    os.makedirs(REPORTS_DIR, exist_ok=True)
    out = os.path.join(REPORTS_DIR, out_name)
    m.save(out)
    print(f"\n[map]  wrote {out}  ({len(latlon)} located nodes, {len(names)} layers)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=15, help="rows per source/sink table (default 15)")
    ap.add_argument("--component", help="restrict to one component (res/commute/retail/school_<lvl>); "
                    "with --sides matches both of that component's sides by prefix")
    ap.add_argument("--sides", action="store_true",
                    help="split each two-leg component into its producer-role and attractor-role "
                         "sides (e.g. commute·worker vs commute·job) — 11 side-signals; ~11 model "
                         "passes; writes reports/imbalance_sides_*")
    ap.add_argument("--exported-only", action="store_true",
                    help="do NOT materialise the intra-zonal self-flow — external gen/con show only "
                         "the exported (inter-zonal) slice (the old behaviour; ratios read wild for "
                         "self-contained coarse zones)")
    ap.add_argument("--doubly", action="store_true",
                    help="apply the doubly_constrained set from tuned_params.json (Furness) when "
                         "computing flows — a sanity check that flagged components' imbalance → ~0 "
                         "(con/gen ratio → 1). SLOW (cold-seeds each doubly leg, minutes). Default is "
                         "the singly diagnostic (the relocatable-fraction estimate of what doubly "
                         "WOULD move). Self-flow materialisation is skipped under --doubly.")
    ap.add_argument("--no-plot", action="store_true", help="skip the scatter PNG")
    ap.add_argument("--no-map", action="store_true", help="skip the interactive folium map")
    ap.add_argument("--no-csv", action="store_true")
    args = ap.parse_args()

    ctx = load_context()
    node_ids = ctx["node_ids"]
    ext_meta = load_ext_meta()
    is_ext = np.array([not str(nid).lstrip("-").isdigit() for nid in node_ids])
    print(f"Loaded {len(node_ids):,} nodes  ({is_ext.sum():,} external / {(~is_ext).sum():,} internal), "
          f"{len(ctx['od_src']):,} OD pairs.")

    if args.doubly:
        if not ctx["doubly_set"]:
            raise SystemExit("--doubly: tuned_params.json has no non-empty doubly_constrained list "
                             "(activate it in tuner_config.json + reset_gravity_params.py first).")
        ctx["doubly"] = True
        print(f"DOUBLY: applying Furness attraction constraint to {sorted(ctx['doubly_set'])} "
              f"(cold-seeded per leg — slow). Flagged components' imbalance should read ~0 (ratio ~1); "
              f"unflagged components (e.g. res) still read their singly imbalance.")

    materialize = not args.exported_only and not args.doubly    # doubly ⇒ gen=con already; skip the singly self-flow
    sf = self_flows(ctx) if materialize else {}       # per-leg implied intra-zonal self-flow
    if sf:
        print("Self-flow: intra-zonal self-term MATERIALISED into external gen & con "
              "(imbalance con−gen unchanged; ratios/relocatable now measured vs full production).")
    elif args.doubly:
        print("Self-flow: not materialised under --doubly (doubly-constrained gen=con on-network "
              "already; raw gen/con used).")
    else:
        print("Self-flow: EXPORTED-ONLY (self-term not materialised; external ratios read wild)."
              if args.exported_only else "Self-flow: none (no self-term/gen-pinning present).")

    results = {}
    if args.sides:
        print("SIDES: each two-leg component split into producer-role & attractor-role.\n"
              "  out = constrained generation of that role (≈ its producer/attractor weight × scale);\n"
              "  in  = emergent return consumption.  imbalance = in − out.")
        sides = build_sides(leg_rowcol(ctx), sf)
        for label, (out, inn) in sides.items():
            if args.component and not label.startswith(args.component):
                continue
            results[label] = report_from_gencon(label, out, inn, node_ids, is_ext, ext_meta,
                                                 args.top, balanced=False)
        suffix, default_show = "_sides", "commute·worker"
    else:
        print("gen = exported (production-constrained input); con = imported (emergent). "
              "imbalance = con − gen.")
        comps = combined_flows(ctx)
        comp_self = {c: sum((sf.get(k, 0.0) for k in COMPONENT_LEGS[c]), np.zeros(ctx["N"]))
                     for c in comps} if sf else {}
        names = [args.component] if args.component else list(comps)
        for name in names:
            if name not in comps:
                raise SystemExit(f"unknown component '{name}' — choose from {list(comps)}")
            results[name] = report_component(name, comps[name], node_ids, ctx["od_src"], ctx["od_dst"],
                                             is_ext, ext_meta, args.top,
                                             self_add=comp_self.get(name))
        suffix, default_show = "", "commute"

    if not args.no_csv:
        write_csv(results, node_ids, is_ext, ext_meta, f"imbalance{suffix}.csv")
    if not args.no_plot:
        make_plot(results, node_ids, is_ext, f"imbalance{suffix}_scatter.png")
    if not args.no_map:
        make_map(results, node_ids, is_ext, ext_meta, f"imbalance{suffix}_map.html", default_show)


if __name__ == "__main__":
    main()
