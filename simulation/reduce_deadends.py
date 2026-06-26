"""
Collapse "residential dead-end" regions in the consolidated routing graph.

A residential dead-end region R (with entrance node E, E not in R) satisfies:
  1. R connects to the rest of the network through exactly one node E (a cut vertex);
  2. R contains no boundary node and no school-demand node (both are *protected* and so
     can never be absorbed — this enforces criterion 2 and the zero-school rule
     structurally);
  3. max directed journey time E->n over n in R is below T_MAX (routing-cost seconds);
  4. total business demand in R is below BIZ_CAP (residential pop is unbounded);
  5. |R| >= 2 (single-node spurs are skipped — collapsing 1->1 saves nothing).

Each accepted region is replaced by one super-node S (summed pop/biz/school of R), joined
to E by a directed link pair (E->S, S->E) whose travel times are population-weighted
averages of the intra-region journey times.

Algorithm: every valid region has a single entrance = a single cut vertex, so it is a
protected-free connected component of H - a for some articulation point a (H = undirected
simple projection). We enumerate all such (entrance, region) candidates, keep the feasible
ones (constraints 2-5 + directed reachability both ways), and select the *maximal feasible*
regions. Because the candidates form a laminar family, maximal-feasible regions are pairwise
disjoint, and selecting maximal-feasible naturally descends into an oversized branch to find
the largest collapsible sub-pockets.

Outputs (opt-in; build_paths.py consumes them only when USE_REDUCED is enabled there):
  simulation/newtownards_reduced.graphml   — reduced routing graph
  simulation/node_weights_reduced.json     — pop/biz/school re-bucketed onto super-nodes
  simulation/deadend_map.json              — provenance: super-node -> absorbed nodes + times
  simulation/deadend_broken_obs.json       — observed/count links whose endpoints were eaten

Run AFTER build_demographics.py (needs pop/biz/school + boundary nodes) and BEFORE
build_paths.py.
"""

import argparse
import json
import os
import sys

import networkx as nx
import osmnx as ox
import pyproj

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from routing_config import HIGHWAY_COST_FACTOR
import model
from demographics_config import PROJECTED_CRS

CONS_GRAPH    = "simulation/newtownards_consolidated.graphml"
WEIGHTS_FILE  = "simulation/node_weights.json"
LINK_AADT     = "data/link_aadt.json"

OUT_GRAPH     = "simulation/newtownards_reduced.graphml"
OUT_WEIGHTS   = "simulation/node_weights_reduced.json"
OUT_MAP       = "simulation/deadend_map.json"
OUT_BROKEN    = "simulation/deadend_broken_obs.json"

T_MAX_DEFAULT   = 60.0    # max entrance->dead-end journey time (routing-cost seconds)
BIZ_CAP_DEFAULT = 100.0   # max business demand a region may absorb
MIN_REGION_SIZE = 2       # skip single-node spurs
SYNTH_HIGHWAY   = "deadend_collapsed"  # not in HIGHWAY_COST_FACTOR -> cost factor 1.0
SYNTH_SPEED_KPH = 30.0    # reference speed for synthetic edge length back-calculation


def hw_factor(highway):
    if isinstance(highway, list):
        highway = highway[0] if highway else "unclassified"
    return HIGHWAY_COST_FACTOR.get(highway, 1.0)


def build_cost_digraph(G):
    """Min routing-cost (travel_time * highway factor) per directed (u, v)."""
    D = nx.DiGraph()
    D.add_nodes_from(G.nodes())
    for u, v, edata in G.edges(data=True):
        cost = float(edata.get("travel_time", 1.0)) * hw_factor(edata.get("highway"))
        if D.has_edge(u, v):
            if cost < D[u][v]["cost"]:
                D[u][v]["cost"] = cost
        else:
            D.add_edge(u, v, cost=cost)
    return D


def main():
    ap = argparse.ArgumentParser(description="Collapse residential dead-ends.")
    ap.add_argument("--t-max", type=float, default=T_MAX_DEFAULT,
                    help=f"max entrance->dead-end journey time, cost-seconds (default {T_MAX_DEFAULT})")
    ap.add_argument("--biz-cap", type=float, default=BIZ_CAP_DEFAULT,
                    help=f"max business demand per region (default {BIZ_CAP_DEFAULT})")
    args = ap.parse_args()
    T_MAX, BIZ_CAP = args.t_max, args.biz_cap

    print(f"T_MAX={T_MAX}s  BIZ_CAP={BIZ_CAP}  MIN_REGION_SIZE={MIN_REGION_SIZE}")

    # ── Load graph + weights ──────────────────────────────────────────────────
    print("Loading graph …")
    G = ox.load_graphml(CONS_GRAPH)
    G = nx.relabel_nodes(G, {n: int(n) for n in G.nodes()})
    G = ox.routing.add_edge_speeds(G)
    G = ox.routing.add_edge_travel_times(G)
    print(f"  {G.number_of_nodes()} nodes  {G.number_of_edges()} edges")

    with open(WEIGHTS_FILE) as f:
        nw = json.load(f)
    pop = nw["node_population"]
    biz = nw["node_business_demand"]
    sch = nw["node_school_demand"]
    park = nw.get("node_parking_equiv", {})

    def w(d, n):
        return float(d.get(str(n), 0.0))

    boundary = set(int(x) for x in nw["boundary_node_ids_cons"])
    school_nodes = set(int(k) for k, v in sch.items()
                       if str(k).lstrip("-").isdigit() and float(v) > 0.0)
    # Official count-site nodes are hardcoded references in build_assignment.py — protect
    # them so a site can never be absorbed (none are at present; this is future-proofing).
    count_site_nodes = set()
    for site in model.COUNT_SITES:
        if site.get("node") is not None:
            count_site_nodes.add(int(site["node"]))
        for (u, v) in (site.get("links") or []):
            count_site_nodes.add(int(u)); count_site_nodes.add(int(v))
    gnodes = set(G.nodes())
    protected = (boundary | school_nodes | count_site_nodes) & gnodes
    print(f"  protected: {len(boundary & gnodes)} boundary + "
          f"{len(school_nodes & gnodes)} school + "
          f"{len(count_site_nodes & gnodes)} count-site = {len(protected)} nodes")

    # ── Undirected simple projection ──────────────────────────────────────────
    H = nx.Graph(G)

    # components with no protected node have no entrance — left intact, reported
    no_entrance = []
    for comp in nx.connected_components(H):
        if not (comp & protected):
            no_entrance.append(sorted(comp))
    n_no_entrance_nodes = sum(len(c) for c in no_entrance)
    print(f"  {len(no_entrance)} component(s) without a protected node "
          f"({n_no_entrance_nodes} nodes) left intact")

    # ── Enumerate candidate regions ───────────────────────────────────────────
    print("Enumerating candidate regions …")
    candidates = {}  # frozenset(region) -> entrance
    for a in sorted(nx.articulation_points(H)):
        incident = list(H.edges(a))
        H.remove_node(a)
        for comp in nx.connected_components(H):
            if comp & protected:
                continue
            region = frozenset(comp)
            candidates.setdefault(region, a)
        H.add_node(a)
        H.add_edges_from(incident)
    print(f"  {len(candidates)} candidate regions")

    # ── Feasibility ───────────────────────────────────────────────────────────
    D = build_cost_digraph(G)
    feasible = []   # (region frozenset, entrance, t_fwd, t_rev, maxT, region_pop, region_biz)
    skipped_unreachable = []

    for region, E in candidates.items():
        if len(region) < MIN_REGION_SIZE:
            continue
        region_biz = sum(w(biz, n) for n in region)
        if region_biz >= BIZ_CAP:
            continue
        # school is 0 by protection; assert defensively
        if any(n in school_nodes for n in region):
            continue

        nodes = region | {E}
        sub = D.subgraph(nodes)
        fwd = nx.single_source_dijkstra_path_length(sub, E, weight="cost")
        rev = nx.single_source_dijkstra_path_length(sub.reverse(copy=False), E, weight="cost")
        if not all(n in fwd and n in rev for n in region):
            skipped_unreachable.append((E, sorted(region)))
            continue
        maxT = max(fwd[n] for n in region)
        if maxT >= T_MAX:
            continue

        region_pop = sum(w(pop, n) for n in region)
        tot = region_pop if region_pop > 0 else 0.0
        if tot > 0:
            t_fwd = sum(w(pop, n) * fwd[n] for n in region) / tot
            t_rev = sum(w(pop, n) * rev[n] for n in region) / tot
        else:
            t_fwd = sum(fwd[n] for n in region) / len(region)
            t_rev = sum(rev[n] for n in region) / len(region)
        feasible.append((region, E, t_fwd, t_rev, maxT, region_pop, region_biz))

    if skipped_unreachable:
        print(f"  WARNING: {len(skipped_unreachable)} region(s) skipped — not all nodes "
              f"directed-reachable from entrance (one-way topology):")
        for E, r in skipped_unreachable[:10]:
            print(f"    entrance {E}: {r}")

    # ── Maximal-feasible selection (laminar -> disjoint) ──────────────────────
    feasible.sort(key=lambda t: len(t[0]), reverse=True)
    feasible_sets = [t[0] for t in feasible]
    selected = []
    for i, (region, E, t_fwd, t_rev, maxT, rpop, rbiz) in enumerate(feasible):
        if any(region < bigger for bigger in feasible_sets if bigger is not region
               and len(bigger) > len(region)):
            continue  # strictly contained in a larger feasible region
        selected.append((region, E, t_fwd, t_rev, maxT, rpop, rbiz))

    # sanity: selected regions are pairwise disjoint
    seen = set()
    for region, *_ in selected:
        assert not (region & seen), "selected regions overlap — laminar assumption broken"
        seen |= set(region)

    print(f"  {len(selected)} regions selected for collapse "
          f"({sum(len(r[0]) for r in selected)} nodes -> {len(selected)} super-nodes)")

    # ── Rewrite graph + weights ───────────────────────────────────────────────
    nw_out = {k: dict(v) if isinstance(v, dict) else list(v) if isinstance(v, list) else v
              for k, v in nw.items()}
    pop_o, biz_o, sch_o, park_o = (nw_out["node_population"], nw_out["node_business_demand"],
                                   nw_out["node_school_demand"], nw_out.get("node_parking_equiv", {}))
    internal = set(int(x) for x in nw_out.get("internal_node_ids", []))

    deadend_map = {}
    absorbed_all = set()
    utm_to_wgs = pyproj.Transformer.from_crs(PROJECTED_CRS, "EPSG:4326", always_xy=True)

    for region, E, t_fwd, t_rev, maxT, rpop, rbiz in sorted(selected, key=lambda t: min(t[0])):
        region_nodes = sorted(region)
        S = region_nodes[0]  # min id — stable, was a removed node so no collision

        # pop-weighted centroid (fallback uniform)
        weights = [(n, w(pop, n)) for n in region_nodes]
        tw = sum(wt for _, wt in weights)
        def wavg(attr):
            vals = [(float(G.nodes[n][attr]), (wt if tw > 0 else 1.0)) for n, wt in weights]
            sw = sum(wt for _, wt in vals)
            return sum(v * wt for v, wt in vals) / sw
        cx, cy = wavg("x"), wavg("y")
        clon, clat = utm_to_wgs.transform(cx, cy)

        absorbed_orig = []
        rbiz_sum = rsch_sum = rpark_sum = 0.0
        rpop_sum = 0.0
        for n in region_nodes:
            o = G.nodes[n].get("osmid_original", str(n))
            absorbed_orig.append(str(o))
            rpop_sum += w(pop, n); rbiz_sum += w(biz, n)
            rsch_sum += w(sch, n); rpark_sum += w(park, n)
            for d in (pop_o, biz_o, sch_o, park_o):
                d.pop(str(n), None)
            internal.discard(n)
            G.remove_node(n)

        # add super-node
        G.add_node(S, x=cx, y=cy, lon=clon, lat=clat, street_count=1,
                   osmid_original=str([int(x) for x in region_nodes]))
        internal.add(S)
        pop_o[str(S)] = rpop_sum
        biz_o[str(S)] = rbiz_sum
        sch_o[str(S)] = rsch_sum
        if park_o is not None:
            park_o[str(S)] = rpark_sum

        # synthetic directed links E<->S reproducing t_fwd / t_rev after osmnx re-augment
        def synth_len(t):
            return max(t, 0.1) * SYNTH_SPEED_KPH / 3.6
        common = dict(highway=SYNTH_HIGHWAY, name="collapsed dead-end",
                      oneway=False, reversed=False, maxspeed=str(int(SYNTH_SPEED_KPH)))
        G.add_edge(E, S, key=0, length=synth_len(t_fwd), **common)
        G.add_edge(S, E, key=0, length=synth_len(t_rev), **common)

        absorbed_all |= set(region_nodes)
        deadend_map[str(S)] = {
            "entrance": int(E), "absorbed": [int(n) for n in region_nodes],
            "absorbed_osmid_original": absorbed_orig,
            "pop": round(rpop_sum, 3), "biz": round(rbiz_sum, 3),
            "school": round(rsch_sum, 3),
            "t_fwd_s": round(t_fwd, 2), "t_rev_s": round(t_rev, 2),
            "maxT_s": round(maxT, 2), "n_nodes": len(region_nodes),
        }

    nw_out["internal_node_ids"] = sorted(internal)

    # ── Broken-observation cross-check ────────────────────────────────────────
    broken = {"link_aadt": [], "count_sites": []}
    if os.path.exists(LINK_AADT):
        with open(LINK_AADT) as f:
            la = json.load(f)
        for key in la.get("links", {}):
            try:
                u, v = (int(x) for x in key.split(","))
            except ValueError:
                continue
            if u in absorbed_all or v in absorbed_all:
                broken["link_aadt"].append(
                    {"link": key, "u_eaten": u in absorbed_all, "v_eaten": v in absorbed_all,
                     "aadt": la["links"][key].get("aadt")})
    for site in model.COUNT_SITES:
        hit = []
        if site.get("node") is not None and int(site["node"]) in absorbed_all:
            hit.append(int(site["node"]))
        for (u, v) in (site.get("links") or []):
            if int(u) in absorbed_all or int(v) in absorbed_all:
                hit.append([int(u), int(v)])
        if hit:
            broken["count_sites"].append({"label": site["label"], "eaten": hit})

    # ── Save ──────────────────────────────────────────────────────────────────
    ox.save_graphml(G, OUT_GRAPH)
    with open(OUT_WEIGHTS, "w") as f:
        json.dump(nw_out, f, indent=2)
    with open(OUT_MAP, "w") as f:
        json.dump(deadend_map, f, indent=2)
    with open(OUT_BROKEN, "w") as f:
        json.dump(broken, f, indent=2)

    print(f"\nReduction summary:")
    print(f"  nodes: {len(gnodes)} -> {G.number_of_nodes()} "
          f"({len(gnodes) - G.number_of_nodes()} removed)")
    print(f"  regions collapsed: {len(selected)}")
    print(f"  pop moved to super-nodes: {sum(deadend_map[s]['pop'] for s in deadend_map):.0f}")
    print(f"  biz moved to super-nodes: {sum(deadend_map[s]['biz'] for s in deadend_map):.0f}")
    n_broken = len(broken['link_aadt']) + len(broken['count_sites'])
    if n_broken:
        print(f"  *** {len(broken['link_aadt'])} observed link(s) and "
              f"{len(broken['count_sites'])} count site(s) hit absorbed nodes — "
              f"see {OUT_BROKEN} (MANUAL REVIEW NEEDED)")
    else:
        print(f"  no observed links or count sites affected")
    print(f"\nWrote: {OUT_GRAPH}\n       {OUT_WEIGHTS}\n       {OUT_MAP}\n       {OUT_BROKEN}")


if __name__ == "__main__":
    main()
