"""
Surgical edits to the road network model.

Deletions requested (all at study-area boundary):
  Cons 91  (OSM 138036003)               — dead-end stub, Portaferry Road 91↔92 link
  Cons 115 (OSM 181111983)               — Andersons Hill terminus
  Cons 118 (OSM 181112908 + 537845355)   — Ballybarnes Road terminus cluster
  Cons 616 (OSM 469938721)               — Ballyalton Road stub beyond boundary node 617
  Cons 732 (OSM 6457073140 + 6457073138) — Bangor Road residential bypass; makes 731
                                           the clean trunk boundary node
  Cons 767 (OSM 548735778)               — Drumhirk Way terminus

Highway type overrides:
  Hardford Link (cons 18↔21, 21↔68)     — reclassified tertiary → primary; it carries
                                           primary-route traffic levels and was being
                                           systematically under-routed as tertiary.

Removes the nodes and all their incident edges from both the consolidated
graph and the raw graph, then overwrites the saved GraphML files.
"""

import osmnx as ox

CONS_PATH = "simulation/newtownards_consolidated.graphml"
RAW_PATH  = "simulation/newtownards_network.graphml"

# ── Consolidated graph ─────────────────────────────────────────────────────────

G_cons = ox.load_graphml(CONS_PATH)
print(f"Consolidated before: {G_cons.number_of_nodes()} nodes, {G_cons.number_of_edges()} edges")

CONS_NODES_TO_REMOVE = [91, 115, 118, 616, 732, 767]
for n in CONS_NODES_TO_REMOVE:
    if n in G_cons.nodes:
        edges_removed = list(G_cons.edges(n)) + list(G_cons.in_edges(n))
        print(f"  Removing cons node {n} ({len(edges_removed)} incident edges)")
        G_cons.remove_node(n)
    else:
        print(f"  Cons node {n}: already absent, skipping")

# Reclassify Hardford Link edges as primary
for u, v in [(18, 21), (21, 18), (21, 68), (68, 21)]:
    if G_cons.has_edge(u, v):
        for key in G_cons[u][v]:
            old = G_cons[u][v][key].get("highway", "?")
            G_cons[u][v][key]["highway"] = "primary"
            print(f"  Hardford Link cons {u}→{v}: {old} → primary")

print(f"Consolidated after:  {G_cons.number_of_nodes()} nodes, {G_cons.number_of_edges()} edges")
ox.save_graphml(G_cons, CONS_PATH)

# ── Raw graph ──────────────────────────────────────────────────────────────────

G_raw = ox.load_graphml(RAW_PATH)
print(f"\nRaw before: {G_raw.number_of_nodes()} nodes, {G_raw.number_of_edges()} edges")

RAW_NODES_TO_REMOVE = [
    138036003,           # cons 91  — Portaferry Road stub
    181111983,           # cons 115 — Andersons Hill / Ballybarnes Road
    181112908,           # cons 118 — Ballybarnes Road cluster
    537845355,           # cons 118 — Ballybarnes Road cluster (merged node)
    469938721,           # cons 616 — Ballyalton Road stub
    6457073140,          # cons 732 — Bangor Road residential bypass (merged node)
    6457073138,          # cons 732 — Bangor Road residential bypass (merged node)
    548735778,           # cons 767 — Drumhirk Way
]
for n in RAW_NODES_TO_REMOVE:
    if n in G_raw.nodes:
        edges_removed = list(G_raw.edges(n)) + list(G_raw.in_edges(n))
        print(f"  Removing raw node {n} ({len(edges_removed)} incident edges)")
        G_raw.remove_node(n)
    else:
        print(f"  Raw node {n}: already absent, skipping")

# Reclassify Hardford Link edges as primary (match by name)
changed = 0
for u, v, key, data in G_raw.edges(data=True, keys=True):
    name = data.get("name", "")
    if isinstance(name, list):
        name = name[0] if name else ""
    if name == "Hardford Link":
        G_raw[u][v][key]["highway"] = "primary"
        changed += 1
print(f"  Hardford Link raw edges reclassified: {changed} (tertiary → primary)")

print(f"Raw after:  {G_raw.number_of_nodes()} nodes, {G_raw.number_of_edges()} edges")
ox.save_graphml(G_raw, RAW_PATH)

print("\nDone. Re-run build_demographics.py to regenerate the map.")

# ── Virtual Dundonald boundary node ───────────────────────────────────────────
# Node 10000: degree-1 stub connected only to node 97. Gives Dundonald its own
# boundary node on the A20 corridor without altering internal routing (no
# interior-to-interior path ever routes through a dead-end stub).

DUNDONALD_NODE_ID = 10000

if DUNDONALD_NODE_ID not in G_cons.nodes:
    x97, y97 = float(G_cons.nodes[97]["x"]), float(G_cons.nodes[97]["y"])
    G_cons.add_node(DUNDONALD_NODE_ID, x=x97, y=y97, street_count=1)
    for u, v in [(97, DUNDONALD_NODE_ID), (DUNDONALD_NODE_ID, 97)]:
        G_cons.add_edge(u, v, key=0,
                        osmid=0, name="Dundonald virtual link",
                        highway="unclassified", oneway=False,
                        length=1.0, lanes=1, ref="")
    print(f"  Added virtual Dundonald node {DUNDONALD_NODE_ID} (stub on node 97)")
    ox.save_graphml(G_cons, CONS_PATH)
    print(f"  Consolidated graph saved: "
          f"{G_cons.number_of_nodes()} nodes, {G_cons.number_of_edges()} edges")
else:
    print(f"  Virtual Dundonald node {DUNDONALD_NODE_ID}: already present, skipping")
