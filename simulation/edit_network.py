"""
Surgical edits to the road network model.

Deletions requested (all at study-area boundary):
  Cons 91  (OSM 138036003)   — dead-end stub, Portaferry Road 91↔92 link
  Cons 115 (OSM 181111983)   — Andersons Hill terminus
  Cons 118 (OSM 181112908 + 537845355) — Ballybarnes Road terminus cluster
  Cons 767 (OSM 548735778)   — Drumhirk Way terminus

Removes the nodes and all their incident edges from both the consolidated
graph and the raw graph, then overwrites the saved GraphML files.
"""

import osmnx as ox

CONS_PATH = "model/newtownards_consolidated.graphml"
RAW_PATH  = "model/newtownards_network.graphml"

# ── Consolidated graph ─────────────────────────────────────────────────────────

G_cons = ox.load_graphml(CONS_PATH)
print(f"Consolidated before: {G_cons.number_of_nodes()} nodes, {G_cons.number_of_edges()} edges")

CONS_NODES_TO_REMOVE = [91, 115, 118, 767]
for n in CONS_NODES_TO_REMOVE:
    edges_removed = list(G_cons.edges(n)) + list(G_cons.in_edges(n))
    print(f"  Removing cons node {n} ({len(edges_removed)} incident edges)")
    G_cons.remove_node(n)

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
    548735778,           # cons 767 — Drumhirk Way
]
for n in RAW_NODES_TO_REMOVE:
    if n in G_raw.nodes:
        edges_removed = list(G_raw.edges(n)) + list(G_raw.in_edges(n))
        print(f"  Removing raw node {n} ({len(edges_removed)} incident edges)")
        G_raw.remove_node(n)
    else:
        print(f"  Raw node {n}: already absent, skipping")

print(f"Raw after:  {G_raw.number_of_nodes()} nodes, {G_raw.number_of_edges()} edges")
ox.save_graphml(G_raw, RAW_PATH)

print("\nDone. Re-run build_demographics.py to regenerate the map.")
