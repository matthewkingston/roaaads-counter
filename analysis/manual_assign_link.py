"""
Manually assign a walking count session to a specific directed link.

Usage:
  python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>

The from_node → to_node direction becomes the "with" direction for the session.
The to_node → from_node direction becomes the "against" direction.

Validates the assignment against the road network, updates
data/manual_link_overrides.json, and patches data/counts_processed.json.
Idempotent: re-running with any arguments overwrites the previous assignment.
If counts_processed.json is wiped and ingest_counts.py re-run, the override
in manual_link_overrides.json is applied automatically.
"""

import json
import os
import sys

import osmnx as ox

CONS_GRAPH            = "simulation/newtownards_consolidated.graphml"
MANUAL_OVERRIDES_FILE = "data/manual_link_overrides.json"
PROCESSED_FILE        = "data/counts_processed.json"


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>")
        sys.exit(1)

    session_id = sys.argv[1]
    try:
        from_node = int(sys.argv[2])
        to_node   = int(sys.argv[3])
    except ValueError:
        print("Error: from_node and to_node must be integers")
        sys.exit(1)

    if from_node == to_node:
        print("Error: from_node and to_node must be different")
        sys.exit(1)

    print(f"Loading graph …")
    G = ox.load_graphml(CONS_GRAPH)

    if from_node not in G.nodes:
        print(f"Error: node {from_node} not in network")
        sys.exit(1)
    if to_node not in G.nodes:
        print(f"Error: node {to_node} not in network")
        sys.exit(1)

    with_edge_ok    = G.has_edge(from_node, to_node)
    against_edge_ok = G.has_edge(to_node, from_node)

    if not with_edge_ok and not against_edge_ok:
        print(f"Error: neither {from_node}→{to_node} nor {to_node}→{from_node} "
              f"exists in the network — check node IDs")
        sys.exit(1)

    link_with    = [from_node, to_node]
    link_against = [to_node, from_node]

    def _road_name():
        for u, v in [link_with, link_against]:
            if G.has_edge(u, v):
                name = G[u][v][0].get("name")
                if name:
                    return "; ".join(name) if isinstance(name, list) else name
        return None

    road_name = _road_name()

    if with_edge_ok:
        print(f"  With direction    {from_node}→{to_node}: edge exists ✓")
    else:
        print(f"  With direction    {from_node}→{to_node}: no edge (walking against one-way traffic)")
    if against_edge_ok:
        print(f"  Against direction {to_node}→{from_node}: edge exists ✓")
    else:
        print(f"  Against direction {to_node}→{from_node}: no edge (one-way road)")

    # Load processed sessions
    if os.path.exists(PROCESSED_FILE):
        with open(PROCESSED_FILE) as f:
            processed = json.load(f)
    else:
        processed = {"sessions": {}}

    rec = processed["sessions"].get(session_id)
    if rec is None:
        print(f"  Session {session_id} not yet in {PROCESSED_FILE} "
              f"— override will apply on next ingest_counts.py run")
    else:
        mode          = rec.get("mode", "unknown")
        with_count    = rec.get("with_count")
        against_count = rec.get("against_count")
        print(f"  Session mode: {mode}, with_count={with_count}, against_count={against_count}")

        # Non-null count must map to a real directed edge (mirrors ingest_counts.py rule)
        if not with_edge_ok and with_count is not None:
            print(
                f"\nError: with_count={with_count} on non-existent edge {from_node}→{to_node}.\n"
                f"If the observer walked against one-way traffic, with_count should be null. "
                f"Swap from/to if the observer was counting the direction cars actually travel."
            )
            sys.exit(1)
        if not against_edge_ok and against_count is not None:
            print(
                f"\nError: against_count={against_count} on non-existent edge {to_node}→{from_node}.\n"
                f"Swap from/to if the observer was counting the direction cars actually travel."
            )
            sys.exit(1)

    # Update overrides file
    if os.path.exists(MANUAL_OVERRIDES_FILE):
        with open(MANUAL_OVERRIDES_FILE) as f:
            overrides = json.load(f)
    else:
        overrides = {}

    prev = overrides.get(session_id)
    overrides[session_id] = {"link_with": link_with, "link_against": link_against}

    with open(MANUAL_OVERRIDES_FILE, "w") as f:
        json.dump(overrides, f, indent=2)

    if prev is None:
        print(f"  Added to {MANUAL_OVERRIDES_FILE}")
    elif prev["link_with"] != link_with:
        print(f"  Updated {MANUAL_OVERRIDES_FILE} "
              f"(was {prev['link_with'][0]}→{prev['link_with'][1]})")
    else:
        print(f"  {MANUAL_OVERRIDES_FILE} unchanged (same assignment)")

    # Patch counts_processed.json directly for immediate effect
    if rec is not None:
        rec["matched_link_with"]    = link_with
        rec["matched_link_against"] = link_against
        rec["match_rmse_m"]         = None
        rec["match_method"]         = "manual"
        rec["matched_link_name"]    = road_name

        with open(PROCESSED_FILE, "w") as f:
            json.dump(processed, f, indent=2)
        print(f"  Patched {PROCESSED_FILE}")

    print(f"\nSession {session_id} → link {from_node}→{to_node} "
          f"(against: {to_node}→{from_node})")
    if rec is not None:
        print("Re-run aggregate_counts.py (and tune_assignment.py) to update the model.")


if __name__ == "__main__":
    main()
