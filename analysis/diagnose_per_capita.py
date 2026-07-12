"""
Fit diagnostic: true trips-per-capita for the core domain.

Broad sanity check on the deployed gravity model.  Sums the modelled daily
car-driver trips whose journey STARTS OR ENDS inside the core network, by trip
category and overall, and divides by the core area's resident POPULATION (a true
per-capita — never by producers/attractors).  Restricting to core-touching trips
keeps the number inside the well-modelled domain; dividing by core population
makes it directly comparable to NI travel-survey headline rates (TSNI car-driver
trips/person/day by purpose).  A per-capita far from survey values flags a
generation/distribution problem before any spatial detail is examined.

Fully portable: core, population, generation rates and AADT weights are all read
from the per-centre artifacts (no hardcoded node IDs, no Newtownards assumptions),
so a moved CENTRE or altered radius flows through unchanged once the pipeline is
re-run.

Per component c the per-OD-pair DAILY trips are  τ^c_ij = K_c · W_c · t^c_ij
(the same "true daily trips" basis build_assignment.py uses for ext_node_trips):
  t^c   pre-K per-OD-pair flow from model.constrained_od_flows
  K_c   tuned component scale (tuned_params.json)
  W_c   daily AADT weight (model.aadt_weights)
Trips are counted one-way (each OD pair is one trip leg), matching NTS/TSNI; a
round trip appears as its two legs (out + return are separate OD pairs).

Read-only — writes nothing.  Run from the repo root:
    python3 analysis/diagnose_per_capita.py
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, "simulation")
from model import (  # noqa: E402
    PATHS_CACHE, WEIGHTS_FILE, TUNED_PARAMS, SCHOOL_LEVELS,
    constrained_od_flows, aadt_weights,
    load_self_terms, load_generation_rates, compute_generation_scales,
    assert_paths_cache_fresh, willingness_keys, willingness_from_flat,
)

# ── Load node weights (mirrors build_assignment.py) ───────────────────────────

with open(WEIGHTS_FILE) as f:
    weights = json.load(f)

_pnid = lambda k: (int(k) if k.lstrip("-").isdigit() else k)
node_population        = {_pnid(k): v for k, v in weights["node_population"].items()}
node_commute_attractor = {_pnid(k): v for k, v in weights.get("node_commute_attractor", {}).items()}
node_retail_spaces     = {_pnid(k): v for k, v in weights.get("node_retail_spaces", {}).items()}
node_commute_producers = {_pnid(k): v for k, v in weights.get("node_commute_producers", {}).items()}
node_school_demand    = {lvl: {_pnid(k): v for k, v in weights.get(f"node_school_demand_{lvl}", {}).items()}
                         for lvl in SCHOOL_LEVELS}
node_school_producers = {lvl: {_pnid(k): v for k, v in weights.get(f"node_school_producers_{lvl}", {}).items()}
                         for lvl in SCHOOL_LEVELS}

# ── Load tuned params (K's, willingness kernels, slot fractions, Furness set) ──

if not os.path.exists(TUNED_PARAMS):
    raise SystemExit(f"tuned_params.json not found ({TUNED_PARAMS}) — run the tuner first.")
with open(TUNED_PARAMS) as f:
    _tp = json.load(f)

if not all(k in _tp for k in willingness_keys()):
    raise SystemExit("tuned_params.json lacks the 18 double-exp willingness params — "
                     "run reset_gravity_params.py (after sync_kernel_anchor.py) then tune_assignment.py")
willingness = willingness_from_flat(_tp)
if not all(k in _tp for k in ("K_res", "K_commute", "K_retail")):
    raise SystemExit("tuned_params.json lacks the multi-component scales (K_res/K_commute/K_retail).")
K_res     = _tp["K_res"]
K_commute = _tp["K_commute"]
K_retail  = _tp["K_retail"]
K_school  = {lvl: _tp.get(f"K_{lvl}", 0.0) for lvl in SCHOOL_LEVELS}

_parse_sf = lambda key: {tuple(int(x) for x in k.split(",")): v for k, v in _tp.get(key, {}).items()}
slot_fracs_res     = _parse_sf("slot_fracs_res")
slot_fracs_commute = _parse_sf("slot_fracs_commute")
slot_fracs_retail  = _parse_sf("slot_fracs_retail")
slot_fracs_school  = {lvl: _parse_sf(f"slot_fracs_school_{lvl}") for lvl in SCHOOL_LEVELS}

_DBLC = None   # doubly-constrained (Furness) components; None ⇒ singly everywhere
_dblc_raw = _tp.get("doubly_constrained")
if _dblc_raw:
    _valid_dblc = {"commute", "retail"} | {f"school_{lvl}" for lvl in SCHOOL_LEVELS}
    _bad = set(_dblc_raw) - _valid_dblc
    if _bad:
        raise SystemExit(f"tuned_params.json doubly_constrained has unknown components {sorted(_bad)}")
    _DBLC = set(_dblc_raw)

# ── Load paths cache ──────────────────────────────────────────────────────────

if not os.path.exists(PATHS_CACHE):
    raise SystemExit(f"paths cache not found ({PATHS_CACHE}) — run build_paths.py first.")
cache = np.load(PATHS_CACHE, allow_pickle=True)
assert_paths_cache_fresh(cache)
node_ids_arr = cache["node_ids"]
od_src       = cache["od_src"]
od_dst       = cache["od_dst"]
od_dist      = cache["od_dist"].astype(np.float64)
N_nodes      = len(node_ids_arr)

# ── Node-indexed weight arrays (mirrors build_assignment.py) ──────────────────

w_pop          = np.array([node_population.get(nid, 0)         for nid in node_ids_arr], dtype=np.float64)
w_commute_attr = np.array([node_commute_attractor.get(nid, 0) for nid in node_ids_arr], dtype=np.float64)
w_retail       = np.array([node_retail_spaces.get(nid, 0)     for nid in node_ids_arr], dtype=np.float64)
w_commute_prod = np.array([node_commute_producers.get(nid, 0) for nid in node_ids_arr], dtype=np.float64)
w_school_levels = {lvl: np.array([node_school_demand[lvl].get(nid, 0) for nid in node_ids_arr],
                                 dtype=np.float64) for lvl in SCHOOL_LEVELS}
w_school_prod_levels = {lvl: np.array([node_school_producers[lvl].get(nid, 0) for nid in node_ids_arr],
                                      dtype=np.float64) for lvl in SCHOOL_LEVELS}

# Generation pinning (per-leg producer scales, vehicle-driver trips/day) + self-term.
_gen_rates = load_generation_rates()
_GEN_SCALE = (compute_generation_scales(weights, _gen_rates) if _gen_rates is not None else None)
self_terms = load_self_terms(list(node_ids_arr))

_active_school = [lvl for lvl in SCHOOL_LEVELS
                  if K_school.get(lvl, 0.0) > 0 and w_school_levels[lvl].sum() > 0]
_use_school = len(_active_school) > 0

# ── Compute pre-K per-OD-pair component flows (same call as build_assignment) ──

print(f"Computing OD flows (doubly_constrained={sorted(_DBLC) if _DBLC else []}) …")
t_res, t_commute, t_retail, t_sch_by_level = constrained_od_flows(
    od_src, od_dst, od_dist, N_nodes,
    w_pop, w_commute_attr, w_retail,
    willingness,
    with_school=_use_school,
    w_school_levels=w_school_levels, w_school_prod_levels=w_school_prod_levels,
    self_terms=self_terms,
    w_commute_prod=w_commute_prod,
    gen_scale=_GEN_SCALE,
    doubly_constrained=_DBLC)

# Daily AADT weights W_c (component daily-trip conversion).
W_res, W_commute, W_retail, W_school = aadt_weights(
    slot_fracs_res, slot_fracs_commute, slot_fracs_retail, slot_fracs_school)

# ── Portable core definition ──────────────────────────────────────────────────
# Core road-node set = the explicit internal_node_ids (rebuilt per-centre by
# build_demographics.py / reduce_deadends.py); this is how build_external_links.py
# defines core, and it correctly excludes pure-numeric RoI external codes.
core = set(int(x) for x in weights["internal_node_ids"])
is_core = np.array([nid in core for nid in node_ids_arr], dtype=bool)
pop_core = float(w_pop[is_core].sum())
if pop_core <= 0:
    raise SystemExit("core population is zero — check internal_node_ids / node_population.")
n_core_nodes = int(is_core.sum())

src_core = is_core[od_src]
dst_core = is_core[od_dst]
touch_mask   = src_core | dst_core            # trip starts OR ends in core (counted once)
out_mask     = src_core                       # origin in core: resident-generated (incl. internal→internal)
inbound_mask = dst_core & ~src_core           # external → core: inbound-only load

_rho = _gen_rates or {}

# component: (label, t_array, K_c, W_c, rho_key)
_components = [
    ("residential",      t_res,     K_res,     W_res,     "res"),
    ("commute",          t_commute, K_commute, W_commute, "commute"),
    ("retail",           t_retail,  K_retail,  W_retail,  "retail"),
]
for lvl in SCHOOL_LEVELS:
    _components.append((f"school_{lvl}", t_sch_by_level[lvl], K_school[lvl],
                        W_school[lvl], f"school_{lvl}"))


def _per_capita(t, K_c, W_c, mask):
    return float((t[mask] * K_c * W_c).sum()) / pop_core


# ── Report ────────────────────────────────────────────────────────────────────

print()
print("Trips-per-capita diagnostic (core domain)")
print(f"  core nodes: {n_core_nodes:,}   core population: {pop_core:,.0f}"
      f"   doubly_constrained: {sorted(_DBLC) if _DBLC else '[]'}")
print("  per-capita = daily car-driver trips (start OR end in core) / core population")
print("  one-way legs; outbound = origin in core (resident-generated); inbound = external→core")
print()
hdr = (f"  {'Component':<20s}  {'per-capita':>10s}  {'outbound':>9s}  {'inbound':>8s}"
       f"  {'ρ (input)':>10s}  {'K_c':>7s}")
print(hdr)
print("  " + "-" * (len(hdr) - 2))

tot_touch = tot_out = tot_in = 0.0
for label, t, K_c, W_c, rho_key in _components:
    pc_touch = _per_capita(t, K_c, W_c, touch_mask)
    pc_out   = _per_capita(t, K_c, W_c, out_mask)
    pc_in    = _per_capita(t, K_c, W_c, inbound_mask)
    tot_touch += pc_touch
    tot_out   += pc_out
    tot_in    += pc_in
    rho = _rho.get(rho_key)
    rho_s = f"{rho:>10.4f}" if rho is not None else f"{'—':>10s}"
    print(f"  {label:<20s}  {pc_touch:>10.4f}  {pc_out:>9.4f}  {pc_in:>8.4f}  {rho_s}  {K_c:>7.3f}")

print("  " + "-" * (len(hdr) - 2))
rho_tot = sum(v for v in _rho.values()) if _rho else None
rho_tot_s = f"{rho_tot:>10.4f}" if rho_tot is not None else f"{'—':>10s}"
print(f"  {'OVERALL':<20s}  {tot_touch:>10.4f}  {tot_out:>9.4f}  {tot_in:>8.4f}  {rho_tot_s}  {'':>7s}")
print()
print("  ρ (input) = generation_rates.json per-capita car-driver trips/person/day (island-wide).")
print("  Compare the outbound column to TSNI car-driver trips/person/day by purpose.")
