"""
Shared constants and functions for the Newtownards gravity model pipeline.

Imported by simulation/build_assignment.py and analysis/tune_assignment.py to
keep their chi²/N calculations, count-site definitions, and gravity kernel
implementations in sync.
"""

import json
import math
import os
import numpy as np

# ── Official AADT count sites ─────────────────────────────────────────────────
# AADT totals retained for the quick site-level sanity check in build_assignment.py.
# The tuner and compute_chi2() use per-hour observations from OFFICIAL_HOURLY instead.

COUNT_SITES = [
    {"label": "site 507, A21 Bangor Road",     "node": None, "links": [(538692601,550205936),(550205936,538692601)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road", "node": 136173611, "links": None, "observed": 10_792},
    {"label": "site 444, A20 Portaferry Road", "node": 449111329, "links": None, "observed":  7_282},
]

# Links present in link_aadt.json but excluded from calibration.
# Directed: (u, v) excludes only that direction.
# The Westmount Park / Old Belfast Road links lie inside dead-end regions collapsed by
# reduce_deadends.py (their endpoints are absorbed into super-nodes and no longer exist in
# the reduced routing graph), so their observations are discarded from calibration. If the
# dead-end reduction parameters change, regenerate this set from deadend_broken_obs.json.
EXCLUDE_LINKS = {
    (181844513, 181839481),
    (448393355, 538253737), (538253737, 448393355),   # Westmount Park
    (538253737, 7085530067), (7085530067, 538253737),  # Westmount Park
    (540663959, 6620711226), (6620711226, 540663959),  # Old Belfast Road
}

# ── File paths ────────────────────────────────────────────────────────────────
# Routing graph + node weights come from reduce_deadends.py (residential dead-ends
# collapsed). build_demographics.py still writes the full node_weights.json and uses the
# full consolidated graph; reduce_deadends.py reads those and writes the *_reduced.* files
# consumed here, by build_paths.py, build_assignment.py and tune_assignment.py.

PATHS_CACHE     = "simulation/newtownards_paths.npz"
WEIGHTS_FILE    = "simulation/node_weights_reduced.json"
ROUTING_GRAPH   = "simulation/newtownards_reduced.graphml"
EXTERNAL_LINKS  = "data/external_links.json"
TUNER_CONFIG    = "simulation/tuner_config.json"
TUNED_PARAMS    = "simulation/tuned_params.json"
LINK_AADT       = "data/link_aadt.json"
OFFICIAL_HOURLY = "data/official_hourly.json"

_DOW_TO_TYPE = {d: (0 if d < 5 else (1 if d == 5 else 2)) for d in range(7)}

# ── Paths-cache freshness guard ─────────────────────────────────────────────────
# build_paths.py stamps a signature of its inputs (routing graph, external links,
# HIGHWAY_COST_FACTOR) into the .npz. tune_assignment.py and build_assignment.py
# re-check it at load time and fail loudly if the cache is stale, rather than
# silently assigning/tuning against an out-of-date cache (a recurring footgun —
# see CLAUDE.md "Paths cache note").

def _file_sha1(path):
    import hashlib
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def paths_cache_signature():
    """Signature of the inputs build_paths.py consumes. Stamped into the npz at
    build time and re-checked at load time. Returns a dict of {field: str}."""
    from routing_config import HIGHWAY_COST_FACTOR, PROBIT_CV, PROBIT_LL_SIGMA
    return {
        "src_graph_sha1":      _file_sha1(ROUTING_GRAPH),
        "src_extlinks_sha1":   _file_sha1(EXTERNAL_LINKS),
        "src_cost_factor":     json.dumps(HIGHWAY_COST_FACTOR, sort_keys=True),
        "src_probit_cv":       repr(float(PROBIT_CV)),
        "src_probit_ll_sigma": repr(float(PROBIT_LL_SIGMA)),
    }

def assert_paths_cache_fresh(cache):
    """Raise SystemExit if the loaded paths cache was built from different inputs
    than the current pipeline state. `cache` is the np.load(...) handle."""
    label = {
        "src_graph_sha1":      "routing graph (newtownards_reduced.graphml)",
        "src_extlinks_sha1":   "external links (data/external_links.json)",
        "src_cost_factor":     "HIGHWAY_COST_FACTOR (simulation/routing_config.py)",
        "src_probit_cv":       "PROBIT_CV (simulation/routing_config.py)",
        "src_probit_ll_sigma": "PROBIT_LL_SIGMA (simulation/routing_config.py)",
    }
    sig = paths_cache_signature()
    stale = []
    for key, current in sig.items():
        stored = cache[key].item() if key in cache else None
        if stored != current:
            stale.append(label.get(key, key))
    if stale:
        msg = ["Paths cache is STALE — it was built from different inputs than the "
               "current pipeline state:"]
        msg += [f"  - {s} changed since the cache was built" for s in stale]
        msg.append("Re-run: python3 simulation/build_paths.py")
        raise SystemExit("\n".join(msg))

# ── Gravity kernel ────────────────────────────────────────────────────────────

def gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
                   W_BIZ, P, ALPHA, w_pop, w_biz,
                   BETA=1.0, THETA=None,
                   P_biz=None, ALPHA_biz=None, BETA_biz=None,
                   W_SCHOOL=None, P_school=None, ALPHA_school=None, w_school=None,
                   od_dist_2=None, pair_idx_2=None, link_idx_2=None,
                   od_dist_3=None, pair_idx_3=None, link_idx_3=None,
                   link_weight=None,
                   return_components=False):
    """
    Generalised rational kernel assignment.

    Kernel: f(d) = (ALPHA+BETA)*u^BETA / (ALPHA + BETA*u^(ALPHA+BETA))  where u = d/P.
    Peak at d=P with f(P)=1; tail ~ 1/d^ALPHA; rise ~ u^BETA near origin.
    BETA=1 (default) recovers the original kernel (ALPHA+1)*u / (ALPHA + u^(ALPHA+1)).

    link_weight: optional float32 array parallel to (pair_idx, link_idx), giving the
      fraction of probit stochastic passes that routed through each (pair, link) entry.
      When provided, each entry's flow contribution is scaled by link_weight[entry]
      instead of 1.0 (binary all-or-nothing).  When None, all-or-nothing on k=1 path.

    THETA / k=2/k=3 arrays: legacy logit stochastic routing retained for backward
      compatibility with old caches.  New caches use link_weight instead.

    P_biz/ALPHA_biz/BETA_biz: optional separate kernel for the business component.
      Only used when return_components=True.

    W_SCHOOL/P_school/ALPHA_school/w_school: optional school-trip component.
      When w_school is provided and return_components=True, returns a third component:
      flow_school = W_SCHOOL·(pop→school + school→pop) using P_school/ALPHA_school kernel.
      w_school uses BETA (shared with residential kernel).

    return_components=False (default): returns combined pre-K flow array (N_links,).
    return_components=True, w_school=None: returns (flow_res, flow_biz).
    return_components=True, w_school provided: returns (flow_res, flow_biz, flow_school).
    """
    _has_stoch = (THETA is not None and od_dist_2 is not None)

    # Resolve effective biz kernel params (fall back to shared params when not set)
    _P_biz    = P_biz    if P_biz    is not None else P
    _A_biz    = ALPHA_biz if ALPHA_biz is not None else ALPHA
    _B_biz    = BETA_biz  if BETA_biz  is not None else BETA

    _has_school = (w_school is not None and W_SCHOOL is not None and W_SCHOOL > 0)
    _P_sch  = P_school    if P_school    is not None else P
    _A_sch  = ALPHA_school if ALPHA_school is not None else ALPHA

    if not _has_stoch:
        u    = od_dist / P
        kern = (ALPHA + BETA) * u**BETA / (ALPHA + BETA * u**(ALPHA + BETA))

        # entry_w: per-entry multiplier (probit fractional weight, or 1.0 for binary)
        entry_w = link_weight if link_weight is not None else 1.0

        if not return_components:
            w_vec = w_pop + W_BIZ * w_biz
            if _has_school:
                w_vec = w_vec + W_SCHOOL * w_school
            t_ij  = w_vec[od_src] * w_vec[od_dst] * kern
            return np.bincount(link_idx, weights=t_ij[pair_idx] * entry_w, minlength=N_links)

        # Separate biz kernel (may equal kern when no biz params supplied)
        if _P_biz == P and _A_biz == ALPHA and _B_biz == BETA:
            kern_biz = kern
        else:
            u_b      = od_dist / _P_biz
            kern_biz = (_A_biz + _B_biz) * u_b**_B_biz / (_A_biz + _B_biz * u_b**(_A_biz + _B_biz))

        pp = w_pop[od_src] * w_pop[od_dst] * kern
        pb = (w_pop[od_src] * w_biz[od_dst] + w_biz[od_src] * w_pop[od_dst]) * kern_biz
        bb = w_biz[od_src] * w_biz[od_dst] * kern_biz
        flow_res = np.bincount(link_idx, weights=pp[pair_idx] * entry_w, minlength=N_links)
        flow_biz = np.bincount(link_idx,
                               weights=((W_BIZ * pb + W_BIZ ** 2 * bb)[pair_idx]) * entry_w,
                               minlength=N_links)
        if not _has_school:
            return flow_res, flow_biz

        # School component: separate kernel, pop×school cross-term only
        if _P_sch == P and _A_sch == ALPHA:
            kern_sch = kern
        else:
            u_s      = od_dist / _P_sch
            kern_sch = (_A_sch + BETA) * u_s**BETA / (_A_sch + BETA * u_s**(_A_sch + BETA))
        ps = (w_pop[od_src] * w_school[od_dst] + w_school[od_src] * w_pop[od_dst]) * kern_sch
        flow_school = np.bincount(link_idx,
                                  weights=(W_SCHOOL * ps[pair_idx]) * entry_w,
                                  minlength=N_links)
        return flow_res, flow_biz, flow_school

    # ── Stochastic logit ──────────────────────────────────────────────────────
    d_mat  = np.stack([od_dist, od_dist_2, od_dist_3], axis=1)
    log_w  = -THETA * d_mat / P   # logit shares always use P (shared routing scale)
    log_w -= log_w.max(axis=1, keepdims=True)
    shares = np.exp(log_w)
    shares /= shares.sum(axis=1, keepdims=True)

    if not return_components:
        w_vec = w_pop + W_BIZ * w_biz
        t_ij  = w_vec[od_src] * w_vec[od_dst]
        flow  = np.zeros(N_links, dtype=np.float64)
        for r, (pidx, lidx, d_r) in enumerate([
            (pair_idx,   link_idx,   od_dist),
            (pair_idx_2, link_idx_2, od_dist_2),
            (pair_idx_3, link_idx_3, od_dist_3),
        ]):
            u_r  = d_r / P
            f_r  = (ALPHA + BETA) * u_r**BETA / (ALPHA + BETA * u_r**(ALPHA + BETA))
            flow += np.bincount(lidx, weights=(t_ij * shares[:, r] * f_r)[pidx],
                                minlength=N_links)
        return flow

    pp_od    = w_pop[od_src] * w_pop[od_dst]
    pb_od    = w_pop[od_src] * w_biz[od_dst] + w_biz[od_src] * w_pop[od_dst]
    bb_od    = w_biz[od_src] * w_biz[od_dst]
    biz_base = W_BIZ * pb_od + W_BIZ ** 2 * bb_od
    flow_res = np.zeros(N_links, dtype=np.float64)
    flow_biz = np.zeros(N_links, dtype=np.float64)
    for r, (pidx, lidx, d_r) in enumerate([
        (pair_idx,   link_idx,   od_dist),
        (pair_idx_2, link_idx_2, od_dist_2),
        (pair_idx_3, link_idx_3, od_dist_3),
    ]):
        u_res = d_r / P
        f_res = (ALPHA + BETA) * u_res**BETA / (ALPHA + BETA * u_res**(ALPHA + BETA))
        u_biz = d_r / _P_biz
        f_biz = (_A_biz + _B_biz) * u_biz**_B_biz / (_A_biz + _B_biz * u_biz**(_A_biz + _B_biz))
        s_r   = shares[:, r]
        flow_res += np.bincount(lidx, weights=(pp_od    * s_r * f_res)[pidx], minlength=N_links)
        flow_biz += np.bincount(lidx, weights=(biz_base * s_r * f_biz)[pidx], minlength=N_links)
    return flow_res, flow_biz

# ── Production-constrained (singly-constrained) assignment ────────────────────
# Replaces the unconstrained T_ij = K·w_i·w_j·F(d) with a singly-constrained form
# applied per component:  T^c_ij = K_c · p^c_i · a^c_j · F_c(d_ij) / D^c_i,
# where D^c_i = Σ_k a^c_k·F_c(d_ik).  Then Σ_j T^c_ij = K_c·p^c_i, i.e. each origin's
# trip production is fixed by its producing weight and is independent of accessibility
# (fixes the generation/distribution conflation — see project memory note
# project_production_constrained_gravity).  Component producer/attractor (locked 2026-06-24):
#   res:  p=pop,  a=pop,    kernel (P,ALPHA,BETA)
#   biz:  symmetric split + per-origin normaliser (pop→biz and biz→pop legs), plus a
#         biz×biz term constrained on biz; weighted by W_BIZ; kernel (P_biz,ALPHA_biz,BETA)
#   sch:  symmetric split + per-origin normaliser (pop→school and school→pop legs);
#         kernel (P_school,ALPHA_school,BETA)
# The per-origin denominators depend on the kernel params, so they are recomputed every
# evaluation (cheap O(N_OD) bincounts).  K_c is applied by the caller (the analytical
# K/φ/f calibration blocks are unchanged — D_i has no K, flow stays linear in K).


def _rational_kernel(d, P, ALPHA, BETA):
    """f(d) = (ALPHA+BETA)·u^BETA / (ALPHA + BETA·u^(ALPHA+BETA)), u = d/P."""
    u = d / P
    return (ALPHA + BETA) * u**BETA / (ALPHA + BETA * u**(ALPHA + BETA))


def constrained_od_flows(od_src, od_dst, od_dist, N_nodes,
                         w_pop, w_biz, w_school,
                         W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz,
                         P_school=None, ALPHA_school=None, with_school=False):
    """Per-OD-pair, pre-K production-constrained component flows.

    Returns (t_res, t_biz, t_sch), each a float64 array parallel to od_src/od_dst.
    These are the per-pair flows BEFORE the K_res/K_biz/K_sch scaling; the caller
    scatters them onto links (all links, or observed rows only) and applies K.

    Denominators are summed over the FULL destination set of each origin (all od
    pairs sharing that od_src — internal, external-routed, and denominator-only
    ext→ext virtual edges), so each origin's production budget is conserved.
    """
    src = od_src
    dst = od_dst
    F_res = _rational_kernel(od_dist, P,     ALPHA,     BETA)
    F_biz = _rational_kernel(od_dist, P_biz, ALPHA_biz, BETA)

    pop_s = w_pop[src]; pop_d = w_pop[dst]
    biz_s = w_biz[src]; biz_d = w_biz[dst]

    # Per-origin denominators D^c_i = Σ_k a^c_k·F_c(d_ik); inverse with 0/0 → 0.
    def _inv_denom(attr_d, F):
        D = np.bincount(src, weights=attr_d * F, minlength=N_nodes)
        return np.where(D > 0, 1.0 / D, 0.0)

    iD_res_pop = _inv_denom(pop_d, F_res)   # res: attraction = pop
    iD_biz_pop = _inv_denom(pop_d, F_biz)   # biz-cross leg biz→pop: attraction = pop
    iD_biz_biz = _inv_denom(biz_d, F_biz)   # biz-cross leg pop→biz + biz×biz: attraction = biz

    # res: pop_i·pop_j·F_res / D^res,pop_i
    t_res = pop_s * pop_d * F_res * iD_res_pop[src]

    # biz: symmetric split (pop→biz, biz→pop) + biz×biz (W_BIZ), each per-origin-normalised
    t_biz = F_biz * (
        pop_s * biz_d * iD_biz_biz[src]            # home i → work j   (attraction biz)
        + biz_s * pop_d * iD_biz_pop[src]          # biz i → pop j     (attraction pop)
        + W_BIZ * biz_s * biz_d * iD_biz_biz[src]  # biz×biz           (attraction biz)
    )

    if with_school and w_school is not None and w_school.sum() > 0:
        F_sch = _rational_kernel(od_dist, P_school, ALPHA_school, BETA)
        sch_s = w_school[src]; sch_d = w_school[dst]
        iD_sch_pop = _inv_denom(pop_d, F_sch)   # school-cross leg school→pop: attraction = pop
        iD_sch_sch = _inv_denom(sch_d, F_sch)   # school-cross leg pop→school: attraction = school
        t_sch = F_sch * (
            pop_s * sch_d * iD_sch_sch[src]       # pop i → school j  (attraction school)
            + sch_s * pop_d * iD_sch_pop[src]     # school i → pop j  (attraction pop)
        )
    else:
        t_sch = np.zeros(len(src), dtype=np.float64)

    return t_res, t_biz, t_sch


def scatter_od_to_links(t_pair, pair_idx, link_idx, link_weight, N_links):
    """Scatter a per-OD-pair flow vector onto links via the probit routing incidence.

    flow[l] = Σ_entries t_pair[pair_idx]·link_weight  bincounted into link_idx.
    Denominator-only OD pairs have no entries in pair_idx/link_idx, so they
    contribute to the denominators (upstream) but carry no link flow here.
    """
    w = t_pair[pair_idx]
    if link_weight is not None:
        w = w * link_weight
    return np.bincount(link_idx, weights=w, minlength=N_links)

# ── Flow extraction ───────────────────────────────────────────────────────────

def site_flow(link_flow_dict, site):
    """Return total modelled flow for a COUNT_SITES entry."""
    if site["links"]:
        return sum(link_flow_dict.get(lnk, 0.0) for lnk in site["links"])
    node = site["node"]
    return sum(f for (u, v), f in link_flow_dict.items() if u == node or v == node)


def _site_flow_2c(flow_res_dict, flow_biz_dict, node, links):
    """Return (m_res, m_biz) for a site defined by node or directed links."""
    if links:
        m_r = sum(flow_res_dict.get(tuple(lnk), 0.0) for lnk in links)
        m_b = sum(flow_biz_dict.get(tuple(lnk), 0.0) for lnk in links)
    else:
        m_r = sum(f for (u, v), f in flow_res_dict.items() if u == node or v == node)
        m_b = sum(f for (u, v), f in flow_biz_dict.items() if u == node or v == node)
    return m_r, m_b


def _site_flow_3c(flow_res_dict, flow_biz_dict, flow_school_dict, node, links):
    """Return (m_res, m_biz, m_school) for a site defined by node or directed links."""
    if links:
        m_r = sum(flow_res_dict.get(tuple(lnk), 0.0)    for lnk in links)
        m_b = sum(flow_biz_dict.get(tuple(lnk), 0.0)    for lnk in links)
        m_s = sum(flow_school_dict.get(tuple(lnk), 0.0) for lnk in links)
    else:
        m_r = sum(f for (u, v), f in flow_res_dict.items()    if u == node or v == node)
        m_b = sum(f for (u, v), f in flow_biz_dict.items()    if u == node or v == node)
        m_s = sum(f for (u, v), f in flow_school_dict.items() if u == node or v == node)
    return m_r, m_b, m_s

# ── Chi²/N ───────────────────────────────────────────────────────────────────

def compute_chi2(link_flow_dict, label_fn=None,
                 link_aadt_file=LINK_AADT, exclude_links=EXCLUDE_LINKS,
                 official_hourly_file=OFFICIAL_HOURLY,
                 link_flow_biz_dict=None, link_flow_school_dict=None,
                 slot_fracs_res=None, slot_fracs_biz=None, slot_fracs_school=None):
    """
    Compute chi²/N matching the tuner's component formulation.

    Three-component mode (link_flow_biz_dict and link_flow_school_dict provided):
      link_flow_dict        — {(u,v): K_res * flow_res}
      link_flow_biz_dict    — {(u,v): K_biz * flow_biz}
      link_flow_school_dict — {(u,v): K_school * flow_school}
      slot_fracs_res/biz/school — {(day_type, hour): f}
      N_eff = N − 3·N_slots.

    Two-component mode (link_flow_biz_dict provided, link_flow_school_dict=None):
      link_flow_dict     — {(u,v): K_res * flow_res}
      link_flow_biz_dict — {(u,v): K_biz * flow_biz}
      slot_fracs_res/biz — {(day_type, hour): f}
      N_eff = N − 2·N_slots.  No coupling penalty (pure data fit).

    Legacy mode (link_flow_biz_dict=None):
      link_flow_dict — {(u,v): K * combined_flow}  (pre-scaled)
      Uses 3 AADT-space obs from COUNT_SITES plus Woodbury-corrected walking obs.
      N_eff = N − N_slots.

    Returns (rows, chi2, n_obs, n_eff).
      rows: list of (kind, label, obs, sig, mod, z) sorted by |z| descending.
    """
    if link_flow_school_dict is not None and link_flow_biz_dict is not None:
        return _compute_chi2_3c(link_flow_dict, link_flow_biz_dict, link_flow_school_dict,
                                slot_fracs_res, slot_fracs_biz, slot_fracs_school,
                                label_fn, link_aadt_file, exclude_links,
                                official_hourly_file)
    if link_flow_biz_dict is not None:
        return _compute_chi2_2c(link_flow_dict, link_flow_biz_dict,
                                slot_fracs_res, slot_fracs_biz,
                                label_fn, link_aadt_file, exclude_links,
                                official_hourly_file)
    return _compute_chi2_legacy(link_flow_dict, label_fn, link_aadt_file, exclude_links)


def _compute_chi2_2c(flow_res_dict, flow_biz_dict,
                     slot_fracs_res, slot_fracs_biz,
                     label_fn, link_aadt_file, exclude_links,
                     official_hourly_file):
    """Two-component chi² matching the tuner's objective (minus coupling penalty)."""
    rows  = []
    chi2  = 0.0
    excl  = exclude_links or set()

    # (day_type, hour) → (f_res, f_biz)
    def _fracs(slot_key):
        f_r = (slot_fracs_res or {}).get(slot_key, 1.0 / 24)
        f_b = (slot_fracs_biz or {}).get(slot_key, 1.0 / 24)
        return f_r, f_b

    # ── Official hourly obs (count-space Gaussian) ────────────────────────────
    n_official = 0
    if official_hourly_file and os.path.exists(official_hourly_file):
        with open(official_hourly_file) as f:
            oh = json.load(f)
        for site_id, site in oh.items():
            node  = site["node"]
            links = [tuple(lnk) for lnk in site["links"]] if site["links"] else None
            m_r, m_b = _site_flow_2c(flow_res_dict, flow_biz_dict, node, links)
            for obs in site["observations"]:
                dt, h    = obs["time_slot"]
                sk       = (dt, h)
                f_r, f_b = _fracs(sk)
                count    = float(obs["count"])
                sigma    = float(obs["sigma"])
                pred     = m_r * f_r + m_b * f_b   # T=3600 → T/3600=1
                z        = (pred - count) / sigma if sigma > 0 else 0.0
                chi2    += z ** 2
                lbl      = f"{site['label']} h{h:02d}"
                rows.append(("official", lbl, count, sigma, pred, z))
                n_official += 1

    # ── Walking obs (count-space Poisson) ─────────────────────────────────────
    n_slots_seen = set()
    if link_aadt_file and os.path.exists(link_aadt_file):
        with open(link_aadt_file) as f:
            link_aadt = json.load(f)["links"]
        for key, entry in sorted(link_aadt.items()):
            u, v = map(int, key.split(","))
            if (u, v) in excl:
                continue
            lbl  = label_fn(u, v) if label_fn else f"{u}→{v}"
            m_r  = flow_res_dict.get((u, v), 0.0)
            m_b  = flow_biz_dict.get((u, v), 0.0)
            for sess in entry.get("observations", []):
                ts   = sess.get("time_slot")
                neff = float(sess.get("n_eff", 0.5))
                dur  = float(sess.get("duration_s", 0.0))
                if ts is None or dur <= 0:
                    continue
                sk       = (_DOW_TO_TYPE[ts[0]], ts[1])
                f_r, f_b = _fracs(sk)
                Th       = dur / 3600.0
                pred     = (m_r * f_r + m_b * f_b) * Th
                z        = (pred - neff) / math.sqrt(neff) if neff > 0 else 0.0
                n_actual  = neff - 0.5
                pred_safe = max(pred, 1e-30)
                if n_actual > 0:
                    chi2 += 2.0 * (pred_safe - n_actual + n_actual * math.log(n_actual / pred_safe))
                else:
                    chi2 += 2.0 * pred_safe
                n_slots_seen.add(sk)
                # display in AADT-space: divide by effective f_s
                m_tot = m_r + m_b + 1e-30
                f_eff = (m_r * f_r + m_b * f_b) / m_tot
                if f_eff > 0 and Th > 0:
                    obs_disp = neff / (Th * f_eff)
                    sig_disp = math.sqrt(neff) / (Th * f_eff)
                    mod_disp = pred / (Th * f_eff)   # = combined AADT (m_r + m_b)
                else:
                    obs_disp = sig_disp = mod_disp = 0.0
                rows.append(("walking", lbl, obs_disp, sig_disp, mod_disp, z))

    n_obs   = len(rows)
    # N_eff: 2 df per slot (f_res and f_biz each).
    if n_official:
        n_slots_seen |= {(dt, h) for dt in range(3) for h in range(24)}
    n_eff   = n_obs - 2 * len(n_slots_seen)

    rows.sort(key=lambda r: abs(r[5]), reverse=True)
    return rows, chi2, n_obs, n_eff


def _compute_chi2_3c(flow_res_dict, flow_biz_dict, flow_school_dict,
                     slot_fracs_res, slot_fracs_biz, slot_fracs_school,
                     label_fn, link_aadt_file, exclude_links,
                     official_hourly_file):
    """Three-component chi² (res + biz + school).  N_eff = N − 3·N_slots."""
    rows  = []
    chi2  = 0.0
    excl  = exclude_links or set()

    def _fracs(slot_key):
        f_r = (slot_fracs_res    or {}).get(slot_key, 1.0 / 24)
        f_b = (slot_fracs_biz    or {}).get(slot_key, 1.0 / 24)
        f_s = (slot_fracs_school or {}).get(slot_key, 0.0)
        return f_r, f_b, f_s

    n_official = 0
    if official_hourly_file and os.path.exists(official_hourly_file):
        with open(official_hourly_file) as f:
            oh = json.load(f)
        for site_id, site in oh.items():
            node  = site["node"]
            links = [tuple(lnk) for lnk in site["links"]] if site["links"] else None
            m_r, m_b, m_s = _site_flow_3c(flow_res_dict, flow_biz_dict,
                                           flow_school_dict, node, links)
            for obs in site["observations"]:
                dt, h    = obs["time_slot"]
                sk       = (dt, h)
                f_r, f_b, f_s = _fracs(sk)
                count    = float(obs["count"])
                sigma    = float(obs["sigma"])
                pred     = m_r * f_r + m_b * f_b + m_s * f_s
                z        = (pred - count) / sigma if sigma > 0 else 0.0
                chi2    += z ** 2
                lbl      = f"{site['label']} h{h:02d}"
                rows.append(("official", lbl, count, sigma, pred, z))
                n_official += 1

    n_slots_seen = set()
    if link_aadt_file and os.path.exists(link_aadt_file):
        with open(link_aadt_file) as f:
            link_aadt = json.load(f)["links"]
        for key, entry in sorted(link_aadt.items()):
            u, v = map(int, key.split(","))
            if (u, v) in excl:
                continue
            lbl  = label_fn(u, v) if label_fn else f"{u}→{v}"
            m_r  = flow_res_dict.get((u, v), 0.0)
            m_b  = flow_biz_dict.get((u, v), 0.0)
            m_s  = flow_school_dict.get((u, v), 0.0)
            for sess in entry.get("observations", []):
                ts   = sess.get("time_slot")
                neff = float(sess.get("n_eff", 0.5))
                dur  = float(sess.get("duration_s", 0.0))
                if ts is None or dur <= 0:
                    continue
                sk       = (_DOW_TO_TYPE[ts[0]], ts[1])
                f_r, f_b, f_s = _fracs(sk)
                Th       = dur / 3600.0
                pred     = (m_r * f_r + m_b * f_b + m_s * f_s) * Th
                z        = (pred - neff) / math.sqrt(neff) if neff > 0 else 0.0
                n_actual  = neff - 0.5
                pred_safe = max(pred, 1e-30)
                if n_actual > 0:
                    chi2 += 2.0 * (pred_safe - n_actual + n_actual * math.log(n_actual / pred_safe))
                else:
                    chi2 += 2.0 * pred_safe
                n_slots_seen.add(sk)
                m_tot = m_r + m_b + m_s + 1e-30
                f_eff = (m_r * f_r + m_b * f_b + m_s * f_s) / m_tot
                if f_eff > 0 and Th > 0:
                    obs_disp = neff / (Th * f_eff)
                    sig_disp = math.sqrt(neff) / (Th * f_eff)
                    mod_disp = pred / (Th * f_eff)   # = combined AADT (m_r + m_b + m_s)
                else:
                    obs_disp = sig_disp = mod_disp = 0.0
                rows.append(("walking", lbl, obs_disp, sig_disp, mod_disp, z))

    n_obs = len(rows)
    if n_official:
        n_slots_seen |= {(dt, h) for dt in range(3) for h in range(24)}
    # N_eff: 3 df per slot (f_res, f_biz, f_school each consume one df)
    n_eff = n_obs - 3 * len(n_slots_seen)

    rows.sort(key=lambda r: abs(r[5]), reverse=True)
    return rows, chi2, n_obs, n_eff


def _compute_chi2_legacy(link_flow_dict, label_fn, link_aadt_file, exclude_links):
    """Legacy single-component chi² with Woodbury correction (backward compat)."""
    obs_data = []

    for s in COUNT_SITES:
        obs_data.append(("official", s["label"],
                         site_flow(link_flow_dict, s),
                         float(s["observed"]), 0.10 * s["observed"],
                         None, None))

    if link_aadt_file and os.path.exists(link_aadt_file):
        with open(link_aadt_file) as f:
            link_aadt = json.load(f)["links"]
        excl = exclude_links or set()
        for key, entry in sorted(link_aadt.items()):
            u, v = map(int, key.split(","))
            if (u, v) in excl:
                continue
            lbl = label_fn(u, v) if label_fn else f"{u}→{v}"
            mod = link_flow_dict.get((u, v), 0.0)
            for sess in entry.get("observations", []):
                ts  = sess.get("time_slot")
                frs = sess.get("frac_rel_std")
                obs_data.append(("walking", lbl, mod,
                                 float(sess["aadt"]), float(sess["aadt_uncertainty"]),
                                 tuple(ts) if ts is not None else None,
                                 float(frs) if frs is not None else None))

    n_obs = len(obs_data)

    sigma_f    = np.array([(d[3] * d[6]) if d[6] is not None else 0.0 for d in obs_data])
    sigma_sq   = np.array([d[4] ** 2 for d in obs_data])
    sigma_c_sq = np.maximum(sigma_sq - sigma_f ** 2, (np.sqrt(sigma_sq) * 1e-6) ** 2)
    r          = np.array([d[2] - d[3] for d in obs_data])

    slot_groups = {}
    for i, d in enumerate(obs_data):
        if d[5] is not None:
            slot_groups.setdefault(d[5], []).append(i)
    slot_list = list(slot_groups.items())
    n_slots   = len(slot_list)
    n_eff     = n_obs - n_slots

    unslotted = [i for i, d in enumerate(obs_data) if d[5] is None]
    chi2 = float(sum((r[i] / obs_data[i][4]) ** 2 for i in unslotted))
    for ts, idxs in slot_list:
        denom = 1.0 + sum(sigma_f[i] ** 2 / sigma_c_sq[i] for i in idxs)
        uf_r  = sum(sigma_f[i] * r[i] / sigma_c_sq[i] for i in idxs)
        chi2 += float(sum(r[i] ** 2 / sigma_c_sq[i] for i in idxs) - uf_r ** 2 / denom)

    rows = [(d[0], d[1], d[3], d[4], d[2], float(r[i] / d[4]))
            for i, d in enumerate(obs_data)]
    rows.sort(key=lambda row: abs(row[5]), reverse=True)
    return rows, chi2, n_obs, n_eff


def print_chi2_table(rows, chi2, n_obs, n_eff=None, header="Goodness of fit"):
    """Print the standard chi²/N table. Returns chi²/N."""
    chi2_per_n = chi2 / n_obs if n_obs else 0.0
    n_eff_str  = f"  χ²/N_eff={chi2 / n_eff:.4f}  N_eff={n_eff}" if n_eff else ""
    LABEL_W = 52
    print(f"\n{header}  χ²={chi2:.2f}  n={n_obs}  χ²/N={chi2_per_n:.4f}{n_eff_str}")
    print(f"  {'':1s}  {'Src':<8}  {'Label':<{LABEL_W}}  {'Obs':>8}  {'σ':>7}  {'Model':>8}  {'z':>6}")
    for kind, lbl, obs, sig, mod, z in rows:
        marker = "*" if abs(z) > 2 else " "
        print(f"  {marker} {kind:<8}  {lbl:<{LABEL_W}}  {obs:>8,.0f}  {sig:>7,.0f}  {mod:>8,.0f}  {z:>+.2f}")
    abs_z  = [abs(row[5]) for row in rows]
    n_out2 = sum(1 for a in abs_z if a > 2)
    n_out3 = sum(1 for a in abs_z if a > 3)
    mean_z = sum(abs_z) / len(abs_z) if abs_z else 0.0
    print(f"\n  n={n_obs}  χ²/N={chi2_per_n:.4f}  mean|z|={mean_z:.2f}"
          f"  |z|>2: {n_out2}  |z|>3: {n_out3}")
    return chi2_per_n
