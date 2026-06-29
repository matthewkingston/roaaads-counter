"""
Shared constants and functions for the Newtownards gravity model pipeline.

Imported by simulation/build_assignment.py and analysis/tune_assignment.py to
keep their chi²/N calculations, count-site definitions, and gravity kernel
implementations in sync.
"""

import json
import math
import os
import sys
import numpy as np

# The kernel's mode-substitution + speed factors live in analysis/ (shared with the
# NTS derivation tooling).  model.py is in simulation/, so put analysis/ on the path.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis"))
from equiv_miles import equiv_miles      # OSRM seconds → equivalent trip length (miles)
from driveshare import driveshare        # vehicle-driver share by trip length (miles)

# ── Official AADT count sites ─────────────────────────────────────────────────
# AADT totals retained for the quick site-level sanity check in build_assignment.py.
# The tuner and compute_chi2() use per-hour observations from OFFICIAL_HOURLY instead.

COUNT_SITES = [
    {"label": "site 507, A21 Bangor Road",     "node": None, "links": [(538692601,550205936),(550205936,538692601)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road", "node": None, "links": [(544419411,136173596),(136173596,544419411)], "observed": 10_792},
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
INTRA_TIMES     = "data/external_intra_times.json"
GENERATION_RATES = "analysis/generation_rates.json"

_DOW_TO_TYPE = {d: (0 if d < 5 else (1 if d == 5 else 2)) for d in range(7)}

# ── Report-table display helpers ────────────────────────────────────────────────
# Shared by print_chi2_table (model.py), the tuner's live fit table + history obs
# dicts (tune_assignment.py), and report_tune.py, so the per-observation
# disagreement table reads the same everywhere.

DT_TIME_NAMES = {0: "Week", 1: "Sat", 2: "Sun"}


def format_slot_time(day_type, hour):
    """(day_type, hour) → 'Week 1400' / 'Sat 0100' / 'Sun 1900'.  '' if either is None."""
    if day_type is None or hour is None:
        return ""
    return f"{DT_TIME_NAMES.get(day_type, '?')} {int(hour):02d}00"


def nice_official(site):
    """
    From an official_hourly.json / COUNT_SITES site dict, return (label, link):
      label — road name annotated as a count site, e.g. 'A21 Bangor Road (count site)'
              (the leading 'site NNN, ' prefix is stripped).
      link  — precise reference: 'u↔v' from site['links'][0] if present, else the node id.
    """
    raw = site.get("label", "") or ""
    # "site 507, A21 Bangor Road" → "A21 Bangor Road"
    name = raw.split(", ", 1)[1] if raw.startswith("site ") and ", " in raw else raw
    label = f"{name} (count site)" if name else "(count site)"
    links = site.get("links")
    if links:
        u, v = links[0]
        link = f"{u}↔{v}"
    else:
        link = str(site.get("node"))
    return label, link

# ── Paths-cache freshness guard ─────────────────────────────────────────────────
# build_paths.py stamps a signature of its inputs (routing graph, external links,
# the calibrated profile + base speeds) into the .npz. tune_assignment.py and
# build_assignment.py re-check it at load time and fail loudly if the cache is
# stale, rather than silently assigning/tuning against an out-of-date cache (a
# recurring footgun — see CLAUDE.md "Paths cache note").

def _file_sha1(path):
    import hashlib
    if not os.path.exists(path):
        return f"MISSING:{path}"
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def paths_cache_signature():
    """Signature of the inputs build_paths.py consumes. Stamped into the npz at
    build time and re-checked at load time. Returns a dict of {field: str}."""
    from routing_config import PROBIT_CV, PROBIT_LL_SIGMA
    from edge_speed import TUNED_PROFILE, BASE_SPEEDS
    return {
        "src_graph_sha1":       _file_sha1(ROUTING_GRAPH),
        "src_extlinks_sha1":    _file_sha1(EXTERNAL_LINKS),
        "src_profile_sha1":     _file_sha1(TUNED_PROFILE),
        "src_base_speeds_sha1": _file_sha1(BASE_SPEEDS),
        "src_probit_cv":        repr(float(PROBIT_CV)),
        "src_probit_ll_sigma":  repr(float(PROBIT_LL_SIGMA)),
    }

def assert_paths_cache_fresh(cache):
    """Raise SystemExit if the loaded paths cache was built from different inputs
    than the current pipeline state. `cache` is the np.load(...) handle."""
    label = {
        "src_graph_sha1":       "routing graph (newtownards_reduced.graphml)",
        "src_extlinks_sha1":    "external links (data/external_links.json)",
        "src_profile_sha1":     "tuned profile (simulation/tuned_profile.json)",
        "src_base_speeds_sha1": "base speeds (data/google_cache/base_speeds.json)",
        "src_probit_cv":        "PROBIT_CV (simulation/routing_config.py)",
        "src_probit_ll_sigma":  "PROBIT_LL_SIGMA (simulation/routing_config.py)",
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
    """f(d) = (ALPHA+BETA)·u^BETA / (ALPHA + BETA·u^(ALPHA+BETA)), u = d/P.

    Evaluated in the algebraically-identical, overflow-safe form
        f = (ALPHA+BETA) / (ALPHA·u^(-BETA) + BETA·u^ALPHA).
    The denominator is a sum of two non-negative powers of u; for u>0 it is always
    > 0 (minimised at u=1, where it equals ALPHA+BETA ⇒ f=1), and at most ONE term
    can overflow (u^(-BETA)→∞ needs u<1, u^ALPHA→∞ needs u>1). So when a power
    overflows the result is (ALPHA+BETA)/∞ = 0 — the correct tail/origin limit —
    instead of the ∞/∞ = NaN the direct form produces at large ALPHA+BETA (e.g. the
    optimizer pushing ALPHA high). errstate silences the benign overflow warning."""
    u = d / P
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        denom = ALPHA * u ** (-BETA) + BETA * u ** ALPHA
        return (ALPHA + BETA) / denom


def _tanner_kernel(d, P, BETA):
    """Tanner deterrence in the same normalised u=d/P form as _rational_kernel, but
    with an EXPONENTIAL tail instead of a power-law one:
        u = d/P;  f(u) = u^BETA · exp(BETA·(1 − u)).
    Peak f(P)=1 (at u=1), rise ~u^BETA near the origin, tail ~exp(−BETA·d/P) — fast
    enough that long trips are cheap to suppress (well-conditioned), unlike the rational
    kernel's heavy 1/d^ALPHA tail.  Numerically safe: exp(−BETA·u) underflows to 0 long
    before u^BETA could overflow; d=0 ⇒ u=0 ⇒ 0 (BETA>0).  γ = BETA/P (decay scale 1/γ)."""
    u = d / P
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return u ** BETA * np.exp(BETA * (1.0 - u))


def _modesub_kernel(d, TAU):
    """The production-constrained deterrence: an empirical mode-substitution rise ×
    an exponential willingness decay (decoupled — see analysis/driveshare.py and the
    project-tanner-kernel-tld memory):
        f(c) = driveshare(equiv_miles(c)) · exp(−c / TAU),   c in OSRM seconds.
    `driveshare(equiv_miles(c))` is the shared, empirical car-mode share by trip
    length (short trips walked ⇒ f→0 at the origin); `exp(−c/TAU)` is the per-component
    willingness, TAU = 1/γ the characteristic willingness time (seconds).  The
    driveshare PLATEAU is a shared constant that cancels in the production constraint.
    f(0)=0 (equiv_miles→0 ⇒ driveshare 0); finite everywhere."""
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        return driveshare(equiv_miles(d)) * np.exp(-d / TAU)


def constrained_od_flows(od_src, od_dst, od_dist, N_nodes,
                         w_pop, w_workplace, w_retail, w_school,
                         TAU, TAU_commute, TAU_retail,
                         TAU_school=None, with_school=False,
                         self_src=None, self_dist=None, self_w=None,
                         w_commute_prod=None, w_school_prod=None,
                         gen_scale=None):
    """Per-OD-pair, pre-K production-constrained component flows.

    Returns (t_res, t_commute, t_retail, t_sch), each a float64 array parallel to
    od_src/od_dst.  These are the per-pair flows BEFORE the K_res/K_commute/
    K_retail/K_sch scaling; the caller scatters them onto links (all links, or
    observed rows only) and applies K.

    The commute and retail components are independent clones of the school
    component: each is a symmetric two-leg, per-origin-normalised pop↔attractor
    split with its OWN kernel — NO weight parameter and NO self/cross term
    (the old single business component's W_BIZ and biz×biz term are gone):
      commute (modesub kernel, willingness TAU_commute): producer = resident commuters
        (w_commute_prod), attractor = workplace jobs (w_workplace).
      retail  (modesub kernel, willingness TAU_retail):   producer = population,
        attractor = retail parking spaces (w_retail).

    Denominators are summed over the FULL destination set of each origin (all od
    pairs sharing that od_src — internal, external-routed, and denominator-only
    ext→ext virtual edges), so each origin's production budget is conserved.

    External intra-zonal self-term (optional, denominator-only):
      self_src  — node indices i of external zones, each repeated M_i× (one entry
                  per sampled intra-zonal OSRM time).
      self_dist — the sampled intra-zonal times (seconds), parallel to self_src.
      self_w    — per-entry multiplicity weight (1/M_i), so the M_i entries of a
                  zone collectively count as one diagonal destination.
    Adds  a^c_i·(1/M_i)·Σ_m F_c(t_im)  to each per-origin denominator D^c_i — i.e.
    restores the k=i diagonal that collapsing a zone to a centroid dropped (see
    project_production_constrained_gravity, "external intra-zonal self-term").  These
    entries touch ONLY the denominators; they contribute no link flow (no pair_idx).
    self_src=None ⇒ exact prior behaviour.
    """
    src = od_src
    dst = od_dst

    # Per-leg generation scales (vehicle-driver-trips/day pinning, see
    # compute_generation_scales).  Multiply the PRODUCER term of each leg only —
    # the same array's attractor use is scale-invariant (cancels in attr_j/D_i) so
    # the denominators are untouched.  gen_scale=None ⇒ all 1.0 ⇒ exact prior flows.
    gs = gen_scale or {}
    gs_res     = gs.get("res",     1.0)
    gs_com_out = gs.get("com_out", 1.0)
    gs_com_ret = gs.get("com_ret", 1.0)
    gs_ret_out = gs.get("ret_out", 1.0)
    gs_ret_ret = gs.get("ret_ret", 1.0)
    gs_sch_out = gs.get("sch_out", 1.0)
    gs_sch_ret = gs.get("sch_ret", 1.0)

    F_res = _modesub_kernel(od_dist, TAU)
    F_com = _modesub_kernel(od_dist, TAU_commute)
    F_ret = _modesub_kernel(od_dist, TAU_retail)

    pop_s  = w_pop[src];       pop_d  = w_pop[dst]
    work_s = w_workplace[src]; work_d = w_workplace[dst]
    ret_s  = w_retail[src];    ret_d  = w_retail[dst]

    _has_self = self_src is not None and len(self_src) > 0
    if _has_self:
        F_res_self = _modesub_kernel(self_dist, TAU)
        F_com_self = _modesub_kernel(self_dist, TAU_commute)
        F_ret_self = _modesub_kernel(self_dist, TAU_retail)
    else:
        F_res_self = F_com_self = F_ret_self = None

    # Per-origin denominators D^c_i = Σ_k a^c_k·F_c(d_ik) (+ intra-zonal self-term);
    # inverse with 0/0 → 0.  attr_d/F = per-pair destination attraction × kernel;
    # attr_full/F_self = per-node own-zone attraction × self-kernel for the diagonal.
    def _inv_denom(attr_d, F, attr_full, F_self):
        D = np.bincount(src, weights=attr_d * F, minlength=N_nodes)
        if _has_self:
            D += np.bincount(self_src,
                             weights=attr_full[self_src] * F_self * self_w,
                             minlength=N_nodes)
        return np.where(D > 0, 1.0 / D, 0.0)

    iD_res_pop  = _inv_denom(pop_d,  F_res, w_pop,       F_res_self)   # res: attraction = pop
    iD_com_work = _inv_denom(work_d, F_com, w_workplace, F_com_self)   # commute leg home→work: attraction = workplace
    iD_com_pop  = _inv_denom(pop_d,  F_com, w_pop,       F_com_self)   # commute leg work→home: attraction = pop
    iD_ret_ret  = _inv_denom(ret_d,  F_ret, w_retail,    F_ret_self)   # retail leg home→shop:  attraction = retail
    iD_ret_pop  = _inv_denom(pop_d,  F_ret, w_pop,       F_ret_self)   # retail leg shop→home:  attraction = pop

    # res: pop_i·pop_j·F_res / D^res,pop_i  (single leg covers both directions)
    t_res = gs_res * pop_s * pop_d * F_res * iD_res_pop[src]

    # commute: symmetric split, each per-origin-normalised (no weight, no cross term).
    # Home→work producer = resident commuters (node_commute_producers) when supplied,
    # else falls back to population.
    commprod_s = (w_commute_prod[src] if w_commute_prod is not None else pop_s)
    t_commute = F_com * (
        gs_com_out * commprod_s * work_d * iD_com_work[src]   # commuters i → work j  (attraction workplace)
        + gs_com_ret * work_s * pop_d * iD_com_pop[src]       # work i → home j        (attraction pop)
    )

    # retail: symmetric pop↔retail split, each per-origin-normalised (no weight, no cross term).
    t_retail = F_ret * (
        gs_ret_out * pop_s * ret_d * iD_ret_ret[src]          # home i → shop j        (attraction retail)
        + gs_ret_ret * ret_s * pop_d * iD_ret_pop[src]        # shop i → home j        (attraction pop)
    )

    if with_school and w_school is not None and w_school.sum() > 0:
        F_sch = _modesub_kernel(od_dist, TAU_school)
        F_sch_self = _modesub_kernel(self_dist, TAU_school) if _has_self else None
        sch_s = w_school[src]; sch_d = w_school[dst]
        iD_sch_pop = _inv_denom(pop_d, F_sch, w_pop,    F_sch_self)   # school-cross leg school→pop: attraction = pop
        iD_sch_sch = _inv_denom(sch_d, F_sch, w_school, F_sch_self)   # school-cross leg pop→school: attraction = school
        # Home→school producer = resident students (node_school_producers) when supplied,
        # else falls back to population (legacy behaviour).
        schprod_s = (w_school_prod[src] if w_school_prod is not None else pop_s)
        t_sch = F_sch * (
            gs_sch_out * schprod_s * sch_d * iD_sch_sch[src]   # students i → school j  (attraction school)
            + gs_sch_ret * sch_s * pop_d * iD_sch_pop[src]     # school i → pop j        (attraction pop)
        )
    else:
        t_sch = np.zeros(len(src), dtype=np.float64)

    return t_res, t_commute, t_retail, t_sch


def load_generation_rates(path=GENERATION_RATES):
    """Per-component vehicle-driver trips/person/day from analysis/generation_rates.json
    (written by analysis/derive_generation_rates.py).  Returns {commute,retail,school,res}
    or None if the file is absent (⇒ caller skips generation pinning, K_c unpinned)."""
    if not os.path.exists(path):
        print(f"  [gen-rates] {path} not found — generation NOT pinned (run "
              f"analysis/derive_generation_rates.py)")
        return None
    with open(path) as f:
        return json.load(f)["rates"]


def compute_generation_scales(node_weights, rates, verbose=False):
    """Per-leg producer-scale coefficients that put each component's production in
    absolute vehicle-driver trips/day, so the tuned K_c should land at ≈ 1.0.

    node_weights : the loaded node_weights(_reduced).json dict.  The producer layers
        are summed island-wide (the external census nodes tile the whole island), so
        the per-capita anchors recompute automatically for any CENTRE.
    rates        : {commute,retail,school,res} per-person/day (load_generation_rates()).

    Returns the gen_scale dict consumed by constrained_od_flows:
        {res, com_out, com_ret, ret_out, ret_ret, sch_out, sch_ret}.
    Two-leg directions carry ρ_c/2; res is single-leg (full ρ_res, both directions are
    separate OD pairs).  Per-producer rate r = (ρ_c·share)/k with island anchor
    k = Σ(producer layer)/Σ(population) over all nodes (k=1 when producer = population).
    """
    def _sum(layer):
        return float(sum(node_weights.get(layer, {}).values()))

    pop_tot = _sum("node_population")
    if pop_tot <= 0:
        raise ValueError("compute_generation_scales: Σ node_population is zero")

    anchors = {
        "k_commuters": _sum("node_commute_producers") / pop_tot,
        "k_jobs":      _sum("node_workplace")         / pop_tot,
        "k_retail":    _sum("node_retail_spaces")     / pop_tot,
        "k_students":  _sum("node_school_producers")  / pop_tot,
        "k_enrolment": _sum("node_school_demand")     / pop_tot,
    }
    for name, k in anchors.items():
        if k <= 0:
            raise ValueError(f"compute_generation_scales: island anchor {name} ≤ 0 "
                             f"(producer layer missing/empty across all nodes)")

    gen_scale = {
        "res":     rates["res"],                              # single leg, k=1
        "com_out": (rates["commute"] / 2) / anchors["k_commuters"],
        "com_ret": (rates["commute"] / 2) / anchors["k_jobs"],
        "ret_out": (rates["retail"]  / 2),                    # producer = population, k=1
        "ret_ret": (rates["retail"]  / 2) / anchors["k_retail"],
        "sch_out": (rates["school"]  / 2) / anchors["k_students"],
        "sch_ret": (rates["school"]  / 2) / anchors["k_enrolment"],
    }
    if verbose:
        print(f"  [gen-scale] island pop {pop_tot:,.0f}; anchors "
              + ", ".join(f"{n}={v:.4f}" for n, v in anchors.items()))
        print("  [gen-scale] producer rates (trips/day per producer): "
              + ", ".join(f"{k}={v:.4g}" for k, v in gen_scale.items()))
    return gen_scale


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


def load_self_terms(node_ids, intra_times_file=INTRA_TIMES):
    """Build (self_src, self_dist, self_w) for constrained_od_flows from the
    intra-zonal OSRM time samples written by build_intra_times.py.

    For each external zone present BOTH in the intra-times file and in node_ids,
    emits one entry per sampled time: self_src = the zone's node index (repeated
    M_i×), self_dist = the sampled times, self_w = 1/M_i (so a zone's M_i samples
    collectively count as one diagonal destination, contributing mean_m F(t_im)).

    Returns (None, None, None) if the file is absent or yields no usable entries
    (⇒ constrained_od_flows reverts to no self-term).  Zones in the file but absent
    from node_ids are skipped (printed); they cannot be indexed into the weight arrays.
    """
    if not os.path.exists(intra_times_file):
        print(f"  [self-term] {intra_times_file} not found — no intra-zonal self-term")
        return None, None, None
    with open(intra_times_file) as f:
        data = json.load(f)
    data.pop("_meta", None)
    node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
    src, dist, wt = [], [], []
    n_zones = 0
    missing = []
    for zid, times in data.items():
        idx = node_to_idx.get(zid)
        if idx is None:
            missing.append(zid)
            continue
        if not times:
            continue
        m = len(times)
        src.extend([idx] * m)
        dist.extend(times)
        wt.extend([1.0 / m] * m)
        n_zones += 1
    if not src:
        print(f"  [self-term] {intra_times_file} has no zones matching the cache nodes — no self-term")
        return None, None, None
    if missing:
        print(f"  [self-term] {len(missing)} intra-times zones absent from cache node_ids (skipped)")
    print(f"  [self-term] intra-zonal self-term active for {n_zones} external zones "
          f"({len(src)} samples)")
    return (np.asarray(src, dtype=np.intp),
            np.asarray(dist, dtype=np.float64),
            np.asarray(wt, dtype=np.float64))

# ── Flow extraction ───────────────────────────────────────────────────────────

def site_flow(link_flow_dict, site):
    """Return total modelled flow for a COUNT_SITES entry."""
    if site["links"]:
        return sum(link_flow_dict.get(lnk, 0.0) for lnk in site["links"])
    node = site["node"]
    return sum(f for (u, v), f in link_flow_dict.items() if u == node or v == node)


# Number of days each model day_type stands for (weekday=Mon–Fri, Sat, Sun).
_AADT_DAY_WEIGHT = {0: 5, 1: 1, 2: 1}


def aadt_weights(slot_fracs_res, slot_fracs_commute, slot_fracs_retail,
                 slot_fracs_school=None):
    """Per-component AADT weights for converting a pre-K component flow to a daily total.

    The annual-average daily contribution of component c is K_c·m_c·W_c with
        W_c = (5·Σ_h f_c[weekday,h] + Σ_h f_c[Sat,h] + Σ_h f_c[Sun,h]) / 7,
    the day-type-weighted sum of the component's hourly fractions.

    The temporal profiles (derive_component_profiles.py) are now per-component SHAPES,
    each normalised so **W_c ≈ 1** (decoupled from magnitude, which generation pins) —
    so K_c·m_c is itself ≈ the component's daily AADT, and the combined daily AADT is
    Σ_c K_c·m_c·W_c.  (Pre-decoupling the components partitioned the aggregate profile,
    so the four W_c summed to ≈1 instead of each being ≈1.)

    Returns (W_res, W_commute, W_retail, W_sch).  Components with no slot_fracs
    return 0.0.
    """
    def _w(sf):
        if not sf:
            return 0.0
        return sum(_AADT_DAY_WEIGHT.get(dt, 0) * f for (dt, h), f in sf.items()) / 7.0
    return (_w(slot_fracs_res), _w(slot_fracs_commute),
            _w(slot_fracs_retail), _w(slot_fracs_school))


def _site_flow_components(flow_dicts, node, links):
    """Return one modelled flow per component dict for a site (node or links).

    flow_dicts: list of {(u,v): flow} dicts.  Returns a list of the same length.
    """
    out = []
    for fd in flow_dicts:
        if links:
            out.append(sum(fd.get(tuple(lnk), 0.0) for lnk in links))
        else:
            out.append(sum(f for (u, v), f in fd.items() if u == node or v == node))
    return out

# ── Chi²/N ───────────────────────────────────────────────────────────────────

def compute_chi2(link_flow_dict, label_fn=None,
                 link_aadt_file=LINK_AADT, exclude_links=EXCLUDE_LINKS,
                 official_hourly_file=OFFICIAL_HOURLY,
                 link_flow_commute_dict=None, link_flow_retail_dict=None,
                 link_flow_school_dict=None,
                 slot_fracs_res=None, slot_fracs_commute=None,
                 slot_fracs_retail=None, slot_fracs_school=None):
    """
    Compute chi²/N matching the tuner's component formulation.

    Four-component mode (link_flow_commute_dict and link_flow_retail_dict provided):
      link_flow_dict         — {(u,v): K_res * flow_res}
      link_flow_commute_dict — {(u,v): K_commute * flow_commute}
      link_flow_retail_dict  — {(u,v): K_retail  * flow_retail}
      link_flow_school_dict  — {(u,v): K_school  * flow_school}  (optional)
      slot_fracs_res/commute/retail/school — {(day_type, hour): f}
      N_eff = N (temporal fractions pinned at NTS — no per-slot df consumed).

    Legacy mode (link_flow_commute_dict=None):
      link_flow_dict — {(u,v): K * combined_flow}  (pre-scaled)
      Uses 3 AADT-space obs from COUNT_SITES plus Woodbury-corrected walking obs.

    Returns (rows, chi2, n_obs, n_eff).
      rows: list of (label, time_str, obs, sig, mod, z, link) sorted by |z| descending.
    """
    if link_flow_commute_dict is not None and link_flow_retail_dict is not None:
        comps = [(link_flow_dict,         slot_fracs_res),
                 (link_flow_commute_dict, slot_fracs_commute),
                 (link_flow_retail_dict,  slot_fracs_retail)]
        if link_flow_school_dict is not None:
            comps.append((link_flow_school_dict, slot_fracs_school))
        return _compute_chi2_components(comps, label_fn, link_aadt_file,
                                        exclude_links, official_hourly_file)
    return _compute_chi2_legacy(link_flow_dict, label_fn, link_aadt_file, exclude_links)


def walking_session_residual(comps, sess):
    """Per-session walking-count residual for the production-constrained model.

    Single source of truth for the per-observation arithmetic used by both
    _compute_chi2_components and the map's residuals layer.

    comps: iterable of (m_c, slot_fracs_c) pairs — one per active gravity
    component (m_c = that component's K-scaled modelled flow on this link,
    slot_fracs_c = {(day_type, hour): f}).  The combined hourly prediction in a
    slot is Σ_c m_c·f_c[slot]; the daily combined AADT is Σ_c m_c.

    Returns None if the session is unusable (no time slot or non-positive
    duration), else a dict:
      sk       — (day_type, hour) slot key
      z        — Poisson Pearson residual (pred_count − n_eff)/√n_eff
      deviance — Poisson deviance contribution to chi²
      obs_disp — observed effective AADT (display)
      sig_disp — observed AADT uncertainty (display)
      mod_disp — model combined AADT (= Σ_c m_c)
    """
    ts   = sess.get("time_slot")
    neff = float(sess.get("n_eff", 0.5))
    dur  = float(sess.get("duration_s", 0.0))
    if ts is None or dur <= 0:
        return None
    sk  = (_DOW_TO_TYPE[ts[0]], ts[1])
    rate  = 0.0   # Σ_c m_c·f_c[slot]  (combined hourly rate)
    m_tot = 0.0   # Σ_c m_c            (combined daily AADT)
    for m_c, sf_c in comps:
        f_c = (sf_c or {}).get(sk, 1.0 / 24)
        rate  += m_c * f_c
        m_tot += m_c
    Th        = dur / 3600.0
    pred      = rate * Th
    z         = (pred - neff) / math.sqrt(neff) if neff > 0 else 0.0
    n_actual  = neff - 0.5
    pred_safe = max(pred, 1e-30)
    if n_actual > 0:
        deviance = 2.0 * (pred_safe - n_actual + n_actual * math.log(n_actual / pred_safe))
    else:
        deviance = 2.0 * pred_safe
    f_eff = rate / (m_tot + 1e-30)
    if f_eff > 0 and Th > 0:
        obs_disp = neff / (Th * f_eff)
        sig_disp = math.sqrt(neff) / (Th * f_eff)
        mod_disp = pred / (Th * f_eff)   # = combined AADT (Σ_c m_c)
    else:
        obs_disp = sig_disp = mod_disp = 0.0
    return {"sk": sk, "z": z, "deviance": deviance,
            "obs_disp": obs_disp, "sig_disp": sig_disp, "mod_disp": mod_disp}


def _compute_chi2_components(comps, label_fn, link_aadt_file, exclude_links,
                             official_hourly_file):
    """Production-constrained chi² over an arbitrary component list.

    comps: list of (flow_dict, slot_fracs) — one per active component
    (res, commute, retail, school).  flow_dict = {(u,v): K_c·flow_c}.

    Official obs use the count-space Gaussian z; walking obs use the Poisson
    deviance via walking_session_residual.  N_eff = N (temporal fractions are
    pinned at NTS, so no per-slot df are consumed).
    """
    rows  = []
    chi2  = 0.0
    excl  = exclude_links or set()
    flow_dicts = [fd for fd, _ in comps]

    # ── Official hourly obs (count-space Gaussian) ────────────────────────────
    n_official = 0
    if official_hourly_file and os.path.exists(official_hourly_file):
        with open(official_hourly_file) as f:
            oh = json.load(f)
        for site_id, site in oh.items():
            node  = site["node"]
            links = [tuple(lnk) for lnk in site["links"]] if site["links"] else None
            m_c   = _site_flow_components(flow_dicts, node, links)
            lbl, link = nice_official(site)
            for obs in site["observations"]:
                dt, h = obs["time_slot"]
                sk    = (dt, h)
                # T=3600 → T/3600=1; pred = Σ_c m_c·f_c[slot]
                pred  = sum(m * ((sf or {}).get(sk, 1.0 / 24))
                            for m, (_, sf) in zip(m_c, comps))
                count = float(obs["count"])
                sigma = float(obs["sigma"])
                z     = (pred - count) / sigma if sigma > 0 else 0.0
                chi2 += z ** 2
                rows.append((lbl, format_slot_time(dt, h), count, sigma, pred, z, link))
                n_official += 1

    # ── Walking obs (count-space Poisson) ─────────────────────────────────────
    if link_aadt_file and os.path.exists(link_aadt_file):
        with open(link_aadt_file) as f:
            link_aadt = json.load(f)["links"]
        for key, entry in sorted(link_aadt.items()):
            u, v = map(int, key.split(","))
            if (u, v) in excl:
                continue
            lbl  = (label_fn(u, v) if label_fn else "") or "(unnamed)"
            link = f"{u}→{v}"
            sess_comps = [(fd.get((u, v), 0.0), sf) for fd, sf in comps]
            for sess in entry.get("observations", []):
                res = walking_session_residual(sess_comps, sess)
                if res is None:
                    continue
                chi2 += res["deviance"]
                rows.append((lbl, format_slot_time(*res["sk"]), res["obs_disp"],
                             res["sig_disp"], res["mod_disp"], res["z"], link))

    n_obs = len(rows)
    n_eff = n_obs   # temporal fractions pinned at NTS — no per-slot df consumed
    rows.sort(key=lambda r: abs(r[5]), reverse=True)
    return rows, chi2, n_obs, n_eff


def _compute_chi2_legacy(link_flow_dict, label_fn, link_aadt_file, exclude_links):
    """Legacy single-component chi² with Woodbury correction (backward compat)."""
    obs_data = []

    for s in COUNT_SITES:
        lbl, link = nice_official(s)
        obs_data.append(("official", lbl,
                         site_flow(link_flow_dict, s),
                         float(s["observed"]), 0.10 * s["observed"],
                         None, None, link))

    if link_aadt_file and os.path.exists(link_aadt_file):
        with open(link_aadt_file) as f:
            link_aadt = json.load(f)["links"]
        excl = exclude_links or set()
        for key, entry in sorted(link_aadt.items()):
            u, v = map(int, key.split(","))
            if (u, v) in excl:
                continue
            lbl  = (label_fn(u, v) if label_fn else "") or "(unnamed)"
            link = f"{u}→{v}"
            mod = link_flow_dict.get((u, v), 0.0)
            for sess in entry.get("observations", []):
                ts  = sess.get("time_slot")
                frs = sess.get("frac_rel_std")
                obs_data.append(("walking", lbl, mod,
                                 float(sess["aadt"]), float(sess["aadt_uncertainty"]),
                                 tuple(ts) if ts is not None else None,
                                 float(frs) if frs is not None else None, link))

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

    def _time(d):
        ts = d[5]
        return format_slot_time(_DOW_TO_TYPE[ts[0]], ts[1]) if ts is not None else ""

    rows = [(d[1], _time(d), d[3], d[4], d[2], float(r[i] / d[4]), d[7])
            for i, d in enumerate(obs_data)]
    rows.sort(key=lambda row: abs(row[5]), reverse=True)
    return rows, chi2, n_obs, n_eff


def print_chi2_table(rows, chi2, n_obs, n_eff=None, header="Goodness of fit"):
    """Print the standard chi²/N table. Returns chi²/N."""
    chi2_per_n = chi2 / n_obs if n_obs else 0.0
    n_eff_str  = f"  χ²/N_eff={chi2 / n_eff:.4f}  N_eff={n_eff}" if n_eff else ""
    LABEL_W = 30
    print(f"\n{header}  χ²={chi2:.2f}  n={n_obs}  χ²/N={chi2_per_n:.4f}{n_eff_str}")
    print(f"  {'':1s}  {'Label':<{LABEL_W}}  {'Time':<9}  {'Obs':>8}  {'σ':>7}  {'Model':>8}  {'z':>6}  {'Link'}")
    for label, time_str, obs, sig, mod, z, link in rows:
        marker = "*" if abs(z) > 2 else " "
        lbl = label if len(label) <= LABEL_W else label[:LABEL_W - 1] + "…"
        print(f"  {marker} {lbl:<{LABEL_W}}  {time_str:<9}  {obs:>8,.0f}  {sig:>7,.0f}  {mod:>8,.0f}  {z:>+.2f}  {link}")
    abs_z  = [abs(row[5]) for row in rows]
    n_out2 = sum(1 for a in abs_z if a > 2)
    n_out3 = sum(1 for a in abs_z if a > 3)
    mean_z = sum(abs_z) / len(abs_z) if abs_z else 0.0
    print(f"\n  n={n_obs}  χ²/N={chi2_per_n:.4f}  mean|z|={mean_z:.2f}"
          f"  |z|>2: {n_out2}  |z|>3: {n_out3}")
    return chi2_per_n
