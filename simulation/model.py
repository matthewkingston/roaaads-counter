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
from collections import defaultdict
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
    # site 444 (A20 Portaferry Road) removed: the ODS count point is artificially misplaced (it sits
    # ~1 km into the core at node 449111329, not its true location much further down the peninsula),
    # so it under-reads relative to where the peninsula traffic actually is — it was dominating and
    # distorting the calibration (it "wanted" K≈0.66, pulling the fit down). Dropped as a bad obs.
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
MOBILISATION_FILE = "analysis/mobilisation.json"

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

# ── Production-constrained gravity kernel + assignment ────────────────────────
# Singly-constrained per component:  T^c_ij = K_c · p^c_i · a^c_j · F_c(d_ij) / D^c_i,
# D^c_i = Σ_k a^c_k·F_c(d_ik), so Σ_j T^c_ij = K_c·p^c_i — each origin's trip production is
# fixed by its producing weight, independent of accessibility (fixes the generation/
# distribution conflation; see project_production_constrained_gravity).  The per-origin
# denominators depend on the kernel, so are recomputed every evaluation (cheap O(N_OD)
# bincounts); K_c is applied by the caller.  Kernel = _modesub_kernel (below).


def _modesub_kernel(d, wparams, component):
    """The production-constrained deterrence: an empirical mode-substitution rise ×
    a double-exponential willingness decay (decoupled — see analysis/driveshare.py,
    analysis/fit_kernel.py and the project-tanner-kernel-tld memory):
        f(c) = driveshare(equiv_miles(c), component) · [ w·exp(−c/τs) + (1−w)·exp(−c/τl) ],
    c in OSRM seconds.  `driveshare(equiv_miles(c), component)` is that component's
    empirical car-mode share by trip length (per-component — the short-range walk↔drive
    substitution differs by purpose; short trips walked ⇒ f→0 at the origin).
    `wparams = (w, τs, τl)` is the per-component willingness: a fast head (weight w, scale
    τs) + a heavier tail (weight 1−w, scale τl>τs), both seconds — the TLD/n_Ire divide
    (fit_kernel) showed a single exponential is too light-tailed.  W(0)=w+(1−w)=1, and the
    willingness amplitude is absorbed by K in the production constraint, so only the shape
    (w, τs, τl) is load-bearing (as the single-exp τ was before).  The driveshare PLATEAU_c
    is likewise a constant that cancels within each component.  `component` is required — one
    of driveshare.CURVES (res/commute/retail, or school_primary/postprimary/tertiary).
    f(0)=0 (equiv_miles→0 ⇒ driveshare 0); finite everywhere."""
    w, tau_s, tau_l = wparams
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        W = w * np.exp(-d / tau_s) + (1.0 - w) * np.exp(-d / tau_l)
        return driveshare(equiv_miles(d), component) * W


SCHOOL_LEVELS = ("primary", "postprimary", "tertiary")   # the three independent school components

# The six gravity components, each carrying its own double-exp willingness kernel
# (res/commute/retail + the three independent school levels).  These strings are both the
# `component` arg to _modesub_kernel/driveshare and the keys of the `willingness` dict.
WILLINGNESS_COMPONENTS = ("res", "commute", "retail",
                          "school_primary", "school_postprimary", "school_tertiary")
_WPARAM_SUFFIXES = ("taus", "taul", "w")   # flat-key order: "<comp>_taus/_taul/_w"


def willingness_keys():
    """The 18 flat willingness param keys in canonical order (tuned_params/tuner_config)."""
    return [f"{c}_{s}" for c in WILLINGNESS_COMPONENTS for s in _WPARAM_SUFFIXES]


def willingness_from_flat(flat):
    """Build the {component: (w, τs, τl)} dict `constrained_od_flows` expects, from a flat
    {"<comp>_taus"/"_taul"/"_w": value} dict (natural units).  Only components whose keys are
    present are built (so a no-school 9-key dict yields 3 components).  **Enforces τl ≥ τs**
    (the fast head is never slower than the tail): τl < τs collapses to a single exponential at
    τs — the sensible limit — which also removes the τs↔τl label-swap degeneracy."""
    out = {}
    for c in WILLINGNESS_COMPONENTS:
        if f"{c}_w" in flat:
            w  = float(flat[f"{c}_w"])
            ts = float(flat[f"{c}_taus"])
            tl = float(flat[f"{c}_taul"])
            out[c] = (w, ts, max(tl, ts))
    return out


def constrained_od_flows(od_src, od_dst, od_dist, N_nodes,
                         w_pop, w_workplace, w_retail,
                         willingness,
                         with_school=False,
                         w_school_levels=None, w_school_prod_levels=None,
                         self_terms=None,
                         w_commute_prod=None,
                         gen_scale=None,
                         doubly_constrained=None,
                         furness_max_sweeps=None,
                         furness_state=None,
                         mu=None,
                         return_legs=False):
    """Per-OD-pair, pre-K production-constrained component flows.

    Returns (t_res, t_commute, t_retail, t_sch_by_level), where the first three are
    float64 arrays parallel to od_src/od_dst and t_sch_by_level is a dict
    {level: array} over SCHOOL_LEVELS (primary/post-primary/tertiary).  These are the
    per-pair flows BEFORE the K_res/K_commute/K_retail/K_<level> scaling; the caller
    scatters them onto links (all links, or observed rows only) and applies K.

    If return_legs=True, additionally returns a dict `legs` of the individual
    producer→attractor legs that make up each summed component (keys: "res",
    "commute_out"/"commute_ret", "retail_out"/"retail_ret", and per active level
    "school_<lvl>_out"/"school_<lvl>_ret").  Each out-leg is producer→attractor
    (home→activity, home = origin); each ret-leg is attractor→producer
    (activity→home, home = destination).  Used by diagnostics that need to anchor
    trips to the home end; the default (return_legs=False) path is unchanged.

    `willingness` is a dict {component: (w, τs, τl)} of per-component double-exp
    willingness params (6 entries — res/commute/retail + school_primary/postprimary/
    tertiary; each school level fully independent, no shared τ_school), consumed by
    _modesub_kernel.  The commute and retail components are independent clones of the
    school component: each is a symmetric two-leg producer↔attractor round-trip, each
    leg per-origin-normalised, with its OWN kernel — NO weight parameter and NO
    self/cross term (the old single business component's W_BIZ and biz×biz term are gone).
    Both legs share the SAME two layers: the outbound is producer→attractor and the
    return is attractor→producer — the return home is attracted by the PRODUCER layer,
    NOT raw population, so returning commuters/students land where their producers live
    (a residential zone's evening inflow is distributed by its resident-commuter /
    resident-student count, not its total population):
      commute (modesub kernel, willingness["commute"]): producer = resident commuters
        (w_commute_prod), attractor = workplace jobs (w_workplace).
      retail  (modesub kernel, willingness["retail"]):  producer = population,
        attractor = retail parking spaces (w_retail).  (Here producer IS population,
        so the pop-attracted return is already the producer-attracted return.)

    Denominators are summed over the FULL destination set of each origin (all od
    pairs sharing that od_src — internal, external-routed, and denominator-only
    ext→ext virtual edges), so each origin's production budget is conserved.

    External intra-zonal self-term (optional, denominator-only):
      self_terms — {component: (self_src, self_dist, self_w)} from model.load_self_terms
                  (mass-weighted per-component intra-zonal time histograms, build_intra_times.py).
                  Per component: self_src = external-zone node indices (one per bin), self_dist =
                  bin centre times (s), self_w = bin weights (Σ=1 per zone).  Both legs of a
                  component share its entry (the p×a interaction is symmetric).
    Adds  a^c_i·Σ_bin w·F_c(t_bin) = a^c_i·S^c_i  to each per-origin denominator D^c_i — i.e.
    restores the k=i diagonal that collapsing a zone to a centroid dropped (see
    project_production_constrained_gravity, "external intra-zonal self-term").  These
    entries touch ONLY the denominators; they contribute no link flow (no pair_idx).
    self_terms=None ⇒ no self-term.

    Doubly-constrained (Furness) option:
      doubly_constrained — a set/list of component names ({"commute","retail",
                  "school_primary","school_postprimary","school_tertiary"}) whose legs
                  are ALSO attraction-constrained.  For a flagged component every leg is
                  balanced (Furness) so BOTH margins hold: Σ_j T_ij (+ self diagonal) =
                  gen-scaled producer_i (production, the absolute magnitude anchor) AND
                  Σ_i T_ij (+ self diagonal) ∝ attractor_j (attraction — the attractor's
                  raw scale is normalised away, only its cross-zone proportions enter).
                  Unflagged components (and residential, always) keep the singly
                  (production) constrained  gs·p_i·a_j·F/D_i.  Flow stays linear in K_c
                  (balancing factors normalise to the raw margins, not K), so the caller's
                  convex K-solve is unchanged.  None/empty ⇒ exact singly-constrained
                  behaviour.
      furness_max_sweeps — approximate-balancing budget for warm-started legs.  Plain IPF
                  converges pathologically slowly on the real short-range kernels (~1000+
                  iters/leg), so for tuning the balancing is run as a FIXED number of
                  warm-started sweeps instead of to tolerance: a leg with a cached warm
                  start (see furness_state) runs exactly `furness_max_sweeps` sweeps and
                  ends on a row-normalisation, so PRODUCTION stays exact and only the
                  attraction margin is approximate (<1% at k≈10, well under count noise).
                  A COLD leg (no warm start) always converges to tolerance to seed the
                  cache.  None ⇒ every leg converges to tolerance (exact; used by
                  build_assignment for the deployed flows).
      furness_state — a mutable {leg_key: b} dict the caller keeps ACROSS evals so each
                  leg's balancing factors warm-start from the previous eval (b drifts only
                  ~1% per tuner step, so k≈10 sweeps keep it current).  None ⇒ no cache
                  (every call cold).
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

    F_res = _modesub_kernel(od_dist, willingness["res"], "res")
    F_com = _modesub_kernel(od_dist, willingness["commute"], "commute")
    F_ret = _modesub_kernel(od_dist, willingness["retail"], "retail")

    # Per-component self-term arrays (self_src, F_self, self_w) or None, from the mass-weighted
    # intra-zonal histograms; the kernel is applied to the stored bin-centre times here.
    def _self(comp, wparams):
        if not self_terms or comp not in self_terms:
            return None
        s_src, s_dist, s_w = self_terms[comp]
        return (s_src, _modesub_kernel(s_dist, wparams, comp), s_w)
    st_res = _self("res", willingness["res"])
    st_com = _self("commute", willingness["commute"])
    st_ret = _self("retail", willingness["retail"])

    # Nodes appearing as an origin / destination in the pair set (the reachable support),
    # used to normalise the doubly-constrained attraction margin to the production total.
    src_present = np.zeros(N_nodes, dtype=bool); src_present[src] = True
    dst_present = np.zeros(N_nodes, dtype=bool); dst_present[dst] = True
    dbl = set(doubly_constrained) if doubly_constrained else set()

    # Per-origin denominators D^c_i = Σ_k a^c_k·F_c(d_ik) (+ intra-zonal self-term);
    # inverse with 0/0 → 0.  attr_d/F = per-pair destination attraction × kernel;
    # st = (self_src, F_self, self_w): the diagonal a^c_i·Σ_bin w·F_c(t_bin).
    def _inv_denom(attr_d, F, attr_full, st):
        D = np.bincount(src, weights=attr_d * F, minlength=N_nodes)
        if st is not None:
            s_src, F_self, s_w = st
            D += np.bincount(s_src, weights=attr_full[s_src] * F_self * s_w, minlength=N_nodes)
        return np.where(D > 0, 1.0 / D, 0.0)

    def _furness(O_full, A_full, F, st, comp, leg_key, max_iter=3000, tol=1e-9):
        """Doubly-constrained (Furness) per-pair flow for one leg.

        O_full = per-node generator (gen-scaled producer; the absolute magnitude anchor).
        A_full = per-node attractor weight, used PROPORTIONALLY — its raw scale is
        normalised away so ΣD == ΣO over the reachable support (both margins consistent).
        The fixed point pins Σ_j T_ij (+ self diagonal) = O_i (production, the hard
        constraint) and Σ_i T_ij (+ self diagonal) = D_j ∝ A_j (attraction).  `st`
        (self_src, F_self, self_w) restores the intra-zonal diagonal to BOTH balancing
        sums (the p·a·f interaction is symmetric, so one histogram serves both margins);
        it is denominator-only (no link flow) exactly as in the singly-constrained path.

        Two modes (see furness_max_sweeps / furness_state on constrained_od_flows):
          • COLD (no cached warm start, or furness_max_sweeps is None): iterate up to
            `max_iter` on the PRODUCTION (row) residual — the load-bearing margin — to seed
            the cache.  Short-kernel legs (school τs≈90–100 s) mix too slowly to reach `tol`;
            when the cap is hit it WARNS (loudly, once) and proceeds with the best-effort b
            — this is the documented approximate-balancing regime, not a failure, since the
            final row-normalisation makes production exact regardless.
          • WARM (cached b + furness_max_sweeps set): run exactly `furness_max_sweeps`
            sweeps from the cached b.  b drifts ~1% per tuner step, so this stays near the
            fixed point; the deliberate approximation lives ENTIRELY in the attraction
            margin (a few tenths of a % at k≈10).
        BOTH modes end on a final row-normalisation, so PRODUCTION is exact either way and
        only the attraction margin is (slightly) approximate in the warm mode.  The final
        b is written back to furness_state[leg_key] for the next eval's warm start."""
        sumO = float(O_full[src_present].sum())
        sumA = float(A_full[dst_present].sum())
        if sumO <= 0.0 or sumA <= 0.0:
            return np.zeros(len(src), dtype=np.float64)
        D_full = A_full * (sumO / sumA)                 # ΣD == ΣO over the reachable support
        has_self = st is not None
        if has_self:
            s_src, F_self, s_w = st

        def _dena(bvec):                                # Σ_j b_j D_j F_ij (+ self diagonal)
            bd = bvec * D_full
            d = np.bincount(src, weights=bd[dst] * F, minlength=N_nodes)
            if has_self:
                d += np.bincount(s_src, weights=bd[s_src] * F_self * s_w, minlength=N_nodes)
            return d

        def _denb(avec):                                # Σ_i a_i O_i F_ij (+ self diagonal)
            aO = avec * O_full
            d = np.bincount(dst, weights=aO[src] * F, minlength=N_nodes)
            if has_self:
                d += np.bincount(s_src, weights=aO[s_src] * F_self * s_w, minlength=N_nodes)
            return d

        prod_mask = O_full > 0
        O_mean = sumO / max(int(prod_mask.sum()), 1)
        b0 = furness_state.get(leg_key) if furness_state is not None else None
        warm = furness_max_sweeps is not None and b0 is not None
        b = (b0.copy() if b0 is not None else np.ones(N_nodes, dtype=np.float64))
        denom_a = _dena(b)
        rel = np.inf
        if warm:                                        # fixed k sweeps from the cached b
            for it in range(1, int(furness_max_sweeps) + 1):
                a = np.where(denom_a > 0, 1.0 / denom_a, 0.0)
                b = np.where((db := _denb(a)) > 0, 1.0 / db, 0.0)
                denom_a = _dena(b)
        else:                                           # cold: iterate up to max_iter to seed the cache
            for it in range(1, max_iter + 1):
                a = np.where(denom_a > 0, 1.0 / denom_a, 0.0)
                b = np.where((db := _denb(a)) > 0, 1.0 / db, 0.0)
                denom_a = _dena(b)                      # for the next iter AND the residual
                rel = float(np.abs(a * O_full * denom_a - O_full)[prod_mask].max() / O_mean) \
                    if prod_mask.any() else 0.0
                if rel < tol:
                    break
            else:
                # Short-kernel legs (school especially) mix so slowly IPF cannot reach `tol`
                # in max_iter — this is the documented approximate-balancing regime, NOT a
                # failure: production is made exact by the final row-normalisation below, and
                # `rel` bounds the residual attraction error.  Warn loudly (once, at the cold
                # seed) rather than crash, and proceed with the best-effort b.
                print(f"  [furness {comp}/{leg_key}] IPF capped at {max_iter} iters — "
                      f"production residual {rel:.1e} (short kernel; production stays exact, "
                      f"attraction ≈{rel:.0e}). Approximate balancing.")
        a = np.where(denom_a > 0, 1.0 / denom_a, 0.0)   # final row-normalisation ⇒ production exact
        if furness_state is not None:
            furness_state[leg_key] = b                  # warm start for the next eval
        if os.environ.get("MODEL_FURNESS_DEBUG"):
            mode = f"warm {it} sweeps" if warm else f"cold {it} iters (resid {rel:.1e})"
            print(f"    [furness {comp}/{leg_key}] {mode}")
        return a[src] * O_full[src] * b[dst] * D_full[dst] * F

    def _leg(prod_full, gs_c, attr_full, F, st, comp, leg_key):
        """One producer→attractor leg: doubly-constrained (Furness) if `comp` is flagged,
        else the singly (production) constrained  gs·p_i·a_j·F/D_i."""
        if comp in dbl:
            return _furness(gs_c * prod_full, attr_full, F, st, comp, leg_key)
        iD = _inv_denom(attr_full[dst], F, attr_full, st)
        return gs_c * prod_full[src] * attr_full[dst] * F * iD[src]

    def _muprod(key, prod):
        """Per-area car-ownership multiplier on a producer array (home-end legs only).
        `mu` is {component: per-node array} normalised (by the caller) to producer-weighted
        mean 1, so it redistributes production spatially without moving the component total.
        `mu=None` (or a missing key) ⇒ prod unchanged ⇒ bit-identical to the pre-μ path."""
        if mu is None:
            return prod
        m = mu.get(key)
        return prod if m is None else m * prod

    legs = {}   # populated for return_legs=True (the individual producer→attractor legs)

    # res: single leg covers both directions (pop↔pop, symmetric).  Held singly-constrained.
    # Producer (home end) carries the car-ownership multiplier; the pop attractor does not.
    t_res = _leg(_muprod("res", w_pop), gs_res, w_pop, F_res, st_res, "res", "res")
    legs["res"] = t_res

    # commute: symmetric producer↔attractor round-trip.  Out leg home→work (producer =
    # resident commuters, attractor = jobs); return leg work→home (producer = jobs,
    # attractor = resident commuters — returning commuters land where commuters live).
    _com_out = _leg(w_commute_prod, gs_com_out, w_workplace,     F_com, st_com, "commute", "commute_out")
    _com_ret = _leg(w_workplace,    gs_com_ret, w_commute_prod, F_com, st_com, "commute", "commute_ret")
    t_commute = _com_out + _com_ret
    legs["commute_out"] = _com_out
    legs["commute_ret"] = _com_ret

    # retail: symmetric pop↔retail round-trip (out home→shop, return shop→home).
    # Only the outbound (home-end) producer carries μ; the return producer is the shop.
    _ret_out = _leg(_muprod("retail", w_pop), gs_ret_out, w_retail, F_ret, st_ret, "retail", "retail_out")
    _ret_ret = _leg(w_retail, gs_ret_ret, w_pop,    F_ret, st_ret, "retail", "retail_ret")
    t_retail = _ret_out + _ret_ret
    legs["retail_out"] = _ret_out
    legs["retail_ret"] = _ret_ret

    # School: three INDEPENDENT components (primary / post-primary / tertiary), each a symmetric
    # two-leg producer↔attractor round-trip with its OWN producer, attractor, generation scale, K,
    # per-level driveshare curve AND its OWN double-exp willingness (willingness["school_<lvl>"]) —
    # the levels are fully independent (no shared τ_school).  Out leg = students→school (attraction
    # = this level's enrolment); return leg = school→home (attraction = this level's resident
    # students, NOT raw population — returning students land where that level's producers live).
    # Each leg is doubly-constrained iff "school_<lvl>" is in doubly_constrained.
    t_sch_by_level = {lvl: np.zeros(len(src), dtype=np.float64) for lvl in SCHOOL_LEVELS}
    if with_school and w_school_levels:
        for lvl in SCHOOL_LEVELS:
            w_sch = w_school_levels.get(lvl)
            if w_sch is None or w_sch.sum() <= 0:
                continue
            comp = f"school_{lvl}"
            w_sch_params = willingness[comp]
            F_sch = _modesub_kernel(od_dist, w_sch_params, comp)
            st_sch = _self(comp, w_sch_params)
            prod = w_school_prod_levels[lvl]                           # resident students of this level
            gs_out = gs.get(f"sch_{lvl}_out", 1.0)
            gs_ret = gs.get(f"sch_{lvl}_ret", 1.0)
            # Only the outbound (home-end, students) producer carries μ; the school does not.
            _sch_out = _leg(_muprod(comp, prod), gs_out, w_sch, F_sch, st_sch, comp, f"{comp}_out")   # students i → school j
            _sch_ret = _leg(w_sch, gs_ret, prod,  F_sch, st_sch, comp, f"{comp}_ret")   # school i → home j
            t_sch_by_level[lvl] = _sch_out + _sch_ret
            legs[f"{comp}_out"] = _sch_out
            legs[f"{comp}_ret"] = _sch_ret

    if return_legs:
        return t_res, t_commute, t_retail, t_sch_by_level, legs
    return t_res, t_commute, t_retail, t_sch_by_level


def load_generation_rates(path=GENERATION_RATES, mobilisation_path=MOBILISATION_FILE):
    """Per-component vehicle-driver trips/person/day from analysis/generation_rates.json
    (written by analysis/derive_generation_rates.py).  Returns {commute, retail, res,
    school_primary, school_postprimary, school_tertiary} or None if the file is absent
    (⇒ caller skips generation pinning, K_c unpinned).

    The England-NTS rates are rescaled to the ISLAND car-driver mobilisation level by a
    single global multiplier m_island (analysis/mobilisation.json, from
    analysis/derive_mobilisation.py), preserving the NTS purpose split, so each K_c
    anchors cleanly at 1 against island-level generation rather than England's.  This is
    the common-mode (level) half of the properly-derived K normalisation; the split and
    spatial-dispersion widths live in the K-prior (analysis/tune_assignment.py).  If
    mobilisation.json is absent, m_island=1 (rates left at the England level) with a note."""
    if not os.path.exists(path):
        print(f"  [gen-rates] {path} not found — generation NOT pinned (run "
              f"analysis/derive_generation_rates.py)")
        return None
    with open(path) as f:
        rates = json.load(f)["rates"]
    m_island = 1.0
    if os.path.exists(mobilisation_path):
        with open(mobilisation_path) as f:
            m_island = float(json.load(f).get("m_island", 1.0))
    else:
        print(f"  [gen-rates] {mobilisation_path} not found — m_island=1 (rates NOT "
              f"rescaled to island level; run analysis/derive_mobilisation.py)")
    return {c: v * m_island for c, v in rates.items()}


def compute_generation_scales(node_weights, rates, verbose=False):
    """Per-leg producer-scale coefficients that put each component's production in
    absolute vehicle-driver trips/day, so the tuned K_c should land at ≈ 1.0.

    node_weights : the loaded node_weights(_reduced).json dict.  The producer layers
        are summed island-wide (the external census nodes tile the whole island), so
        the per-capita anchors recompute automatically for any CENTRE.
    rates        : {commute, retail, res, school_primary, school_postprimary,
        school_tertiary} per-person/day (load_generation_rates()).  **All rates are
        per-capita** — each encodes the island TOTAL journeys of its type (rate × pop);
        the producer/attractor layer only distributes them spatially (no generation
        scaling).

    Returns the gen_scale dict consumed by constrained_od_flows:
        {res, com_out, com_ret, ret_out, ret_ret,
         sch_primary_out/ret, sch_postprimary_out/ret, sch_tertiary_out/ret}.
    Two-leg directions carry ρ_c/2; res is single-leg (full ρ_res, both directions are
    separate OD pairs).  Per-producer rate r = (ρ_c/2)/k with island anchor
    k = Σ(producer layer)/Σ(population) over all nodes (k=1 when producer = population),
    so production = producer_share × (ρ_c/2) × pop — the producer contributes only its
    spatial share.  The school per-capita rates were built as (per-student behaviour ×
    island student/pop) in derive_generation_rates, so ρ_school/k_students recovers the
    per-student rate exactly; each of the three school levels is fully independent.
    """
    def _sum(layer):
        return float(sum(node_weights.get(layer, {}).values()))

    pop_tot = _sum("node_population")
    if pop_tot <= 0:
        raise ValueError("compute_generation_scales: Σ node_population is zero")

    anchors = {
        "k_commuters": _sum("node_commute_producers")  / pop_tot,
        "k_jobs":      _sum("node_commute_attractor")   / pop_tot,
        "k_retail":    _sum("node_retail_spaces")     / pop_tot,
    }
    for lvl in SCHOOL_LEVELS:
        anchors[f"k_students_{lvl}"]  = _sum(f"node_school_producers_{lvl}") / pop_tot
        anchors[f"k_enrolment_{lvl}"] = _sum(f"node_school_demand_{lvl}")    / pop_tot
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
    }
    for lvl in SCHOOL_LEVELS:
        rho = rates[f"school_{lvl}"]                          # per-capita (encodes island total)
        gen_scale[f"sch_{lvl}_out"] = (rho / 2) / anchors[f"k_students_{lvl}"]
        gen_scale[f"sch_{lvl}_ret"] = (rho / 2) / anchors[f"k_enrolment_{lvl}"]
    if verbose:
        print(f"  [gen-scale] island pop {pop_tot:,.0f}; anchors "
              + ", ".join(f"{n}={v:.4f}" for n, v in anchors.items()))
        print("  [gen-scale] producer rates (trips/day per producer): "
              + ", ".join(f"{k}={v:.4g}" for k, v in gen_scale.items()))
    return gen_scale


def build_mu_arrays(weights, node_ids, w_pop, w_school_prod_levels, verbose=False):
    """Per-node car-ownership multiplier arrays for constrained_od_flows(mu=…) (M3).

    Reads the raw `node_mu_<component>` layers (written by build_demographics /
    reduce_deadends) and normalises each to **producer-weighted mean 1** over the model's
    node set — res/retail by population (`w_pop`), each school level by that level's students
    (`w_school_prod_levels[lvl]`).  That is the exact level-preserving condition
    (Σ μ_c·prod_c = Σ prod_c), so μ redistributes production spatially without moving any
    component's total.  Returns `None` when no `node_mu_*` layers are present (⇒ μ off,
    pre-M3 behaviour, and constrained_od_flows(mu=None) is the identical old path)."""
    if not any(str(k).startswith("node_mu_") for k in weights):
        return None
    _pnid = lambda k: (int(k) if str(k).lstrip("-").isdigit() else k)
    def _load(layer):
        d = {_pnid(k): v for k, v in weights.get(layer, {}).items()}
        return np.array([d.get(nid, 1.0) for nid in node_ids], dtype=np.float64)  # 1.0 = neutral
    def _norm(arr, wt):
        m = float((arr * wt).sum())
        return arr * (float(wt.sum()) / m) if m > 0 else arr
    mu = {"res":    _norm(_load("node_mu_res"),    w_pop),
          "retail": _norm(_load("node_mu_retail"), w_pop)}
    for lvl in SCHOOL_LEVELS:
        mu[f"school_{lvl}"] = _norm(_load(f"node_mu_school_{lvl}"), w_school_prod_levels[lvl])
    if verbose:
        print("  [μ] car-ownership multiplier active (producer-weighted mean 1): "
              + ", ".join(f"{k}∈[{a.min():.2f},{a.max():.2f}]" for k, a in mu.items()))
    return mu


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
    """Build the PER-COMPONENT intra-zonal self-term arrays for constrained_od_flows from
    the mass-weighted time histograms written by build_intra_times.py.

    File format: {"<zone census code>": {"<component>": {"t": [bin centres s], "w": [weights,
    Σ=1]}}} — one weighted histogram per external zone per component (res/commute/retail +
    the three school levels), the producer×attractor mass-weighted intra-zonal time
    distribution.  For each component this returns (self_src, self_dist, self_w) where
    self_src = the zone's node index (one per histogram bin), self_dist = bin centres,
    self_w = bin weights (so the zone contributes  a^c_i · Σ_bin w·F_c(t_bin) = a^c_i·S^c_i
    to its denominator).  Both legs of a component share its self-term (symmetric interaction).

    Returns {component: (self_src, self_dist, self_w)} (only components with ≥1 entry), or
    None if the file is absent / yields nothing (⇒ constrained_od_flows reverts to no
    self-term).  Zones in the file but absent from node_ids are skipped (printed).
    """
    if not os.path.exists(intra_times_file):
        print(f"  [self-term] {intra_times_file} not found — no intra-zonal self-term")
        return None
    with open(intra_times_file) as f:
        data = json.load(f)
    data.pop("_meta", None)
    node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
    per = defaultdict(lambda: ([], [], []))
    n_zones = 0
    missing = []
    for zid, by_comp in data.items():
        idx = node_to_idx.get(zid)
        if idx is None:
            missing.append(zid)
            continue
        used = False
        for comp, hist in by_comp.items():
            t, w = hist.get("t", []), hist.get("w", [])
            if not t:
                continue
            s_src, s_dist, s_w = per[comp]
            s_src.extend([idx] * len(t))
            s_dist.extend(t)
            s_w.extend(w)
            used = True
        if used:
            n_zones += 1
    out = {c: (np.asarray(s, dtype=np.intp),
               np.asarray(d, dtype=np.float64),
               np.asarray(w, dtype=np.float64))
           for c, (s, d, w) in per.items() if s}
    if not out:
        print(f"  [self-term] {intra_times_file} has no zones matching the cache nodes — no self-term")
        return None
    if missing:
        print(f"  [self-term] {len(missing)} intra-times zones absent from cache node_ids (skipped)")
    print(f"  [self-term] mass-weighted intra-zonal self-term active for {n_zones} external zones, "
          f"components " + ", ".join(f"{c}({len(v[0])})" for c, v in out.items()))
    return out

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
                 slot_fracs_school_levels=None):
    """Per-component AADT weights for converting a pre-K component flow to a daily total.

    The annual-average daily contribution of component c is K_c·m_c·W_c with
        W_c = (5·Σ_h f_c[weekday,h] + Σ_h f_c[Sat,h] + Σ_h f_c[Sun,h]) / 7,
    the day-type-weighted sum of the component's hourly fractions.

    The temporal profiles (derive_component_profiles.py) are per-component SHAPES,
    each normalised so **W_c ≈ 1** (decoupled from magnitude, which generation pins) —
    so K_c·m_c is itself ≈ the component's daily AADT, and the combined daily AADT is
    Σ_c K_c·m_c·W_c.

    Returns (W_res, W_commute, W_retail, W_school_by_level) where W_school_by_level is a
    dict {level: W} over SCHOOL_LEVELS.  Components with no slot_fracs return 0.0.
    """
    def _w(sf):
        if not sf:
            return 0.0
        return sum(_AADT_DAY_WEIGHT.get(dt, 0) * f for (dt, h), f in sf.items()) / 7.0
    sfl = slot_fracs_school_levels or {}
    W_school = {lvl: _w(sfl.get(lvl)) for lvl in SCHOOL_LEVELS}
    return (_w(slot_fracs_res), _w(slot_fracs_commute), _w(slot_fracs_retail), W_school)


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
                 link_flow_school_dicts=None,
                 slot_fracs_res=None, slot_fracs_commute=None,
                 slot_fracs_retail=None, slot_fracs_school_levels=None):
    """
    Compute chi²/N matching the tuner's component formulation.

    Four-component mode (link_flow_commute_dict and link_flow_retail_dict provided):
      link_flow_dict         — {(u,v): K_res * flow_res}
      link_flow_commute_dict — {(u,v): K_commute * flow_commute}
      link_flow_retail_dict  — {(u,v): K_retail  * flow_retail}
      link_flow_school_dicts — {level: {(u,v): K_<level> * flow_<level>}}  (optional, per school level)
      slot_fracs_res/commute/retail — {(day_type, hour): f};
      slot_fracs_school_levels — {level: {(day_type, hour): f}}
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
        if link_flow_school_dicts:                      # one comp per active school level
            _sfl = slot_fracs_school_levels or {}
            for lvl in SCHOOL_LEVELS:
                if link_flow_school_dicts.get(lvl) is not None:
                    comps.append((link_flow_school_dicts[lvl], _sfl.get(lvl)))
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
