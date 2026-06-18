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
    {"label": "site 507, A21 Bangor Road",     "node": 731, "links": [(731, 730), (730, 731)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road", "node":  47, "links": None,                      "observed": 10_792},
    {"label": "site 444, A20 Portaferry Road", "node":  92, "links": None,                      "observed":  7_282},
]

# Links present in link_aadt.json but excluded from calibration.
# Directed: (u, v) excludes only that direction.
EXCLUDE_LINKS = {
    (161, 160),  # too-short count, likely distorted by traffic light timing
}

# ── File paths ────────────────────────────────────────────────────────────────

PATHS_CACHE     = "simulation/newtownards_paths.npz"
WEIGHTS_FILE    = "simulation/node_weights.json"
TUNER_CONFIG    = "simulation/tuner_config.json"
TUNED_PARAMS    = "simulation/tuned_params.json"
LINK_AADT       = "data/link_aadt.json"
OFFICIAL_HOURLY = "data/official_hourly.json"

_DOW_TO_TYPE = {d: (0 if d < 5 else (1 if d == 5 else 2)) for d in range(7)}

# ── Gravity kernel ────────────────────────────────────────────────────────────

def gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
                   W_BIZ, P, ALPHA, w_pop, w_biz,
                   BETA=1.0, THETA=None,
                   od_dist_2=None, pair_idx_2=None, link_idx_2=None,
                   od_dist_3=None, pair_idx_3=None, link_idx_3=None,
                   return_components=False):
    """
    Generalised rational kernel assignment.

    Kernel: f(d) = (ALPHA+BETA)*u^BETA / (ALPHA + BETA*u^(ALPHA+BETA))  where u = d/P.
    Peak at d=P with f(P)=1; tail ~ 1/d^ALPHA; rise ~ u^BETA near origin.
    BETA=1 (default) recovers the original kernel (ALPHA+1)*u / (ALPHA + u^(ALPHA+1)).

    When THETA is None or k=2/k=3 arrays absent: all-or-nothing on k=1 path.
    When THETA is given with k=2/k=3 arrays: logit spread across 3 paths.
      share(r) ∝ exp(−THETA · d_r / P); THETA→∞ collapses to all-or-nothing.

    return_components=False (default): returns combined pre-K flow array (N_links,).
    return_components=True: returns (flow_res, flow_biz) where
      flow_res = pop→pop component, flow_biz = W_BIZ·pb + W_BIZ²·bb component.
    """
    _has_stoch = (THETA is not None and od_dist_2 is not None)

    if not _has_stoch:
        u    = od_dist / P
        kern = (ALPHA + BETA) * u**BETA / (ALPHA + BETA * u**(ALPHA + BETA))

        if not return_components:
            w_vec = w_pop + W_BIZ * w_biz
            t_ij  = w_vec[od_src] * w_vec[od_dst] * kern
            return np.bincount(link_idx, weights=t_ij[pair_idx], minlength=N_links)

        pp = w_pop[od_src] * w_pop[od_dst] * kern
        pb = (w_pop[od_src] * w_biz[od_dst] + w_biz[od_src] * w_pop[od_dst]) * kern
        bb = w_biz[od_src] * w_biz[od_dst] * kern
        flow_res = np.bincount(link_idx, weights=pp[pair_idx], minlength=N_links)
        flow_biz = np.bincount(link_idx,
                               weights=((W_BIZ * pb + W_BIZ ** 2 * bb)[pair_idx]),
                               minlength=N_links)
        return flow_res, flow_biz

    # ── Stochastic logit ──────────────────────────────────────────────────────
    d_mat  = np.stack([od_dist, od_dist_2, od_dist_3], axis=1)
    log_w  = -THETA * d_mat / P
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

    pp_od = w_pop[od_src] * w_pop[od_dst]
    pb_od = w_pop[od_src] * w_biz[od_dst] + w_biz[od_src] * w_pop[od_dst]
    bb_od = w_biz[od_src] * w_biz[od_dst]
    flow_res = np.zeros(N_links, dtype=np.float64)
    flow_biz = np.zeros(N_links, dtype=np.float64)
    for r, (pidx, lidx, d_r) in enumerate([
        (pair_idx,   link_idx,   od_dist),
        (pair_idx_2, link_idx_2, od_dist_2),
        (pair_idx_3, link_idx_3, od_dist_3),
    ]):
        u_r  = d_r / P
        f_r  = (ALPHA + BETA) * u_r**BETA / (ALPHA + BETA * u_r**(ALPHA + BETA))
        s_r  = shares[:, r]
        flow_res += np.bincount(lidx, weights=(pp_od * s_r * f_r)[pidx], minlength=N_links)
        flow_biz += np.bincount(lidx,
                                weights=((W_BIZ * pb_od + W_BIZ ** 2 * bb_od) * s_r * f_r)[pidx],
                                minlength=N_links)
    return flow_res, flow_biz

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

# ── Chi²/N ───────────────────────────────────────────────────────────────────

def compute_chi2(link_flow_dict, label_fn=None,
                 link_aadt_file=LINK_AADT, exclude_links=EXCLUDE_LINKS,
                 official_hourly_file=OFFICIAL_HOURLY,
                 link_flow_biz_dict=None,
                 slot_fracs_res=None, slot_fracs_biz=None):
    """
    Compute chi²/N matching the tuner's two-component formulation.

    Two-component mode (link_flow_biz_dict provided):
      link_flow_dict     — {(u,v): K_res * flow_res}  (residential component, scaled)
      link_flow_biz_dict — {(u,v): K_biz * flow_biz}  (business component, scaled)
      slot_fracs_res/biz — {(day_type, hour): f}  from tuned_params.json
      Uses 216 official hourly obs (count-space Gaussian) from official_hourly_file
      plus per-session walking obs (count-space Poisson) from link_aadt_file.
      N_eff = N − 2·N_slots.  No coupling penalty (pure data fit).

    Legacy mode (link_flow_biz_dict=None):
      link_flow_dict — {(u,v): K * combined_flow}  (pre-scaled)
      Uses 3 AADT-space obs from COUNT_SITES plus Woodbury-corrected walking obs.
      N_eff = N − N_slots.

    Returns (rows, chi2, n_obs, n_eff).
      rows: list of (kind, label, obs, sig, mod, z) sorted by |z| descending.
    """
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
                chi2    += z ** 2
                n_slots_seen.add(sk)
                # display in AADT-space: divide by effective f_s
                m_tot = m_r + m_b + 1e-30
                f_eff = (m_r * f_r + m_b * f_b) / m_tot
                if f_eff > 0 and Th > 0:
                    obs_disp = neff / (Th * f_eff)
                    sig_disp = math.sqrt(neff) / (Th * f_eff)
                    mod_disp = pred / f_eff
                else:
                    obs_disp = sig_disp = mod_disp = 0.0
                rows.append(("walking", lbl, obs_disp, sig_disp, mod_disp, z))

    n_obs   = len(rows)
    # N_eff: 2 df per slot (f_res and f_biz each).
    # Official obs cover all 72 possible (day_type, hour) slots; walking slots are a subset.
    # Take the union so overlapping slots aren't counted twice.
    if n_official:
        n_slots_seen |= {(dt, h) for dt in range(3) for h in range(24)}
    n_eff   = n_obs - 2 * len(n_slots_seen)

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
