"""
Shared constants and functions for the Newtownards gravity model pipeline.

Imported by simulation/build_assignment.py and analysis/tune_assignment.py to
keep their chi²/N calculations, count-site definitions, and gravity kernel
implementations in sync.
"""

import json
import os
import numpy as np

# ── Official AADT count sites ─────────────────────────────────────────────────

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

PATHS_CACHE  = "simulation/newtownards_paths.npz"
WEIGHTS_FILE = "simulation/node_weights.json"
TUNER_CONFIG = "simulation/tuner_config.json"
TUNED_PARAMS = "simulation/tuned_params.json"
LINK_AADT    = "data/link_aadt.json"

# ── Gravity kernel ────────────────────────────────────────────────────────────

def gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx, N_links,
                   W_BIZ, P, ALPHA, w_pop, w_biz):
    """
    Rational kernel all-or-nothing assignment.

    Kernel: f(d) = (ALPHA+1)*u / (ALPHA + u^(ALPHA+1))  where u = d/P.
    Peak at d=P with f(P)=1; tail ~ 1/d^ALPHA.

    Returns pre-K raw link flow array (length N_links).
    """
    w_vec = w_pop + W_BIZ * w_biz
    u     = od_dist / P
    t_ij  = w_vec[od_src] * w_vec[od_dst] * (ALPHA + 1) * u / (ALPHA + u ** (ALPHA + 1))
    return np.bincount(link_idx, weights=t_ij[pair_idx], minlength=N_links)

# ── Flow extraction ───────────────────────────────────────────────────────────

def site_flow(link_flow_dict, site):
    """Return total modelled flow for a COUNT_SITES entry."""
    if site["links"]:
        return sum(link_flow_dict.get(lnk, 0.0) for lnk in site["links"])
    node = site["node"]
    return sum(f for (u, v), f in link_flow_dict.items() if u == node or v == node)

# ── Chi²/N ───────────────────────────────────────────────────────────────────

def compute_chi2(link_flow_dict, label_fn=None,
                 link_aadt_file=LINK_AADT, exclude_links=EXCLUDE_LINKS):
    """
    Compute chi²/N using per-session walking count observations with Woodbury correction.

    Walking observations in the same (weekday, hour) time slot share a correlated
    fractional AADT uncertainty (same hourly-fraction draw). The Woodbury rank-1
    correction removes this double-counting, matching the tuner's objective exactly.

    label_fn: optional callable(u, v) -> str for road-name labels.
    Returns (rows, chi2, n_obs, n_eff).
      rows: list of (kind, label, obs, sig, mod, z) sorted by |z| descending.
      chi2: Woodbury-corrected chi² sum.
      n_obs: total observation count.
      n_eff: effective df = n_obs - n_slots (one df lost per correlated slot).
    """
    # ── Build raw observation list ────────────────────────────────────────────
    # Each entry: (kind, label, mod, obs, sigma, time_slot_key, frac_rel_std)
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

    # ── Woodbury decomposition ────────────────────────────────────────────────
    sigma_f = np.array([
        (d[3] * d[6]) if d[6] is not None else 0.0   # obs * frac_rel_std
        for d in obs_data
    ])
    sigma_sq   = np.array([d[4] ** 2 for d in obs_data])
    sigma_c_sq = np.maximum(sigma_sq - sigma_f ** 2, (np.sqrt(sigma_sq) * 1e-6) ** 2)
    r          = np.array([d[2] - d[3] for d in obs_data])  # mod - obs

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

    # ── Build display rows (z uses total sigma, same as tuner's fit table) ────
    rows = [
        (d[0], d[1], d[3], d[4], d[2], float(r[i] / d[4]))
        for i, d in enumerate(obs_data)
    ]
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
