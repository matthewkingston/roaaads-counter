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
    Compute chi²/N using per-session walking count observations.

    Uses the same observation data as tune_assignment.py. No Woodbury
    correction is applied, so the returned chi²/N will be slightly higher
    than the tuner's Woodbury-corrected value.

    label_fn: optional callable(u, v) -> str for road-name labels on walking
              rows. Defaults to "u→v" if None.

    Returns (rows, chi2, n) where rows are sorted by |z| descending.
    Each row: (kind, label, obs, sig, mod, z).
    """
    rows = []
    chi2 = 0.0

    for s in COUNT_SITES:
        mod = site_flow(link_flow_dict, s)
        obs = float(s["observed"])
        sig = 0.10 * obs
        z   = (mod - obs) / sig
        chi2 += z * z
        rows.append(("official", s["label"], obs, sig, mod, z))

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
            for sess_obs in entry.get("observations", []):
                obs = float(sess_obs["aadt"])
                sig = float(sess_obs["aadt_uncertainty"])
                z   = (mod - obs) / sig
                chi2 += z * z
                rows.append(("walking", lbl, obs, sig, mod, z))

    rows.sort(key=lambda r: abs(r[5]), reverse=True)
    return rows, chi2, len(rows)


def print_chi2_table(rows, chi2, n, header="Goodness of fit"):
    """Print the standard chi²/N table. Returns chi²/N."""
    chi2_per_n = chi2 / n if n else 0.0
    LABEL_W = 52
    print(f"\n{header}  χ²={chi2:.2f}  n={n}  χ²/N={chi2_per_n:.4f}")
    print(f"  {'':1s}  {'Src':<8}  {'Label':<{LABEL_W}}  {'Obs':>8}  {'σ':>7}  {'Model':>8}  {'z':>6}")
    for kind, lbl, obs, sig, mod, z in rows:
        marker = "*" if abs(z) > 2 else " "
        print(f"  {marker} {kind:<8}  {lbl:<{LABEL_W}}  {obs:>8,.0f}  {sig:>7,.0f}  {mod:>8,.0f}  {z:>+.2f}")
    abs_z  = [abs(r[5]) for r in rows]
    n_out2 = sum(1 for a in abs_z if a > 2)
    n_out3 = sum(1 for a in abs_z if a > 3)
    mean_z = sum(abs_z) / len(abs_z) if abs_z else 0.0
    print(f"\n  n={n}  χ²/N={chi2_per_n:.4f}  mean|z|={mean_z:.2f}"
          f"  |z|>2: {n_out2}  |z|>3: {n_out3}")
    return chi2_per_n
