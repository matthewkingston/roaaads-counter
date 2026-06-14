"""
Powell's-method tuning of gravity model parameters against combined chi-squared
objective across all traffic count observations.

K is analytically calibrated at each evaluation (not in the optimizer):
  K = sum(model_i * obs_i / sigma_i^2) / sum(model_i^2 / sigma_i^2)

All optimizer parameters are stored in log-space to enforce positivity.

Stage 1 (--gravity, default): tune W_BIZ, MU, SIGMA, ALPHA (4 params).
Stage 2 (--full): also tune city-level pop/wp and sub-1 dampings (27 params)
  with L2 regularisation relative to simulation/tuner_config.json references.

Node 180 (pop=50, wp=0) is excluded from stage 2 — too small to tune.

Results:
  simulation/tuned_params.json   best params from this run (read by build_assignment.py)
  simulation/tuning_history.jsonl  appended record of every run

Usage:
  python3 analysis/tune_assignment.py
  python3 analysis/tune_assignment.py --full
  python3 analysis/tune_assignment.py --tag "added-june-counts"
"""

import json, math, os, subprocess, sys, time
from datetime import datetime, timezone

import numpy as np
import scipy.optimize

# ── Paths ─────────────────────────────────────────────────────────────────────

PATHS_CACHE  = "simulation/newtownards_paths.npz"
WEIGHTS_FILE = "simulation/node_weights.json"
TUNER_CONFIG = "simulation/tuner_config.json"
LINK_AADT    = "data/link_aadt.json"
TUNED_PARAMS = "simulation/tuned_params.json"
HISTORY_FILE = "simulation/tuning_history.jsonl"

COUNT_SITES = [
    {"label": "site 507, A21 Bangor Road",     "node": 731, "links": [(731,730),(730,731)], "observed": 21_202},
    {"label": "site 508, A48 Donaghadee Road", "node":  47, "links": None,                  "observed": 10_792},
    {"label": "site 444, A20 Portaferry Road", "node":  92, "links": None,                  "observed":  7_282},
]

# ── CLI args ──────────────────────────────────────────────────────────────────

stage = "gravity"
tag   = None
argv  = sys.argv[1:]
i = 0
while i < len(argv):
    if argv[i] == "--full":
        stage = "full"
    elif argv[i] == "--tag" and i + 1 < len(argv):
        i += 1
        tag = argv[i]
    i += 1

if tag is None:
    try:
        tag = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        tag = "untagged"

print(f"Stage: {stage}  tag: {tag}")

# ── Load paths cache ──────────────────────────────────────────────────────────

if not os.path.exists(PATHS_CACHE):
    print(f"ERROR: {PATHS_CACHE} not found — run build_paths.py first")
    raise SystemExit(1)

print("Loading paths cache …")
cache        = np.load(PATHS_CACHE)
node_ids_arr = cache["node_ids"]
od_src       = cache["od_src"]
od_dst       = cache["od_dst"]
od_dist      = cache["od_dist"].astype(np.float64)
pair_idx     = cache["pair_idx"]
link_idx_arr = cache["link_idx"]
link_u       = cache["link_u"]
link_v       = cache["link_v"]
N_links      = len(link_u)
node_ids     = [int(nid) for nid in node_ids_arr]
N_nodes      = len(node_ids)

ln_od_dist = np.log(od_dist)
link_index = {(int(link_u[k]), int(link_v[k])): k for k in range(N_links)}

print(f"  {N_nodes} nodes  {N_links} links  {len(od_src):,} OD pairs")

# ── Load node weights ─────────────────────────────────────────────────────────

print("Loading node weights …")
with open(WEIGHTS_FILE) as f:
    wdata = json.load(f)

node_pop_full = {int(k): v for k, v in wdata["node_population"].items()}
node_biz_full = {int(k): v for k, v in wdata["node_business_demand"].items()}

# Precomputed base weight arrays (used in stage 1 and as fallback in stage 2)
base_w_pop = np.array([node_pop_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)
base_w_biz = np.array([node_biz_full.get(nid, 0.0) for nid in node_ids], dtype=np.float64)

# ── Load tuner config ─────────────────────────────────────────────────────────

with open(TUNER_CONFIG) as f:
    config = json.load(f)

lam       = config["lambda"]
city_list = list(config["cities"].items())

# External nodes covered by tuner_config (node 180 excluded)
external_nodes = set()
for _, city_cfg in city_list:
    external_nodes.update(city_cfg["nodes"])

# Node array indices for external nodes (for fast overriding in stage 2)
ext_indices = [(i, nid) for i, nid in enumerate(node_ids) if nid in external_nodes]

# Tunable dampings: those with value < 1.0, in deterministic order
tunable_dampings = []  # [(city_name, node_id, ref_damping), ...]
for city_name, city_cfg in city_list:
    for node_str in sorted(city_cfg["dampings"], key=int):
        damp = city_cfg["dampings"][node_str]
        if damp < 1.0:
            tunable_dampings.append((city_name, int(node_str), damp))

n_gravity = 4
n_city    = len(city_list) * 2
n_damp    = len(tunable_dampings)
n_ext     = n_city + n_damp  # external params in stage 2

print(f"  {len(city_list)} cities  {len(tunable_dampings)} tunable dampings")

# ── Build observation list ────────────────────────────────────────────────────

observations = []
for s in COUNT_SITES:
    observations.append((
        "official", s["node"], s["links"],
        float(s["observed"]), 0.10 * s["observed"]
    ))

if os.path.exists(LINK_AADT):
    with open(LINK_AADT) as f:
        link_aadt_data = json.load(f)["links"]
    for key, entry in sorted(link_aadt_data.items()):
        u, v = map(int, key.split(","))
        observations.append((
            "walking", (u, v), None,
            float(entry["aadt"]), float(entry["aadt_uncertainty"])
        ))

n_obs   = len(observations)
obs_arr = np.array([o[3] for o in observations], dtype=np.float64)
sig_arr = np.array([o[4] for o in observations], dtype=np.float64)

# Precompute link index sets per observation for fast model flow extraction
obs_link_idxs = []
for kind, target, links, _, _ in observations:
    if kind == "official":
        if links:
            idxs = [link_index[lnk] for lnk in links if lnk in link_index]
        else:
            idxs = [k for k in range(N_links)
                    if int(link_u[k]) == target or int(link_v[k]) == target]
    else:
        k = link_index.get(target, -1)
        idxs = [k] if k >= 0 else []
    obs_link_idxs.append(idxs)

print(f"  {n_obs} observations ({len(COUNT_SITES)} official, {n_obs - len(COUNT_SITES)} walking)")

# ── Assignment and chi-squared helpers ───────────────────────────────────────

def run_assignment(W_BIZ, MU, SIGMA, ALPHA, w_pop, w_biz):
    w_vec = w_pop + W_BIZ * w_biz
    t_ij  = (w_vec[od_src] * w_vec[od_dst]
             * np.exp(-0.5 * ((ln_od_dist - MU) / SIGMA) ** 2)
             / od_dist ** ALPHA)
    return np.bincount(link_idx_arr, weights=t_ij[pair_idx], minlength=N_links)


def model_obs(raw_flow):
    return np.array([raw_flow[idxs].sum() if idxs else 0.0
                     for idxs in obs_link_idxs])


def calibrate_K(m_arr):
    w2  = 1.0 / sig_arr ** 2
    num = float(np.dot(w2 * m_arr, obs_arr))
    den = float(np.dot(w2 * m_arr, m_arr))
    return num / den if den > 0 else 1.0


# ── Objective function ────────────────────────────────────────────────────────

eval_count = [0]
best       = {"chi2": float("inf"), "log_params": None, "K": 1.0}
t0         = time.time()


def objective(log_params, log_ref=None):
    W_BIZ = math.exp(log_params[0])
    MU    = math.exp(log_params[1])
    SIGMA = math.exp(log_params[2])
    ALPHA = math.exp(log_params[3])

    if stage == "full":
        # Build per-node weight arrays from city-level params and dampings
        w_pop = base_w_pop.copy()
        w_biz = base_w_biz.copy()

        idx = n_gravity
        city_pops = {}
        city_wps  = {}
        for city_name, _ in city_list:
            city_pops[city_name] = math.exp(log_params[idx])
            city_wps[city_name]  = math.exp(log_params[idx + 1])
            idx += 2

        # Start from reference dampings, then override tunable ones
        curr_dampings = {}
        for city_name, city_cfg in city_list:
            for node_str, damp in city_cfg["dampings"].items():
                curr_dampings[int(node_str)] = damp
        for i_td, (_, node_id, _) in enumerate(tunable_dampings):
            curr_dampings[node_id] = math.exp(log_params[n_gravity + n_city + i_td])

        # Apply to node arrays
        for city_name, city_cfg in city_list:
            for node_id in city_cfg["nodes"]:
                damp = curr_dampings[node_id]
                # Find this node's index in node_ids
                pass  # done below via ext_indices

        # Use ext_indices for efficient override
        ext_pop_map = {}
        ext_biz_map = {}
        for city_name, city_cfg in city_list:
            for node_id in city_cfg["nodes"]:
                damp = curr_dampings[node_id]
                ext_pop_map[node_id] = city_pops[city_name] * damp
                ext_biz_map[node_id] = city_wps[city_name]  * damp
        for arr_i, nid in ext_indices:
            w_pop[arr_i] = ext_pop_map.get(nid, base_w_pop[arr_i])
            w_biz[arr_i] = ext_biz_map.get(nid, base_w_biz[arr_i])
    else:
        w_pop = base_w_pop
        w_biz = base_w_biz

    raw_flow = run_assignment(W_BIZ, MU, SIGMA, ALPHA, w_pop, w_biz)
    m_arr    = model_obs(raw_flow)
    K        = calibrate_K(m_arr)
    resid    = (K * m_arr - obs_arr) / sig_arr
    chi2     = float(np.dot(resid, resid))

    if stage == "full" and log_ref is not None:
        chi2 += lam * float(np.sum((log_params[n_gravity:] - log_ref[n_gravity:]) ** 2))

    eval_count[0] += 1
    if chi2 < best["chi2"]:
        best["chi2"]       = chi2
        best["log_params"] = log_params.copy()
        best["K"]          = K

    if eval_count[0] % 20 == 0:
        elapsed = time.time() - t0
        print(f"  {eval_count[0]:4d}  χ²={chi2:.2f}  χ²/N={chi2/n_obs:.2f}"
              f"  K={K:.4f}  best={best['chi2']:.2f}  ({elapsed:.0f}s)")

    return chi2


# ── Build initial parameter vector ────────────────────────────────────────────

# Gravity start: from tuned_params.json if available, else hardcoded defaults
grav_start = {"W_BIZ": 1.0, "MU": 7.5, "SIGMA": 1.0, "ALPHA": 2.0}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        tp = json.load(f)
    for k in ("W_BIZ", "MU", "SIGMA", "ALPHA"):
        if k in tp:
            grav_start[k] = tp[k]
    print(f"Starting gravity params from {TUNED_PARAMS}")

log_p0 = np.array([
    math.log(grav_start["W_BIZ"]),
    math.log(grav_start["MU"]),
    math.log(grav_start["SIGMA"]),
    math.log(grav_start["ALPHA"]),
], dtype=np.float64)

log_ref = None

if stage == "full":
    # City pop/wp: reference values from config
    for city_name, city_cfg in city_list:
        log_p0 = np.append(log_p0, [
            math.log(city_cfg["ref_pop"]),
            math.log(max(city_cfg["ref_wp"], 1)),  # guard against 0 (not present in config)
        ])
    # Tunable dampings: reference values from config
    for _, _, ref_damp in tunable_dampings:
        log_p0 = np.append(log_p0, math.log(ref_damp))

    log_ref = log_p0.copy()
    print(f"Full stage: {len(log_p0)} params "
          f"({n_gravity} gravity + {n_city} city pop/wp + {n_damp} dampings)")
else:
    print(f"Gravity stage: {len(log_p0)} params")

# ── Run optimization ──────────────────────────────────────────────────────────

print(f"\nRunning Powell's method (λ={lam}) …")
print(f"  {'eval':>4s}  χ²        χ²/N    K        best     elapsed")

# Evaluate initial point
objective(log_p0, log_ref)

result = scipy.optimize.minimize(
    lambda p: objective(p, log_ref),
    log_p0,
    method="powell",
    options={"maxiter": 5000, "ftol": 1e-5 if stage == "gravity" else 1e-4,
             "xtol":  1e-5 if stage == "gravity" else 1e-4},
)

# Use best params seen (Powell may backtrack at convergence)
log_best = best["log_params"]

# ── Unpack best params ────────────────────────────────────────────────────────

W_BIZ = math.exp(log_best[0])
MU    = math.exp(log_best[1])
SIGMA = math.exp(log_best[2])
ALPHA = math.exp(log_best[3])

ext_pop_map  = {}
ext_biz_map  = {}
city_pops_out = {}
city_wps_out  = {}
dampings_out  = {}

if stage == "full":
    idx = n_gravity
    for city_name, city_cfg in city_list:
        city_pops_out[city_name] = math.exp(log_best[idx])
        city_wps_out[city_name]  = math.exp(log_best[idx + 1])
        idx += 2

    curr_dampings = {}
    for city_name, city_cfg in city_list:
        for node_str, damp in city_cfg["dampings"].items():
            curr_dampings[int(node_str)] = damp
    for i_td, (_, node_id, _) in enumerate(tunable_dampings):
        d = math.exp(log_best[n_gravity + n_city + i_td])
        curr_dampings[node_id] = d
        dampings_out[str(node_id)] = d

    for city_name, city_cfg in city_list:
        for node_id in city_cfg["nodes"]:
            damp = curr_dampings[node_id]
            ext_pop_map[node_id] = city_pops_out[city_name] * damp
            ext_biz_map[node_id] = city_wps_out[city_name]  * damp

# Final evaluation for clean chi2 and K (no L2 term)
w_pop_f = base_w_pop.copy() if ext_pop_map else base_w_pop
w_biz_f = base_w_biz.copy() if ext_biz_map else base_w_biz
for arr_i, nid in ext_indices:
    if nid in ext_pop_map:
        w_pop_f[arr_i] = ext_pop_map[nid]
        w_biz_f[arr_i] = ext_biz_map[nid]

raw_flow  = run_assignment(W_BIZ, MU, SIGMA, ALPHA, w_pop_f, w_biz_f)
m_arr     = model_obs(raw_flow)
K         = calibrate_K(m_arr)
resid     = (K * m_arr - obs_arr) / sig_arr
chi2      = float(np.dot(resid, resid))
chi2_per_n = chi2 / n_obs

elapsed = time.time() - t0
print(f"\nResult  ({eval_count[0]} evals, {elapsed:.0f}s)")
print(f"  K={K:.4f}  W_BIZ={W_BIZ:.4f}  MU={MU:.4f}  SIGMA={SIGMA:.4f}  ALPHA={ALPHA:.4f}")
print(f"  χ²={chi2:.2f}  χ²/N={chi2_per_n:.3f}  (target ~1.0)")

# ── Print goodness-of-fit table ───────────────────────────────────────────────

print(f"\n  {'Source':<10s}  {'Label':<42s}  {'Obs':>7s}  {'σ':>6s}  {'Model':>7s}  {'z':>6s}")
for i_obs, (kind, target, links, obs, sig) in enumerate(observations):
    if kind == "official":
        lbl = next(s["label"] for s in COUNT_SITES if s["node"] == target)
    else:
        lbl = f"{target[0]}→{target[1]}"
    mod = K * m_arr[i_obs]
    z   = (mod - obs) / sig
    print(f"  {kind:<10s}  {lbl:<42s}  {obs:>7,.0f}  {sig:>6,.0f}  {mod:>7,.0f}  {z:>+.2f}")

# ── Save tuned_params.json ────────────────────────────────────────────────────

tuned = {
    "K":     round(K, 6),
    "W_BIZ": round(W_BIZ, 6),
    "MU":    round(MU, 6),
    "SIGMA": round(SIGMA, 6),
    "ALPHA": round(ALPHA, 6),
    "external_node_pop": {str(k): round(v) for k, v in ext_pop_map.items()},
    "external_node_biz": {str(k): round(v) for k, v in ext_biz_map.items()},
    "chi2":       round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "n_obs":      n_obs,
    "stage":      stage,
}
if stage == "full":
    tuned["external_city_pop"] = {k: round(v) for k, v in city_pops_out.items()}
    tuned["external_city_wp"]  = {k: round(v) for k, v in city_wps_out.items()}
    tuned["external_dampings"] = {k: round(v, 4) for k, v in dampings_out.items()}

with open(TUNED_PARAMS, "w") as f:
    json.dump(tuned, f, indent=2)
print(f"\nSaved → {TUNED_PARAMS}")

# ── Append to tuning history ──────────────────────────────────────────────────

history_entry = {
    "timestamp":  datetime.now(timezone.utc).isoformat(),
    "tag":        tag,
    "stage":      stage,
    "n_evals":    eval_count[0],
    "n_obs":      n_obs,
    "n_params":   len(log_p0),
    "chi2":       round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "K":     round(K, 6),
    "W_BIZ": round(W_BIZ, 6),
    "MU":    round(MU, 6),
    "SIGMA": round(SIGMA, 6),
    "ALPHA": round(ALPHA, 6),
    "observations": [
        {
            "kind":     obs[0],
            "target":   list(obs[1]) if isinstance(obs[1], tuple) else obs[1],
            "observed": obs[3],
            "sigma":    round(obs[4], 1),
        }
        for obs in observations
    ],
}

with open(HISTORY_FILE, "a") as f:
    f.write(json.dumps(history_entry) + "\n")
print(f"Appended → {HISTORY_FILE}")
