"""
Powell's-method tuning of gravity model parameters against combined chi-squared
objective across all traffic count observations.

K is analytically calibrated at each evaluation (not in the optimizer):
  K = sum(model_i * obs_i / sigma_i^2) / sum(model_i^2 / sigma_i^2)

All optimizer parameters are stored in log-space to enforce positivity.

Stage 1 (--gravity, default): tune W_BIZ, P, ALPHA (3 params).
Stage 2 (--full): also tune city-level pop/wp and sub-1 dampings (26 params)
  with L2 regularisation relative to simulation/tuner_config.json references.

Node 180 (pop=50, wp=0) is excluded from stage 2 — too small to tune.

Results:
  simulation/tuned_params.json   best params from this run (read by build_assignment.py)
  simulation/tuning_history.jsonl  appended record of every run
  simulation/gravity_model_curve.png  kernel shape plot

Usage:
  python3 analysis/tune_assignment.py
  python3 analysis/tune_assignment.py --full
  python3 analysis/tune_assignment.py --note "added-june-counts"
"""

import json, math, os, secrets, subprocess, sys, time, xml.etree.ElementTree as ET
from datetime import datetime, timezone

import numpy as np
import scipy.optimize

sys.path.insert(0, "simulation")
from model import (COUNT_SITES, EXCLUDE_LINKS, PATHS_CACHE, WEIGHTS_FILE,
                   TUNER_CONFIG, LINK_AADT, TUNED_PARAMS,
                   gravity_assign, print_chi2_table)

# ── Paths ─────────────────────────────────────────────────────────────────────

CONS_GRAPH   = "simulation/newtownards_consolidated.graphml"
HISTORY_FILE = "simulation/tuning_history.jsonl"
CURVE_PNG    = "simulation/gravity_model_curve.png"

# ── CLI args ──────────────────────────────────────────────────────────────────

stage = "gravity"
note  = None
argv  = sys.argv[1:]
i = 0
while i < len(argv):
    if argv[i] == "--full":
        stage = "full"
    elif argv[i] == "--note" and i + 1 < len(argv):
        i += 1
        note = argv[i]
    i += 1

run_id = secrets.token_hex(4)

try:
    git_hash = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
    ).decode().strip()
except Exception:
    git_hash = "unknown"

print(f"Run ID: {run_id}  stage: {stage}  git: {git_hash}" +
      (f"  note: {note}" if note else ""))

# ── Load previous best from history ──────────────────────────────────────────

prev_chi2_per_n = None
prev_id         = None
if os.path.exists(HISTORY_FILE):
    with open(HISTORY_FILE) as _hf:
        _lines = [ln for ln in _hf if ln.strip()]
    if _lines:
        _last = json.loads(_lines[-1])
        prev_chi2_per_n = _last.get("chi2_per_n")
        prev_id         = _last.get("id", "")
        print(f"Previous run:  χ²/N={prev_chi2_per_n:.4f}  (id={prev_id})")

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

# ── Load street names from consolidated GraphML (sequential node IDs) ─────────

link_name = {}
if os.path.exists(CONS_GRAPH):
    _ns  = "http://graphml.graphstruct.org/graphml"
    _ns2 = "http://graphml.graphstruct.org/graphml"
    try:
        _tree = ET.parse(CONS_GRAPH)
        _root = _tree.getroot()
        _nsmap = {"g": _root.tag.split("}")[0].lstrip("{")} if "}" in _root.tag else {"g": ""}
        _pfx = "{" + _nsmap["g"] + "}" if _nsmap["g"] else ""
        for _edge in _root.iter(f"{_pfx}edge"):
            _u = int(_edge.get("source"))
            _v = int(_edge.get("target"))
            for _data in _edge:
                if _data.get("key") == "d14" and _data.text:
                    link_name[(_u, _v)] = _data.text
                    break
    except Exception as _e:
        print(f"Warning: could not load street names ({_e})")


def _link_label(u, v):
    name = link_name.get((u, v), "")
    return f"{u}→{v}  {name}" if name else f"{u}→{v}"

# ── Load tuner config ─────────────────────────────────────────────────────────

with open(TUNER_CONFIG) as f:
    config = json.load(f)

lam       = config["lambda"]
city_list = list(config["cities"].items())

grav_ref = config.get("gravity_ref", {})
grav_lam = config.get("gravity_lambda", 0.0)
log_grav_ref = np.array([
    math.log(max(grav_ref.get("W_BIZ", 1.0),   1e-4)),
    math.log(max(grav_ref.get("P",     300.0),  1e-4)),
    math.log(max(grav_ref.get("ALPHA", 2.0),    1e-4)),
])

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

n_gravity = 3
n_city    = len(city_list) * 2
n_damp    = len(tunable_dampings)
n_ext     = n_city + n_damp  # external params in stage 2

print(f"  {len(city_list)} cities  {len(tunable_dampings)} tunable dampings")

# Override external node weights from tuner_config so gravity-only tuning
# always uses the same city values as the --full stage, regardless of what
# node_weights.json contains.
_ext_pop_cfg = {}
_ext_biz_cfg = {}
for _city_name, _city_cfg in city_list:
    for _nid in _city_cfg["nodes"]:
        _damp = _city_cfg["dampings"][str(_nid)]
        _ext_pop_cfg[_nid] = _city_cfg["ref_pop"] * _damp
        _ext_biz_cfg[_nid] = _city_cfg["ref_wp"]  * _damp
for _arr_i, _nid in ext_indices:
    base_w_pop[_arr_i] = _ext_pop_cfg[_nid]
    base_w_biz[_arr_i] = _ext_biz_cfg[_nid]

# ── Build observation list ────────────────────────────────────────────────────

observations     = []
obs_time_slots   = []
obs_frac_rel_std = []

for s in COUNT_SITES:
    observations.append(("official", s["node"], s["links"],
                         float(s["observed"]), 0.10 * s["observed"]))
    obs_time_slots.append(None)
    obs_frac_rel_std.append(None)

if os.path.exists(LINK_AADT):
    with open(LINK_AADT) as f:
        link_aadt_data = json.load(f)["links"]
    for key, entry in sorted(link_aadt_data.items()):
        u, v = map(int, key.split(","))
        if (u, v) in EXCLUDE_LINKS:
            continue
        for sess_obs in entry["observations"]:
            observations.append(("walking", (u, v), None,
                float(sess_obs["aadt"]), float(sess_obs["aadt_uncertainty"])))
            ts  = sess_obs.get("time_slot")
            frs = sess_obs.get("frac_rel_std")
            obs_time_slots.append(tuple(ts) if ts is not None else None)
            obs_frac_rel_std.append(float(frs) if frs is not None else None)

n_obs   = len(observations)
obs_arr = np.array([o[3] for o in observations], dtype=np.float64)
sig_arr = np.array([o[4] for o in observations], dtype=np.float64)

# Woodbury: decompose total uncertainty into correlated (σ_f) and independent (σ_c)
sigma_f_arr = np.zeros(n_obs)
for i, frs in enumerate(obs_frac_rel_std):
    if frs is not None:
        sigma_f_arr[i] = obs_arr[i] * frs

sigma_c_sq = np.maximum(sig_arr**2 - sigma_f_arr**2, (sig_arr * 1e-6)**2)

# Group slotted observations by (weekday, hour) — each slot shares one correlated mode
slot_groups = {}
for i, ts in enumerate(obs_time_slots):
    if ts is not None:
        slot_groups.setdefault(ts, []).append(i)
slot_list = list(slot_groups.items())

# Precompute per-slot Woodbury denominators and observation terms (obs-dependent, constant)
slot_denom  = {}
slot_uf_obs = {}
for ts, idxs in slot_list:
    slot_denom[ts]  = 1.0 + sum(sigma_f_arr[i]**2 / sigma_c_sq[i] for i in idxs)
    slot_uf_obs[ts] = sum(sigma_f_arr[i] * obs_arr[i] / sigma_c_sq[i] for i in idxs)

unslotted_idxs = [i for i, ts in enumerate(obs_time_slots) if ts is None]

n_slots   = len(slot_list)
n_slotted = sum(len(idxs) for _, idxs in slot_list)
n_eff     = n_obs - n_slots
n_official = len(COUNT_SITES)
print(f"  {n_obs} observations ({n_official} official, {n_obs - n_official} walking"
      f" in {n_slots} time slot(s))")
if n_slotted < n_obs - n_official:
    print(f"  Warning: {n_obs - n_official - n_slotted} walking obs have no time-slot"
          f" data — treated as independent")
for ts, idxs in sorted(slot_list):
    if len(idxs) > 2:
        print(f"  Slot {ts}: {len(idxs)} correlated observations")

# Precompute link index sets per observation for fast model flow extraction
obs_link_idxs = []
for kind, target, links, *_ in observations:
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

# ── Assignment and chi-squared helpers ───────────────────────────────────────

def run_assignment(W_BIZ, P, ALPHA, w_pop, w_biz):
    return gravity_assign(od_src, od_dst, od_dist, pair_idx, link_idx_arr, N_links,
                          W_BIZ, P, ALPHA, w_pop, w_biz)


def model_obs(raw_flow):
    return np.array([raw_flow[idxs].sum() if idxs else 0.0
                     for idxs in obs_link_idxs])


def calibrate_K(m_arr):
    A = 0.0
    B = 0.0
    for i in unslotted_idxs:
        w = 1.0 / sig_arr[i] ** 2
        A += w * m_arr[i] ** 2
        B += w * m_arr[i] * obs_arr[i]
    for ts, idxs in slot_list:
        denom = slot_denom[ts]
        uf_m  = sum(sigma_f_arr[i] * m_arr[i] / sigma_c_sq[i] for i in idxs)
        A += sum(m_arr[i] ** 2 / sigma_c_sq[i] for i in idxs) - uf_m ** 2 / denom
        B += (sum(obs_arr[i] * m_arr[i] / sigma_c_sq[i] for i in idxs)
              - slot_uf_obs[ts] * uf_m / denom)
    return B / A if A > 0 else 1.0


# ── Objective function ────────────────────────────────────────────────────────

eval_count = [0]
best       = {"chi2": float("inf"), "log_params": None, "K": 1.0}
t0         = time.time()


def objective(log_params, log_ref=None):
    log_params = np.clip(log_params, -100, 100)
    W_BIZ = math.exp(log_params[0])
    P     = math.exp(log_params[1])
    ALPHA = math.exp(log_params[2])

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

    raw_flow = run_assignment(W_BIZ, P, ALPHA, w_pop, w_biz)
    m_arr    = model_obs(raw_flow)
    K        = calibrate_K(m_arr)
    r        = K * m_arr - obs_arr
    chi2     = float(sum((r[i] / sig_arr[i]) ** 2 for i in unslotted_idxs))
    for ts, idxs in slot_list:
        denom = slot_denom[ts]
        uf_r  = sum(sigma_f_arr[i] * r[i] / sigma_c_sq[i] for i in idxs)
        chi2 += float(sum(r[i] ** 2 / sigma_c_sq[i] for i in idxs) - uf_r ** 2 / denom)

    if stage == "full" and log_ref is not None:
        chi2 += lam * float(np.sum((log_params[n_gravity:] - log_ref[n_gravity:]) ** 2))
    if grav_lam > 0:
        chi2 += grav_lam * float(np.sum((log_params[:n_gravity] - log_grav_ref) ** 2))

    eval_count[0] += 1
    if chi2 < best["chi2"]:
        best["chi2"]       = chi2
        best["log_params"] = log_params.copy()
        best["K"]          = K

    if eval_count[0] % 20 == 0:
        elapsed = time.time() - t0
        print(f"  {eval_count[0]:4d}  χ²={chi2:.2f}  χ²/N={chi2/n_obs:.3f}"
              f"  best={best['chi2']/n_obs:.3f}  ({elapsed:.0f}s)")

    return chi2


# ── Build initial parameter vector ────────────────────────────────────────────

# Gravity start: from tuned_params.json if available, else hardcoded defaults
grav_start = {"W_BIZ": 1.0, "P": 300.0, "ALPHA": 2.0}
if os.path.exists(TUNED_PARAMS):
    with open(TUNED_PARAMS) as f:
        tp = json.load(f)
    for k in ("W_BIZ", "P", "ALPHA"):
        if k in tp:
            grav_start[k] = tp[k]
    print(f"Starting gravity params from {TUNED_PARAMS}")

# Clamp to a safe minimum before log-transform (guards against degenerate prior runs)
_LOG_MIN = 1e-4
log_p0 = np.array([
    math.log(max(grav_start["W_BIZ"], _LOG_MIN)),
    math.log(max(grav_start["P"],     _LOG_MIN)),
    math.log(max(grav_start["ALPHA"], _LOG_MIN)),
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
print(f"  {'eval':>4s}  χ²/N(curr)  χ²/N(best)  elapsed")

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
P     = math.exp(log_best[1])
ALPHA = math.exp(log_best[2])

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

raw_flow  = run_assignment(W_BIZ, P, ALPHA, w_pop_f, w_biz_f)
m_arr     = model_obs(raw_flow)
K         = calibrate_K(m_arr)
r         = K * m_arr - obs_arr
chi2      = float(sum((r[i] / sig_arr[i]) ** 2 for i in unslotted_idxs))
for ts, idxs in slot_list:
    denom = slot_denom[ts]
    uf_r  = sum(sigma_f_arr[i] * r[i] / sigma_c_sq[i] for i in idxs)
    chi2 += float(sum(r[i] ** 2 / sigma_c_sq[i] for i in idxs) - uf_r ** 2 / denom)
chi2_per_n = chi2 / n_obs
resid      = r / sig_arr  # for fit table and history entry z-scores

elapsed = time.time() - t0
print(f"\nResult  ({eval_count[0]} evals, {elapsed:.0f}s)")
print(f"  K={K:.4e}  W_BIZ={W_BIZ:.4f}  P={P:.2f}s  ALPHA={ALPHA:.4f}")
print(f"  χ²={chi2:.2f}  χ²/N={chi2_per_n:.4f}  χ²/N_eff={chi2/n_eff:.3f}  (N={n_obs}, N_eff={n_eff})")
if prev_chi2_per_n is not None:
    delta = chi2_per_n - prev_chi2_per_n
    direction = "improvement" if delta < 0 else "regression"
    print(f"  vs previous ({prev_id}):  Δχ²/N={delta:+.4f}  ({direction})")

# ── City parameter delta table (full stage) ───────────────────────────────────

if stage == "full":
    print(f"\n  {'City':<12}  {'ref_pop':>10}  {'tuned_pop':>10}  {'Δpop%':>7}  "
          f"{'ref_wp':>8}  {'tuned_wp':>8}  {'Δwp%':>7}")
    for city_name, city_cfg in city_list:
        rp = city_cfg["ref_pop"]
        tp = city_pops_out[city_name]
        rw = city_cfg["ref_wp"]
        tw = city_wps_out[city_name]
        dp = 100.0 * (tp - rp) / rp
        dw = 100.0 * (tw - rw) / rw if rw > 0 else float("nan")
        print(f"  {city_name:<12}  {rp:>10,.0f}  {tp:>10,.0f}  {dp:>+7.1f}%  "
              f"{rw:>8,.0f}  {tw:>8,.0f}  {dw:>+7.1f}%")

# ── Goodness-of-fit table (sorted by |z|) ────────────────────────────────────

fit_rows = []
for i_obs, (kind, target, links, obs, sig) in enumerate(observations):
    if kind == "official":
        lbl = next(s["label"] for s in COUNT_SITES if s["node"] == target)
    else:
        lbl = _link_label(target[0], target[1])
    mod = K * m_arr[i_obs]
    z   = resid[i_obs]
    fit_rows.append((kind, lbl, obs, sig, mod, z))

fit_rows.sort(key=lambda r: abs(r[5]), reverse=True)
print_chi2_table(fit_rows, chi2, len(fit_rows))

# ── Save tuned_params.json ────────────────────────────────────────────────────

tuned = {
    "kernel": "rational",
    "K":      round(K, 6),
    "W_BIZ":  round(W_BIZ, 6),
    "P":      round(P, 4),
    "ALPHA":  round(ALPHA, 6),
    "external_node_pop": {str(k): round(v) for k, v in ext_pop_map.items()},
    "external_node_biz": {str(k): round(v) for k, v in ext_biz_map.items()},
    "chi2":       round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "n_obs":      n_obs,
    "n_slots":    n_slots,
    "n_eff":      n_eff,
    "stage":      stage,
}
if stage == "full":
    tuned["external_city_pop"] = {k: round(v) for k, v in city_pops_out.items()}
    tuned["external_city_wp"]  = {k: round(v) for k, v in city_wps_out.items()}
    tuned["external_dampings"] = {k: round(v, 4) for k, v in dampings_out.items()}

with open(TUNED_PARAMS, "w") as f:
    json.dump(tuned, f, indent=2)
print(f"\nSaved → {TUNED_PARAMS}")

# ── Gravity model kernel plot ─────────────────────────────────────────────────

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d_sec  = np.logspace(np.log10(30), np.log10(7200), 500)  # 30s – 120 min
    d_min  = d_sec / 60.0
    u      = d_sec / P
    kernel = K * (ALPHA + 1) * u / (ALPHA + u ** (ALPHA + 1))

    d_peak_sec = P          # peak is exactly at d = P
    k_peak     = K          # f(P) = 1, so kernel peak = K

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(d_min, kernel, linewidth=1.8)
    ax.axvline(d_peak_sec / 60.0, color="r", linestyle="--", alpha=0.7,
               label=f"peak  {d_peak_sec/60.0:.1f} min  (K={k_peak:.3e})")
    ax.set_xlabel("Travel time (minutes)")
    ax.set_ylabel("K · (ALPHA+1) · u / (ALPHA + u^(ALPHA+1))  where u = d/P")
    ax.set_title(
        f"Gravity kernel (rational)  "
        f"P={P:.1f}s  ALPHA={ALPHA:.3f}  K={K:.3e}\n"
        f"stage={stage}  χ²/N={chi2_per_n:.4f}  id={run_id}  git={git_hash}"
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(CURVE_PNG, dpi=150)
    plt.close(fig)
    print(f"Saved → {CURVE_PNG}")
except Exception as _e:
    print(f"Warning: could not save gravity curve plot ({_e})")

# ── Append to tuning history ──────────────────────────────────────────────────

params = {
    "kernel": "rational",
    "K":      round(K, 6),
    "W_BIZ":  round(W_BIZ, 6),
    "P":      round(P, 4),
    "ALPHA":  round(ALPHA, 6),
}
if stage == "full":
    params["external_node_pop"] = {str(k): round(v) for k, v in ext_pop_map.items()}
    params["external_node_biz"] = {str(k): round(v) for k, v in ext_biz_map.items()}
    params["external_city_pop"] = {k: round(v) for k, v in city_pops_out.items()}
    params["external_city_wp"]  = {k: round(v) for k, v in city_wps_out.items()}
    params["external_dampings"] = {k: round(v, 4) for k, v in dampings_out.items()}

history_entry = {
    "id":        run_id,
    "git_hash":  git_hash,
    "timestamp": datetime.now(timezone.utc).isoformat(),
    "stage":     stage,
    "n_evals":   eval_count[0],
    "n_obs":     n_obs,
    "n_slots":   n_slots,
    "n_eff":     n_eff,
    "n_params":  len(log_p0),
    "chi2":      round(chi2, 3),
    "chi2_per_n": round(chi2_per_n, 4),
    "params":    params,
    "observations": [
        {
            "kind":     observations[i_obs][0],
            "target":   (list(observations[i_obs][1])
                         if isinstance(observations[i_obs][1], tuple)
                         else observations[i_obs][1]),
            "observed": observations[i_obs][3],
            "sigma":    round(observations[i_obs][4], 1),
            "model":    round(float(K * m_arr[i_obs]), 1),
            "z":        round(float(resid[i_obs]), 3),
        }
        for i_obs in range(n_obs)
    ],
}
if note:
    history_entry["note"] = note

with open(HISTORY_FILE, "a") as f:
    f.write(json.dumps(history_entry) + "\n")
print(f"Appended → {HISTORY_FILE}")
