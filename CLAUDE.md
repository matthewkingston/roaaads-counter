# Newtownards Traffic Model — Project Overview

A gravity-model traffic assignment pipeline for Newtownards, calibrated against
walking count data and official AADT figures. The pipeline is fully reproducible:
running the scripts in order regenerates all outputs from raw data.

**Agent instruction:** Keep this file up to date. After any tuning run, count data
ingest, model change, or reference value update, edit the relevant sections before
committing. This file is the authoritative record of model state.

---

## Pipeline (run in this order)

```
python3 simulation/build_network.py          # build road network from OSM
python3 simulation/build_demographics.py     # node weights + map scaffold
python3 simulation/build_paths.py            # precompute k=3 shortest paths (~30 min)

python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/tune_assignment.py                        # Stage 1: tune gravity params (4 params incl. THETA, ~30 s)
python3 analysis/tune_assignment.py --full                 # Stage 2: + external zones (24 params, ~2-3 min)
python3 analysis/tune_assignment.py --note "description"   # optional human label in history

python3 simulation/build_assignment.py       # apply tuned params, write flows
python3 simulation/build_demographics.py --map-only   # rebuild map HTML only (fast)
```

After adding new count data, re-run from `ingest_counts.py` onward. The tuner
reads `simulation/tuned_params.json` as its starting point, so repeated runs
refine rather than restart.

---

## Key Files

| File | Role |
|------|------|
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, assigns external zone weights, builds map. `--map-only` skips demographic recomputation and rebuilds only the HTML. `--zones-only` patches boundary node weights without rebuilding. External zone pop/wp/damping values are read from `tuner_config.json` (lat/lon centroids remain hardcoded in the script). **Population distribution:** for each DZ, OSM building centroids (residential tags + `addr:housenumber`, cached in `data/cache_osm_buildings.geojson`) are snapped to the nearest road edge via Shapely STRtree; population is split between the two endpoint nodes proportional to the building's position along that edge. DZs with <3 buildings or pop/building ≥ 10 fall back to road-length weighting. A per-DZ quality table is printed on each run. **Business demand:** NISRA workplace population is allocated by POI count within each DZ, then augmented with OSM car park polygon area (cached in `data/cache_osm_parking.geojson`): public car parks add `area/25` equivalent persons per node; `access=private` car parks add `area/50` (half weight). This gives W_BIZ a spatially concentrated signal from major retail car parks (Ards Shopping Centre, Castlebawn, supermarkets) rather than relying solely on workplace headcount. The map includes a parking layer (blue = public, red = private) for realism checking. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Re-run if road network changes or `HIGHWAY_COST_FACTOR` values change. Edge costs are travel time × a road-class multiplier (trunk/primary: ×0.67, residential/unclassified/living_street: ×1.2, others: ×1.0) to bias routing toward major roads. Also reads `tuner_config.json` to filter which external→external OD pairs to include as through routes. Produces **k=3 alternative paths** per OD pair (k=2/k=3 via progressive edge penalisation ×10) for stochastic logit routing. Build time ~30 min (3 Dijkstra passes). |
| `simulation/model.py` | **Shared constants and functions** imported by both `build_assignment.py` and `tune_assignment.py`: `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (rational kernel, used by `build_assignment.py` only), `site_flow()`, `compute_chi2()` (Woodbury-corrected, used by `build_assignment.py` only), `print_chi2_table()` (used by both). **Note:** `tune_assignment.py` has its own parallel `calibrate_K()` and vectorized chi² loop (in `objective()`) that implement the same Woodbury algebra as `model.py`'s `compute_chi2()` — but using precomputed arrays for speed. If you change the Woodbury formula or observation-weighting logic in one place, change it in both. |
| `simulation/build_assignment.py` | Gravity model assignment. Loads `tuned_params.json` automatically if present. Uses stochastic logit routing (k=3 paths, THETA parameter) if the paths cache and tuned params both support it; otherwise falls back to all-or-nothing. Prints Woodbury-corrected χ²/N goodness-of-fit table via `model.py`'s `compute_chi2()`. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Reference values for L2 regularization, city→node groupings, `through_route_pairs` whitelist, and gravity param regularization. `lambda` regularises external zones; `gravity_lambda` + `gravity_ref` regularise P/ALPHA/W_BIZ/THETA toward physically plausible values (prevents K-drift pathology). Default P=300 s sets the peak travel time; ALPHA=2 gives 1/d² tail decay; THETA=1.0 is the logit dispersion anchor. Edit to change external zone priors, allowed through routes, or gravity anchors. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry now carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s` so the tuner can work in count space. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. When the paths cache has k=3 alternative paths, Stage 1 tunes 4 gravity params (W_BIZ, P, ALPHA, THETA); otherwise 3. `--full` adds 14 city pop/wp + 6 dampings. Uses per-session observations from `link_aadt.json`. K and per-slot hourly fractions {f_s} are jointly calibrated at each evaluation via alternating minimisation; walking obs in count space, official sites in AADT space. Per-slot fraction priors derived from `hourly_fractions.csv` grouped by weekday/Saturday/Sunday. **Performance:** all-or-nothing (THETA=None/old cache): ~0.12 ms/eval stage 1, ~5.5 ms stage 2 via precomputed bin-matmul path. Stochastic (THETA given + k=3 paths): ~150 ms/eval via CSR SpMV scatter; ~1 min stage 1, ~8 min stage 2. Startup builds three sparse (N_links×N_OD) matrices (~0.6 s one-off). |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. |
| `simulation/reset_gravity_params.py` | Reset only the gravity params (K, W_BIZ, P, ALPHA, THETA) in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. External zone params are preserved. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile (fraction of AADT per hour×day-of-week). Used for AADT estimation from short-duration counts. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_flows.json`, `simulation/newtownards_map.html`,
`simulation/tuned_params.json` — all regenerated by the pipeline.

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`simulation/tuner_config.json` — committed as source config (gitignore exception).
`analysis/hourly_fractions.csv` — committed as source data (single authoritative version).

### Large reference data (gitignored, kept locally only)
`data/*.ods`, `data/*.xlsx` — too large to commit; keep a local copy for reference.
Currently present: `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods`
(the 3 AADT values already in use are in `simulation/model.py` `COUNT_SITES` — no
further data from this file is needed) and `data/census-2021-apwp001.xlsx`.

---

## Model Design

### Gravity model
OD flow: `T_ij = K × w_i × w_j × f(d_ij)`

Rational kernel: `f(d) = (ALPHA+1) × P^ALPHA × d / (ALPHA × P^(ALPHA+1) + d^(ALPHA+1))`

Equivalently (numerically stable): `u = d/P; f(d) = (ALPHA+1) × u / (ALPHA + u^(ALPHA+1))`

Properties: f(P) = 1 (peak always at d = P seconds), f(0) = 0, tail ~ 1/d^ALPHA for large d.

Node weight: `w = population + W_BIZ × business_demand`

Distances are least-time shortest paths (seconds), with an off-network leg added
for boundary nodes using their real-world centroid position.

### Stochastic route choice (logit)
When the paths cache contains k=3 alternatives and THETA is a tuned parameter,
demand is split across 3 paths per OD pair:

`share(path r) ∝ exp(−THETA × d_r / P)`

THETA → ∞: collapses to all-or-nothing (k=1 path only).
THETA = 0: equal split across all 3 paths.

Alternative paths k=2/k=3 are found by penalising k=1 (and k=1+k=2) edges ×10
in the Dijkstra adjacency matrix. Pairs with no alternative fall back to k=1
(which is equivalent to all-or-nothing for those pairs under any THETA).

### External zones
14 boundary nodes grouped into 7 cities in `tuner_config.json`. Each city shares
one pop and one wp value; individual nodes scale by a fixed damping ratio.
Nodes with damping=1.0 are fixed; nodes with damping<1.0 are tunable (but
L2-penalised toward the config reference values).

### Through routes
External→external OD pairs are included for a whitelisted set of city pairs
stored in `tuner_config.json` under `through_route_pairs`. All other
boundary→boundary pairs are excluded. The whitelist captures trips that genuinely
traverse Newtownards (e.g. Comber↔Donaghadee via A22→A48) while excluding trips
that use an out-of-network bypass (e.g. Belfast↔Bangor via the A2 coast road).
Changing the whitelist requires rebuilding the paths cache (`build_paths.py`).

Current whitelist: Comber↔Donaghadee, Comber↔LowerArds, Comber↔Millisle,
Comber↔Bangor, Bangor↔LowerArds, Belfast↔LowerArds.

### K and slot fractions: joint analytical calibration
At each optimizer evaluation, K (global scale) and per-slot hourly fractions {f_s} are
jointly calibrated via alternating minimisation (10 iterations, converges in 3–5).

**K-step** (for fixed {f_s}): K = B / A where A and B accumulate official AADT terms
(AADT space) and walking count terms (count space, weight m·T·f_s).

**f_s-step** (for fixed K): each f_s = (Σ C_i + mean_f/std_f²) / (Σ C_i²/n_eff_i + 1/std_f²)
where C_i = K · m_i · T_i/3600 and the second term is the Gaussian prior.

Slot key is (day_type, hour): day_type = 0 (weekday, Mon–Fri), 1 (Saturday), 2 (Sunday).
Prior hyperparameters from `hourly_fractions.csv` via law of total variance:
  std_f = sqrt(between_day_var + mean(within_day_var across dows in slot)).

### Goodness of fit
`χ²/N` (mean squared z-score, target ~1.0). Three official AADT sites (±10%) plus
all per-session walking-count observations from `link_aadt.json` observations lists.

Walking obs compared in count space: chi²_walk = Σ (K·m·T·f_s/3600 − n_eff)² / n_eff,
where n_eff = n + 0.5 (Jeffreys count). Shared f_s per (day_type, hour) slot is the
correlation mechanism (all obs in a slot scale together), anchored by a Gaussian prior.
Official sites compared in AADT space: chi²_off = Σ (K·m − y)² / σ².
Total objective: chi² = chi²_off + chi²_walk + Σ (f_s − mean_f)² / std_f².

`N_eff = N − N_slots` is printed as a diagnostic (each slot loses one effective df).
With weekday/Sat/Sun grouping, 15 original (dow, hour) slots typically collapse to ~10
(day_type, hour) slots, yielding N_eff = N − N_slots.

**Note:** `model.py` / `build_assignment.py` still use the older Woodbury formulation
for `compute_chi2()` — their reported χ²/N will differ slightly from the tuner's.

---

## Count Data

**Official AADT sites (fixed):**
- Site 507: A21 Bangor Road — 21,202 AADT
- Site 508: A48 Donaghadee Road — 10,792 AADT
- Site 444: A20 Portaferry Road — 7,282 AADT

**Walking counts:** 5 CSV files, 130 sessions, 256 individual session-direction observations across 136 directed links. Day-type grouping (weekday/Sat/Sun) of the 15 original (dow, hour) time slots yields fewer (day_type, hour) slots (around 10), increasing N_eff. Combined with 3 official sites: N=258 total tuner observations. The tuner uses per-session observations directly (not per-link aggregates); per-link aggregates are retained in `link_aadt.json` for reference.

---

## Tuning History

| Date | Stage | N obs | N params | χ²/N | Notes |
|------|-------|-------|----------|------|-------|
| 2026-06-14 | gravity | 15 | 4 | 5.72 | |
| 2026-06-14 | gravity | 25 | 4 | 3.97 | |
| 2026-06-14 | full | 25 | 24 | 0.98 | |
| 2026-06-14 | full | 25 | 24 | 0.90 | |
| 2026-06-14 | full | 62 | 24 | 0.956 | road-class routing, Hardford Link primary, excl 161→160 |
| 2026-06-15 | full | 62 | 24 | **0.895** | + through routes (6 city pairs); refs updated |
| 2026-06-15 | gravity | 109 | 4 | 2.346 | + 4th count session (107 directed links); per-link agg, no Woodbury |
| 2026-06-15 | gravity | 161 | 4 | 2.00 | Woodbury correction; per-session obs (N_eff=151, 10 slots) |
| 2026-06-15 | full | 161 | 24 | 1.1754 | Jeffreys v3 reprocess; paths cache stale (no through-routes) |
| 2026-06-15 | gravity | 161 | 4 | 1.1687 | rebuilt paths cache with through-routes (+56 OD pairs) |
| 2026-06-15 | full | 161 | 24 | 1.1207 | through-routes active; LowerArds resolved (+92% not +514%) |
| 2026-06-16 | gravity | 161 | 3 | 1.1565 | rational kernel (P/ALPHA replaces MU/SIGMA/ALPHA) |
| 2026-06-16 | full | 161 | 23 | **1.0833** | rational kernel full tune; P=190s, ALPHA=4.88 |

Current best full-tune: χ²/N = 1.0833 (161 obs, N_eff=151, 10 slots; rational kernel).
Official sites: 507 z=−1.38, 508 z=−0.46, 444 z=+0.50.
Persistent structural outliers: `719↔325` Messines Road (z=−4.36/−3.05), `328→326` Comber Road (z=−3.56), `18→21` Hardford Link (z=−3.15). All are internal-traffic corridors; external zone tuning cannot resolve them.
All cities tuning well above tuner_config refs (Donaghadee +184%, Belfast +284%, LowerArds +226%, Bangor +30%) — refs need updating before next run.
`22→159` (model=0) was a data error (snap direction bug, fixed 2026-06-15): now recorded as 159→22.
Belfast Road `20→18` zero-count obs gives z=−2.80 (obs=628, model=2,277) — persistent underprediction of internal Circular Road traffic.

### Paths cache note
The paths cache (`newtownards_paths.npz`) must be rebuilt with `build_paths.py` whenever
`through_route_pairs` changes or whenever stochastic routing (THETA) is to be used.

**Stochastic routing requires a fresh cache rebuild.** The cache from 2026-06-15 contains
only k=1 paths — it predates the k=3 alternative path feature. Until `build_paths.py` is
re-run, `_has_stoch = False` and THETA is not included in the parameter space; tuning
proceeds as all-or-nothing. After rebuilding, the cache will contain `pair_idx_2`/`pair_idx_3`
keys and THETA becomes a 4th gravity parameter automatically.

The cache previously lacked all through-routes despite the whitelist being correct — the
cache predated the through-route feature. This caused the LowerArds pop to blow up to
+514% as the tuner compensated for missing through-traffic. After rebuilding (2026-06-15)
LowerArds settled at +92%.

### Known model behaviour
- `W_BIZ` was converging to ~0 when `node_business_demand` was based solely on NISRA workplace population (spatially correlated with residential density). As of 2026-06-17, `node_business_demand` now includes OSM car park polygon area (public: area/25, private: area/50 equivalent persons), which adds a strongly spatially concentrated signal from major retail car parks. W_BIZ convergence after next tuning run should confirm whether this resolves the identifiability issue.
- `K` is analytically set at each optimizer step to rescale the raw flow field to
  match observed AADT. It absorbs the overall scale of unnormalized gravity flows,
  which shifts by many orders of magnitude as ALPHA and P change (e.g. ALPHA 2→5
  changes the tail by ~10^9 for a typical d=1000s path). The degeneracy is
  between K and the gravity parameters, not K and the external zone populations
  (which only vary O(100%) under L2 regularization and can only contribute O(4×)
  to raw flows). χ²/N is reliable; K is not interpretable in isolation.
- After a structural model change (e.g. adding through routes or new count data), a gravity-only
  stage 1 run with fixed external zone weights will show inflated χ²/N and
  spurious outliers at boundary sites (esp. site 508). A full 24-param re-tune
  is needed to restore fit quality.
- **Snap direction bug (fixed 2026-06-15):** `ingest_counts.py` stored canonical
  `(min(u,v), max(u,v))` indices in `edge_geoms` instead of the actual directed
  `(u, v)`. For one-way roads where u > v (e.g. 159→22), the dot-product sign was
  correct but `link_with` was flipped. Fix: store `(u, v)` not `(pair[0], pair[1])`.
  Only `f56b2ce4` was materially affected — re-snapped from 22→159 to 159→22.
- Slot fractions f_s are inferred per (day_type, hour) — grouping weekday / Saturday / Sunday.
  They are anchored by a Gaussian prior derived from `hourly_fractions.csv` via the law of
  total variance (between-day variance + mean within-day variance). The prior prevents runaway
  values in slots with sparse data. All obs in a slot share the same f_s, preserving the
  within-slot correlation that Woodbury previously handled. `N_eff = N − N_slots` is the
  effective degrees of freedom. With day_type grouping, N_slots is typically ~10 vs 15 with
  the old (dow, hour) formulation, gaining ~5 effective degrees of freedom.

---

## External Zone Reference Values (`tuner_config.json`)

**These values must not be changed without explicit user approval.** After a full tuning run,
updating them is something to *consider and discuss*, not an automatic step — the refs anchor
L2 regularization and changing them shifts the penalty basin for all future runs.

Last updated: 2026-06-16 full tuning run (χ²/N=1.2546, 258 obs, stochastic k=3).

| City | Nodes | ref_pop | ref_wp | Tunable dampings |
|------|-------|---------|--------|-----------------|
| Donaghadee | 47 | 190,201 | 7,018 | — |
| Comber | 65, 617, 618, 620 | 53,571 | 2,996 | 617 (×0.38), 618 (×0.35), 620 (×0.43) |
| LowerArds | 92 | 84,500 | 5,024 | — |
| Belfast | 97, 119 | 1,034,719 | 183,661 | 119 (×0.31) |
| Bangor | 98, 731 | 95,426 | 21,246 | 98 (×0.39) |
| Holywood | 99 | 3,652 | 1,203 | — |
| Millisle | 748, 749 | 2,570 | 498 | 749 (×0.47) |

---

## Worktree Convention

Background/parallel Claude work is done in `.claude/worktrees/` and
cherry-picked to `main` after review. All work is on `main`.
