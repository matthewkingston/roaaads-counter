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
python3 simulation/build_paths.py            # probit stochastic paths (N_PASSES=25, CV=0.25; ~30-60 min)

python3 analysis/parse_official_hourly.py    # parse ODS hourly counts → data/official_hourly.json (one-off)
python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/tune_assignment.py                        # Stage 1: tune gravity params (6 params, no THETA; ~14s after probit cache rebuild)
python3 analysis/tune_assignment.py --full                 # Stage 2: + external zones (26 params, ~40s after probit cache rebuild)
python3 analysis/tune_assignment.py --fast                 # looser tolerances + fewer alt-min iters (~2× faster, minimal precision loss)
python3 analysis/tune_assignment.py --note "description"   # optional human label in history

python3 simulation/build_assignment.py       # apply tuned params, write flows
python3 simulation/build_demographics.py --map-only   # rebuild map HTML only (fast)

python3 analysis/report_tune.py              # generate reports/ text + pull plot for last run
```

After adding new count data, re-run from `ingest_counts.py` onward. Re-run
`parse_official_hourly.py` only if the ODS source file changes. The tuner reads
`simulation/tuned_params.json` as its starting point, so repeated runs refine rather
than restart.

---

## Key Files

| File | Role |
|------|------|
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, assigns external zone weights, builds map. `--map-only` skips demographic recomputation and rebuilds only the HTML. `--zones-only` patches boundary node weights without rebuilding. **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Business demand:** workplace population × POI count, augmented with OSM car park polygon area (public: area/25, private: area/50 equivalent persons). **Flow map layers:** reads `newtownards_flows.json` and adds: (1) combined AADT layer (shown by default, blue→yellow→red, tooltip includes res/biz breakdown if component flows present); (2) residential layer (teal, off by default); (3) business-adjacent layer (amber→red, off by default). Component layers only appear if `flows_res`/`flows_biz` keys exist in the flows file. Map also includes parking (blue/red) and POI layers. **Adding a new external zone:** add an entry to `_EXT_GEO` (the dict at the top of the file) and to `tuner_config.json` cities — `BOUNDARY_NODE_IDS` is now derived automatically from `_EXT_GEO.keys()` so no separate update is needed there. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Re-run if road network changes, `HIGHWAY_COST_FACTOR` values change, or `N_PASSES`/`PROBIT_CV` change. Edge costs are travel time × a road-class multiplier (trunk/primary: ×0.67, residential/unclassified/living_street: ×1.2, others: ×1.0). Also reads `tuner_config.json` to filter external→external OD pairs. **Probit stochastic loading:** runs `N_PASSES=25` Dijkstra passes with log-normal edge-cost noise (CV=0.25), accumulates fractional link-assignment weights (`link_weight` per entry = fraction of passes using that link for that OD pair). `od_dist` is the mean path distance across passes. Build time ~30–60 min. |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (rational kernel, `return_components=True` returns `(flow_res, flow_biz)` tuple), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `gravity_assign()` accepts optional `link_weight` array (float, parallel to `pair_idx`/`link_idx`): when provided, each entry's flow contribution is scaled by its fractional weight (probit loading); when `None`, falls back to binary all-or-nothing. `compute_chi2()` has two modes: **two-component** (when `link_flow_biz_dict` provided) uses 216 official hourly obs + count-space formula matching the tuner; **legacy** (single flow dict) uses 3 AADT obs + Woodbury correction. |
| `simulation/build_assignment.py` | Gravity model assignment. Requires `simulation/newtownards_paths.npz` (exits with an informative error if absent). Calls `gravity_assign(return_components=True)` and `compute_chi2` two-component mode when K_res/K_biz are in `tuned_params.json`. Saves `flows_res` and `flows_biz` alongside `flows` in `newtownards_flows.json` when two-component. Falls back to legacy single-K mode for old param files. **Two-component mode does not require `official_hourly.json` to be present** — `compute_chi2` handles a missing file gracefully; the res/biz assignment and map layers work from K_res/K_biz alone. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Reference values for L2 regularization, city→node groupings, `through_route_pairs` whitelist, gravity param regularization, and temporal coupling. `lambda` regularises external zones; `gravity_lambda` + `gravity_ref` regularise P/ALPHA/BETA/W_BIZ/P_biz/ALPHA_biz; `gamma_coupling_scale` controls the per-slot aggregate coupling (γ = scale/std_f²); `phi_prior` + `phi_std` set the Gaussian prior on the business flow fraction (see Model Design). |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. Run once (or when the ODS file changes). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. **`MANUAL_LINK_OVERRIDES`** dict at the top maps session IDs to forced directed links; use when the observer was on a parallel carriageway and GPS snap would land on the wrong road (e.g. A20 Kempe Stones sessions `e644eae2`/`760b0c8e` → link 8→7). After every new link assignment (manual or auto), validates that each non-null count direction maps to a real directed edge in G; raises `ValueError` if not, preventing counts from silently hitting a zero-flow phantom edge. Edges without geometry (virtual stub nodes such as Dundonald 10000) are skipped during snap candidate construction. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry now carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s` so the tuner can work in count space. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Two-component model:** gravity flows split into residential (pop→pop, `flow_res`) and business-adjacent (pb+bb, `flow_biz`) components, each with its own temporal profile (f_s_res, f_s_biz) and scale (K_res, K_biz). Stage 1 tunes 6 gravity params (W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz); `--full` adds 14 city pop/wp + 6 dampings (26 params total). Uses 216 official hourly obs from `official_hourly.json` plus per-session walking obs from `link_aadt.json`. **Alternating minimisation (4 blocks, up to 10 iters, early-exit on K convergence):** K-step (1D, total scale); phi-step (1D, K_biz/K ratio, Gaussian prior N(0.35, 0.15²) prevents K_biz→0); f_res-step (per-slot analytical); f_biz-step (symmetric) + aggregate coupling γ·(f_res+f_biz−f_agg)² per slot. **Performance:** With probit cache (build_paths.py rebuilt): distance-bin matmul path active (~100× faster run_assignment) + vectorised calibrate_Ks_and_fracs (28× faster, replaces 72-slot Python loops with bincount over 559 obs) → ~10–12 ms/eval; ~14 s stage 1, ~40 s stage 2. With old k=3 cache (current): ~250 ms/eval. **`--fast` mode:** ftol/xtol loosened 5×, alt-min capped at 5 iters for intermediate evals; ~2× fewer optimizer iterations with negligible change to final χ²/N. Recorded in history under `tuner_hyperparams.fast`. |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` (summary, chi² by measurement, chi² by link, gravity params with K_res/K_biz/phi, external zones, slot fractions table showing f_res and f_biz side by side with component pulls) and `reports/slot_pulls_{id}.png` (two side-by-side heatmaps: 24 h × 3 day-types per component, colour = pull from prior). History `slot_prior` entries carry 4 values: `[mean_f_agg, std_f, mean_f_res, mean_f_biz]`. History `tuner_hyperparams` carries phi_prior, phi_std, gamma_coupling_scale, gravity_lambda, lambda at run time. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. Writes all keys from `params` dict (K_res, K_biz, slot_fracs_res, slot_fracs_biz, external zone params etc.); clears stale `slot_fracs` legacy key. |
| `simulation/reset_gravity_params.py` | Reset only the gravity params (K, W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz) in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. External zone params are preserved. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile (fraction of AADT per hour×day-of-week). Includes `mean_fraction_res` and `mean_fraction_biz` columns: temporal profile priors for the two-component model. **Derived from NTS Table NTS0502a** via `analysis/derive_component_profiles.py` (see below). Re-running `hourly_fractions.py` preserves these columns automatically; do not edit the aggregate columns by hand. **Summation convention:** `mean_fraction[D,H] = count[D,H] / AADT` where `AADT = weekly_total / 7`. Consequently the 168 rows sum to **7.0** (not 1.0): each day sums to that day's traffic relative to AADT (Mon ≈ 1.00, Fri ≈ 1.12, Sun ≈ 0.77). This is intentional — the day-of-week volume differences are encoded in the magnitude. The AADT-weighted business share across the whole week is `Σ(mean_fraction_biz) / Σ(mean_fraction) = 3.07 / 7.0 ≈ 44%`. |
| `analysis/derive_component_profiles.py` | Derives `mean_fraction_res` and `mean_fraction_biz` from DfT NTS data (2023–2024 rolling average). **Weekdays (NTS0502a, `data/nts0502.ods`):** per-hour `biz_share(h)` = (commuting + employer's business + education×⅕ + escort education + shopping) / adjusted_all. **Education is downweighted by ÷5**: the NTS records all modes; the standalone "education" trip is predominantly the child travelling, with far lower car trip generation than commuting — the adult car driver is already counted under "escort education". The denominator shrinks by the same amount so all other shares are proportionally increased. **Weekends (NTS0504b, `data/nts0504.ods`):** Saturday and Sunday each get a distinct flat `biz_share` from the actual day-of-week trip rates (2023–2024), with "Just walk" trips removed from the denominator (vehicle counts contain no pedestrian trips) and the same education ÷5 applied. Results: Saturday 36.8% business, Sunday 31.2%. The temporal shape within each weekend day comes from the aggregate `mean_fraction` profile. Re-run whenever the NTS files or purpose classification change. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_map.html`, `simulation/tuned_params.json` — all regenerated by the pipeline.
`simulation/newtownards_flows.json` — combined flows plus optional `flows_res`/`flows_biz` keys when two-component params active.
`reports/` — generated by `report_tune.py` and `tune_assignment.py`; not tracked.

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`, `data/official_hourly.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`simulation/tuner_config.json` — committed as source config (gitignore exception).
`analysis/hourly_fractions.csv` — committed as source data (single authoritative version).

### Large reference data (gitignored, kept locally only)
`data/*.ods`, `data/*.xlsx` — too large to commit; keep a local copy for reference.
Currently present:
- `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods` — used by `parse_official_hourly.py` to extract hourly counts for sites 444/507/508; annual AADT values in `model.py` `COUNT_SITES` are no longer used by the tuner but retained for `build_assignment.py`.
- `data/nts0502.ods` — DfT NTS Table NTS0502a, "Trip start time by trip purpose (Mon–Fri): England, 2002 onwards". Used by `derive_component_profiles.py` for weekday hourly business shares. Download from GOV.UK NTS data tables (NTS05 Trips).
- `data/nts0504.ods` — DfT NTS Table NTS0504b, "Average trips by day of the week and purpose: England, 2002 onwards". Used by `derive_component_profiles.py` for Saturday/Sunday business shares. Same download page.
- `data/census-2021-apwp001.xlsx` — census workplace population data.

---

## Model Design

### Gravity model
OD flow: `T_ij = K × w_i × w_j × f(d_ij)`

Generalised rational kernel: `u = d/P; f(u) = (ALPHA+BETA) × u^BETA / (ALPHA + BETA × u^(ALPHA+BETA))`

Properties: f(P) = 1 (peak always at d = P seconds, for any positive ALPHA, BETA), f(0) = 0,
tail ~ 1/d^ALPHA for large d, rise ~ u^BETA near origin.
BETA=1 recovers the original kernel `(ALPHA+1) × u / (ALPHA + u^(ALPHA+1))`.
ALPHA controls the right-side tail decay; BETA controls the left-side approach to the peak.

Node weight: `w = population + W_BIZ × business_demand`

Distances are least-time shortest paths (seconds), with an off-network leg added
for boundary nodes using their real-world centroid position.

### Stochastic route choice (probit loading)
The paths cache stores fractional link-assignment weights computed from `N_PASSES=25`
Dijkstra runs, each with log-normal edge-cost noise (CV=0.25). For each OD pair,
`link_weight[entry]` is the fraction of passes that routed through that link. Pairs
with no topological route diversity (degree-1 stubs, single-access nodes) converge to
weight=1.0 on their forced route. `od_dist` is the mean path distance across passes.

This replaces the previous k=2/k=3 global-penalisation scheme, which was ineffective
in a dense network (global penalisation of all k=1-used links preserves relative ordering
and produces identical alternative paths for most OD pairs). THETA is no longer tuned.

### External zones
15 boundary nodes grouped into 8 cities in `tuner_config.json`. Each city shares
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
Comber↔Bangor, Bangor↔LowerArds, Belfast↔LowerArds, Dundonald↔LowerArds,
Dundonald↔Donaghadee, Dundonald↔Millisle.

### Two-component flow decomposition
The gravity OD flows are split into two spatial components at each tuner evaluation:

- **Residential** (`flow_res`): `all_bin_pp @ f_b` — purely pop×pop trips
- **Business-adjacent** (`flow_biz`): `W_BIZ·(all_bin_pb @ f_b) + W_BIZ²·(all_bin_bb @ f_b)` — home↔work/retail trips

Each component has its own temporal profile (f_s_res, f_s_biz) and scale (K_res, K_biz).
Predicted count for observation i in slot s:
`pred_i = K_res · flow_res[link_i] · (T_i/3600) · f_s_res[s] + K_biz · flow_biz[link_i] · (T_i/3600) · f_s_biz[s]`

### Four-block analytical calibration
At each optimizer evaluation, (K, phi, f_s_res, f_s_biz) are calibrated via alternating
minimisation (up to 10 iterations, early exit on K convergence; typically 3–5 suffice),
where K_res = K·(1−phi), K_biz = K·phi. `--fast` caps at 5 iterations for intermediate evals.

**K-step:** 1D solve (same structure as old single-K step), using combined coefficient
`(1−phi)·c_r·f_res + phi·c_b·f_biz` per observation.

**phi-step:** 1D solve for business fraction, with Gaussian prior phi ~ N(0.35, 0.15²).
This prior prevents K_biz collapsing to 0 (the K_biz/W_BIZ degeneracy otherwise exploited
by the optimizer). phi = 0 means all-residential; phi = 1 means all-business.

**f_res-step / f_biz-step:** per-slot analytical update, anchored by component-specific
priors from `hourly_fractions.csv` columns `mean_fraction_res` / `mean_fraction_biz`
(synthetic profiles: business peaks weekday AM/PM; residential is complement).

**Aggregate coupling:** each slot also carries a penalty γ·(f_res + f_biz − f_agg)²
where γ = `gamma_coupling_scale` / std_f_agg² (per-slot, from `tuner_config.json`).
This prevents the two profiles collectively drifting from the known NI aggregate profile.
`gamma_coupling_scale = 0.0` disables coupling; `1.0` gives coupling as strong as the
individual component priors.

Slot key: (day_type, hour), day_type = 0 (weekday), 1 (Saturday), 2 (Sunday).
Prior std derived from `hourly_fractions.csv` via law of total variance.

### Observations
All 559 observations are in count space with per-obs weights:
- **Official hourly** (216 obs, 24 h × 3 day-types × 3 sites): from `data/official_hourly.json`;
  Gaussian error (sigma from between-weekday std, 10% floor); weight = 1/sigma².
- **Walking** (343 obs): from `data/link_aadt.json`; Poisson error; weight = 1/n_eff.

### Goodness of fit
`χ²/N` (mean squared z-score; N=559 obs, N_eff = N − 2·N_slots = 559 − 144 = 415).
Two df lost per slot (one each for f_s_res and f_s_biz). With coupling enabled,
chi²/N includes the coupling penalty terms; pure data-fit chi²/N is lower.

`build_assignment.py` uses the two-component `compute_chi2()` when K_res/K_biz are present in tuned_params.json. This gives a **data-only** chi²/N (pure sum of squared z-scores) — it excludes the f-prior penalties `(f_r−mfr)²/std_f²` and the aggregate coupling penalty that the tuner includes in its chi²/N. Expect the build_assignment chi²/N to be somewhat lower than the tuner's; the two are directionally comparable but not numerically equal. The legacy Woodbury path is used only for old single-K param files.

---

## Count Data

**Official hourly obs (from ODS, parsed by `parse_official_hourly.py`):**
- Site 507: A21 Bangor Road — 72 obs (24 h × 3 day-types); Gaussian sigma, 10%/15% floors
- Site 508: A48 Donaghadee Road — 72 obs
- Site 444: A20 Portaferry Road — 72 obs

Annual AADT values (507: 21,202; 508: 10,792; 444: 7,282) are retained in `model.py`
`COUNT_SITES` for `build_assignment.py` backward compatibility but are no longer used
by the tuner directly.

**Walking counts:** 7 CSV files, 177 sessions, 343 per-session observations (after EXCLUDE_LINKS). Two sessions manually assigned: `e644eae2` and `760b0c8e` (A20 Kempe Stones eastbound, link 8→7; observer was on the westbound carriageway). The tuner uses per-session observations directly; per-link aggregates are retained in `link_aadt.json`.
New sessions added 2026-06-18 (7 sessions): Saratoga Avenue (333↔335), Glenford Road (332→331, 67→21), Hardford Link (21→18), Belfast Road (18→20).
**Total: 559 observations (216 official hourly + 343 walking) in 72 time slots. N_eff = 559 − 2×72 = 415.**

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
| 2026-06-17 | gravity | 374 | 4 | 1.3064 | two-component model; 216 ODS hourly obs; K_res=1.47e-05, K_biz=4.43e-06 (phi≈23%), γ=0 |
| 2026-06-17 | gravity | 374 | 4 | 1.4286 | + aggregate coupling γ=1/std_f²; K_res=8.96e-05, K_biz=1.29e-05 (phi≈13%) |
| 2026-06-17 | gravity | 545 | 4 | 1.9582 | + new count data (329 walking obs); sqrt(count) sigma floor active |
| 2026-06-17 | full | 545 | 26 | 1.6640 | first two-component full tune; phi=16.5%; LowerArds wp +1303% flag |
| 2026-06-18 | full | 545 | 26 | 1.6432 | NTS-derived component priors; LowerArds wp +645% (improved); Belfast wp +1083% new flag |
| 2026-06-19 | full | 559 | 28 | **1.3742** | first probit-cache tune; P=125s, ALPHA=4.10; phi=30.4%; city refs updated |

**Note on comparability:** runs from 2026-06-17 onward use the two-component model with coupling penalty terms in chi²/N; not directly comparable to earlier single-component runs. From 2026-06-18 count ingest onward: 559 observations (216 official hourly + 343 walking, 72 time slots, N_eff=415).

Current best full-tune: chi²/N = 1.3742 (559 obs, N_eff=415; two-component with probit cache, run 8a0fe24b).
W_BIZ=3.123, P=125.1s, ALPHA=4.10. phi=30.4% business fraction.
mean|z|=0.81  |z|>2: 42  |z|>3: 13.
Probit cache (N_PASSES=25, CV=0.25) active from 2026-06-19.

**Confirmed working:**
- Temporal profiles separating meaningfully: business peaks weekday h06 earlier than residential (Δ/σ_biz=+1.54 vs −1.39); overnight business fraction higher (deliveries/early commuters).
- Site 444 overnight z-scores improved: previously z≈−5 at h04; now worst official-hourly is h06 at z=−3.04.
- Map layers (residential/business) confirmed — `build_assignment.py` populated `flows_res`/`flows_biz`.

**Outstanding concerns:**
- **Belfast wp = 1,822,862 (ref 183,661, +891%)** — persistent external-zone inflation. P improved from 57s → 125s with gravity_lambda=0.5 but Belfast wp still strongly inflated. Further tightening or a Belfast-specific lambda may be warranted.
- Structural outliers: `22→12 Regent Street` (z=+4.36), `23→295 Frances Street` (z=+3.98), `296→297 Nursery Road` (z=−3.61), `139→137 Portaferry Road` (z=−3.64) — not necessarily model failures.
- `73→70` Mill Street severe underprediction (z=−3.33; obs 23,644 vs model 2,222).
- `719→325` / `325→719` Messines Road persistent (z=−3.13/−2.83).
- `18→21` / `68→21` / `21→68` Hardford Link persistent (z=−3.28/−3.28/−2.78).

### Paths cache note
The paths cache (`newtownards_paths.npz`) must be rebuilt with `build_paths.py` whenever
`through_route_pairs` changes, the road network changes, or `N_PASSES`/`PROBIT_CV` change.

**Current cache format** (probit): contains `link_weight` (float32, per entry), `od_dist`
(mean across passes), `probit_n_passes`, `probit_cv`. No `pair_idx_2/3` keys —
`_has_stoch = False`, THETA not tuned. Rebuild with `build_paths.py` if the road network
or `through_route_pairs` changes.

### Known model behaviour
- **Two-component K_biz/W_BIZ degeneracy:** Without the phi prior, the optimizer exploits
  `K_biz × W_BIZ → 0 × ∞` to collapse K_biz to zero while using large W_BIZ to reshape
  the spatial flow. The phi prior phi ~ N(0.35, 0.15²) prevents this by anchoring the
  business fraction. phi ≈ 0.13–0.23 in current runs.
- `W_BIZ` was converging to ~0 when `node_business_demand` was based solely on NISRA
  workplace population. As of 2026-06-17, it includes OSM car park polygon area
  (public: area/25, private: area/50 equivalent persons). W_BIZ=1.20 in current runs
  (previously blowing up to ~10 when K_biz was unconstrained).
- `K` (total scale) is analytically calibrated at each optimizer step, absorbing the
  overall magnitude of unnormalised gravity flows (shifts by many orders of magnitude as
  ALPHA/P/BETA change). K_res and K_biz are derived from K × (1−phi) and K × phi respectively.
  K is not interpretable in isolation; chi²/N is reliable.
- After a structural model change (e.g. adding through routes or new count data), a gravity-only
  stage 1 run will show inflated chi²/N. A full `--full` re-tune is needed to restore fit quality.
- **Dundonald virtual node (added 2026-06-17):** Node 10000 is a degree-1 stub connected
  only to node 97. With probit loading, its `link_weight` values will be 1.0 on the forced
  stub edge for every pass (no route diversity possible for a degree-1 node). Route diversity
  for Dundonald OD pairs is zero, which is topologically correct.
- **Manual link overrides:** `MANUAL_LINK_OVERRIDES` in `ingest_counts.py` hard-assigns specific sessions to a directed link, bypassing GPS snap. Use when the observer stood on a parallel carriageway (e.g. a dual one-way road) and the snap would land on the wrong physical road. The override is idempotent and takes effect even if `counts_processed.json` is wiped and rebuilt. After assignment (manual or auto), the script validates each non-null count direction against the directed graph and raises `ValueError` if the edge doesn't exist.
- **Snap direction bug (fixed 2026-06-15):** `ingest_counts.py` previously stored canonical
  `(min(u,v), max(u,v))` — fixed to store actual directed `(u, v)`. Only session `f56b2ce4`
  was materially affected (re-snapped from 22→159 to 159→22).
- Two temporal profiles (f_s_res, f_s_biz) are inferred per (day_type, hour) slot, each
  anchored by component-specific priors from `hourly_fractions.csv`. The aggregate coupling
  (gamma_coupling_scale / std_f²) per slot keeps their sum near f_agg. With 72 slots
  and 2 df each, N_eff = N − 2×N_slots = 559 − 144 = 415.
- **Dead-end street absorption (ghost edges, fixed 2026-06-18):** OSMnx `simplify_graph`
  treats bidirectional dead-end terminus nodes as degree-2 (in=1, out=1 in the directed
  graph) and removes them, causing the dead-end edge to vanish from the consolidated graph.
  Without correction, buildings on absorbed dead-end streets would snap to the nearest
  surviving consolidated edge — often the main road but not reliably so for longer stubs
  in dense areas. `build_demographics.py` now detects these absorbed termini by comparing
  raw and consolidated network nodes, reconstructs their UTM geometry from the raw network,
  and adds ~761 "ghost" edges to the STRtree. Buildings snapping to a ghost edge have all
  their demand attributed to the surviving junction consolidated node (the only network
  entry point for that street). No change to `build_paths.py`, `model.py`, or the paths
  cache. Running `build_demographics.py` now prints "Added N ghost dead-end edges to
  STRtree (absorbed termini)".
- **`tuned_params.json` structure:** contains `slot_fracs_res` and `slot_fracs_biz` (dicts keyed `"dt,h"`); does **not** contain a `slot_fracs` key (legacy combined average, removed). `restore_params.py` also strips `slot_fracs` if found in old history entries.

---

## External Zone Reference Values (`tuner_config.json`)

**These values must not be changed without explicit user approval.** After a full tuning run,
updating them is something to *consider and discuss*, not an automatic step — the refs anchor
L2 regularization and changing them shifts the penalty basin for all future runs.

Last updated: 2026-06-19 (run 8a0fe24b, chi²/N=1.3742, 559 obs).
- gravity_lambda raised 0.05→0.5, gravity_ref P updated 300→600s (TSNI average journey ≈10 min) to prevent P drifting to unrealistic sub-minute values.

| City | Nodes | ref_pop | ref_wp | Tunable dampings |
|------|-------|---------|--------|-----------------|
| Donaghadee | 47 | 386,202 | 11,093 | — |
| Comber | 65, 617, 618, 620 | 71,927 | 2,208 | 617 (×0.29), 618 (×0.23), 620 (×0.36) |
| LowerArds | 92 | 136,840 | 5,850 | — |
| Belfast | 97, 119 | 517,496 | 1,822,862 | 119 (×0.39) |
| Dundonald | 10000 | 7,370 | 6,820 | — |
| Bangor | 98, 731 | 150,387 | 6,862 | 98 (×0.44) |
| Holywood | 99 | 3,785 | 1,204 | — |
| Millisle | 748, 749 | 2,525 | 496 | 749 (×0.46) |

---

## Worktree Convention

Background/parallel Claude work is done in `.claude/worktrees/` and
cherry-picked to `main` after review. All work is on `main`.
