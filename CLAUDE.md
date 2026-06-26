# Newtownards Traffic Model — Project Overview

A gravity-model traffic assignment pipeline for Newtownards, calibrated against
walking count data and official AADT figures. The pipeline is fully reproducible:
running the scripts in order regenerates all outputs from raw data.

**Agent instruction:** Keep this file up to date. After any tuning run, count data
ingest, model change, or reference value update, edit the relevant sections before
committing. This file is the authoritative record of model state.

**Dependencies:** Python package requirements are pinned in `requirements.txt`
(`python3 -m pip install -r requirements.txt`). It also documents the non-pip
dependencies (Docker + local OSRM backend) and the gitignored reference-data
downloads. **Agents: keep `requirements.txt` current** — add a pinned entry
whenever a new third-party import is introduced, and remove ones no longer used.

---

## Pipeline (run in this order)

```
python3 simulation/build_wz_apportionment.py # WZ→SA workplace apportionment via POI-weighted geometric intersection → data/ireland_data/cache_sa_workplace.csv (one-off; needs Docker + WZ boundary shapefile; re-run only when WZ/SA boundaries or WZ SAPS change)
python3 simulation/build_census_zones.py     # classify NI+RoI census areas → data/census_zones.json (one-off; needs SDZ/DEA boundary files + cache_sa_workplace.csv)
python3 simulation/build_network.py          # build road network from local NI .osm.pbf via osmium (core polygon + 5km bbox; needs Docker)
python3 simulation/build_demographics.py     # node weights + boundary detection + external weights (no longer builds the map — see build_map.py)
python3 simulation/build_external_links.py   # OSRM queries → external↔boundary links + through-route allowlist (needs local OSRM)
python3 simulation/build_intra_times.py      # OSRM-sample intra-zonal times per external zone → data/external_intra_times.json (denominator self-term; needs local OSRM; independent of build_paths — no cache rebuild)
python3 simulation/reduce_deadends.py        # collapse residential dead-ends → newtownards_reduced.graphml + node_weights_reduced.json (consumed by build_paths/build_assignment/tune; see reduce_deadends.py row)
python3 simulation/build_paths.py            # probit stochastic paths incl. external nodes (N_PASSES=25, CV=0.25, N_WORKERS=1; build time depends on hardware)

python3 analysis/parse_official_hourly.py    # parse ODS hourly counts → data/official_hourly.json (one-off)
python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/derive_component_profiles.py              # regenerate hourly_fractions.csv school column (re-run when NTS files change)

python3 analysis/tune_assignment.py                        # tune gravity params (8 params, production-constrained; external zones fixed from census)
python3 analysis/tune_assignment.py --fast                 # looser tolerances + fewer alt-min iters (~2× faster, minimal precision loss)
python3 analysis/tune_assignment.py --note "description"   # optional human label in history
python3 analysis/tune_assignment.py --f-frozen             # pin temporal fractions at the NTS profile (f-steps skipped); residuals become purely spatial. Writes params/history (f_frozen marker)

python3 simulation/build_assignment.py       # apply tuned params, write flows
python3 simulation/build_map.py              # build interactive map HTML (run after build_assignment.py to refresh flow layers)

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
| `simulation/zones_config.py` | Single source of truth for study-area geography: `CENTRE` (lat/lon), `CORE_RADIUS` (3 km), `SDZ_ZONE_RADIUS` (10 km). `CENTRE` is defined here only. Changing radii requires re-running `build_census_zones.py`. |
| `simulation/build_census_zones.py` | Classifies all of NI into a three-level hierarchy (DEA → SDZ → DZ) centred on `CENTRE`. Outputs `data/census_zones.json` (core polygon, external node list with census codes, centroids, population, workplace_pop). Needs `sdz2021/SDZ2021.geojson` and `dea2021/DEA2021.geojson` from NISRA/OpenDataNI. |
| `simulation/build_network.py` | Builds road graph from the local NI `.osm.pbf` via osmctools Docker image + osmnx. Outputs `newtownards_network.graphml` (raw, 5 km bbox) + `newtownards_consolidated.graphml` (junction-consolidated, clipped to core polygon). Needs Docker + pbf on disk. |
| `simulation/build_external_links.py` | OSRM-derives all external zone connectivity: X→B and B→X directed links (first/last boundary node in route), boundary↔boundary exterior shortcuts, and through-route allowlist. Outputs `data/external_links.json`. |
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, detects boundary nodes from core polygon, adds external node weights, writes `node_weights.json` + `newtownards_demographics.geojson`. `--zones-only` patches only external entries. |
| `simulation/build_map.py` | Builds the interactive folium map (`newtownards_map.html`). Run after `build_demographics.py`, and again after `build_assignment.py` to refresh flow layers. |
| `simulation/demographics_config.py` | Shared constants imported by `build_demographics.py`, `build_map.py`, and `build_network.py`: file paths, OSM tag handling, map styling, `PBF_PATH`, `BOUNDARY_BBOX_MARGIN_M`. Defines `PROJECTED_CRS = "EPSG:2157"` (ITM) — the single source of truth for all projected spatial ops. **Do not hardcode `EPSG:32630` anywhere.** |
| `simulation/build_intra_times.py` | OSRM-samples intra-zonal travel times per external census zone for the production-constrained self-term. Outputs `data/external_intra_times.json`. Run after `build_census_zones.py`. Independent of `build_paths.py` — no paths rebuild needed. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths (road + external nodes) with probit stochastic loading (N_PASSES=25, CV=0.25, length-scaled noise `PROBIT_LL_SIGMA=120s`). Result cached in `newtownards_paths.npz`. Re-run if road network, external links, `HIGHWAY_COST_FACTOR`, or probit noise params change. |
| `simulation/model.py` | Shared constants and core functions: `COUNT_SITES`, `EXCLUDE_LINKS`, `constrained_od_flows()`, `scatter_od_to_links()`, `load_self_terms()`, `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `WEIGHTS_FILE`/`ROUTING_GRAPH` point at the dead-end-reduced artifacts. COUNT_SITES: site 507 links 538692601↔549139252; site 508 node 136173611; site 444 node 449111329. EXCLUDE_LINKS: `{(181844513, 181839481)}` plus the Westmount Park and Old Belfast Road links absorbed by `reduce_deadends.py` (regenerate from `deadend_broken_obs.json` if reduction params change). |
| `simulation/build_assignment.py` | Production-constrained gravity assignment. Three-component mode when K_sch > 0 and school params/weights present; falls back to two-component or legacy single-K for old param files. |
| `simulation/routing_config.py` | Single source of truth for `HIGHWAY_COST_FACTOR` and probit noise params (`PROBIT_CV`, `PROBIT_LL_SIGMA`). Imported by `build_paths.py` and `build_osrm_profile.py`. |
| `simulation/build_osrm_profile.py` | Generates `car_roaaads.lua` from `HIGHWAY_COST_FACTOR`. Re-run when factors change, then re-preprocess OSRM. Output: `/home/matthew/Documents/CodingFun/osrm/car_roaaads.lua`. |
| `simulation/reduce_deadends.py` | Collapses residential dead-end regions (single cut-vertex, no boundary/school nodes, ≤60s max travel time, <100 biz demand) to super-nodes, shrinking the routing graph for faster tuning. Run after `build_demographics.py`, before `build_paths.py`. Outputs `newtownards_reduced.graphml`, `node_weights_reduced.json`, `deadend_map.json`, `deadend_broken_obs.json`. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Gravity param regularization and temporal coupling priors. `ext_biz_scale` compensates for the OSM car park area contribution external census nodes don't receive (recompute with `build_demographics.py` ratio if core coverage changes). `ext_school_per_pop` is the pupils-per-person ratio from the core, applied uniformly to external nodes (recompute via `simulation/compute_ext_school_scale.py` if core school POIs change significantly). |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. Run once, or when the ODS file changes. |
| `analysis/ingest_counts.py` | Reads CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT. Idempotent. Loads manual overrides from `data/manual_link_overrides.json`. Validates each non-null count direction against the directed graph. |
| `analysis/manual_assign_link.py` | CLI: `python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>`. Bypasses GPS snap; writes to `data/manual_link_overrides.json`. After correcting, re-run `aggregate_counts.py` then `tune_assignment.py`. |
| `analysis/aggregate_counts.py` | Combines per-session AADT into per-link estimates (inverse-variance weighting). Always regenerates from scratch. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method tuning of 8 gravity params (W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, P_school, ALPHA_school). External zone values fixed from census. Five-block alternating calibration (K, phi_biz, phi_sch, f_res/f_biz/f_school) at each eval. `CALIBRATE_PROBE=1` env var runs a diagnostic reporting residual global scale λ. |
| `analysis/report_tune.py` | Generates `reports/tune_report_{id}.txt` and `reports/slot_pulls_{id}.png` from a history entry. |
| `simulation/restore_params.py` | Restores `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial prefix matching supported. |
| `simulation/reset_gravity_params.py` | Resets only gravity params in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile (res/biz/school columns). Derived from NTS via `derive_component_profiles.py`. Constraint: res + biz + school = agg; rows sum to 7.0. Re-run `derive_component_profiles.py` when NTS files change. |
| `analysis/derive_component_profiles.py` | Derives component hourly fractions from DfT NTS data. Re-run when NTS files or purpose classification change. |
| `analysis/google_routing_common.py` | Shared stdlib helpers for Google routing calibration tooling. See **Google Routing-Time Calibration** section. |
| `analysis/google_feasibility.py` | One-shot feasibility pilot (superseded by manifest+runner for production; kept as reference). |
| `analysis/build_od_manifest.py` | Writes the fixed, deterministic OD sample → `data/google_cache/od_manifest.json`. No API/OSRM calls. |
| `analysis/google_query_routes.py` | Crash-safe, resumable Google Routes runner. Phase A (spendy): queries each uncached OD, writes atomically to `data/google_cache/raw/`. Phase B (free): rebuilds `results.jsonl` from cache via OSRM `/match`. Requires explicit per-run user approval for Phase A. |
| `simulation/profile_spec.py` | Single source of truth for a calibrated OSRM time profile: per-`(highway_class × speed_band)` multiplicative factors + global turn params. `factor=1.0` = stock OSRM speed; `factor>1 ⇒ slower`. Holds bucketisation logic + `ProfileSpec` (JSON load/save). Maxspeed resolution mirrors OSRM's `WayHandlers.maxspeed` exactly (key precedence, symbolic/nsl values). |
| `simulation/osrm_lua.py` | OSRM car.lua/Docker plumbing: `pull_base_lua`, `copy_lib`, injection helpers, Lua emitters (`emit_factor_block(spec)`, `apply_turn_overrides`). |
| `simulation/build_edge_index.py` | Builds the raw OSM edge cache. `--match`: single `/match` pass → `match_cache.jsonl` (resumable). `--extract`: osmctools+iterparse → `osm_ways.jsonl`/`osm_nodes.jsonl`. |
| `simulation/build_skeleton_index.py` | Rebuilds `skeletons.jsonl` (profile-independent route skeletons) from `match_cache.jsonl` + edge index. No OSRM calls; free to re-run. `--base-speeds` derives empirical per-bucket base speeds → `base_speeds.json`. |
| `analysis/skeleton_model.py` | Fast offline model: `predict_duration(skel, spec)`, `evaluate(skeletons, spec)`. Milliseconds for the whole cache. |
| `analysis/eval_profile.py` | Scores a `ProfileSpec` against `skeletons.jsonl`: loss, ratio distribution, per-leg/per-cell breakdown, per-bucket coverage. No spend. |
| `simulation/compile_profile.py` | `tuned_profile.json` → deployable `car_roaaads.lua`. Prints re-extract/partition/customize commands. |
| `analysis/verify_profile.py` | Fidelity gate after deploying a compiled profile: per-leg median `predict/real` within ±0.03 (default). Exits non-zero on fail. |
| `analysis/tune_profile.py` | Fits per-bucket speed factors to minimise squared-log-ratio error vs Google (external-focused; INT down-weighted). Linear-in-factor, scipy L-BFGS-B. Writes `simulation/tuned_profile.json` + `profile_tuning_history.jsonl`. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_map.html`, `simulation/tuned_params.json` — all regenerated by the pipeline.

**Node ID scheme:** road nodes use stable OSM IDs (integers; a consolidated junction gets `min(osmid_original)`). External census nodes use census-area-code string IDs (e.g. `"N21000219"`). `node_to_idx` mixes both. `node_weights.json` keys: `node_population`, `node_business_demand`, `node_school_demand`, `node_parking_equiv`, `boundary_node_ids`.

`simulation/newtownards_flows.json` — combined flows; `flows_res`/`flows_biz`/`flows_school` present when three-component params active.

`simulation/newtownards_reduced.graphml`, `simulation/node_weights_reduced.json` — dead-end-reduced graph + weights; consumed by `build_paths.py`, `build_assignment.py`, `tune_assignment.py`. `deadend_map.json` (provenance) and `deadend_broken_obs.json` (absorbed observations — review before adoption) also written.

`reports/` — generated by `report_tune.py` and `tune_assignment.py`; not tracked.

`data/google_cache/` — **gitignored** (Google ToS). Holds `od_manifest.json`, `raw/<od_id>.json`, `results.jsonl`, `match_cache.jsonl`, `osm_ways.jsonl`, `osm_nodes.jsonl`, `edge_index_meta.json`, `skeletons.jsonl`, `base_speeds.json`, `profile_tuning_history.jsonl`. Survives worktree removal; at risk from `git clean -xfd`.

`simulation/tuned_profile.json` (candidate `ProfileSpec`) — gitignored.

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`, `data/official_hourly.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`data/census_zones.json` — output of `build_census_zones.py`. Re-run only if NISRA boundary files or census data change.
`data/external_links.json` — output of `build_external_links.py`. Re-run when boundary nodes change or OSRM data is updated.
`data/external_intra_times.json` — output of `build_intra_times.py`. Re-run only when external zones change.
`data/manual_link_overrides.json` — committed so manual assignments survive a wipe of `counts_processed.json`.
`simulation/tuner_config.json` — committed as source config (gitignore exception).
`analysis/hourly_fractions.csv` — committed as source data.

### Large reference data (gitignored, kept locally only)
- `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods` — used by `parse_official_hourly.py`.
- `data/nts0502.ods` — DfT NTS Table NTS0502a. Used by `derive_component_profiles.py`.
- `data/nts0504.ods` — DfT NTS Table NTS0504b. Used by `derive_component_profiles.py`.
- `data/census-2021-apwp001.xlsx` — DZ-level workplace population. Used by `build_demographics.py` and `build_census_zones.py`.

Boundary files needed by `build_census_zones.py` (download from NISRA / OpenDataNI):
- `simulation/dz2021/DZ2021.geojson` — DZ polygon boundaries (present, gitignored).
- `simulation/sdz2021/SDZ2021.geojson` — SDZ polygon boundaries (present, gitignored).
- `simulation/dea2021/DEA2021.geojson` — DEA polygon boundaries (present, gitignored).

RoI data files for `build_wz_apportionment.py` + `build_census_zones.py` (in `data/ireland_data/`):
- `Small_Area_National_Statistical_Boundaries_2022_Ungeneralised_view_*.geojson` — 2022 SA boundaries (~410 MB).
- `Complete_set_of_Census_2022_SAPs/SAPS_2022_Small_Area_UR_171024.csv` — SA population (`T1_1AGETT`).
- `Workplace_Zones_ITM/Workplace_Zones_ITM.shp` — 2016 WZ boundaries in EPSG:2157 with headcount (`T1_T`); used only by `build_wz_apportionment.py`.
- `cache_sa_workplace.csv` — **generated** by `build_wz_apportionment.py`; committed once computed.

---

## Model Design

### Gravity model (production-constrained)
Each component is **singly (production) constrained**:
`T^c_ij = K_c · p^c_i · a^c_j · f_c(d_ij) / D^c_i`,  `D^c_i = Σ_k a^c_k·f_c(d_ik)`.
Each origin's trip production is fixed by its producing weight and independent of accessibility.

**Three components:**
- **Residential** (`flow_res`): `pop_i·pop_j·f_res(d)/D^res_i`
- **Business** (`flow_biz`): symmetric pop↔biz + `W_BIZ`·biz×biz, each term constrained on its origin type.
- **School** (`flow_school`): symmetric pop↔school, magnitude `K_sch`.

**Kernel:** `u = d/P; f(u) = (ALPHA+BETA) / (ALPHA·u^(-BETA) + BETA·u^ALPHA)` (overflow-safe form).
f(P)=1, f(0)=0, tail ~1/d^ALPHA, rise ~u^BETA near origin. BETA=1 recovers the original single-parameter kernel.

**Known-broken school component:** external nodes have `school_demand=0`, so all external school budgets dump into the core. The self-term machinery is wired and will activate automatically once external school demand is populated — no code change needed.

**Distances:** least-time shortest paths (seconds). External→internal paths traverse an OSRM-derived X→B edge then the internal network; Dijkstra selects the optimal boundary entry.

### External intra-zonal self-term
Each external origin's denominator `D^c_i` is augmented with `a^c_i·E[f_c(t)]` over OSRM-sampled intra-zonal times from `data/external_intra_times.json` — denominator-only, no link flow. Prevents the fixed external budget from over-allocating to the core (worst for large, isolated DEAs). Absent file ⇒ no self-term (prior behaviour). Re-tune after enabling to realize the fit gain.

### Stochastic route choice (probit loading)
N_PASSES=25 Dijkstra runs, each with log-normal edge-cost noise (CV=0.25, length-scaled: gain `w = σ_ll/(σ_ll + CV·cost)`, `PROBIT_LL_SIGMA=120s`). The length-scaling caps long external legs' perturbation at an absolute sigma rather than a multiplicative one, so boundary entry is driven by real time differences. `link_weight[entry]` = fraction of passes using that link. `od_dist` = mean path distance across passes.

### External zones (big-world network)
Three-level hierarchy centred on Newtownards:
- **Core area** (DZ level): union of DZs whose parent SDZ intersects CORE_RADIUS (3 km). Irregular polygon boundary.
- **SDZ nodes**: SDZs within SDZ_ZONE_RADIUS (10 km) not in the core.
- **DEA nodes**: DEAs entirely outside SDZ_ZONE_RADIUS.

External node IDs are census-area-code strings (e.g. `"N21000219"`). Demand from Census 2021; not tuned by the optimizer. Business demand scaled by `ext_biz_scale`; school demand by `population × ext_school_per_pop` (both from `tuner_config.json`).

### Through routes
External→external pairs allowed only if their OSRM route transits a boundary node (auto-generated allowlist in `data/external_links.json`). Changing the allowlist requires re-running `build_external_links.py` and `build_paths.py`.

### Five-block analytical calibration
At each optimizer evaluation, (K, phi_biz, phi_sch, f_res, f_biz, f_school) are calibrated via alternating minimisation (max_iter 40; `--fast` 20). K_res = K·(1−phi_b−phi_s), K_biz = K·phi_b, K_sch = K·phi_s.

The alternation is non-monotonic (single Poisson Newton steps). `calibrate_Ks_and_fracs` returns the **best-seen iterate** by the full regularized objective, not the last. `CALIBRATE_PROBE=1` reports the residual global scale λ (≈1 ⇒ K is at its optimum).

**f-steps:** per-slot analytical update anchored by NTS priors from `hourly_fractions.csv`. Aggregate coupling γ = `gamma_coupling_scale`/std_f² keeps res+biz+school sum near f_agg. **`--f-frozen`** skips the f-steps entirely, making residuals purely spatial (diagnostic mode).

**phi priors:** `phi_biz ~ N(phi_biz_prior, phi_biz_std²)`, `phi_sch ~ N(phi_school_prior, phi_school_std²)` from `tuner_config.json`. Prevent K_biz/K_sch → 0 degeneracy.

### Observations and goodness of fit
- **Official hourly** (216 obs, 24h × 3 day-types × 3 sites): Gaussian error, sigma from between-weekday std with 10%/15% relative floors.
- **Walking** (from `data/link_aadt.json`): Poisson error, weight = 1/n_eff.

`χ²/N`: mean squared z-score. N_eff = N_obs − 3·N_slots (three df per slot for three-component model). The tuner's χ²/N includes f-prior and coupling penalty terms; `build_assignment.py`'s is data-only (lower, but directionally comparable).

**Reading "modelled flow" across reports:**
- `build_assignment.py` / `newtownards_flows.json`: directed daily AADT (`K_res·flow_res + K_biz·flow_biz + K_sch·flow_school`).
- Map combined layer: two-way sum (`flow(u,v)+flow(v,u)`) — roughly 2× a single directed link.
- Tuner fit table: per-observation count-space (vehicles/hour for official, combined AADT for walking).

### `tuned_params.json` structure
Contains: `K_res`, `K_biz`, `K_sch`, `P_school`, `ALPHA_school`, `slot_fracs_res`, `slot_fracs_biz`, `slot_fracs_school` (dicts keyed `"dt,h"`). Does **not** contain `W_SCHOOL` (removed), `slot_fracs` (legacy), or external zone keys. Old param files without school keys fall back to two-component or legacy mode.

---

## Count Data

**Official hourly obs:**
- Site 507: A21 Bangor Road — 72 obs (24h × 3 day-types)
- Site 508: A48 Donaghadee Road — 72 obs
- Site 444: A20 Portaferry Road — 72 obs

Annual AADT values in `model.py` `COUNT_SITES` are for `build_assignment.py` display only; not used by the tuner.

**Walking counts:** manual link overrides (OSM node IDs): sessions e644eae2, 760b0c8e, bb934ba7 → 86604223↔86604221 (A20 Kempe Stones); e66989f4 → 150995265↔6622295361; b0043fd1, 32d425d6 → 181844516↔4688250384; 29d5f5f0 → 181844516↔538692566.

Current observation totals and N_eff are printed by `tune_assignment.py` at run start — treat the live output as authoritative over any number here.

---

## Tuning History

See `simulation/tuning_history.jsonl` (committed). Use `simulation/restore_params.py --list` to browse runs and restore any entry to `tuned_params.json`.

**Outstanding TBDs:**
- **phi_sch** unidentifiable without school-peak count sessions.
- External school demand: `population × ext_school_per_pop` is a coarse approximation; census school-age population per zone would be better.
- External business demand: `workplace_pop × ext_biz_scale` uses a different unit basis from internal OSM POI demand. May need separate `W_BIZ_ext` scaling.

---

## Paths Cache

The cache (`newtownards_paths.npz`) must be rebuilt with `build_paths.py` whenever:
- Road network changes (`newtownards_consolidated.graphml`)
- External links change (`data/external_links.json`)
- `HIGHWAY_COST_FACTOR`, `PROBIT_CV`, or `PROBIT_LL_SIGMA` change

**Staleness guard:** `build_paths.py` stamps a SHA-1 signature of its inputs. `tune_assignment.py` and `build_assignment.py` call `model.assert_paths_cache_fresh()` on load and raise `SystemExit` if stale, naming what changed.

**Cache format:** `node_ids` = road nodes (OSM ints) + external nodes (census strings); `link_u`/`link_v` = road links only (external edges excluded from flow accumulation); `link_weight` = float32 probit fractions; `od_dist` = mean path distance. `n_routed_pairs` marks flow-carrying pairs (0..n_routed_pairs−1); remainder are denominator-only non-through ext→ext pairs.

---

## External Zone Configuration

External zone values are fully data-driven (Census 2021 via `data/census_zones.json`; OSRM via `data/external_links.json`). No hand-crafted reference values, dampings, or city groupings.

**Gravity param refs** in `tuner_config.json` (`gravity_ref`): anchored for L2 regularization. Do not change without explicit approval.

**To update external zone coverage** (e.g. NISRA boundary update):
1. Re-run `build_census_zones.py`
2. Re-run `build_demographics.py`
3. Re-run `build_external_links.py`
4. Re-run `build_paths.py`
5. Re-tune, then `build_assignment.py` and `build_map.py`

---

## Google Routing-Time Calibration (offline, optional — NOT part of the main pipeline)

**Purpose:** OSRM (current profile) is too fast, especially on external approach corridors. This workflow calibrates a per-`(road-class × speed-band)` OSRM time profile against Google Maps as a source-of-truth for journey times, decoupling impedance (realistic travel time) from route preference (generalised cost).

**⚠️ Paid external API.** Never run a live Google query without explicit per-run user approval. `--dry-run`, cached re-runs, and all OSRM calls are free.

**Workflow:**
```
# 1. Build fixed OD sample (no API/OSRM calls):
python3 analysis/build_od_manifest.py            # → data/google_cache/od_manifest.json

# 2. With explicit approval:
python3 analysis/google_query_routes.py --dry-run            # counts + cost, no spend
GOOGLE_MAPS_API_KEY=... python3 analysis/google_query_routes.py --limit 100
python3 analysis/google_query_routes.py --reprocess-only     # rebuild results.jsonl, no spend

# 3. One-time skeleton cache (slow: one cached /match pass at ~1.7s/route):
python3 simulation/build_edge_index.py --match
python3 simulation/build_edge_index.py --extract
python3 simulation/build_skeleton_index.py
python3 simulation/build_skeleton_index.py --base-speeds     # empirical base speeds

# 4. Fast offline benchmark (no OSRM/spend):
python3 analysis/eval_profile.py                             # stock profile
python3 analysis/eval_profile.py --spec simulation/tuned_profile.json

# 5. Tune factors:
python3 analysis/tune_profile.py                             # → tuned_profile.json

# 6. Deploy + gate (once per accepted profile):
python3 simulation/compile_profile.py --spec simulation/tuned_profile.json
#   ... rebuild :5000 OSRM with the printed commands ...
python3 analysis/verify_profile.py --spec simulation/tuned_profile.json
#   then: build_external_links → reduce_deadends → build_paths → tune_assignment
```

`verify_profile.py` gates per-leg median `predict/real` within ±0.03 over external legs (X2B/B2X/X2X). INT is reported but not gated (offline turn model under-counts in-town junctions). The verify gate must pass before any `car_roaaads.lua` is adopted.

---

## Worktree Convention

Background/parallel Claude work is done in `.claude/worktrees/` and
cherry-picked to `main` after review. All work is on `main`.

---

## Agent Behaviour

Even in auto mode, when a user reports a bug or asks a question, agents must not make functional code changes without explicit user approval. Investigate and propose — do not implement unless the user has agreed to the specific change. "Defensive" workarounds (e.g. silent clamps, fallback defaults) are especially suspect: they tend to mask future bugs rather than surface them.

**Paid external APIs (Google Maps Routes) — never run without explicit approval.** Every Google query costs real money. Before any live run, state the planned query count + estimated cost, then stop and wait for an explicit go. Prefer having the user supply `GOOGLE_MAPS_API_KEY` and/or run the command themselves so the key never enters agent tool calls.
