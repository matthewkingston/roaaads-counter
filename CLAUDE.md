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
python3 simulation/build_census_zones.py     # classify NI census areas → data/census_zones.json (one-off; needs SDZ/DEA boundary files)
python3 simulation/build_network.py          # build road network from OSM (RADIUS_M=5000m)
python3 simulation/build_demographics.py     # node weights + boundary detection + external weights + map scaffold
python3 simulation/build_external_links.py   # OSRM queries → external↔boundary links + through-route allowlist (needs local OSRM)
python3 simulation/build_paths.py            # probit stochastic paths incl. external nodes (N_PASSES=25, CV=0.25, N_WORKERS=1; build time depends on hardware)

python3 analysis/parse_official_hourly.py    # parse ODS hourly counts → data/official_hourly.json (one-off)
python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/derive_component_profiles.py              # regenerate hourly_fractions.csv school column (re-run when NTS files change)

python3 analysis/tune_assignment.py                        # tune gravity params (9 params; external zones fixed from census)
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
| `simulation/build_census_zones.py` | **NEW.** Classifies all of NI into a three-level hierarchy (DEA → SDZ → DZ) centred on CENTRE. DEAs intersecting `SDZ_ZONE_RADIUS` (10 km) are broken into constituent SDZs; SDZs intersecting `CORE_RADIUS` (3 km) contribute their DZs to the core area. Population-weighted centroids are computed for each external (DEA/SDZ) node from DZ-level census population. Outputs `data/census_zones.json` containing: core area polygon (union of core DZs, in WGS84), list of external nodes with integer IDs (1-based sequential), census area codes, centroid lat/lon/UTM, population, and workplace_pop. Requires `simulation/sdz2021/SDZ2021.geojson` and `simulation/dea2021/DEA2021.geojson` (download from NISRA/OpenDataNI). School and business demand for external nodes are **TBD** — currently `business_demand = workplace_pop` and `school_demand = 0`. |
| `simulation/build_external_links.py` | **NEW.** Queries a local OSRM instance (NI extract, **biased `car_roaaads.lua` profile** — see `build_osrm_profile.py`, `http://localhost:5000`) to derive all external zone connectivity. **X→B links:** for each (external node, boundary node) ordered pair, keeps the link only if that boundary node is the first boundary node encountered in the OSRM route (i.e., it is the natural entry point into the core). **B→X links:** symmetric with X→B — keeps B→X only if no other boundary node appears in the OSRM route sequence (i.e., B is the last boundary node departed on the way to X). If another boundary node B' appears, the journey is already covered by B→B' + B'→X. **Boundary→boundary exterior shortcuts:** for each ordered boundary pair, if the route exits the core first, adds a directed shortcut with duration summed from OSRM annotations up to the first boundary node re-encountered. **Through-route allowlist:** for each ordered external-external pair, checks if any OSRM route node is a boundary node; if so adds to `allowed_through_pairs`. Outputs `data/external_links.json`. ~28,000 OSRM queries; under a minute on a local instance. |
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, detects boundary nodes, adds external node weights, builds map. `--map-only` skips demographic recomputation. `--zones-only` re-reads `data/census_zones.json` and patches only the external node entries in `node_weights.json`. **Boundary node detection:** loads core polygon from `census_zones.json`, identifies internal nodes (within core polygon), then boundary nodes = internal nodes with at least one edge going outside. Writes `boundary_node_ids` to `node_weights.json` (replaces the previous hand-specified list). **External node weights:** reads external node list from `census_zones.json` and writes population/workplace_pop to `node_weights.json`. **Study area = core polygon (not a circle):** DZ selection and all OSM downloads (buildings/POIs/parking) are bounded by the core polygon from `census_zones.json` (extent `max_core_vertex_dist_m`, ~10.2 km), matching the road graph built by `build_network.py`. Core DZs are selected by centroid-within the polygon (recovers exactly the `n_core_dzs` core DZs) and use **full** DZ population/workplace_pop (no area-fraction clipping — the legacy 3 km `RADIUS_M` circle is removed). OSM downloads use a circle sized to the polygon (+1 km margin); POIs and parking are then filtered to within the core polygon to avoid margin leakage (buildings are already DZ-bounded by sjoin). **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Business demand:** workplace population × POI count (school/college/university excluded), augmented with OSM car park polygon area (public: area/25, private: area/50 equivalent persons). **School demand:** separate `node_school_demand` layer from OSM school POIs. Enrollment fallbacks: school→300, secondary_school→900, college→2000, university→3000 pupils. **Flow map layers:** combined AADT (default), residential (teal), business (amber→red), school (violet→purple). |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Now covers both internal road nodes and external census-area nodes. **Graph augmentation:** loads external nodes and edges from `data/external_links.json`; adds them to the routing graph before Dijkstra. External edges (X↔B, boundary shortcuts) are included in the adjacency matrix but NOT in `link_list` — they contribute to path distance but not to flow accumulation. **Probit loading:** all edges (road and external) perturbed with the same log-normal noise (CV=0.25, N_PASSES=25), giving stochastic spread in boundary node selection for external→internal OD pairs with similarly-weighted entry options. **OD pair filter:** external→external pairs only allowed if in `allowed_through_pairs` (from `external_links.json`). No offscreen leg calculation. `HIGHWAY_COST_FACTOR` is imported from `simulation/routing_config.py`. Re-run if road network, external links, `HIGHWAY_COST_FACTOR`, `N_PASSES`, or `PROBIT_CV` change. **Performance constants:** `N_WORKERS` (default 1) controls parallel pass workers via `multiprocessing.Pool` — increase on machines with sufficient RAM (each worker uses ~100–150 MB extra); `MAX_HOPS` (default 120) caps per-pair path-trace iterations. Inner path-tracing loop is vectorised (numpy). |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (rational kernel; `return_components=True` with `w_school` provided returns `(flow_res, flow_biz, flow_school)` tuple), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `gravity_assign()` accepts optional `link_weight` array (probit fractional weights). `compute_chi2()` dispatches to three-component mode when `link_flow_school_dict` provided (N_eff = N − 3·N_slots); two-component when only `link_flow_biz_dict` provided; legacy single-flow otherwise. All node IDs are OSM IDs (stable). COUNT_SITES: site 507 links 538692601↔549139252; site 508 node 136173611; site 444 node 449111329. EXCLUDE_LINKS: `{(181844513, 181839481)}`. |
| `simulation/build_assignment.py` | Gravity model assignment. Requires `simulation/newtownards_paths.npz`. Three-component mode activated when K_sch > 0 and W_SCHOOL are in `tuned_params.json` and `node_school_demand` is in `node_weights.json`. Saves `flows_res`, `flows_biz`, `flows_school` in `newtownards_flows.json`. Falls back to two-component (K_res/K_biz only) or legacy single-K for old param files. External node weights come from `node_weights.json` directly (no override from tuned params). |
| `simulation/routing_config.py` | **Single source of truth for `HIGHWAY_COST_FACTOR`.** Imported by `build_paths.py` (internal Dijkstra) and `build_osrm_profile.py` (OSRM profile generator). Edit factors here; both systems pick up the change. |
| `simulation/build_osrm_profile.py` | Generates `car_roaaads.lua` — the road-class-biased OSRM car profile. Pulls the default `car.lua` from the `osrm/osrm-backend` Docker image, injects a block after the `forward_rate` assignment that divides `forward_speed`/`forward_rate` by `HIGHWAY_COST_FACTOR` (matching internal Dijkstra biasing). Re-run whenever `HIGHWAY_COST_FACTOR` changes, then re-preprocess OSRM (`osrm-extract -p car_roaaads.lua`, `osrm-partition`, `osrm-customize`). Output: `/home/matthew/Documents/CodingFun/osrm/car_roaaads.lua`. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Gravity param regularization and temporal coupling priors. `gravity_lambda` + `gravity_ref` regularise P/ALPHA/BETA/W_BIZ/P_biz/ALPHA_biz; `gamma_coupling_scale` controls the per-slot aggregate coupling (γ = scale/std_f²); `phi_biz_prior`/`phi_biz_std` and `phi_school_prior`/`phi_school_std` set Gaussian priors on business and school flow fractions. `ext_biz_scale` (default 1.0; current value 1.702) scales external node business demand in `build_demographics.py` to compensate for the OSM car park area contribution that internal nodes receive but external census nodes don't (ratio of total business demand to total NISRA workplace pop across all internal nodes, including boundary). Recomputed 2026-06-22 after the core-polygon demographics fix widened internal coverage: total internal biz 22,517 / workplace 13,230 = 1.702. **Removed:** `cities` block (replaced by `census_zones.json`) and `through_route_pairs` whitelist (replaced by OSRM-derived `external_links.json`). `lambda` is retained but no longer used (external zone params are not tuned). |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. Run once (or when the ODS file changes). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. Loads manual link overrides from `data/manual_link_overrides.json`. After every new link assignment, validates each non-null count direction against the directed graph; raises `ValueError` if the edge doesn't exist. |
| `analysis/manual_assign_link.py` | CLI tool to manually assign a session to a specific directed link, bypassing GPS snap. Usage: `python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>`. Validates both nodes exist and checks count-edge consistency. Writes to `data/manual_link_overrides.json` and patches `counts_processed.json` directly. After correcting an assignment, re-run `aggregate_counts.py` then `tune_assignment.py`. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s`. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Three-component model:** gravity flows split into residential (`flow_res`), business-adjacent (`flow_biz`), and school (`flow_school`) components. Tunes 9 gravity params (W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, W_SCHOOL, P_school, ALPHA_school). External zone values are fixed from census data and are not tuned. **Alternating minimisation (5 blocks, up to 10 iters):** K-step (1D); phi_biz-step; phi_sch-step; f_res/f_biz/f_school steps (per-slot analytical) + aggregate coupling γ per slot. All OD pairs (including external-involved) use a single unified bin-matrix path; no internal/external split. **Performance estimate:** ~40–50 s per run. `--fast` mode caps alt-min at 5 iters. |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` and `reports/slot_pulls_{id}.png`. History `slot_prior` entries carry 5 values: `[mean_f_agg, std_f, mean_f_res, mean_f_biz, mean_f_school]`. Old entries with 4 values handled gracefully. Note: report still attempts to print external city delta table — will be a no-op for new-format runs. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. New-format param dicts no longer contain external zone keys. |
| `simulation/reset_gravity_params.py` | Reset only the gravity params (K, W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz) in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile. Includes `mean_fraction_res`, `mean_fraction_biz`, and `mean_fraction_school` columns. **Derived from NTS** via `analysis/derive_component_profiles.py`. Constraint: res + biz + school = agg for all 168 rows. Rows sum to 7.0 (AADT normalisation). Re-run `derive_component_profiles.py` whenever NTS files change. |
| `analysis/derive_component_profiles.py` | Derives `mean_fraction_res`, `mean_fraction_biz`, and `mean_fraction_school` from DfT NTS data (2023–2024 rolling average). Purpose classification: biz = commuting + employer's business + shopping; school = education×⅕ + escort education; res = remainder. Re-run whenever NTS files or purpose classification change. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_map.html`, `simulation/tuned_params.json` — all regenerated by the pipeline.
`simulation/node_weights.json` keys: `node_population`, `node_business_demand`, `node_school_demand`, `node_parking_equiv`, `boundary_node_ids` (auto-detected from core polygon). External node entries (integer IDs 1, 2, …) are included alongside internal OSM node IDs. `node_effective_utm` is removed.
**Node ID scheme (as of 2026-06-20):** `build_network.py` relabels consolidated graph nodes to stable OSM IDs after junction consolidation. A consolidated junction (multiple OSM nodes merged) gets `min(osmid_original)` as its ID; a non-merged node gets `int(osmid_original)`. All road node IDs are therefore genuine OSM node IDs (in the hundreds of millions) and stable across graph regenerations. External census nodes retain their small integer IDs (1, 2, …). Map tooltips and all downstream files use these OSM IDs.
`simulation/newtownards_flows.json` — combined flows plus optional `flows_res`/`flows_biz`/`flows_school` keys when three-component params active.
`reports/` — generated by `report_tune.py` and `tune_assignment.py`; not tracked.

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`, `data/official_hourly.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`data/census_zones.json` — committed; output of `build_census_zones.py`. Contains core polygon, external node list with IDs/codes/centroids/census demand. Re-run `build_census_zones.py` only if NISRA boundary files or census data change.
`data/external_links.json` — committed; output of `build_external_links.py`. Contains OSRM-derived X↔B links, boundary shortcuts, and through-route allowlist. Re-run `build_external_links.py` when boundary nodes change or OSRM data is updated.
`data/manual_link_overrides.json` — committed so manual assignments survive a wipe of `counts_processed.json`.
`simulation/tuner_config.json` — committed as source config (gitignore exception).
`analysis/hourly_fractions.csv` — committed as source data (single authoritative version).

### Large reference data (gitignored, kept locally only)
`data/*.ods`, `data/*.xlsx`, boundary GeoJSON files — too large to commit; keep local copies.
Currently present:
- `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods` — used by `parse_official_hourly.py`.
- `data/nts0502.ods` — DfT NTS Table NTS0502a, weekday trip start times. Used by `derive_component_profiles.py`.
- `data/nts0504.ods` — DfT NTS Table NTS0504b, average trips by day/purpose. Used by `derive_component_profiles.py`.
- `data/census-2021-apwp001.xlsx` — DZ-level workplace population. Used by `build_demographics.py` (internal nodes) and `build_census_zones.py` (external nodes).

Boundary files needed by `build_census_zones.py` (download from NISRA / OpenDataNI):
- `simulation/dz2021/DZ2021.geojson` — DZ polygon boundaries (already present as .qmd stub).
- `simulation/sdz2021/SDZ2021.geojson` — SDZ polygon boundaries (**not yet downloaded**).
- `simulation/dea2021/DEA2021.geojson` — DEA polygon boundaries (**not yet downloaded**).

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

Distances are least-time shortest paths (seconds). For external→internal OD pairs, the path traverses an OSRM-derived external edge (X→B, fixed weight) then the internal road network (B→J). Dijkstra selects the optimal boundary entry node for each destination.

### Stochastic route choice (probit loading)
The paths cache stores fractional link-assignment weights computed from `N_PASSES=25`
Dijkstra runs, each with log-normal edge-cost noise (CV=0.25). For each OD pair,
`link_weight[entry]` is the fraction of passes that routed through that link. Pairs
with no topological route diversity (degree-1 stubs, single-access nodes) converge to
weight=1.0 on their forced route. `od_dist` is the mean path distance across passes.

This replaces the previous k=2/k=3 global-penalisation scheme, which was ineffective
in a dense network (global penalisation of all k=1-used links preserves relative ordering
and produces identical alternative paths for most OD pairs). THETA is no longer tuned.

### External zones (big-world network)
NI is represented as a three-level hierarchy centred on Newtownards (CENTRE):

- **Core area** (DZ level): union of all DZs whose parent SDZ intersects `CORE_RADIUS` (3 km). Boundary is irregular (follows census polygon edges, not a circle).
- **SDZ external nodes**: SDZs within `SDZ_ZONE_RADIUS` (10 km) that are not in the core — one centroid node per SDZ.
- **DEA external nodes**: DEAs entirely outside `SDZ_ZONE_RADIUS` — one centroid node per DEA.

Once a DEA is broken into SDZs, all its constituent SDZs become nodes (even those beyond `SDZ_ZONE_RADIUS`). External node integer IDs are sequential (1, 2, …), well below any OSM node ID. Node identifiers in census area code form (e.g. `N09000001`) are stored in `census_zones.json` for traceability.

**Demand:** population and workplace_pop from Census 2021 (DZ-level data aggregated to SDZ/DEA). Business demand for external nodes is currently set equal to workplace_pop (TBD refinement). School demand = 0 for external nodes (TBD).

**Connectivity:** boundary nodes are all internal nodes with at least one road edge crossing the core polygon boundary. OSRM-derived directed edges connect each external node to its valid boundary nodes; the gravity model path distance for any external→internal pair is the OSRM leg plus the internal shortest path, computed automatically by Dijkstra on the augmented graph.

**No tuning of external zone values.** Pop/wp come directly from census data and are not adjusted by the optimizer. This replaces the previous hand-crafted city configs, ref_pop/ref_wp values, and damping factors.

### Through routes
External→external OD pairs are allowed only for pairs whose OSRM route passes through at least one boundary node (i.e. genuinely transits the core area). This allowlist is auto-generated by `build_external_links.py` and stored in `data/external_links.json`. Changing the allowlist requires re-running `build_external_links.py` and `build_paths.py` (the paths cache must be rebuilt).

### Three-component flow decomposition
The gravity OD flows are split into three spatial components at each tuner evaluation:

- **Residential** (`flow_res`): `all_bin_pp @ f_b` — purely pop×pop trips
- **Business-adjacent** (`flow_biz`): `W_BIZ·(all_bin_pb @ f_b) + W_BIZ²·(all_bin_bb @ f_b)` — home↔work/retail trips
- **School** (`flow_school`): `W_SCHOOL·(all_bin_ps @ f_b)` — pop→school cross-term, using (P_school, ALPHA_school, BETA) kernel

Each component has its own temporal profile and scale (K_res, K_biz, K_sch).
Predicted count for observation i in slot s:
`pred_i = K_res·flow_res·(T/3600)·f_res[s] + K_biz·flow_biz·(T/3600)·f_biz[s] + K_sch·flow_school·(T/3600)·f_school[s]`

### Five-block analytical calibration
At each optimizer evaluation, (K, phi_biz, phi_sch, f_res, f_biz, f_school) are calibrated via
alternating minimisation (up to 10 iterations; `--fast` caps at 5).
K_res = K·(1−phi_biz−phi_sch), K_biz = K·phi_biz, K_sch = K·phi_sch.

**K-step:** 1D solve, using combined coefficient `(1−phi_b−phi_s)·c_r·f_r + phi_b·c_b·f_b + phi_s·c_s·f_s`.

**phi_biz-step / phi_sch-step:** sequential 1D solves with Gaussian priors.
phi_biz ~ N(phi_biz_prior, phi_biz_std²); phi_sch ~ N(phi_school_prior, phi_school_std²) from `tuner_config.json`.
These priors prevent degeneracy K_biz→0 or K_sch→0. Applied sequentially (fix one, solve for the other).

**f_res / f_biz / f_school steps:** per-slot analytical update, anchored by NTS-derived priors from
`hourly_fractions.csv` columns `mean_fraction_res` / `mean_fraction_biz` / `mean_fraction_school`.
School profile: sharp weekday double-peak (h08/h15), near-zero weekends.

**Aggregate coupling:** each slot carries γ·(f_res + f_biz + f_school − f_agg)² where
γ = `gamma_coupling_scale` / std_f_agg². Updated in all three f-steps.

Slot key: (day_type, hour), day_type = 0 (weekday), 1 (Saturday), 2 (Sunday).
Prior std from `hourly_fractions.csv` via law of total variance.

### Observations
All 559 observations are in count space with per-obs weights:
- **Official hourly** (216 obs, 24 h × 3 day-types × 3 sites): from `data/official_hourly.json`;
  Gaussian error (sigma from between-weekday std, 10% floor); weight = 1/sigma².
- **Walking** (343 obs): from `data/link_aadt.json`; Poisson error; weight = 1/n_eff.

### Goodness of fit
`χ²/N` (mean squared z-score; N=559 obs, N_eff = N − 3·N_slots = 559 − 216 = 343).
Three df lost per slot (one each for f_s_res, f_s_biz, f_s_school). With coupling enabled,
chi²/N includes coupling penalty terms; pure data-fit chi²/N is lower.

`build_assignment.py` uses the two-component `compute_chi2()` when K_res/K_biz are present in tuned_params.json. This gives a **data-only** chi²/N (pure sum of squared z-scores) — it excludes the f-prior penalties `(f_r−mfr)²/std_f²` and the aggregate coupling penalty that the tuner includes in its chi²/N. Expect the build_assignment chi²/N to be somewhat lower than the tuner's; the two are directionally comparable but not numerically equal. The legacy Woodbury path is used only for old single-K param files.

**Reading "modelled flow" across reports.** The three reporting surfaces print *different
projections* of the same tuned model — they are not directly comparable line-for-line:
- `build_assignment.py` "Official count sites" block and the `"flows"` values in
  `newtownards_flows.json` → **directed daily AADT** = `K_res·flow_res + K_biz·flow_biz + K_sch·flow_school`.
  This is the canonical modelled link flow. Node-based sites (508/444) sum every directed link at the node.
- `newtownards_map.html` combined layer → the same AADT but **summed over both directions** of each
  edge (`flow(u,v)+flow(v,u)`), i.e. a two-way total (~2× a single directed link).
- The tuner / `report_tune.py` fit table → **per-observation, count-space**: official rows are
  *vehicles/hour* in one (day_type, hour) slot (≈ AADT × hourly fraction), walking rows are
  reconstructed to combined AADT. Correct for goodness-of-fit; not a table of link AADTs.
- Walking "Model" column convention (display only, chi²/N unaffected): both `compute_chi2()`
  (`model.py`) and the tuner fit table show **combined directed AADT** for walking links. (Fixed
  2026-06-21: `model.py` previously divided `pred` by `f_eff` only, omitting the `Th` session-duration
  factor, so it showed AADT×Th; the tuner used a K-weighted reconstruction. Both now use `m_r+m_b+m_s`.)
- Fit-table street names come from the consolidated GraphML edge `name` attribute. `tune_assignment.py`
  resolves the GraphML data-key id dynamically from the `<key>` header (fixed 2026-06-21: it had
  hardcoded `d14`, which is unstable across network regenerations and had become `oneway`).
  `report_tune.py` echoes the labels the tuner stored in history, so names appear only after a fresh
  tune run regenerates `tuning_history.jsonl`.

---

## Count Data

**Official hourly obs (from ODS, parsed by `parse_official_hourly.py`):**
- Site 507: A21 Bangor Road — 72 obs (24 h × 3 day-types); Gaussian sigma, 10%/15% floors
- Site 508: A48 Donaghadee Road — 72 obs
- Site 444: A20 Portaferry Road — 72 obs

Annual AADT values (507: 21,202; 508: 10,792; 444: 7,282) are retained in `model.py`
`COUNT_SITES` for `build_assignment.py` backward compatibility but are no longer used
by the tuner directly.

**Walking counts:** 7 CSV files, 177 sessions, 343 per-session observations (after EXCLUDE_LINKS). Manual overrides re-entered with OSM node IDs: sessions e644eae2, 760b0c8e, bb934ba7 → link 86604223↔86604221 (A20 Kempe Stones); e66989f4 → 150995265↔6622295361; b0043fd1, 32d425d6 → 181844516↔4688250384; 29d5f5f0 → 181844516↔538692566.
New sessions added 2026-06-18 (7 sessions): Saratoga Avenue, Glenford Road, Hardford Link, Belfast Road.
**Total: 559 observations (216 official hourly + 343 walking) in 72 time slots. N_eff = 559 − 3×72 = 343.**

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
| 2026-06-19 | gravity | 559 | 9 | 1.3292 | three-component model (school added); phi_biz=27%, phi_sch=1.2%; school at ref |
| 2026-06-19 | full | 559 | 31 | **1.3146** | three-component full tune; phi_biz=25.6%, phi_sch=1.4%; school params at ref |

**Note on comparability:** runs from 2026-06-17 onward use the two-component model with coupling penalty terms in chi²/N; not directly comparable to earlier single-component runs. From 2026-06-19 three-component model: N_eff = 559 − 3×72 = 343 (one extra df per slot for f_school). **Runs from the big-world architecture are not directly comparable to earlier runs** — external zone representation has fundamentally changed (census-derived vs hand-crafted; many more external nodes; OSRM-based connectivity vs offscreen Euclidean leg).

Last pre-big-world best: chi²/N = **1.3146** (559 obs, N_eff=343; three-component with probit cache, run f09a003e).
K_res=1.23e-04, K_biz=4.31e-05, K_sch=2.40e-06. phi_biz=25.6%, phi_sch=1.4%.
W_BIZ=3.82, P=117.6s, ALPHA=4.02, BETA=7.67. P_biz=83.4s, ALPHA_biz=3.66.
W_SCHOOL=1.00, P_school=600s, ALPHA_school=2.00 (at ref).

**First big-world tune:** not yet run — requires SDZ/DEA boundary files and local OSRM instance.

**Outstanding concerns (carry-forward):**
- **phi_sch=1.4%** — school component unidentifiable without school-peak count sessions.
- Structural outliers (sequential node IDs — will need re-identifying in new OSM-ID graph): `22→12 Regent Street` (z=+4.03), `23→295 Frances Street` (z=+3.94), `296→297 Nursery Road` (z=−3.70), `139→137 Portaferry Road` (z=−3.70).
- `73→70` Mill Street severe underprediction (old IDs; z=−3.30; obs 26,377 vs model 2,682).
- `719→325` / `325→719` Messines Road persistent (old IDs; z=−3.33/−2.53).
- Hardford Link persistent (old IDs; z=−3.25/−3.19/−2.73).
- Business demand units mismatch: external nodes use census workplace_pop, internal nodes use OSM POI proxy. TBD whether `W_BIZ` needs separate scaling for external vs internal demand.
- **Paths cache stale** — must rebuild with `build_paths.py` (existing cache has old sequential node IDs) before re-tuning.

### Paths cache note
The paths cache (`newtownards_paths.npz`) must be rebuilt with `build_paths.py` whenever:
- The road network changes (`newtownards_consolidated.graphml`)
- External links change (`data/external_links.json` — re-run `build_external_links.py` first)
- `HIGHWAY_COST_FACTOR` values change
- `N_PASSES` or `PROBIT_CV` change

**Current cache format** (probit): `node_ids` covers road nodes (OSM IDs) + external nodes (small integer IDs 1, 2, …); `link_u`/`link_v` are road-link endpoints only (external edges are not in `link_list`); `link_weight` (float32, fraction of passes using that link for each OD pair); `od_dist` (mean path distance across passes including external legs); `probit_n_passes`, `probit_cv`. No `pair_idx_2/3` keys — `_has_stoch = False`, THETA not tuned. **The existing cache is stale** (built with old sequential node IDs) and must be rebuilt with `build_paths.py` after the OSM-ID migration.

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
- After a structural model change (e.g. new count data or external link regeneration), a fresh tune is needed to restore fit quality.
- **External node probit loading:** all edges (road and external) receive the same log-normal noise each pass (CV=0.25). Route diversity for external-internal OD pairs comes from both the X→B external leg and the internal B→J portion, giving stochastic spread across similarly-weighted boundary entry points. **Known issue (TBD after first big-world tune):** the noise is multiplicative, so long external legs (e.g. distant DEAs with 90-min legs) receive ±22-minute perturbations — large enough to swamp 5-minute differences between boundary nodes. Gravity suppression limits the practical impact (distant zones have similar total distances regardless of boundary choice), but a separate `PROBIT_CV_EXT` (e.g. 0.05) for external edges vs road edges (0.25) should be evaluated.
- **Dundonald virtual node (10000) is removed** in the big-world system. Dundonald is now represented by an SDZ or DEA external centroid node with a proper census-derived population.
- **Manual link overrides:** Use `analysis/manual_assign_link.py <session_id> <from_node> <to_node>` to assign a session to a specific directed link, bypassing GPS snap. Use when the observer stood on a parallel carriageway and the snap would land on the wrong physical road. The override is stored in `data/manual_link_overrides.json` and takes effect even if `counts_processed.json` is wiped and rebuilt. After assignment (manual or auto), `ingest_counts.py` validates each non-null count direction against the directed graph and raises `ValueError` if the edge doesn't exist.
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
- **`tuned_params.json` structure:** contains `K_res`, `K_biz`, `K_sch`, `W_SCHOOL`, `P_school`, `ALPHA_school`, `slot_fracs_res`, `slot_fracs_biz`, `slot_fracs_school` (dicts keyed `"dt,h"`); does **not** contain a `slot_fracs` key (legacy) or `external_node_pop/biz/city_pop/wp/dampings` keys (removed). Old param files without school keys fall back to two-component or legacy mode in `build_assignment.py`.

---

## External Zone Configuration

External zone values are now fully data-driven from Census 2021 (via `data/census_zones.json`) and OSRM routing (via `data/external_links.json`). There are no hand-crafted reference values, dampings, or city groupings to maintain.

**Gravity param refs** (`tuner_config.json` `gravity_ref`): still anchored for L2 regularization. Last updated 2026-06-19 — `gravity_ref P = 600s`, `gravity_lambda P = 2.0`. These must not be changed without explicit approval.

**To update external zone coverage** (e.g. after a NISRA boundary update):
1. Re-run `build_census_zones.py` (updates `data/census_zones.json`)
2. Re-run `build_demographics.py` (updates `node_weights.json`)
3. Re-run `build_external_links.py` (updates `data/external_links.json`)
4. Re-run `build_paths.py` (rebuilds paths cache with new external nodes)
5. Re-tune

**Outstanding TBDs for external nodes:**
- School demand: currently 0 for all external nodes. Should use census school-age population or similar.
- Business demand: currently set equal to workplace_pop. Units differ from internal nodes (which use OSM POI proxy). A separate `W_BIZ_ext` or normalisation may be needed.

---

## Worktree Convention

Background/parallel Claude work is done in `.claude/worktrees/` and
cherry-picked to `main` after review. All work is on `main`.

---

## Agent Behaviour

Even in auto mode, when a user reports a bug or asks a question, agents must not make functional code changes without explicit user approval. Investigate and propose — do not implement unless the user has agreed to the specific change. "Defensive" workarounds (e.g. silent clamps, fallback defaults) are especially suspect: they tend to mask future bugs rather than surface them, and should be proposed and justified before being applied.
