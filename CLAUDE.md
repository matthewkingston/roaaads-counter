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
python3 simulation/build_census_zones.py     # classify NI census areas → data/census_zones.json (one-off; needs SDZ/DEA boundary files)
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
| `simulation/zones_config.py` | **NEW.** Single source of truth for the study-area geographic knobs: `CENTRE` (lat/lon), `CORE_RADIUS` (3 km), `SDZ_ZONE_RADIUS` (10 km). Imported by `build_census_zones.py` (uses the radii), `build_network.py` (uses `CENTRE`), and `demographics_config.py` (re-exports `CENTRE`). `CENTRE` is defined here and nowhere else. Editing the radii takes effect only after re-running `build_census_zones.py`. |
| `simulation/build_census_zones.py` | **NEW.** Classifies all of NI into a three-level hierarchy (DEA → SDZ → DZ) centred on `CENTRE`. `CENTRE`/`CORE_RADIUS`/`SDZ_ZONE_RADIUS` are imported from `zones_config.py`. DEAs intersecting `SDZ_ZONE_RADIUS` (10 km) are broken into constituent SDZs; SDZs intersecting `CORE_RADIUS` (3 km) contribute their DZs to the core area. Population-weighted centroids are computed for each external (DEA/SDZ) node from DZ-level census population. Outputs `data/census_zones.json` containing: core area polygon (union of core DZs, in WGS84), list of external nodes whose `id` **is** the census-area code (SDZ/DEA 2021 code string, e.g. `"N21000219"` — no separate integer ID), with `level` (SDZ/DEA), centroid lat/lon/UTM, population, and workplace_pop. Requires `simulation/sdz2021/SDZ2021.geojson` and `simulation/dea2021/DEA2021.geojson` (download from NISRA/OpenDataNI). Business demand for external nodes is `workplace_pop × ext_biz_scale` (from `tuner_config.json`). School demand for external nodes is `population × ext_school_per_pop` (from `tuner_config.json`; see `simulation/compute_ext_school_scale.py`). |
| `simulation/build_network.py` | Builds the road graph from the **local NI `.osm.pbf`** (the same Geofabrik snapshot OSRM is built from — `demographics_config.PBF_PATH`), so road/boundary/internal node IDs share one OSM snapshot with OSRM's route node IDs. The full ~400 MB island pbf OOMs an in-process parse, so a small extract is streamed out with **osmctools** (`osmconvert` + `osmfilter`; Docker image `osmctools-roaaads`, auto-built from `simulation/osmctools.Dockerfile`, ~0.5 GB peak RAM): `osmconvert -b=<bbox> --complete-ways` (bbox = core polygon buffered by `demographics_config.BOUNDARY_BBOX_MARGIN_M` = 5 km — supersedes the old 1 km Overpass `dist` margin) then `osmfilter --keep="highway=<drive set>"` (positive form of osmnx's `drive` filter), written to `simulation/_pbf_drive_extract.osm`. (osmctools is used rather than osmium-tool, whose referenced-node id-set is sized by OSM's max node id and needs several GB regardless of extract area.) `ox.graph_from_xml` reads it — identical graph semantics to the old `graph_from_point("drive")` path, **except** `graph_from_xml` omits the `street_count` node attribute, which `consolidate_intersections` needs; it is re-added via `ox.stats.count_streets_per_node` (without it the core under-merges, ≈1416 vs ≈1004 nodes). Raw graph extends 5 km beyond the core (for boundary nodes' external neighbours + `build_external_links.py` positions); the consolidated routing graph is still clipped to the core polygon, then junction-consolidated (tol 15 m) and relabelled to OSM IDs. Outputs `newtownards_network.graphml` (raw) + `newtownards_consolidated.graphml`. **Needs Docker + the pbf on disk.** |
| `simulation/build_external_links.py` | **NEW.** Queries a local OSRM instance (NI extract, **biased `car_roaaads.lua` profile** — see `build_osrm_profile.py`, `http://localhost:5000`) to derive all external zone connectivity. **X→B links:** for each (external node, boundary node) ordered pair, keeps the link only if that boundary node is the first boundary node encountered in the OSRM route (i.e., it is the natural entry point into the core). **B→X links:** symmetric with X→B — keeps B→X only if no other boundary node appears in the OSRM route sequence (i.e., B is the last boundary node departed on the way to X). If another boundary node B' appears, the journey is already covered by B→B' + B'→X. **Boundary→boundary exterior shortcuts:** for each ordered boundary pair, if the route exits the core first, adds a directed shortcut with duration summed from OSRM annotations up to the first boundary node re-encountered. **Through-route allowlist:** for each ordered external-external pair, checks if any OSRM route node is a boundary node; if so adds to `allowed_through_pairs`. Outputs `data/external_links.json`. ~28,000 OSRM queries; under a minute on a local instance. |
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, detects boundary nodes, adds external node weights, writes `node_weights.json` + `newtownards_demographics.geojson`. **Does not build the map** — that moved to `build_map.py` (run it afterwards). `--zones-only` re-reads `data/census_zones.json` and patches only the external node entries in `node_weights.json`. Shared constants (paths, OSM tag handling, map styling) live in `simulation/demographics_config.py`. **Boundary node detection:** loads core polygon from `census_zones.json` and the **pbf-sourced** raw graph (`newtownards_network.graphml` from `build_network.py`), identifies internal nodes (within core polygon), then boundary nodes = internal nodes with at least one edge going outside. Writes `boundary_node_ids`/`internal_node_ids` to `node_weights.json` (replaces the previous hand-specified list). Because the raw graph now comes from the same OSM snapshot as OSRM, these IDs match OSRM's route node IDs exactly, so `build_external_links.py`'s boundary/internal route-sequence matching is no longer subject to Overpass-vs-pbf snapshot drift. **External node weights:** reads external node list from `census_zones.json` and writes population/workplace_pop to `node_weights.json`. **Study area = core polygon (not a circle):** DZ selection and all OSM downloads (buildings/POIs/parking) are bounded by the core polygon from `census_zones.json` (extent `max_core_vertex_dist_m`, ~10.2 km), matching the road graph built by `build_network.py`. Core DZs are selected by centroid-within the polygon (recovers exactly the `n_core_dzs` core DZs) and use **full** DZ population/workplace_pop (no area-fraction clipping — the legacy 3 km `RADIUS_M` circle is removed). OSM downloads use a circle sized to the polygon (+1 km margin); POIs and parking are then filtered to within the core polygon to avoid margin leakage (buildings are already DZ-bounded by sjoin). **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Business demand:** workplace population × POI count (school/college/university excluded), augmented with OSM car park polygon area (public: area/25, private: area/50 equivalent persons). **School demand:** separate `node_school_demand` layer from OSM school POIs. Enrollment fallbacks: school→300, secondary_school→900, college→2000, university→3000 pupils. External node school demand is `population × ext_school_per_pop` from `tuner_config.json` (fails loud if absent — run `simulation/compute_ext_school_scale.py` to derive the value). |
| `simulation/build_map.py` | **NEW.** Builds the interactive folium map (`newtownards_map.html`) from artifacts written by `build_demographics.py` (`node_weights.json`, `newtownards_demographics.geojson`), the road graphs, the cached OSM POI/parking layers, and — if present — `newtownards_flows.json`. This was the old `build_demographics.py --map-only` path, now a standalone step (it always reloads POI/parking from cache). Run after `build_demographics.py`, and again after `build_assignment.py` to refresh flow layers. **Flow map layers:** combined AADT (default), residential (teal), business (amber→red), school (violet→purple). No args (`--help` only). |
| `simulation/demographics_config.py` | **NEW.** Shared pure-constant config imported by `build_demographics.py`, `build_map.py` **and `build_network.py`** (file paths, OSM tag handling — `EXCLUDE_AMENITY`/`POI_WEIGHTS`/`SCHOOL_ENROLL_FALLBACK` — and map styling). `CENTRE` is re-exported from `zones_config.py` (not defined here). Also holds the road-network source knobs: **`PBF_PATH`** (absolute path to the NI `.osm.pbf` OSRM is built from) and **`BOUNDARY_BBOX_MARGIN_M`** (5 km buffer around the core polygon for `build_network.py`'s extract). The separate `NETWORK_MARGIN_M` (1 km) here sizes the OSM POI/building/parking download circle only — unrelated to the road graph. Also defines **`PROJECTED_CRS = "EPSG:2157"`** (Irish Transverse Mercator / ITM) — the single source of truth for all projected spatial operations in the pipeline. ITM covers the whole island of Ireland with uniform accuracy, avoiding UTM Zone 30N's distortion for Republic of Ireland towns west of ~6°W. All simulation and analysis scripts import this constant; `build_network.py` passes it explicitly to `ox.project_graph`. **Do not hardcode `EPSG:32630` anywhere.** Single source of truth so the split scripts don't drift. |
| `simulation/build_intra_times.py` | **NEW.** OSRM-samples intra-zonal travel times per external census zone for the production-suppression **self-term**. For each of the 166 external nodes, recovers its census polygon (SDZ→`SDZ2021_cd`, DZ→`DZ2021_cd`, DEA→`FinalR_DEA`, reprojected to WGS84), rejection-samples `M`=30 uniform point-pairs inside it, and routes each on the local OSRM (`localhost:5000`, same profile as `build_external_links.py`) → `data/external_intra_times.json` (`{census_code: [t1..tM seconds]}` + `_meta`). ~5,000 routes, seconds. Loud on any zone with a missing polygon or `<M` successful routes. `--m N` overrides the pair count. **Run after `build_census_zones.py`, OSRM up. Independent of `build_paths.py`** — the self-term lives in the model layer, so re-running needs no paths-cache rebuild. Re-run only when external zones change. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Now covers both internal road nodes and external census-area nodes. **Graph augmentation:** loads external nodes and edges from `data/external_links.json`; adds them to the routing graph before Dijkstra. External edges (X↔B, boundary shortcuts) are included in the adjacency matrix but NOT in `link_list` — they contribute to path distance but not to flow accumulation. **Probit loading:** all edges (road and external) perturbed each pass with log-normal noise `exp(eps·w)`, `eps ~ N(0, CV=0.25)`, N_PASSES=25, giving stochastic spread in boundary node selection for external→internal OD pairs with similarly-weighted entry options. **Length-scaled noise (`PROBIT_LL_SIGMA`, default 120 s ≈ 2 min):** the per-edge gain `w = σ_ll/(σ_ll + CV·cost)` ∈ (0,1] keeps the noise multiplicative for short legs (`w→1`) but saturates it to a fixed *absolute* sigma `σ_ll` for long legs (`w→σ_ll/(CV·cost)`), so a long single-edge external↔boundary leg's perturbation no longer swamps the few-minute differences between competing boundary entries. The adjusted perturbation never exceeds the pure multiplicative one, and `eps=0 ⇒ no bias`. Crossover at cost ≈ σ_ll/CV (~8 min). `PROBIT_CV` and `PROBIT_LL_SIGMA` are imported from `simulation/routing_config.py` (the gain vector is precomputed once, constant across passes). **OD pair filter:** through-routed external→external pairs (in `allowed_through_pairs`) are routed flow pairs (Dijkstra path through the core). Non-through external→external pairs (from `external_external_times`) are appended as **denominator-only** pairs — entries in `od_src/od_dst/od_dist` (distance = direct OSRM time) but NOT in `pair_idx/link_idx` and excluded from `src_groups`/probit passes, so they carry no flow; they complete each external origin's production-constrained denominator. The cache stamps `n_routed_pairs` (flow-carrying pairs occupy `0..n_routed_pairs-1`). No offscreen leg calculation. `HIGHWAY_COST_FACTOR` is imported from `simulation/routing_config.py`. Re-run if road network, external links, `HIGHWAY_COST_FACTOR`, `N_PASSES`, `PROBIT_CV`, or `PROBIT_LL_SIGMA` change. **Performance constants:** `N_WORKERS` (default 1) controls parallel pass workers via `multiprocessing.Pool` — increase on machines with sufficient RAM (each worker uses ~100–150 MB extra); `MAX_HOPS` (default 120) caps per-pair path-trace iterations. Inner path-tracing loop is vectorised (numpy). |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (legacy unconstrained rational kernel, kept for old param files), **`constrained_od_flows()`** (production-constrained per-component per-pair pre-K flows + per-origin denominators; optional `self_src`/`self_dist`/`self_w` add the **external intra-zonal self-term** to each denominator — denominator-only, `None` ⇒ exact prior behaviour) and **`scatter_od_to_links()`** (the production-constrained assignment core, used by `build_assignment.py` and `tune_assignment.py`), **`load_self_terms(node_ids)`** (builds the self-term arrays from `data/external_intra_times.json`; emits one entry per sampled time with weight `1/M_i`; skips zones absent from `node_ids`; returns `(None,None,None)` if the file is missing), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `gravity_assign()` accepts optional `link_weight` array (probit fractional weights). `compute_chi2()` dispatches to three-component mode when `link_flow_school_dict` provided (N_eff = N − 3·N_slots); two-component when only `link_flow_biz_dict` provided; legacy single-flow otherwise. Road node IDs are OSM integers (stable); external census node IDs are census-area-code strings (e.g. `"N21000219"`) — not OSM IDs. COUNT_SITES: site 507 links 538692601↔549139252; site 508 node 136173611; site 444 node 449111329. **`WEIGHTS_FILE` and `ROUTING_GRAPH` point at the dead-end-reduced artifacts** (`node_weights_reduced.json`, `newtownards_reduced.graphml`) from `reduce_deadends.py`. EXCLUDE_LINKS: `{(181844513, 181839481)}` plus the Westmount Park and Old Belfast Road directed links (both directions) whose endpoints are absorbed by `reduce_deadends.py` and no longer exist in the reduced graph — their walking observations are discarded from calibration (regenerate this set from `deadend_broken_obs.json` if the reduction params change). |
| `simulation/build_assignment.py` | **Production-constrained** gravity assignment (via `model.constrained_od_flows` + `scatter_od_to_links`). Requires `simulation/newtownards_paths.npz`. Three-component mode activated when K_sch > 0, `P_school` is in `tuned_params.json`, and `node_school_demand` is in `node_weights.json` (no longer requires `W_SCHOOL` — removed). Saves `flows_res`, `flows_biz`, `flows_school` in `newtownards_flows.json`. Falls back to two-component (K_res/K_biz only); legacy single-K unconstrained `gravity_assign` path kept for old param files. External node weights come from `node_weights.json` directly (no override from tuned params). |
| `simulation/routing_config.py` | **Single source of truth for `HIGHWAY_COST_FACTOR`.** Imported by `build_paths.py` (internal Dijkstra) and `build_osrm_profile.py` (OSRM profile generator). Edit factors here; both systems pick up the change. |
| `simulation/build_osrm_profile.py` | Generates `car_roaaads.lua` — the road-class-biased OSRM car profile. Pulls the default `car.lua` from the `osrm/osrm-backend` Docker image, injects a block after the `forward_rate` assignment that divides `forward_speed`/`forward_rate` by `HIGHWAY_COST_FACTOR` (matching internal Dijkstra biasing). Re-run whenever `HIGHWAY_COST_FACTOR` changes, then re-preprocess OSRM (`osrm-extract -p car_roaaads.lua`, `osrm-partition`, `osrm-customize`). Output: `/home/matthew/Documents/CodingFun/osrm/car_roaaads.lua`. |
| `simulation/reduce_deadends.py` | **NEW.** Collapses "residential dead-end" regions in the consolidated routing graph to shrink node count (speeds up `build_paths.py`/tuning, enables a larger core area). A region R (entrance E ∉ R) qualifies iff: (1) R connects to the rest of the network through exactly one cut vertex E; (2) R contains no boundary node and no school-demand node (both *protected* — never absorbed — which enforces the no-boundary and zero-school rules structurally); (3) max directed journey time E→n over n∈R < `T_MAX` (default 60 routing-cost seconds); (4) total business demand < `BIZ_CAP` (default 100; residential pop unbounded); (5) `|R| ≥ 2` (single-node spurs skipped — 1→1 saves nothing). **Algorithm:** every valid region is a protected-free connected component of H−a (H = undirected simple projection) for some articulation point a, so it enumerates all such (entrance, region) candidates, filters by constraints 2–5 + directed reachability both ways, and selects the *maximal feasible* regions (laminar family ⇒ disjoint; naturally descends into an oversized branch to find the largest collapsible sub-pockets — catches cyclic closes that leaf-pruning would miss). Each region → one super-node S (=min id, summed pop/biz/school, pop-weighted UTM centroid) joined to E by directed links E→S, S→E whose travel times are population-weighted means of the intra-region directed times. Synthetic edges use `highway="deadend_collapsed"` (factor 1.0 in `HIGHWAY_COST_FACTOR`) with `maxspeed`+`length` set so osmnx's `add_edge_speeds`/`add_edge_travel_times` (re-run by `build_paths.py`) reproduce the target `travel_time`. **Run after `build_demographics.py` (needs pop/biz/school + boundary) and before `build_paths.py`.** Outputs (gitignored): `newtownards_reduced.graphml`, `node_weights_reduced.json`, `deadend_map.json` (provenance: super-node→absorbed nodes + times), `deadend_broken_obs.json` (observed/count links whose endpoints were eaten — **manual review before adoption**; observed-link endpoints are deliberately *not* protected). Params: `--t-max`, `--biz-cap`. Current run: 1002→727 nodes (275 absorbed, 64 regions); build_paths.py runs ~5.3 s/probit-pass on the reduced graph with 0 fallbacks. **Wired into the pipeline:** `build_paths.py` (`CONS_GRAPH`), `build_assignment.py` (`CONS_GRAPH`), `tune_assignment.py` (`CONS_GRAPH`) read `newtownards_reduced.graphml`, and `model.WEIGHTS_FILE`/`ROUTING_GRAPH` point at the reduced files — so this step must run after `build_demographics.py`. The 6 absorbed walking observations (Westmount Park, Old Belfast Road) are discarded via `EXCLUDE_LINKS` in `model.py`. **Map caveat:** `build_map.py` still draws the *full* consolidated graph, so flow on collapsed interior streets is not shown on the map (demand layers are unaffected; main-road flows and the fit are unaffected). Re-mapping collapsed regions via their super-nodes is a possible follow-up. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Gravity param regularization and temporal coupling priors. `gravity_lambda` + `gravity_ref` regularise P/ALPHA/BETA/W_BIZ/P_biz/ALPHA_biz; `gamma_coupling_scale` controls the per-slot aggregate coupling (γ = scale/std_f²); `phi_biz_prior`/`phi_biz_std` and `phi_school_prior`/`phi_school_std` set Gaussian priors on business and school flow fractions. `ext_biz_scale` (current value 1.702) scales external node business demand in `build_demographics.py` to compensate for the OSM car park area contribution that internal nodes receive but external census nodes don't (ratio of total business demand to total NISRA workplace pop across all internal nodes, including boundary). Recomputed 2026-06-22 after the core-polygon demographics fix widened internal coverage: total internal biz 22,517 / workplace 13,230 = 1.702. `ext_school_per_pop` (current value 0.159913) is the pupils-per-person ratio from the core area, applied to external nodes as `population × ext_school_per_pop` — a uniform approximation that lets the school self-term and school component activate for external zones. Computed by `simulation/compute_ext_school_scale.py` (core internal nodes: 33,143 pop, 5,300 pupils). **Required** — `build_demographics.py` fails loud if absent; re-run the script and update this value if core OSM school POIs change significantly. **Removed:** `cities` block (replaced by `census_zones.json`) and `through_route_pairs` whitelist (replaced by OSRM-derived `external_links.json`). `lambda` is retained but no longer used (external zone params are not tuned). |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. Run once (or when the ODS file changes). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. Loads manual link overrides from `data/manual_link_overrides.json`. After every new link assignment, validates each non-null count direction against the directed graph; raises `ValueError` if the edge doesn't exist. |
| `analysis/manual_assign_link.py` | CLI tool to manually assign a session to a specific directed link, bypassing GPS snap. Usage: `python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>`. Validates both nodes exist and checks count-edge consistency. Writes to `data/manual_link_overrides.json` and patches `counts_processed.json` directly. After correcting an assignment, re-run `aggregate_counts.py` then `tune_assignment.py`. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s`. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Three-component, production-constrained model:** gravity flows split into residential (`flow_res`), business-adjacent (`flow_biz`), and school (`flow_school`) components, each singly (production) constrained. Tunes **8 gravity params** (W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, P_school, ALPHA_school) — `W_SCHOOL` removed (redundant with K_sch under the constraint). External zone values are fixed from census data and are not tuned. **Alternating minimisation (5 blocks) SURVIVES the constraint** (`D_i` has no `K`, so flow stays linear in K and the f's): K-step (1D); phi_biz-step; phi_sch-step; f_res/f_biz/f_school steps (per-slot analytical) + aggregate coupling γ per slot. **The alternation is non-monotonic (single Poisson Newton steps) and can collapse K to the 1e-30 floor, so `calibrate_Ks_and_fracs` returns the BEST-SEEN (K,φ,f) iterate by the full regularized objective, max_iter 40 (`--fast` 20) — see "Five-block analytical calibration" (fix 2026-06-25, commit af90de7).** `run_assignment` now calls `model.constrained_od_flows` (per-pair flows + per-origin denominator bincounts) and scatters via the probit routing incidence — the old distance-bin-matrix path is gone. **Performance:** ~2.5–3.5 s/eval (run_assignment-dominated) ⇒ ~20–30 min per run. `CALIBRATE_PROBE=1` is an env-gated diagnostic that reports the post-calibrate residual global scale λ. |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` and `reports/slot_pulls_{id}.png`. History `slot_prior` entries carry 5 values: `[mean_f_agg, std_f, mean_f_res, mean_f_biz, mean_f_school]`. Old entries with 4 values handled gracefully. Note: report still attempts to print external city delta table — will be a no-op for new-format runs. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. New-format param dicts no longer contain external zone keys. |
| `simulation/reset_gravity_params.py` | Reset only the gravity params (K, W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz) in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile. Includes `mean_fraction_res`, `mean_fraction_biz`, and `mean_fraction_school` columns. **Derived from NTS** via `analysis/derive_component_profiles.py`. Constraint: res + biz + school = agg for all 168 rows. Rows sum to 7.0 (AADT normalisation). Re-run `derive_component_profiles.py` whenever NTS files change. |
| `analysis/derive_component_profiles.py` | Derives `mean_fraction_res`, `mean_fraction_biz`, and `mean_fraction_school` from DfT NTS data (2023–2024 rolling average). Purpose classification: biz = commuting + employer's business + shopping; school = education×⅕ + escort education; res = remainder. Re-run whenever NTS files or purpose classification change. |
| `analysis/google_routing_common.py` | **NEW (Google calibration).** Shared pure-stdlib helpers for the Google routing-time calibration tooling: encoded-polyline decode, polyline downsampling, OSRM `/route` and `/match` calls, Google Routes API v2 `computeRoutes` call, `CONF_MIN` (0.5 /match-confidence floor). No third-party deps beyond `networkx` (used only by the manifest builder). See **Google Routing-Time Calibration** section below. |
| `analysis/google_feasibility.py` | **NEW (Google calibration — pilot).** One-shot feasibility experiment: a small hardcoded data-driven OD sample → Google Routes, decode polyline, OSRM `/match` (time error on Google's geometry) vs OSRM `/route` (route-choice divergence). Matches best route + all alternatives; caches raw responses in `data/google_cache/`. `--dry-run` makes no API calls; key only needed on a live cache-miss. Superseded for production sampling by `build_od_manifest.py` + `google_query_routes.py`, kept as the reference pilot. |
| `analysis/build_od_manifest.py` | **NEW (Google calibration).** Writes the fixed, deterministic (seed 20260622), model-aligned, length-skewed OD sample → `data/google_cache/od_manifest.json`. **Makes NO API/OSRM calls** — reads only `census_zones.json`, `external_links.json`, `node_weights.json`, and the raw graph. Leg types + default quotas (of `--n`, default 1000): X2B 45% (external centroid→boundary entry), X2X 25% (allowlisted through-routes), B2X 15% (boundary→external), INT 15% (internal→internal, for in-town junction realism). Within each leg type, 4 length quartile-bands allocated with a long-skew (`BAND_WEIGHTS` 0.15/0.20/0.30/0.35); X2B/B2X length = model `duration_s`, X2X/INT length = haversine. Road-class×speed-band coverage is **emergent in v1** (not stratified on realised road mix — check after a first batch; possible refinement = pre-route candidates with local OSRM). |
| `analysis/google_query_routes.py` | **NEW (Google calibration — runner).** Crash-safe, resumable runner over `od_manifest.json`. **Phase A (spendy, resumable):** queries each uncached OD and writes its raw Google response to `data/google_cache/raw/<od_id>.json` **atomically and immediately** (temp+rename), so a crash loses at most the in-flight query; re-running skips cached ODs; `--limit N` caps queries per run to batch spend. **Phase B (free, idempotent):** rebuilds `data/google_cache/results.jsonl` by running OSRM `/match` over all cached raw responses (best + alternatives) + OSRM `/route` per OD — no API calls, re-run any time via `--reprocess-only`. `--dry-run` reports counts/cost with no spend. **Refuses to start Phase A without `GOOGLE_MAPS_API_KEY`, and a live run requires explicit per-run user approval (see Agent Behaviour).** |
| `simulation/profile_spec.py` | **NEW (profile calibration — single source of truth).** Pure-stdlib definition of a calibrated OSRM time profile: a grid of multiplicative speed **factors** per `(highway_class × speed_band)` bucket (full `DRIVE_HIGHWAYS` classification × NI mph bands `{untagged,20,30,40,50,60,70,other}`) + the four global turn params (`turn_penalty`, `traffic_light_penalty`, `u_turn_penalty`, `turn_bias`). `factor=1.0` = stock-OSRM base speed; `factor>1 ⇒ slower` (OSRM is currently too fast). Holds the bucketisation (`norm_class`/`parse_band`/`band_from_tags`/`bucket_of`/`bucket_index`), the stock base-speed table, `base_speed_for`, and `ProfileSpec` (JSON load/save). **maxspeed resolution mirrors OSRM's `WayHandlers.maxspeed` exactly:** `bucket_of(tags)` takes a way's full tag dict and `band_from_tags` honours OSRM's key precedence (`maxspeed:advisory` > `maxspeed` > `source:maxspeed` > `maxspeed:type`); `osrm_maxspeed_kmh` resolves numeric *and* symbolic/national-speed-limit values (`gb:nsl_single`→60, `gb:nsl_dual`/`gb:motorway`→70, `none`→140, plus `maxspeed_table_default` urban/rural/trunk/motorway) — so nsl roads land in their real speed band instead of collapsing into `untagged`. **Replaces `routing_config.py`/`HIGHWAY_COST_FACTOR` for the calibration work** (the old module still feeds `build_paths.py`/`build_osrm_profile.py` until a calibrated profile is adopted). Imported by both the stdlib `analysis/` tooling and the simulation-side Lua generators, so the offline model and the emitted Lua key on the *same* buckets. |
| `simulation/osrm_lua.py` | **NEW (profile calibration).** Shared OSRM car.lua/Docker plumbing: `pull_base_lua`, `copy_lib`, the 3-strategy injection-point `find_injection_point`/`inject` (refactored out of `build_osrm_profile.py`), and the Lua emitters `emit_probe_block` (legacy probe — now unused; encodes a way's bucket id as its speed) + `emit_factor_block(spec)` (divides speed by the tuned per-bucket factor) + `apply_turn_overrides`. The bucket-index Lua (`_bucket_index_lua`) is generated from `profile_spec`'s `CLASSES`/`BANDS`/`MAXSPEED_KEYS`/symbolic tables and **replicates the full nsl-aware band resolution** (key precedence + symbolic/national-speed-limit lookup), cross-checked against the Python `bucket_of` via a Lua harness, so the compiled profile and the Python offline model can never disagree about a way's bucket. |
| `simulation/build_edge_index.py` | **NEW (profile calibration — raw OSM cache; replaces the probe).** The rounding-free replacement for the OSRM probe. `--match`: the **single** `/match` pass over the route set on the deployed OSRM (:5000) — caches full match detail per route (node sequence + per-segment `distance` + step maneuvers + Google duration) to `data/google_cache/match_cache.jsonl` (resumable; `--limit N` batches the slow ~1.7 s/match pass). `--extract`: streams the NI pbf via the `osmctools-roaaads` image (osmconvert→o5m, `osmfilter --keep="highway="`) into an all-highway `.osm`, then a low-RAM stdlib `xml.etree.iterparse` pass (root-cleared) writes **the complete tag dict of every way and node touching the route set** + geometry → `osm_ways.jsonl` / `osm_nodes.jsonl` (+ `edge_index_meta.json`). **Deliberately raw — no bucketing/filtering here;** the cache stores everything (lanes, surface, oneway, lit, junction, crossing, signals, …) so a future model can use it without re-querying. No new pip dep (pyosmium has no Python-3.8 wheel; streaming iterparse + a matched-node filter is the low-RAM substitute). |
| `simulation/build_skeleton_index.py` | **NEW (profile calibration — skeleton builder).** Rebuilds the profile-independent `data/google_cache/skeletons.jsonl` from `match_cache.jsonl` + the edge index — **no OSRM calls**, pure recompute, so it is free to re-run after any `profile_spec`/bucket change. Each matched segment `(node_u,node_v)` is resolved to its way's tags via the edge index and bucketed with `profile_spec.bucket_of` (exact node-id lookup, not the old probe `annotation.speed` readout that corrupted short urban edges) → `length_by_bucket`; `turns` from cached step maneuvers, `n_signals` from cached node tags (`highway=traffic_signals`), `coverage`/`valid` from geometry. `--base-speeds`: the one remaining `/match` step — samples ~800 routes on a factor-free speed source (stock OSRM with `--no-defactor`, or `:5000` defactored), labels each segment via the edge index, and writes **length-weighted harmonic-mean** per-bucket base speeds → `base_speeds.json`. |
| `analysis/skeleton_model.py` | **NEW (profile calibration — fast offline model).** Pure-stdlib `predict_duration(skel, spec)` = `Σ_bucket factor·length·3.6/base_speed` (edge) + OSRM-style turn sigmoid (gated on degree>2 / u-turn, NI left-hand bias) + `n_signals·traffic_light_penalty`. `evaluate(skeletons, spec)` scores against Google with squared-log-ratio loss, **equal weight per valid route** (conf≥0.5 & coverage in-band), and returns per-leg-type / per-(leg×band) diagnostics. `bucket_coverage` (factor identifiability) and `legacy_spec_from_highway_cost_factor` (deployed-profile reference) included. Milliseconds for the whole cache — no OSRM/Docker. |
| `analysis/eval_profile.py` | **NEW (profile calibration — benchmark entrypoint).** Scores a `ProfileSpec` JSON (default all-1.0 stock; `--spec`; `--legacy-factors`) against `skeletons.jsonl`: aggregate loss, predicted/Google ratio distribution, per-leg-type + per-cell breakdown, per-bucket coverage table, turn-time fraction. No spend. `--legacy-factors` is the faithfulness sanity check (should track the deployed `te_matched`). |
| `simulation/compile_profile.py` | **NEW (profile calibration — compiler).** `tuned_profile.json` (a `ProfileSpec`) → deployable `car_roaaads.lua`: applies tuned turn params in `setup()`, injects `emit_factor_block`, copies `lib/`, prints the re-extract/partition/customize commands for the deployed :5000 instance. |
| `analysis/verify_profile.py` | **NEW (profile calibration — fidelity gate).** After the deployed OSRM is rebuilt from a compiled profile (or against the live deployed instance with `--legacy-factors`), `/match`es a validation subset through real OSRM (:5000) and compares real `match_dur` to `predict_duration(skel, spec)`. **Gate (median-based):** per-leg median `predict/real` within ±`--gate-median-tol` (default 0.03) over `--gate-legs` (default external `X2B,B2X,X2X`); per-route scatter (med/p90 |resid|) is reported but not gated (inherent to probe-vs-deployed re-match). INT reported, not gated. Exits non-zero on fail. Read-only, no Google calls. |
| `analysis/tune_profile.py` | **NEW (profile calibration — external-focused tuner).** Fits per-`(class×band)` speed **factors** to minimise the weighted squared-log-ratio time error vs Google over `skeletons.jsonl`, leg-weighted (default `X2B/B2X/X2X=1, INT=0`). With turn params + base speeds fixed, predicted time is **linear in the factor vector**, so it's a vectorised (numpy) scipy `L-BFGS-B` fit with L2 reg toward 1.0 and `[0.2,5]` bounds; only buckets above `--min-km` weighted coverage are tuned (rest stay 1.0). Writes `simulation/tuned_profile.json` + appends `profile_tuning_history.jsonl`; reports before/after per-leg medians + top factor moves. **Tunes factors only** — global turn params are held at defaults (external turn fraction is small; INT excluded) until the in-town turn model is improved. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_map.html`, `simulation/tuned_params.json` — all regenerated by the pipeline.
`simulation/node_weights.json` keys: `node_population`, `node_business_demand`, `node_school_demand`, `node_parking_equiv`, `boundary_node_ids` (auto-detected from core polygon). External node entries (census-area-code string IDs, e.g. `"N21000219"`) are included alongside internal OSM node IDs. `node_effective_utm` is removed.
**Node ID scheme (as of 2026-06-20):** `build_network.py` relabels consolidated graph nodes to stable OSM IDs after junction consolidation. A consolidated junction (multiple OSM nodes merged) gets `min(osmid_original)` as its ID; a non-merged node gets `int(osmid_original)`. All road node IDs are therefore genuine OSM node IDs (in the hundreds of millions) and stable across graph regenerations. External census nodes use their **census-area-code string IDs** (SDZ/DEA 2021 codes, e.g. `"N21000219"`) — these are the `id` values in `census_zones.json`, *not* small integers (downstream code that consumes `external_links.json`/`census_zones.json` must treat external IDs as strings — `build_paths.py` does). Road node IDs are ints, external node IDs are strings; `node_to_idx` mixes both.
`simulation/newtownards_flows.json` — combined flows plus optional `flows_res`/`flows_biz`/`flows_school` keys when three-component params active.
`simulation/newtownards_reduced.graphml`, `simulation/node_weights_reduced.json` — dead-end-reduced routing graph + weights from `reduce_deadends.py`; consumed by `build_paths.py`, `build_assignment.py`, `tune_assignment.py`. `simulation/deadend_map.json` (super-node→absorbed-nodes provenance + link times) and `simulation/deadend_broken_obs.json` (observed links eaten by collapse) are also written here.
`reports/` — generated by `report_tune.py` and `tune_assignment.py`; not tracked.
`data/google_cache/` — **gitignored** (Google ToS: cached responses kept local, never
committed/redistributed). Holds `od_manifest.json` (the fixed OD sample), `raw/<od_id>.json`
(one raw Google response per OD — the resumable cache and re-processing source of truth),
`results.jsonl` (derived OSRM-match metrics, rebuilt for free from `raw/`), `match_cache.jsonl`
(the single cached `/match` pass per route — node sequence + per-segment distances + maneuvers +
Google duration, from `build_edge_index.py --match`; the slow ~1.7 s/route artifact everything
else derives from for free), `osm_ways.jsonl` / `osm_nodes.jsonl` / `edge_index_meta.json` (the
raw OSM edge index from `build_edge_index.py --extract` — **every tag of every way/node along the
route set**, plus geometry; bucketed downstream), `skeletons.jsonl` (profile-independent route
skeletons from `build_skeleton_index.py` — the fast-benchmark cache), and `base_speeds.json`
(empirical realised per-`(class×band)` base speeds from `--base-speeds`; auto-loaded by the
offline model, overrides the analytical estimate), and `profile_tuning_history.jsonl` (one line
per `tune_profile.py` run). Survives worktree removal (lives in the main checkout); only at risk
from `git clean -xfd` or manual `rm`. `simulation/tuned_profile.json` (a candidate `ProfileSpec`,
gitignored) is also generated. (The legacy probe OSRM instance / `car_probe.lua` / `:5001` and
`signal_nodes.json` are no longer used — superseded by the edge index.)

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`, `data/official_hourly.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`data/census_zones.json` — committed; output of `build_census_zones.py`. Contains core polygon, external node list with IDs/codes/centroids/census demand. Re-run `build_census_zones.py` only if NISRA boundary files or census data change.
`data/external_links.json` — committed; output of `build_external_links.py`. Contains OSRM-derived X↔B links, boundary shortcuts, and through-route allowlist. Re-run `build_external_links.py` when boundary nodes change or OSRM data is updated.
`data/external_intra_times.json` — committed; output of `build_intra_times.py`. Per external zone, `M`=30 sampled intra-zonal OSRM times (s) for the production-suppression self-term (`model.load_self_terms` → `constrained_od_flows`). Committed so the model runs without re-querying OSRM. Re-run `build_intra_times.py` only when external zones change.
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
> **As of the production-constrained change (committed on `main`, code in `model.py`
> `constrained_od_flows`; pending end-to-end re-tune/validation):** the assignment is no longer
> the unconstrained product below. Each component is **singly (production) constrained**:
> `T^c_ij = K_c · p^c_i · a^c_j · f_c(d_ij) / D^c_i`, `D^c_i = Σ_k a^c_k·f_c(d_ik)`, so each origin's
> trip production is fixed by its producing weight `p^c_i` and is independent of accessibility
> (fixes the generation/distribution conflation). Per component: res `p=a=pop`; biz = symmetric
> pop↔biz split (per-origin normalised) + `W_BIZ`·(biz×biz, constrained on biz); school = symmetric
> pop↔school split (magnitude `K_sch`; `W_SCHOOL` removed as redundant). `W_BIZ` no longer appears
> in a combined node weight — components are separate. The kernel `f` and its params below are
> unchanged; only the OD-flow structure changed. See the agent memory note
> `project_production_constrained_gravity`. The lines below describe the *prior* unconstrained model.

OD flow (prior, unconstrained): `T_ij = K × w_i × w_j × f(d_ij)`

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

Once a DEA is broken into SDZs, all its constituent SDZs become nodes (even those beyond `SDZ_ZONE_RADIUS`). Each external node's ID **is** its census-area code (the SDZ/DEA 2021 code, a string such as `"N21000219"`) — the `id` field in `census_zones.json`; there is no separate small-integer ID. (Road node IDs are OSM integers; external node IDs are these strings.)

**Demand:** population and workplace_pop from Census 2021 (DZ-level data aggregated to SDZ/DEA). Business demand for external nodes is `workplace_pop × ext_biz_scale` (TBD refinement — units differ from internal OSM POI proxy). School demand for external nodes is `population × ext_school_per_pop` (a uniform approximation: same pupils/person ratio as the core area; computed by `simulation/compute_ext_school_scale.py`).

**Connectivity:** boundary nodes are all internal nodes with at least one road edge crossing the core polygon boundary. OSRM-derived directed edges connect each external node to its valid boundary nodes; the gravity model path distance for any external→internal pair is the OSRM leg plus the internal shortest path, computed automatically by Dijkstra on the augmented graph.

**No tuning of external zone values.** Pop/wp come directly from census data and are not adjusted by the optimizer. This replaces the previous hand-crafted city configs, ref_pop/ref_wp values, and damping factors.

### Through routes
External→external OD pairs are allowed only for pairs whose OSRM route passes through at least one boundary node (i.e. genuinely transits the core area). This allowlist is auto-generated by `build_external_links.py` and stored in `data/external_links.json`. Changing the allowlist requires re-running `build_external_links.py` and `build_paths.py` (the paths cache must be rebuilt).

### Three-component flow decomposition
The gravity OD flows are split into three **production-constrained** spatial components at each
tuner evaluation (per-pair pre-K flows from `model.constrained_od_flows`, scattered onto links;
`p^c_i`/`a^c_j` = producing/attracting weight, `D^c_i = Σ_k a^c_k·f_c(d_ik)`):

- **Residential** (`flow_res`): `T^res_ij = pop_i·pop_j·f_res(d_ij)/D^res,pop_i` — pop×pop trips.
- **Business-adjacent** (`flow_biz`): symmetric pop↔biz split, each leg normalised by its own origin
  denominator, plus a `W_BIZ`-weighted biz×biz term constrained on biz —
  `f_biz·( pop_i·biz_j/D^biz,biz_i + biz_i·pop_j/D^biz,pop_i + W_BIZ·biz_i·biz_j/D^biz,biz_i )`.
- **School** (`flow_school`): symmetric pop↔school split, per-origin normalised, (P_school,
  ALPHA_school, BETA) kernel — `f_sch·( pop_i·school_j/D^sch,sch_i + school_i·pop_j/D^sch,pop_i )`.
  Magnitude is `K_sch` (no `W_SCHOOL`). **KNOWN-BROKEN (data gap, not code):** external nodes have
  `school_demand=0`, so every external origin's school budget dumps into the core ("the whole world goes
  to school in Newtownards"). The intra-zonal self-term machinery (below) that would let external zones
  retain school trips is fully wired but inert until external `school_demand` is populated — see memory note.

**External intra-zonal self-term (denominator-only).** Each per-origin denominator
`D^c_i = Σ_k a^c_k·f_c(d_ik)` runs over *other* zones; collapsing an external zone to one centroid
drops its `k=i` diagonal (its intra-zonal trips), so `D^c_i` is too small and the external origin's
fixed budget over-allocates to the observed core (worst for large, isolated, far DEAs). `build_intra_times.py`
OSRM-samples `M=30` uniform intra-zonal point-pairs per external zone → `data/external_intra_times.json`;
`model.constrained_od_flows` then adds `a^c_i·(1/M)·Σ_m f_c(t_im)` to each denominator (the mean kernel
over the sample, `E[f]`, not `f(mean)`). It is **denominator-only** — no link flow — and applies to
**external zones only** (internal road nodes have no zone area). Direct OSRM sampling avoids any
characteristic-distance constant, speed assumption, or zone-shape model (real per-zone times carry local
speed + network detour). Effect (measured): exported external→core budget −0.8% overall, concentrated on
the far DEAs (which collapse to ~0 core-bound flow, as they should — bounded because coarse ⇒ distant ⇒ tiny
`f(d→core)`); near-town SDZs barely move. Wired into `build_assignment.py` and `tune_assignment.py` via
`model.load_self_terms`; absent file ⇒ no self-term (exact prior behaviour). **Re-tune required** to realize
the benefit — at the pre-self-term params the removed flow makes the fixed-param χ²/N worse; the tuner
absorbed the bias (likely into the `ALPHA` tail), so re-running `tune_assignment.py` is what converts the
structural correction into a fit gain. Independent of the paths cache (model-layer, not a routing input) —
no `build_paths` rebuild needed.

The **school** component's self-term is fully wired (`f_sch_self` feeds both school denominators with the
correct pop/school attraction) but **currently inert** because external `school_demand=0`. The instant an
external school_demand is populated (in `build_census_zones.py`/`build_demographics.py`), the school
self-term activates with no code change — external zones retain school trips intra-zonally — directly
mitigating the known-broken "whole world goes to school in Newtownards" dump (verified on the real cache
with a synthetic external school_demand).

Each component has its own temporal profile and scale (K_res, K_biz, K_sch).
Predicted count for observation i in slot s:
`pred_i = K_res·flow_res·(T/3600)·f_res[s] + K_biz·flow_biz·(T/3600)·f_biz[s] + K_sch·flow_school·(T/3600)·f_school[s]`

### Five-block analytical calibration
At each optimizer evaluation, (K, phi_biz, phi_sch, f_res, f_biz, f_school) are calibrated via
alternating minimisation (max_iter 40; `--fast` 20).
K_res = K·(1−phi_biz−phi_sch), K_biz = K·phi_biz, K_sch = K·phi_sch.

**Non-monotonic — returns the best-seen iterate (fix 2026-06-25, commit af90de7).** The K-step's
Poisson (walking) Newton correction and the per-slot f-steps' Poisson corrections are *single*
Newton steps that do **not** guarantee descent, so the alternation is non-monotonic and under the
production-constrained flow magnitudes can **collapse K toward the 1e-30 floor** at some iteration
counts. `calibrate_Ks_and_fracs` therefore evaluates the full regularized objective
(Gaussian + Poisson deviance + f-prior/coupling penalty) each iteration and **returns the
lowest-objective (K, φ, f) state seen**, not the last iterate; it runs to `max_iter` with no early
break (the objective oscillates, so a premature stop can miss the good basin). Before this fix the
prior run stopped at a partially-collapsed point, leaving every modelled flow ~6× too small (a
uniform global scale λ≈6.2 halved χ²). `CALIBRATE_PROBE=1` runs calibrate at the start params for
several `max_iter` values and reports the residual global scale λ (≈1 ⇒ K at its optimum).
**TODO (cleaner follow-up):** harden the Poisson Newton steps (line-search / clamp) so the
alternation is monotone and best-iterate tracking becomes belt-and-braces rather than load-bearing.

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
All 890 observations are in count space with per-obs weights:
- **Official hourly** (216 obs, 24 h × 3 day-types × 3 sites): from `data/official_hourly.json`;
  Gaussian error (sigma from between-weekday std, 10% floor); weight = 1/sigma².
- **Walking** (343 obs): from `data/link_aadt.json`; Poisson error; weight = 1/n_eff.

### Goodness of fit
`χ²/N` (mean squared z-score; N=890 obs, N_eff = N − 3·N_slots = 890 − 216 = 674).
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

**Walking counts:** 10 CSV files, 357 sessions, 674 per-session observations (after EXCLUDE_LINKS). Manual overrides re-entered with OSM node IDs: sessions e644eae2, 760b0c8e, bb934ba7 → link 86604223↔86604221 (A20 Kempe Stones); e66989f4 → 150995265↔6622295361; b0043fd1, 32d425d6 → 181844516↔4688250384; 29d5f5f0 → 181844516↔538692566.
New sessions added 2026-06-18 (7 sessions): Saratoga Avenue, Glenford Road, Hardford Link, Belfast Road.
New sessions added 2026-06-23 (97 sessions, recorded 2026-06-20 to 2026-06-23).
**Total: 890 observations (216 official hourly + 674 walking) in 72 time slots. N_eff = 890 − 3×72 = 674.**
(Figures above reflect the last pre-reduction tune. Dead-end reduction additionally excludes
18 walking session-obs across 6 directed links — Westmount Park and Old Belfast Road, absorbed
into super-nodes — via `EXCLUDE_LINKS`. The live observation count is printed by
`tune_assignment.py` at run start; re-tune to refresh these figures.)

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

**Note on comparability:** runs from 2026-06-17 onward use the two-component model with coupling penalty terms in chi²/N; not directly comparable to earlier single-component runs. From 2026-06-19 three-component model: N_eff = 559 − 3×72 = 343 (one extra df per slot for f_school). After 2026-06-23 data addition: N_eff = 890 − 3×72 = 674. **Runs from the big-world architecture are not directly comparable to earlier runs** — external zone representation has fundamentally changed (census-derived vs hand-crafted; many more external nodes; OSRM-based connectivity vs offscreen Euclidean leg).

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
- `N_PASSES`, `PROBIT_CV`, or `PROBIT_LL_SIGMA` change

**Staleness guard (loud failure):** `build_paths.py` stamps a signature of its inputs into the
`.npz` — SHA-1 of the reduced routing graph (`newtownards_reduced.graphml`), SHA-1 of
`data/external_links.json`, the serialised `HIGHWAY_COST_FACTOR`, and the probit noise params
`PROBIT_CV` / `PROBIT_LL_SIGMA` (all from `simulation/routing_config.py`). `tune_assignment.py` and
`build_assignment.py` call `model.assert_paths_cache_fresh(cache)` right after loading the cache
and **raise `SystemExit`** if any input changed (or the cache predates the guard), naming what
changed and telling you to re-run `build_paths.py`. This replaces the previous silent-stale
footgun. Helper lives in `simulation/model.py` (`paths_cache_signature`, `assert_paths_cache_fresh`).

**Current cache format** (probit): `node_ids` covers road nodes (OSM integer IDs) + external nodes (census-area-code strings, e.g. `"N21000219"`); `link_u`/`link_v` are road-link endpoints only (external edges are not in `link_list`); `link_weight` (float32, fraction of passes using that link for each OD pair); `od_dist` (mean path distance across passes including external legs); `probit_n_passes`, `probit_cv`, `probit_ll_sigma`. `n_routed_pairs` marks the flow-carrying OD pairs occupying indices `0..n_routed_pairs-1`; the remainder are **denominator-only non-through ext→ext virtual edges** (entries in `od_src/od_dst/od_dist` but NOT in `pair_idx/link_idx` — they complete each external origin's production-constrained denominator and carry no link flow). No `pair_idx_2/3` keys — `_has_stoch = False`, THETA not tuned.

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
- **External node probit loading:** all edges (road and external) receive log-normal noise each pass (CV=0.25), length-scaled by the per-edge gain `w = σ_ll/(σ_ll + CV·cost)` (`PROBIT_LL_SIGMA`, default 120 s). Route diversity for external-internal OD pairs comes from both the X→B external leg and the internal B→J portion, giving stochastic spread across similarly-weighted boundary entry points. **Resolved (2026-06-23, length-scaled noise):** the noise was previously purely multiplicative, so a long external leg (e.g. a distant DEA's 90-min X→B edge) received a ±~22-min perturbation that swamped the few-minute differences between boundary nodes, making boundary entry effectively random. The length-scaled gain now caps a long leg's perturbation at an *absolute* sigma of `σ_ll` (~2 min) while leaving short internal edges' multiplicative noise essentially unchanged (`w≈1`), so boundary selection is driven by real time differences. This supersedes the earlier proposed `PROBIT_CV_EXT` (separate CV for external edges) — a single length scale handles both regimes smoothly, with no mean bias and the adjusted perturbation never exceeding the multiplicative one. Default `σ_ll=120 s` is anchored to `CV × a typical in-town journey (~8 min)`, so external legs get the same absolute route-choice jitter an internal journey already experiences; the knob lives in `simulation/routing_config.py` and is part of the paths-cache staleness signature.
- **Dundonald virtual node (10000) is removed** in the big-world system. Dundonald is now represented by an SDZ or DEA external centroid node with a proper census-derived population.
- **Manual link overrides:** Use `analysis/manual_assign_link.py <session_id> <from_node> <to_node>` to assign a session to a specific directed link, bypassing GPS snap. Use when the observer stood on a parallel carriageway and the snap would land on the wrong physical road. The override is stored in `data/manual_link_overrides.json` and takes effect even if `counts_processed.json` is wiped and rebuilt. After assignment (manual or auto), `ingest_counts.py` validates each non-null count direction against the directed graph and raises `ValueError` if the edge doesn't exist.
- **Snap direction bug (fixed 2026-06-15):** `ingest_counts.py` previously stored canonical
  `(min(u,v), max(u,v))` — fixed to store actual directed `(u, v)`. Only session `f56b2ce4`
  was materially affected (re-snapped from 22→159 to 159→22).
- Two temporal profiles (f_s_res, f_s_biz) are inferred per (day_type, hour) slot, each
  anchored by component-specific priors from `hourly_fractions.csv`. The aggregate coupling
  (gamma_coupling_scale / std_f²) per slot keeps their sum near f_agg. With 72 slots
  and 2 df each, N_eff = N − 2×N_slots = 890 − 144 = 746.
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
- **`tuned_params.json` structure:** contains `K_res`, `K_biz`, `K_sch`, `P_school`, `ALPHA_school`, `slot_fracs_res`, `slot_fracs_biz`, `slot_fracs_school` (dicts keyed `"dt,h"`); **no longer contains `W_SCHOOL`** (removed — redundant with K_sch under the production constraint); does **not** contain a `slot_fracs` key (legacy) or `external_node_pop/biz/city_pop/wp/dampings` keys (removed). Old param files without school keys fall back to two-component or legacy mode in `build_assignment.py` (a stale `W_SCHOOL` key in an old file is read but ignored).

---

## External Zone Configuration

External zone values are now fully data-driven from Census 2021 (via `data/census_zones.json`) and OSRM routing (via `data/external_links.json`). There are no hand-crafted reference values, dampings, or city groupings to maintain.

**Gravity param refs** (`tuner_config.json` `gravity_ref`): still anchored for L2 regularization. Last updated 2026-06-19 — `gravity_ref P = 600s`, `gravity_lambda P = 2.0`. These must not be changed without explicit approval.

**To update external zone coverage** (e.g. after a NISRA boundary update):
1. Re-run `build_census_zones.py` (updates `data/census_zones.json`)
2. Re-run `build_demographics.py` (updates `node_weights.json`)
3. Re-run `build_external_links.py` (updates `data/external_links.json`)
4. Re-run `build_paths.py` (rebuilds paths cache with new external nodes)
5. Re-tune, then `build_assignment.py` and `build_map.py`

**Outstanding TBDs for external nodes:**
- School demand: now set as `population × ext_school_per_pop` (uniform core ratio, ~0.160 pupils/person). A better approach would use census school-age population per zone; this approximation is intentionally ad-hoc and flagged in `tuner_config.json`.
- Business demand: currently set equal to workplace_pop. Units differ from internal nodes (which use OSM POI proxy). A separate `W_BIZ_ext` or normalisation may be needed.

---

## Google Routing-Time Calibration (offline, optional — NOT part of the main pipeline)

**Purpose.** OSRM (current profile) is systematically *too fast*, badly so on the external
approach corridors (e.g. Ballyrainey), which inflates external→core flow and hurts the fit.
This workflow uses **Google Maps Routes API as a source-of-truth for journey times** to
calibrate a more realistic OSRM time profile. The guiding design is to **decouple impedance
(realistic travel time) from route preference (generalised cost)** — currently conflated in
the single `HIGHWAY_COST_FACTOR`. Feasibility confirmed OSRM is ~26% too fast on average
(median matched-geometry time ratio ≈0.79; up to ~0.55 on short/urban + Ballyrainey), with a
length-structured error implying a **turn/junction penalty** is needed, not just per-class
speed factors. This is a multi-block research project; only the data-collection tooling exists
so far (no calibrated profile yet).

**⚠️ Paid external API.** Every Google query costs money (Routes API ~$5/1000 requests on a
pay-as-you-go account). **Never run a live Google query without explicit, per-run user
approval** (see Agent Behaviour). Building/editing scripts, `--dry-run`, fully-cached
re-runs, and all OSRM calls are free and need no approval.

**Workflow.**
```
# 1. Build the fixed OD sample (no API/OSRM calls, deterministic, safe to re-run):
python3 analysis/build_od_manifest.py            # → data/google_cache/od_manifest.json (~1000 ODs)
python3 analysis/build_od_manifest.py --n 400    # smaller sample

# 2. Inspect the manifest, get explicit approval, then run the resumable query runner.
#    Cleanest: caller supplies the key so it stays out of the agent's tool calls.
python3 analysis/google_query_routes.py --dry-run            # counts + cost, zero spend
GOOGLE_MAPS_API_KEY=... python3 analysis/google_query_routes.py --limit 100   # batch spend
GOOGLE_MAPS_API_KEY=... python3 analysis/google_query_routes.py               # all remaining
python3 analysis/google_query_routes.py --reprocess-only     # rebuild results.jsonl, no spend
```
Phase A is crash-safe and resumable (each raw response written immediately; re-runs skip
cached ODs; `--limit` batches spend). Phase B (OSRM `/match` over the cache) is free and
idempotent. Time basis is **free-flow** (`TRAFFIC_UNAWARE`) to match the daily-average AADT
model. Outputs live in the gitignored `data/google_cache/` (see Generated/gitignored outputs).
Decisions, feasibility numbers, and design rationale are tracked in the agent memory note
`project_google_routing_calibration`.

**Profile benchmark + compile pipeline (decouples impedance from route preference).**
The calibrated profile is a grid of per-`(road-class × speed-band)` multiplicative speed
factors + global turn costs (`simulation/profile_spec.py`). Benchmarking a candidate profile
must be **fast** (a real OSRM re-extract is ~15-25 min, far too slow per tuning step), so a
one-time *skeleton* pass decouples map-matching from scoring: each Google route is matched once
on the deployed OSRM and its segments labelled by an **exact OSM node-id → tag lookup** (the
raw `build_edge_index.py` cache), yielding a profile-independent skeleton (per-bucket metres +
turn features + signal count). A pure-Python model then re-scores any profile against the whole
cache in milliseconds. (This replaced the earlier **probe profile** approach, whose
`annotation.speed = distance/round(duration,0.1s)` readout corrupted short urban-edge buckets and
poisoned the empirical base speeds — the INT diagnosis in `project_google_routing_calibration`.)
```
# One-time edge index + skeletons (the slow part is ONE cached /match pass):
python3 simulation/build_edge_index.py --match     # /match -> match_cache.jsonl (~1.7s/route; --limit N to batch)
python3 simulation/build_edge_index.py --extract   # osmctools+iterparse -> osm_ways/osm_nodes.jsonl
python3 simulation/build_skeleton_index.py         # match_cache + edge index -> skeletons.jsonl (no OSRM; free to re-run)

# Empirical base speeds (closes the offline↔real gap): one /match per sampled route on a
# factor-free speed source, segments labelled exactly via the edge index, aggregated
# length-weighted (harmonic mean) per bucket. Point --speed-url at a factor-free stock
# instance with --no-defactor, or the deployed :5000 with the default ×factor defactor.
python3 simulation/build_skeleton_index.py --base-speeds  # samples ~800 routes; eval/verify auto-load

# Fast offline benchmark (no OSRM/Docker/spend) — score any candidate profile:
python3 analysis/eval_profile.py                         # stock (all factors 1.0)
python3 analysis/eval_profile.py --legacy-factors        # faithfulness check vs deployed profile
python3 analysis/eval_profile.py --spec simulation/tuned_profile.json

# Tune the bucket factors (external-focused; INT down-weighted to 0 by default,
# its offline turn model under-counts in-town junctions). Linear-in-factor fit,
# scipy, vectorised. → simulation/tuned_profile.json + profile_tuning_history.jsonl
python3 analysis/tune_profile.py                         # default external weights
python3 analysis/tune_profile.py --dry-run               # report without writing
python3 analysis/tune_profile.py --leg-weights X2B=1,B2X=1,X2X=1,INT=0.2 --min-km 100

# Deploy + fidelity gate (once per accepted profile):
python3 simulation/compile_profile.py --spec simulation/tuned_profile.json   # → car_roaaads.lua
#   ... rebuild the DEPLOYED :5000 OSRM with the printed commands ...
python3 analysis/verify_profile.py --spec simulation/tuned_profile.json      # gate before adopting
```
Scoring uses a squared-log-ratio loss, equal weight per valid route (no `1/n_alts`); per-leg
and per-bucket breakdowns are diagnostics only. `verify_profile.py` gates **per-leg median**
`predict/real` within ±`--gate-median-tol` (default 0.03) over `--gate-legs` (default external
`X2B,B2X,X2X`); per-route scatter is inherent (probe-matched skeleton vs deployed re-match) and
is reported but not gated. **INT is reported but not gated/tuned** — the offline turn model
under-counts in-town junctions (verify: offline ≈ 0.72× real on INT), so in-town accuracy waits
on a better turn model. The verify gate is the contract that lets the fast loop be trusted before
any tuned `car_roaaads.lua` is adopted (then re-run the downstream chain:
`build_external_links → reduce_deadends → build_paths → tune_assignment`).

**Empirical-calibration status (2026-06-23):** empirical base speeds make the offline model a
faithful proxy for real OSRM on external corridors (verify per-leg medians ≈ 1.00–1.03); the
external-focused factor tune lands X2B/B2X/X2X medians ≈ 0.97–1.00 vs Google with physically
sensible factors (motorway sped up to ~Google free-flow, urban A/B-roads slowed). Outstanding:
in-town (INT) turn model; tuning the global turn params; the route-preference (stage-2) layer.

---

## Worktree Convention

Background/parallel Claude work is done in `.claude/worktrees/` and
cherry-picked to `main` after review. All work is on `main`.

---

## Agent Behaviour

Even in auto mode, when a user reports a bug or asks a question, agents must not make functional code changes without explicit user approval. Investigate and propose — do not implement unless the user has agreed to the specific change. "Defensive" workarounds (e.g. silent clamps, fallback defaults) are especially suspect: they tend to mask future bugs rather than surface them, and should be proposed and justified before being applied.

**Paid external APIs (Google Maps Routes) — never run without explicit approval.** Every
Google query costs real money on the user's pay-as-you-go account. Agents must **not** make a
live Google API call without explicit, per-run user approval — do not infer standing approval
from an earlier "let's run it" or from a prior approved batch. Building/editing the calibration
scripts, `--dry-run`, fully-cached re-runs, and all (local, free) OSRM calls are fine without
asking. Before any live run, state the planned query count + estimated cost, then stop and wait
for an explicit go. Prefer having the user supply `GOOGLE_MAPS_API_KEY` and/or run the command
themselves so the key never enters agent tool calls. (Mirrored in agent memory
`feedback_no_google_api_without_approval`.)
