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
python3 simulation/build_parking.py          # island-wide OSM parking → data/cache_osm_parking_island.geojson (one-off; osmctools, RAM-light, reuses ni.o5m; feeds retail demand)
python3 simulation/build_schools.py          # island-wide OSM schools + per-POI enrolment → data/cache_osm_schools_island.geojson (one-off; osmctools + school_demand estimator; needs DEA boundary for NI/RoI tagging; feeds school demand)
python3 simulation/build_census_zones.py     # classify NI+RoI census areas → data/census_zones.json incl. per-zone retail_spaces + school_demand (one-off; needs SDZ/DEA boundary files + cache_sa_workplace.csv + island parking/school caches)
python3 simulation/build_network.py          # build road network from local all-island .osm.pbf via osmium (core polygon + 5km bbox; needs Docker)
python3 simulation/build_demographics.py     # node weights + boundary detection + external weights (map is built separately by build_map.py)
python3 simulation/build_external_links.py   # OSRM queries → external↔boundary links + through-route allowlist (needs local OSRM)
python3 simulation/build_intra_times.py      # OSRM-sample intra-zonal times per external zone → data/external_intra_times.json (denominator self-term; needs local OSRM; independent of build_paths — no cache rebuild)
python3 simulation/reduce_deadends.py        # collapse residential dead-ends → newtownards_reduced.graphml + node_weights_reduced.json (consumed by build_paths/build_assignment/tune; see reduce_deadends.py row)
python3 simulation/build_paths.py            # probit stochastic paths incl. external nodes (N_PASSES=25, CV=0.25, N_WORKERS=1; build time depends on hardware)

python3 analysis/parse_official_hourly.py    # parse ODS hourly counts → data/official_hourly.json (one-off)
python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/derive_generation_rates.py                # regenerate generation_rates.json — NTS0409a vehicle-driver trips/person/day + per-purpose ρ_p (re-run when NTS0409 or the purpose mapping changes)
python3 analysis/derive_component_profiles.py              # regenerate hourly_fractions.csv per-component shape columns (reads generation_rates.json ρ_p; run AFTER derive_generation_rates.py)

python3 analysis/tune_assignment.py                        # tune gravity params (8 Tanner shape params P_c/BETA_c, 4 production-constrained scales; external zones fixed from census)
python3 analysis/tune_assignment.py --fast                 # looser tolerances (~2× faster, minimal precision loss)
python3 analysis/tune_assignment.py --note "description"   # optional human label in history

python3 simulation/build_assignment.py       # apply tuned params, write flows
python3 simulation/build_map.py              # build interactive map HTML (run after build_assignment.py to refresh flow layers)

python3 analysis/report_tune.py              # generate reports/ text + pull plot for last run
```

After adding new count data, re-run from `ingest_counts.py` onward. Re-run
`parse_official_hourly.py` if the ODS source file **or** a `model.COUNT_SITES` site
geometry (`node`/`links`) changes — it imports COUNT_SITES as the single source of
truth for site geometry, so `data/official_hourly.json` goes stale until regenerated. The tuner reads
`simulation/tuned_params.json` as its starting point, so repeated runs refine rather
than restart.

---

## Key Files

| File | Role |
|------|------|
| `simulation/zones_config.py` | Single source of truth for the study-area geographic knobs: `CENTRE` (lat/lon), `CORE_RADIUS` (3 km), `SDZ_ZONE_RADIUS` (10 km). Imported by `build_census_zones.py` (uses the radii), `build_network.py` (uses `CENTRE`), and `demographics_config.py` (re-exports `CENTRE`). `CENTRE` is defined here and nowhere else. Editing the radii takes effect only after re-running `build_census_zones.py`. |
| `simulation/build_wz_apportionment.py` | **(RoI data prep — one-off).** Pre-computes the WZ→SA workplace apportionment for all of RoI and writes `data/ireland_data/cache_sa_workplace.csv` (columns: `sa_code`, `workplace_pop`). CSO 2016 Workplace Zone (WZ) boundaries do not align with 2022 SA boundaries; this script intersects them geometrically via `gpd.overlay(wz, sa, how="intersection")`, bypassing 2016→2022 SA boundary change codes entirely (geometry is the ground truth). Each WZ's `T11_C1` headcount ("total workers in workplace zone" — place-of-work jobs; **not** `T1_T`, which is total daytime population) is split across the intersection pieces weighted by the sum of POI weights (`EXCLUDE_AMENITY`/`POI_WEIGHTS` from `demographics_config`) falling within each piece; area-proportional fallback for zero-POI pieces. POIs are extracted from the local PBF via the `osmctools-roaaads` Docker image (`osmfilter --keep-nodes="amenity= shop= office=" --drop-ways --drop-relations`) and cached to `data/ireland_data/cache_roi_pois.geojson`. Will reuse `osrm/edge_index/ni.o5m` if present to skip the slow PBF→o5m conversion step. Re-run only when WZ or SA boundaries change or OSM POI data is significantly stale. **Needs Docker + local PBF.** |
| `simulation/ingest_ni_census.py` | Loads NI DZ/SDZ/DEA boundaries + NISRA population + workplace into standardised GeoDataFrames consumed by `build_census_zones.py`. Public API: `load_ni_census() → (dz_gdf, sdz_gdf, dea_gdf)`. Standardised columns: `area_code`, `parent_code`, `level`, `population`, `workplace_pop`, `geometry` (in `PROJECTED_CRS`). Handles DZ→SDZ parent lookup via column or spatial join fallback; SDZ→DEA similarly. Population fetched from NISRA API (cached to `data/cache_nisra_population.csv`). Workplace from `data/census-2021-apwp001.xlsx`. |
| `simulation/ingest_roi_census.py` | Loads RoI SA/ED/LEA boundaries + CSO 2022 population + pre-computed WZ workplace into standardised GeoDataFrames consumed by `build_census_zones.py`. Public API: `load_roi_census() → (sa_gdf, ed_gdf, lea_gdf)`. Fails loud if `data/ireland_data/cache_sa_workplace.csv` is missing (run `build_wz_apportionment.py` first). ED and LEA GeoDataFrames are derived by dissolving SAs — no separate boundary file needed. Standardised columns match `ingest_ni_census.py`. (The SAPS `T1_1AGETT` *column* sums to ~2× the national population because the file carries a "State" aggregate row equal to the sum of all SAs; the per-SA `SA_PUB2022` join excludes it, so the loaded per-SA population is the correct 1×.) |
| `simulation/census_supply.py` | `load_supply() → {area_code: {commute, school}}` — per-small-area trip *producers* harmonised NI (DZ) + RoI (SA): **commute** = resident workers who physically commute (RoI `SAP2022 T11T1` travel-to-work Total − WFH; NI `distance_to_work` "Work within/outside" bands + "No fixed place", excl. WFH/"No code required"); **school** = student headcount (RoI `T11T1` travel-to-school/college Total; NI `in_full_time_education` "Full-time student or schoolchild"). The producing weights for the school component (now) and commute (business-split). RoI key = JSON-stat SA-code label (the GUID index is ignored); a national "State" aggregate row (~2M) is dropped so per-SA values are clean 1×. Source data (gitignored): `data/ireland_census/` (RoI CSO SAP JSON-stat), `data/ni_census/` (NISRA DZ CSVs). |
| `simulation/build_census_zones.py` | Classifies the full island of Ireland into a three-level census hierarchy centred on `CENTRE` — works for any CENTRE on the island. Calls `load_ni_census()` (NI DZ/SDZ/DEA) and `load_roi_census()` (RoI SA/ED/LEA), concatenates the two hierarchies, then runs unified classification: small areas intersecting `CORE_RADIUS` (3 km) → core; intermediate zones in broken outer zones → SDZ/ED external nodes; non-core small areas in partially-core intermediate zones → orphan DZ/SA external nodes; outer zones outside `SDZ_ZONE_RADIUS` (10 km) → single DEA/LEA centroid nodes. Population-weighted centroids computed from constituent small areas. Outputs `data/census_zones.json`: core polygon (WGS84), external node list with `id` = census-area code (`"N21000219"` for NI, `"017001001"` for RoI SA), `level`, centroid, population, workplace_pop, `retail_spaces`, `school_demand`. NI codes start with `'N'`; RoI codes are pure-numeric. `retail_spaces` = sum of `parking_demand.parking_spaces` over island-parking polygons within each zone (sjoin; workplace-derived fallback for zones with no mapped parking). `school_demand` = sum of per-POI enrolment from the island school cache within each zone (0 for zones with no mapped school). |
| `simulation/build_parking.py` | Builds the island-wide OSM parking cache → `data/cache_osm_parking_island.geojson` (gitignored), the single parking source for `build_census_zones.py` (external zones) and `build_demographics.py` (internal core). Streams parking ways from the all-island pbf via **osmctools** (reuses `ni.o5m`, then `osmfilter --keep="amenity=parking landuse=parking"`), assembles closed-way polygons (RAM-light ~0.5 GB). Saves each polygon with the tags the estimator reads (`access`, `parking`, `building`, `building:levels`, `parking:levels`, `capacity`, `fee`, `amenity`, `landuse`, `name`). **Needs Docker + the pbf/ni.o5m.** |
| `simulation/parking_demand.py` | Pure-stdlib `parking_spaces(tags, area_m2)` → estimated retail parking **spaces** for one OSM parking polygon. Recipe: exclude `access ∈ {private,no,permit}`; decks (`parking ∈ {multi-storey,underground,rooftop}` or `building=parking`) trust `capacity` (else `area×levels/30`), gate-exempt; else `capacity` only if implied `area/capacity ∈ [8,80] m²/space`, else area fallback `÷13` on-street (`street_side`/`lane`) or `÷30` otherwise. Constants in `demographics_config.py`. Destination car parks land at ~29 m²/space in both NI and RoI. Tests: `simulation/test_parking_demand.py`. |
| `simulation/build_schools.py` | Builds the island-wide OSM school cache → `data/cache_osm_schools_island.geojson` (gitignored), the single school source for `build_census_zones.py` and `build_demographics.py`. Streams `amenity=school/college/university/kindergarten` from the pbf via **osmctools** (reuses `ni.o5m`), tags each POI's jurisdiction (NI vs RoI via the DEA boundary), applies `school_demand.assign_enrolments` globally, and saves one point per kept POI with `enrolment`, `amenity`, `name`. **Needs Docker + pbf/ni.o5m + DEA boundary.** |
| `simulation/school_demand.py` | `assign_enrolments(features)` → per-POI school enrolment (pupils), operating on the full island set (cross-feature clustering + institution splitting). **Schools/kindergartens:** light dedup (drop same-name node+way dupes and unnamed sub-buildings within a campus cluster; keep distinct co-located schools); primary/secondary classified by `school=` tag then name, valued by jurisdiction-aware sourced averages (`SCHOOL_ENROLL`: NI 210/820, RoI 170/575), SEN 80, kindergarten 40. **Third-level:** name-matched to a curated table of sourced institution totals (`INSTITUTIONS`, HEA/HESA), split equally across each institution's POIs; unmatched university 300 / college 700; obvious non-teaching (research stations, accommodation, …) dropped. Tests: `simulation/test_school_demand.py`. |
| `simulation/build_network.py` | Builds the road graph from the **local NI `.osm.pbf`** (the same Geofabrik snapshot OSRM is built from — `demographics_config.PBF_PATH`), so road/boundary/internal node IDs share one OSM snapshot with OSRM's route node IDs. The full ~400 MB island pbf OOMs an in-process parse, so a small extract is streamed out with **osmctools** (`osmconvert` + `osmfilter`; Docker image `osmctools-roaaads`, auto-built from `simulation/osmctools.Dockerfile`, ~0.5 GB peak RAM): `osmconvert -b=<bbox> --complete-ways` (bbox = core polygon buffered by `demographics_config.BOUNDARY_BBOX_MARGIN_M` = 5 km) then `osmfilter --keep="highway=<drive set>"` (positive form of osmnx's `drive` filter), written to `simulation/_pbf_drive_extract.osm`. (osmctools is used rather than osmium-tool, whose referenced-node id-set is sized by OSM's max node id and needs several GB regardless of extract area.) `ox.graph_from_xml` reads it; `graph_from_xml` omits the `street_count` node attribute, which `consolidate_intersections` needs, so it is re-added via `ox.stats.count_streets_per_node` (without it the core under-merges). Raw graph extends 5 km beyond the core (for boundary nodes' external neighbours + `build_external_links.py` positions); the consolidated routing graph is still clipped to the core polygon, then junction-consolidated (tol 15 m) and relabelled to OSM IDs. Outputs `newtownards_network.graphml` (raw) + `newtownards_consolidated.graphml`. **Needs Docker + the pbf on disk.** |
| `simulation/build_external_links.py` | Queries a local OSRM instance (all-island extract, **biased `car_roaaads.lua` profile** — see `build_osrm_profile.py`, `http://localhost:5000`) to derive all external zone connectivity. **X→B links:** for each (external node, boundary node) ordered pair, keeps the link only if that boundary node is the first boundary node encountered in the OSRM route (i.e., it is the natural entry point into the core). **B→X links:** symmetric with X→B — keeps B→X only if no other boundary node appears in the OSRM route sequence (i.e., B is the last boundary node departed on the way to X). If another boundary node B' appears, the journey is already covered by B→B' + B'→X. **Boundary→boundary exterior shortcuts:** for each ordered boundary pair, if the route exits the core first, adds a directed shortcut with duration summed from OSRM annotations up to the first boundary node re-encountered. **Through-route allowlist:** for each ordered external-external pair, checks if any OSRM route node is a boundary node; if so adds to `allowed_through_pairs`. Outputs `data/external_links.json`. ~28,000 OSRM queries; under a minute on a local instance. |
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, detects boundary nodes, adds external node weights, writes `node_weights.json` + `newtownards_demographics.geojson`. The map is built separately by `build_map.py` (run it afterwards). `--zones-only` re-reads `data/census_zones.json` and patches only the external node entries in `node_weights.json`. Shared constants (paths, OSM tag handling, map styling) live in `simulation/demographics_config.py`. **Boundary node detection:** loads core polygon from `census_zones.json` and the **pbf-sourced** raw graph (`newtownards_network.graphml` from `build_network.py`), identifies internal nodes (within core polygon), then boundary nodes = internal nodes with at least one edge going outside. Writes `boundary_node_ids`/`internal_node_ids` to `node_weights.json`. Because the raw graph comes from the same OSM snapshot as OSRM, these IDs match OSRM's route node IDs exactly (so `build_external_links.py`'s boundary/internal route-sequence matching is exact). **External node weights:** reads the external node list from `census_zones.json` and writes population + the demand layers (workplace, retail_spaces, commute/school producers, school_demand) to `node_weights.json`. **Study area = core polygon (not a circle):** DZ selection and all OSM downloads (buildings/POIs/parking) are bounded by the core polygon from `census_zones.json` (extent `max_core_vertex_dist_m`, ~10.2 km), matching the road graph built by `build_network.py`. Core DZs are selected by centroid-within the polygon (recovers exactly the `n_core_dzs` core DZs) and use **full** DZ population/workplace_pop (no area-fraction clipping). OSM downloads use a circle sized to the polygon (+1 km margin); POIs and parking are then filtered to within the core polygon to avoid margin leakage (buildings are already DZ-bounded by sjoin). **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Demand layers (separate, never summed):** `node_workplace` (workplace jobs distributed within each DZ by POI count — commute attractor) and `node_retail_spaces` (retail parking spaces via `parking_demand.parking_spaces` on the island parking cache clipped to the core, snapped to road edges — retail attractor) are written as independent layers (no `node_business_demand`); `node_commute_producers`/`node_school_producers` come from census (`census_supply`). External nodes take each per-zone layer from `census_zones.json`. **School demand:** per-POI enrolment from the island school cache → `node_school_demand` (internal POIs snapped to core road edges; external nodes take per-zone `school_demand` from `census_zones.json`). |
| `simulation/build_map.py` | Builds the interactive folium map (`newtownards_map.html`) from artifacts written by `build_demographics.py` (`node_weights.json`, `newtownards_demographics.geojson`), the road graphs, the cached OSM POI/parking layers, and — if present — `newtownards_flows.json`. This was the old `build_demographics.py --map-only` path, now a standalone step (it always reloads POI/parking from cache). Run after `build_demographics.py`, and again after `build_assignment.py` to refresh flow layers. **Flow map layers:** combined AADT (default), plus per-component residential / commute / retail / school layers. No args (`--help` only). |
| `simulation/demographics_config.py` | Shared pure-constant config imported by `build_demographics.py`, `build_map.py` **and `build_network.py`** (file paths, OSM tag handling — `EXCLUDE_AMENITY`/`POI_WEIGHTS` — and map styling). `CENTRE` is re-exported from `zones_config.py` (not defined here). Also holds the road-network source knobs: **`PBF_PATH`** (absolute path to the NI `.osm.pbf` OSRM is built from) and **`BOUNDARY_BBOX_MARGIN_M`** (5 km buffer around the core polygon for `build_network.py`'s extract). The separate `NETWORK_MARGIN_M` (1 km) here sizes the OSM POI/building/parking download circle only — unrelated to the road graph. Also defines **`PROJECTED_CRS = "EPSG:2157"`** (Irish Transverse Mercator / ITM) — the single source of truth for all projected spatial operations in the pipeline. ITM covers the whole island of Ireland with uniform accuracy, avoiding UTM Zone 30N's distortion for Republic of Ireland towns west of ~6°W. All simulation and analysis scripts import this constant; `build_network.py` passes it explicitly to `ox.project_graph`. **Do not hardcode `EPSG:32630` anywhere.** Single source of truth so the split scripts don't drift. Also holds the **parking→retail-spaces estimator constants** (`PARKING_M2_PER_SPACE_OFFSTREET=30`/`_ONSTREET=13`, `PARKING_GATE_LO=8`/`PARKING_GATE_HI=80`, `PARKING_EXCLUDE_ACCESS`, `PARKING_DECK_TYPES`, `PARKING_ONSTREET_TYPES`) consumed by `parking_demand.py`, plus `PARKING_ISLAND_CACHE` and `SCHOOL_ISLAND_CACHE` (island parking/school cache paths). |
| `simulation/build_intra_times.py` | OSRM-samples intra-zonal travel times per external census zone for the production-suppression **self-term**. For each external node, recovers its census polygon (NI: SDZ/DZ/DEA from boundary files; RoI: SA directly, ED/LEA dissolved from SA boundaries), rejection-samples `M`=30 uniform point-pairs inside it, and routes each on the local OSRM (`localhost:5000`, same profile as `build_external_links.py`) → `data/external_intra_times.json` (`{census_code: [t1..tM seconds]}` + `_meta`). Loud on any zone with a missing polygon or `<M` successful routes. `--m N` overrides the pair count. **Run after `build_census_zones.py`, OSRM up. Independent of `build_paths.py`** — the self-term lives in the model layer, so re-running needs no paths-cache rebuild. Re-run only when external zones change. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Covers both internal road nodes and external census-area nodes. **Graph augmentation:** loads external nodes and edges from `data/external_links.json`; adds them to the routing graph before Dijkstra. External edges (X↔B, boundary shortcuts) are included in the adjacency matrix but NOT in `link_list` — they contribute to path distance but not to flow accumulation. **Probit loading:** all edges (road and external) perturbed each pass with log-normal noise `exp(eps·w)`, `eps ~ N(0, CV=0.25)`, N_PASSES=25, giving stochastic spread in boundary node selection for external→internal OD pairs with similarly-weighted entry options. **Length-scaled noise (`PROBIT_LL_SIGMA`, default 120 s ≈ 2 min):** the per-edge gain `w = σ_ll/(σ_ll + CV·cost)` ∈ (0,1] keeps the noise multiplicative for short legs (`w→1`) but saturates it to a fixed *absolute* sigma `σ_ll` for long legs (`w→σ_ll/(CV·cost)`), so a long single-edge external↔boundary leg's perturbation no longer swamps the few-minute differences between competing boundary entries. The adjusted perturbation never exceeds the pure multiplicative one, and `eps=0 ⇒ no bias`. Crossover at cost ≈ σ_ll/CV (~8 min). `PROBIT_CV` and `PROBIT_LL_SIGMA` are imported from `simulation/routing_config.py` (the gain vector is precomputed once, constant across passes). **OD pair filter:** through-routed external→external pairs (in `allowed_through_pairs`) are routed flow pairs (Dijkstra path through the core). Non-through external→external pairs (from `external_external_times`) are appended as **denominator-only** pairs — entries in `od_src/od_dst/od_dist` (distance = direct OSRM time) but NOT in `pair_idx/link_idx` and excluded from `src_groups`/probit passes, so they carry no flow; they complete each external origin's production-constrained denominator. The cache stamps `n_routed_pairs` (flow-carrying pairs occupy `0..n_routed_pairs-1`). No offscreen leg calculation. **Internal edge costs:** each road edge's routing cost is the OSRM-equivalent travel time `factor(class,band)·length·3.6/base_speed` from the Google-calibrated profile (`simulation/tuned_profile.json` + empirical `base_speeds.json`), via `simulation/edge_speed.py`. `add_edge_speeds`/`add_edge_travel_times` are still run because `deadend_collapsed` synthetic edges keep their osmnx-encoded intra-region time (factor 1.0). Internal routes are chosen on realistic time alone (no route-preference biasing); the node-based Dijkstra applies no turn penalties. Re-run if road network, external links, the tuned profile or base speeds, `N_PASSES`, `PROBIT_CV`, or `PROBIT_LL_SIGMA` change. **Performance constants:** `N_WORKERS` (default 1) controls parallel pass workers via `multiprocessing.Pool` — increase on machines with sufficient RAM (each worker uses ~100–150 MB extra); `MAX_HOPS` (default 120) caps per-pair path-trace iterations. Inner path-tracing loop is vectorised (numpy). |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (legacy unconstrained rational kernel, kept for old param files), **`constrained_od_flows()`** (production-constrained per-component per-pair pre-K flows + per-origin denominators; optional `self_src`/`self_dist`/`self_w` add the **external intra-zonal self-term** to each denominator — denominator-only, `None` ⇒ exact prior behaviour) and **`scatter_od_to_links()`** (the production-constrained assignment core, used by `build_assignment.py` and `tune_assignment.py`), **`load_self_terms(node_ids)`** (builds the self-term arrays from `data/external_intra_times.json`; emits one entry per sampled time with weight `1/M_i`; skips zones absent from `node_ids`; returns `(None,None,None)` if the file is missing), **`load_generation_rates()`** + **`compute_generation_scales(node_weights, rates)`** (the per-leg producer scales that pin generation to vehicle-driver trips/day — island per-capita anchors from the node-weight layer sums; feed `constrained_od_flows(..., gen_scale=…)`; see "Generation pinning"), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `gravity_assign()` accepts optional `link_weight` array (probit fractional weights). `compute_chi2()` takes the commute/retail/school link-flow dicts + per-component `slot_fracs_*`; in four-component mode (the production-constrained components) N_eff = N (no per-slot df subtracted); legacy single-flow path otherwise. Road node IDs are OSM integers (stable); external census node IDs are census-area-code strings (e.g. `"N21000219"`) — not OSM IDs. COUNT_SITES: site 507 links 538692601↔549139252; site 508 node 136173611; site 444 node 449111329. **`WEIGHTS_FILE` and `ROUTING_GRAPH` point at the dead-end-reduced artifacts** (`node_weights_reduced.json`, `newtownards_reduced.graphml`) from `reduce_deadends.py`. EXCLUDE_LINKS: `{(181844513, 181839481)}` plus the Westmount Park and Old Belfast Road directed links (both directions) whose endpoints are absorbed by `reduce_deadends.py` and no longer exist in the reduced graph — their walking observations are discarded from calibration (regenerate this set from `deadend_broken_obs.json` if the reduction params change). |
| `simulation/build_assignment.py` | **Production-constrained** gravity assignment (via `model.constrained_od_flows` + `scatter_od_to_links`). Requires `simulation/newtownards_paths.npz`. Four-component mode activated when `K_res`/`K_commute`/`K_retail` + `P_commute`/`P_retail` are present in `tuned_params.json` (school sub-component when `K_sch > 0` and `P_school`/`node_school_demand` present). Applies the per-leg `gen_scale` (`model.compute_generation_scales`) so component magnitudes are vehicle-driver trips/day. Saves `flows_res`, `flows_commute`, `flows_retail`, `flows_school` in `newtownards_flows.json`. Legacy single-K unconstrained `gravity_assign` path kept for old pre-split param files. External node weights come from `node_weights.json` directly (no override from tuned params). |
| `simulation/edge_speed.py` | OSRM-equivalent internal edge-time model shared by `build_paths.py` and `reduce_deadends.py` (paths used in `model.paths_cache_signature`). `load_profile()` loads the tuned `ProfileSpec` (`simulation/tuned_profile.json`) + empirical base speeds (`data/google_cache/base_speeds.json`) — **fails loud** if the tuned profile is missing (warns if base speeds absent → analytical fallback). `edge_time_seconds(tags, length_m, spec) = factor(class,band)·length·3.6/base_speed(class,band)` reuses `profile_spec` for all bucketisation/base-speed/factor logic (same `(class×band)` buckets the deployed `car_roaaads.lua` keys on, so internal routing matches the tuned OSRM instance). Robust to `highway`/`maxspeed` stored as a list. Edge impedance only — no turn penalties (the internal Dijkstra is node-based). |
| `simulation/routing_config.py` | Holds `HIGHWAY_COST_FACTOR` + the probit noise params (`PROBIT_CV`/`PROBIT_LL_SIGMA`). `HIGHWAY_COST_FACTOR` does not drive internal routing (`build_paths.py`) or the dead-end reducer (`reduce_deadends.py`) — those use the calibrated `(class×band)` profile via `simulation/edge_speed.py`. It is used only by the legacy tooling that references it: `build_osrm_profile.py`, `build_skeleton_index.py --base-speeds`, and `skeleton_model.legacy_spec_from_highway_cost_factor`. |
| `simulation/build_osrm_profile.py` | Generates `car_roaaads.lua` — the road-class-biased OSRM car profile. Pulls the default `car.lua` from the `osrm/osrm-backend` Docker image, injects a block after the `forward_rate` assignment that divides `forward_speed`/`forward_rate` by `HIGHWAY_COST_FACTOR` (matching internal Dijkstra biasing). Re-run whenever `HIGHWAY_COST_FACTOR` changes, then re-preprocess OSRM (`osrm-extract -p car_roaaads.lua`, `osrm-partition`, `osrm-customize`). Output: `/home/matthew/Documents/CodingFun/osrm/car_roaaads.lua`. |
| `simulation/reduce_deadends.py` | Collapses "residential dead-end" regions in the consolidated routing graph to shrink node count (speeds up `build_paths.py`/tuning, enables a larger core area). A region R (entrance E ∉ R) qualifies iff: (1) R connects to the rest of the network through exactly one cut vertex E; (2) R contains no boundary node and no school-demand node (both *protected* — never absorbed — which enforces the no-boundary and zero-school rules structurally); (3) max directed journey time E→n over n∈R < `T_MAX` (default 60 routing-cost seconds); (4) total workplace+retail demand < `BIZ_CAP` (default 100; residential pop unbounded); (5) `|R| ≥ 2` (single-node spurs skipped — 1→1 saves nothing). **Algorithm:** every valid region is a protected-free connected component of H−a (H = undirected simple projection) for some articulation point a, so it enumerates all such (entrance, region) candidates, filters by constraints 2–5 + directed reachability both ways, and selects the *maximal feasible* regions (laminar family ⇒ disjoint; naturally descends into an oversized branch to find the largest collapsible sub-pockets — catches cyclic closes that leaf-pruning would miss). Each region → one super-node S (=min id, summed pop/workplace/retail/school, pop-weighted UTM centroid) joined to E by directed links E→S, S→E whose travel times are population-weighted means of the intra-region directed times. **Intra-region times use the same OSRM-equivalent `(class×band)` edge model as `build_paths.py`** (`build_cost_digraph(G, spec)` via `simulation/edge_speed.py`), so the collapse/`T_MAX` decisions match the speeds the reduced graph is later routed on. Synthetic edges use `highway="deadend_collapsed"` (factor 1.0 in `build_paths.py`) with `maxspeed`+`length` set so osmnx's `add_edge_speeds`/`add_edge_travel_times` (re-run by `build_paths.py`) reproduce that target time. **Run after `build_demographics.py` (needs pop/biz/school + boundary) and before `build_paths.py`.** Outputs (gitignored): `newtownards_reduced.graphml`, `node_weights_reduced.json`, `deadend_map.json` (provenance: super-node→absorbed nodes + times), `deadend_broken_obs.json` (observed/count links whose endpoints were eaten — **manual review before adoption**; observed-link endpoints are deliberately *not* protected). Params: `--t-max`, `--biz-cap`. **Wired into the pipeline:** `build_paths.py` (`CONS_GRAPH`), `build_assignment.py` (`CONS_GRAPH`), `tune_assignment.py` (`CONS_GRAPH`) read `newtownards_reduced.graphml`, and `model.WEIGHTS_FILE`/`ROUTING_GRAPH` point at the reduced files — so this step must run after `build_demographics.py`. The 6 absorbed walking observations (Westmount Park, Old Belfast Road) are discarded via `EXCLUDE_LINKS` in `model.py`. **Map caveat:** `build_map.py` still draws the *full* consolidated graph, so flow on collapsed interior streets is not shown on the map (demand layers are unaffected; main-road flows and the fit are unaffected). Re-mapping collapsed regions via their super-nodes is a possible follow-up. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Gravity param regularization and scale-share priors. `gravity_lambda` + `gravity_ref` regularise the 8 Tanner shape params (per-component peak+shape: P/BETA/P_commute/BETA_commute/P_retail/BETA_retail/P_school/BETA_school); `phi_commute_prior`/`phi_commute_std`, `phi_retail_prior`/`phi_retail_std`, `phi_school_prior`/`phi_school_std` set Gaussian priors on the commute/retail/school **scale shares** (`φ = K_c/ΣK`) inside `solve_scales` (degeneracy break). `gamma_coupling_scale` and `lambda` are unused dead keys left in the file. External demand layers (workplace, retail_spaces, school_demand, commute/school producers) are measured per zone in `census_zones.json` — no external scale factors. |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. **Imports `model.COUNT_SITES` as the single source of truth for site geometry** (`SITE_MAP`) — it stamps each site's `node`/`links` from COUNT_SITES into the output. **Re-run when the ODS file OR a COUNT_SITES site location changes** (otherwise `official_hourly.json`, which the tuner reads, drifts stale from COUNT_SITES). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. Loads manual link overrides from `data/manual_link_overrides.json`. After every new link assignment, validates each non-null count direction against the directed graph; raises `ValueError` if the edge doesn't exist. |
| `analysis/manual_assign_link.py` | CLI tool to manually assign a session to a specific directed link, bypassing GPS snap. Usage: `python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>`. Validates both nodes exist and checks count-edge consistency. Writes to `data/manual_link_overrides.json` and patches `counts_processed.json` directly. After correcting an assignment, re-run `aggregate_counts.py` then `tune_assignment.py`. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s`. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Four-component, production-constrained model:** gravity flows split into residential (`flow_res`), commute (`flow_commute`), retail (`flow_retail`), and school (`flow_school`) components, each singly (production) constrained. Tunes **8 gravity params** — per-component Tanner kernel `f(u)=u^BETA·exp(BETA·(1−u))`, u=d/P_c (peak P_c + shape BETA_c each): P, BETA, P_commute, BETA_commute, P_retail, BETA_retail, P_school, BETA_school — no `W_BIZ`/`W_SCHOOL`, no shared BETA, no ALPHA. External zone values are fixed from census data and are not tuned. Producer weights are scaled to vehicle-driver trips/day via `model.compute_generation_scales` (generation pinning ⇒ each `K_c ≈ 1`). **Inner calibration = direct-K convex scale solve (`solve_scales`):** the temporal fractions `f_res/f_commute/f_retail/f_school` are **pinned at the NTS profile** (never tuned), so with `f` fixed each prediction is linear in `(K_res, K_commute, K_retail, K_sch)` and the inner objective (Gaussian WLS + Poisson identity-link deviance + scale-share prior) is **convex**, solved by a damped-Newton + line-search step — **monotone, no K-collapse, no best-iterate hack**. `run_assignment` calls `model.constrained_od_flows` and scatters via the probit routing incidence. **Observed-link scatter restriction (tuner-only):** the objective reads modelled flow on only the ~230 observed links, so `run_assignment` scatters just the incidence entries landing on those links (≈32% of the ~62M), precomputed once into a compact observed-link space — bit-identical results, ~3× faster per eval (`build_assignment.py` keeps the full scatter for the map). Once the compact arrays are built the tuner frees the full-incidence cache arrays (`pair_idx`/`link_idx`/`link_weight`) to keep steady-state memory low. **Performance:** ~1.85 s/eval (run_assignment-dominated; more under memory pressure). A full Powell run is a few thousand evals (e.g. run `868d9604` = 4217 evals ⇒ **~2 hours**), so it is a heavy, long-running pass — not a quick verify. `CALIBRATE_PROBE=1` is an env-gated diagnostic that reports the post-calibrate residual global scale λ. |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` and `reports/slot_pulls_{id}.png`. Echoes the labels the tuner stored in history, so street names appear only after a fresh tune run regenerates `tuning_history.jsonl`. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. |
| `simulation/reset_gravity_params.py` | Reset the gravity params in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`: every `gravity_ref` shape param (Tanner: P, BETA, P_commute, BETA_commute, P_retail, BETA_retail, THETA, P_school, BETA_school — iterates `gravity_ref`, so no rename drift) plus the four scales `K_res/K_commute/K_retail/K_sch` → 1.0. Strips dead keys (`K`, `K_biz`, `W_BIZ`, `W_SCHOOL`, `P_biz`, `ALPHA_biz`, the rational-kernel tail exponents `ALPHA`/`ALPHA_commute`/`ALPHA_retail`/`ALPHA_school`, `MU`, `SIGMA`). External params and `slot_fracs_*` are preserved. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** Per-component temporal-**shape** profiles (168 rows = 7 days × 24 h): `mean_fraction_res`, `mean_fraction_commute`, `mean_fraction_retail`, `mean_fraction_school`, plus the aggregate `mean_fraction` (input). **Derived from NTS** via `analysis/derive_component_profiles.py`. Each component column is an **independent shape** normalised so its day-weighted daily sum `W_c = 1` (⇒ each column sums to 7.0 over the 168 rows) — **no** aggregate-partition constraint (magnitude/split is generation's job, see "Generation pinning"). Re-run `derive_component_profiles.py` when the NTS files, the purpose mapping, or `generation_rates.json` change. |
| `analysis/derive_component_profiles.py` | Derives the per-component hourly-**shape** columns of `hourly_fractions.csv`. Each component's weekly profile `f_c(dow,h) = V_c(dow)·H_c(daytype,h)`: weekday hourly shape `H_c` = ρ-weighted blend (ρ_p from `generation_rates.json`) of its purposes' NTS0502a hourly distributions; day-of-week volume `V_c` from NTS0504b (Mon–Sun × purpose), normalised so `Σ_7 V_c = 7` (⇒ `W_c = 1`); weekend hourly shape from the aggregate (no per-purpose weekend data), shared. Purpose→component map shared with `derive_generation_rates.py` via `analysis/purpose_mapping.py`; school uses the **escort-education** shape (the car school-run). Re-run when NTS files / the mapping / generation rates change. |
| `analysis/purpose_mapping.py` | **NEW.** Single source of truth for the NTS purpose→gravity-component mapping (canonical purposes + weights + `LEISURE_RETAIL_FRAC`), imported by **both** `derive_generation_rates.py` (magnitudes, NTS0409a) and `derive_component_profiles.py` (shapes, NTS0502a/0504b) so generation and temporal cannot drift. Carries the judgment allocations (Business/Personal business→retail; Leisure ½/½; Education/escort→school). |
| `analysis/derive_generation_rates.py` | **Tracked output `analysis/generation_rates.json`.** Derives per-component vehicle-driver trips/person/day from England **NTS0409a** (`data/nts0409.ods`, 2023/24 avg, modes *Car or van driver + Motorcycle + Taxi or minicab*) → the generation-pinning rates `ρ_c` (commute 0.178, retail 0.495, school 0.082, res 0.262). Encodes the **purpose→component map** and `LEISURE_RETAIL_FRAC = 0.5` as named, commented module constants (the judgment allocations — see "Generation pinning"). Re-run when the NTS file or the mapping changes. Consumed via `model.load_generation_rates` / `model.compute_generation_scales`. |
| `analysis/google_routing_common.py` | **(Google calibration).** Shared pure-stdlib helpers for the Google routing-time calibration tooling: encoded-polyline decode, polyline downsampling, OSRM `/route` and `/match` calls, Google Routes API v2 `computeRoutes` call, `CONF_MIN` (0.5 /match-confidence floor). No third-party deps beyond `networkx` (used only by the manifest builder). See **Google Routing-Time Calibration** section below. |
| `analysis/google_feasibility.py` | **(Google calibration — pilot).** One-shot feasibility experiment: a small hardcoded data-driven OD sample → Google Routes, decode polyline, OSRM `/match` (time error on Google's geometry) vs OSRM `/route` (route-choice divergence). Matches best route + all alternatives; caches raw responses in `data/google_cache/`. `--dry-run` makes no API calls; key only needed on a live cache-miss. Superseded for production sampling by `build_od_manifest.py` + `google_query_routes.py`, kept as the reference pilot. |
| `analysis/build_od_manifest.py` | **(Google calibration — v1 batch).** Writes the fixed, deterministic (seed 20260622), model-aligned, length-skewed OD sample → `data/google_cache/od_manifest.json`. **Makes NO API/OSRM calls** — reads only `census_zones.json`, `external_links.json`, `node_weights.json`, and the raw graph. Leg types + default quotas (of `--n`, default 1000): X2B 45% (external centroid→boundary entry), X2X 25% (allowlisted through-routes), B2X 15% (boundary→external), INT 15% (internal→internal, for in-town junction realism). Within each leg type, 4 length quartile-bands allocated with a long-skew (`BAND_WEIGHTS` 0.15/0.20/0.30/0.35); X2B/B2X length = model `duration_s`, X2X/INT length = haversine. See also `build_od_manifest_v2.py` for the second batch. |
| `analysis/build_od_manifest_v2.py` | **(Google calibration — v2 batch).** Builds a second fixed 1000-OD sample (seed 20260626) with zero `(origin_label, dest_label)` overlap with v1 → `data/google_cache/od_manifest_v2.json`. All od\_ids prefixed `v2_` to avoid filename collisions in `raw/`. Leg-type distribution re-weighted toward B2X (25% vs 15% in v1, observed higher violation rate). Both v1 and v2 append to the shared `match_cache.jsonl` — `google_query_routes.py --manifest od_manifest_v2.json` followed by `build_skeleton_index.py` picks up the combined set. Combined (v1+v2): ~3932 routed routes, 4471 skeletons. |
| `analysis/google_query_routes.py` | **(Google calibration — runner).** Crash-safe, resumable runner over a manifest. `--manifest` overrides the default `od_manifest.json` (pass `od_manifest_v2.json` for the second batch). **Phase A (spendy, resumable):** queries each uncached OD and writes its raw Google response to `data/google_cache/raw/<od_id>.json` **atomically and immediately** (temp+rename); re-runs skip cached ODs; `--limit N` caps queries per run. **Phase B (information-greedy, single `/match` pass per route):** rebuilds `data/google_cache/results.jsonl` AND simultaneously appends full match detail (node sequence + per-segment distances + maneuvers + `match_dur`) to `match_cache.jsonl` — one `osrm_match_detail` call per route, never two; routes already in `match_cache.jsonl` read from cache with no second call. Also runs OSRM `/route` per OD (free). `--reprocess-only` skips Phase A; safe to re-run any time. **For new batches fully processed via Phase B, `build_edge_index.py --match` is not needed** — Phase B writes the identical `match_cache.jsonl` format. **Refuses to start Phase A without `GOOGLE_MAPS_API_KEY`, and a live run requires explicit per-run user approval (see Agent Behaviour).** |
| `simulation/profile_spec.py` | **(profile calibration — single source of truth).** Pure-stdlib definition of a calibrated OSRM time profile: a grid of multiplicative speed **factors** per `(highway_class × speed_band)` bucket (full `DRIVE_HIGHWAYS` classification × NI mph bands `{untagged,20,30,40,50,60,70,other}`) + the four global turn params (`turn_penalty`, `traffic_light_penalty`, `u_turn_penalty`, `turn_bias`). `factor=1.0` = stock-OSRM base speed; `factor>1 ⇒ slower` (OSRM is currently too fast). Holds the bucketisation (`norm_class`/`parse_band`/`band_from_tags`/`bucket_of`/`bucket_index`), the stock base-speed table, `base_speed_for`, and `ProfileSpec` (JSON load/save). **maxspeed resolution mirrors OSRM's `WayHandlers.maxspeed` exactly:** `bucket_of(tags)` takes a way's full tag dict and `band_from_tags` honours OSRM's key precedence (`maxspeed:advisory` > `maxspeed` > `source:maxspeed` > `maxspeed:type`); `osrm_maxspeed_kmh` resolves numeric *and* symbolic/national-speed-limit values (`gb:nsl_single`→60, `gb:nsl_dual`/`gb:motorway`→70, `none`→140, plus `maxspeed_table_default` urban/rural/trunk/motorway) — so nsl roads land in their real speed band instead of collapsing into `untagged`. **Replaces `routing_config.py`/`HIGHWAY_COST_FACTOR` for the calibration work** (the old module still feeds `build_paths.py`/`build_osrm_profile.py` until a calibrated profile is adopted). Imported by both the stdlib `analysis/` tooling and the simulation-side Lua generators, so the offline model and the emitted Lua key on the *same* buckets. |
| `simulation/osrm_lua.py` | **(profile calibration).** Shared OSRM car.lua/Docker plumbing: `pull_base_lua`, `copy_lib`, the 3-strategy injection-point `find_injection_point`/`inject` (refactored out of `build_osrm_profile.py`), and the Lua emitters `emit_probe_block` (legacy probe — now unused) + `emit_factor_block(spec, pref_dict=None)` (divides `forward_speed` by the tuned per-bucket `_FAC` factor, then optionally divides `forward_rate`/`backward_rate` by a per-class `_PREF` preference multiplier — both in one `do...end` block sharing the highway-tag lookup) + `apply_turn_overrides`. The `_PREF` block resolves link classes to their parent via a `_LPAR` table, then splits trunk/primary/secondary/tertiary into urban (`≤30 mph`) and rural (`>30 mph`) sub-keys using the same speed-source logic as the Python `_pref_key` function (tagged maxspeed → OSRM class default for untagged). The bucket-index Lua replicates the full nsl-aware band resolution, cross-checked against Python `bucket_of`. |
| `simulation/build_edge_index.py` | **(profile calibration — raw OSM cache; replaces the probe).** `--match`: the **single** `/match` pass over the route set on the deployed OSRM (:5000) — caches full match detail per route (node sequence + per-segment `distance` + step maneuvers + `match_dur`) to `data/google_cache/match_cache.jsonl` (resumable; `--limit N` batches the slow ~1.7 s/match pass; `--manifest` processes a second batch into the same shared cache). **Use only for the initial v1 batch or re-matching an existing manifest**; for new batches queried via `google_query_routes.py`, Phase B already appends to `match_cache.jsonl` in the same format. `--extract`: streams the NI pbf via `osmctools-roaaads` → `osm_ways.jsonl` / `osm_nodes.jsonl` (complete raw tag dict for every way/node in the route set). |
| `simulation/build_skeleton_index.py` | **(profile calibration — skeleton builder).** Rebuilds the profile-independent `data/google_cache/skeletons.jsonl` from `match_cache.jsonl` + the edge index — **no OSRM calls**, pure recompute, so it is free to re-run after any `profile_spec`/bucket change. Each matched segment `(node_u,node_v)` is resolved to its way's tags via the edge index and bucketed with `profile_spec.bucket_of` (exact node-id lookup, not the old probe `annotation.speed` readout that corrupted short urban edges) → `length_by_bucket`; `turns` from cached step maneuvers, `n_signals` from cached node tags (`highway=traffic_signals`), `coverage`/`valid` from geometry. `--base-speeds`: the one remaining `/match` step — samples ~800 routes on a factor-free speed source (stock OSRM with `--no-defactor`, or `:5000` defactored), labels each segment via the edge index, and writes **length-weighted harmonic-mean** per-bucket base speeds → `base_speeds.json`. |
| `analysis/skeleton_model.py` | **(profile calibration — fast offline model).** Pure-stdlib `predict_duration(skel, spec)` = `Σ_bucket factor·length·3.6/base_speed` (edge) + OSRM-style turn sigmoid (gated on degree>2 / u-turn, NI left-hand bias) + `n_signals·traffic_light_penalty`. `evaluate(skeletons, spec)` scores against Google with squared-log-ratio loss, **equal weight per valid route** (conf≥0.5 & coverage in-band), and returns per-leg-type / per-(leg×band) diagnostics. `bucket_coverage` (factor identifiability) and `legacy_spec_from_highway_cost_factor` (deployed-profile reference) included. Milliseconds for the whole cache — no OSRM/Docker. |
| `analysis/eval_profile.py` | **(profile calibration — benchmark entrypoint).** Scores a `ProfileSpec` JSON (default all-1.0 stock; `--spec`; `--legacy-factors`) against `skeletons.jsonl`: aggregate loss, predicted/Google ratio distribution, per-leg-type + per-cell breakdown, per-bucket coverage table, turn-time fraction. No spend. `--legacy-factors` is the faithfulness sanity check (should track the deployed `te_matched`). |
| `simulation/compile_profile.py` | **(profile calibration — compiler).** `tuned_profile.json` (a `ProfileSpec`) → deployable `car_roaaads.lua`: applies tuned turn params in `setup()`, injects `emit_factor_block(spec, pref_dict)`, copies `lib/`, prints the re-extract/partition/customize commands. `--pref simulation/tuned_preference.json` (default: auto-loads if the file exists) injects the `_PREF` preference block alongside the `_FAC` timing block; `--no-pref` skips it. **The `_FAC` (timing) and `_PREF` (preference) tables are separate** — `_FAC` divides `forward_speed` (→ duration); `_PREF` divides `forward_rate`/`backward_rate` only (→ routing cost). Timing is never changed by preference factors. |
| `analysis/tune_preference.py` | **(route-preference calibration).** Fits per-highway-class preference multipliers `p_c` (applied to `forward_rate` only, not `forward_speed`) so OSRM routes toward Google's preferred road hierarchy. Uses a **scale-invariant log-ratio ranking loss**: `max(0, log(cost(r0)/cost(rk)) + log_margin)²` over true preference violation pairs (both Google and offline say `rk` is faster, but Google chose `r0`). 13 classes: motorway, trunk/trunk_rural, primary/primary_rural, secondary/secondary_rural, tertiary/tertiary_rural (urban/rural split at 30 mph — tagged maxspeed or OSRM class default for untagged), plus unclassified/residential/living_street/service. `p_c < 1` = preferred; `p_c > 1` = avoided. L2 reg (lam) toward `p_c=1`; bounds `[0.33, 3.0]`; scipy L-BFGS-B. Writes `simulation/tuned_preference.json` + appends `preference_tuning_history.jsonl`. **Status: deferred — class-only factors cannot achieve a net ranking improvement** with the current data (67 violations vs 1033 external concordant pairs; net is negative at every lam 0.001–0.5). |
| `analysis/eval_preference.py` | **(route-preference benchmark).** Scores a `tuned_preference.json` offline against all multi-route skeletons: correctly-ranked true violations (calibration target), timing divergences (OSRM already routes r0 on timing alone — no preference fix needed), concordant regressions (pairs flipped wrong by the preference factors, broken down by leg type and timing-error category). `--unit` scores the `p_c=1` baseline. No OSRM/Google calls. Key insight from combined v1+v2 run: 67 violations vs 1165 concordant pairs (1033 external); flips 4–5× more concordant pairs than violations resolved at every lam, confirming class-level granularity is too coarse for deployment. |
| `analysis/verify_profile.py` | **(profile calibration — fidelity gate).** After the deployed OSRM is rebuilt from a compiled profile (or against the live deployed instance with `--legacy-factors`), `/match`es a validation subset through real OSRM (:5000) and compares real `match_dur` to `predict_duration(skel, spec)`. **Gate (median-based):** per-leg median `predict/real` within ±`--gate-median-tol` (default 0.03) over `--gate-legs` (default external `X2B,B2X,X2X`); per-route scatter (med/p90 |resid|) is reported but not gated (inherent to probe-vs-deployed re-match). INT reported, not gated. Exits non-zero on fail. Read-only, no Google calls. |
| `analysis/tune_profile.py` | **(profile calibration — external-focused tuner).** Fits per-`(class×band)` speed **factors** to minimise the weighted squared-log-ratio time error vs Google over `skeletons.jsonl`, leg-weighted (default `X2B/B2X/X2X=1, INT=0`). With turn params + base speeds fixed, predicted time is **linear in the factor vector**, so it's a vectorised (numpy) scipy `L-BFGS-B` fit with L2 reg toward 1.0 and `[0.2,5]` bounds; only buckets above `--min-km` weighted coverage are tuned (rest stay 1.0). Writes `simulation/tuned_profile.json` + appends `profile_tuning_history.jsonl`; reports before/after per-leg medians + top factor moves. **Tunes factors only** — global turn params are held at defaults (external turn fraction is small; INT excluded) until the in-town turn model is improved. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_map.html`, `simulation/tuned_params.json` — all regenerated by the pipeline.
`simulation/node_weights.json` keys: `node_population`, `node_workplace` (place-of-work jobs = commute attractor, internal + external), `node_retail_spaces` (estimated retail parking spaces = retail attractor, internal + external), `node_school_demand` (school places = school attractor), `node_commute_producers` / `node_school_producers` (census resident commuters / students = trip producers, internal + external; from `census_supply.py`), `boundary_node_ids` (auto-detected from core polygon). There is no combined `node_business_demand` layer — the commute (`node_workplace`) and retail (`node_retail_spaces`) attractors are separate. External node entries (census-area-code string IDs, e.g. `"N21000219"`) are included alongside internal OSM node IDs.
**Node ID scheme:** `build_network.py` relabels consolidated graph nodes to stable OSM IDs after junction consolidation. A consolidated junction (multiple OSM nodes merged) gets `min(osmid_original)` as its ID; a non-merged node gets `int(osmid_original)`. All road node IDs are therefore genuine OSM node IDs (in the hundreds of millions) and stable across graph regenerations. External census nodes use their **census-area-code string IDs** (SDZ/DEA 2021 codes, e.g. `"N21000219"`) — these are the `id` values in `census_zones.json`, *not* small integers (downstream code that consumes `external_links.json`/`census_zones.json` must treat external IDs as strings — `build_paths.py` does). Road node IDs are ints, external node IDs are strings; `node_to_idx` mixes both.
`simulation/newtownards_flows.json` — combined flows plus `flows_res`/`flows_commute`/`flows_retail`/`flows_school` keys (W-weighted directed AADT) when four-component params active, plus an `aadt_weights` block.
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
offline model, overrides the analytical estimate), `profile_tuning_history.jsonl` (one line per
`tune_profile.py` run), `od_manifest_v2.json` (second batch, seed 20260626; `raw/v2_*.json` for
its responses), and `preference_tuning_history.jsonl` (one line per `tune_preference.py` run).
Survives worktree removal (lives in the main checkout); only at risk from `git clean -xfd` or
manual `rm`. `simulation/tuned_profile.json` (a candidate `ProfileSpec`, gitignored) and
`simulation/tuned_preference.json` (a candidate preference dict, gitignored) are also generated.

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`, `data/official_hourly.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`data/census_zones.json` — committed; output of `build_census_zones.py`. Contains core polygon, external node list with IDs/codes/centroids/census demand. Re-run `build_census_zones.py` only if NISRA boundary files or census data change.
`data/external_links.json` — committed; output of `build_external_links.py`. Contains OSRM-derived X↔B links, boundary shortcuts, and through-route allowlist. Re-run `build_external_links.py` when boundary nodes change or OSRM data is updated.
`data/external_intra_times.json` — committed; output of `build_intra_times.py`. Per external zone, `M`=30 sampled intra-zonal OSRM times (s) for the production-suppression self-term (`model.load_self_terms` → `constrained_od_flows`). Committed so the model runs without re-querying OSRM. Re-run `build_intra_times.py` only when external zones change.
`data/manual_link_overrides.json` — committed so manual assignments survive a wipe of `counts_processed.json`.
`simulation/tuner_config.json` — committed as source config (gitignore exception).
`analysis/hourly_fractions.csv` — committed as source data (single authoritative version).
`analysis/generation_rates.json` — committed; per-component vehicle-driver trips/person/day (output of `derive_generation_rates.py`). Source data for generation pinning (`model.compute_generation_scales`).

### Large reference data (gitignored, kept locally only)
`data/*.ods`, `data/*.xlsx`, boundary GeoJSON files — too large to commit; keep local copies.
Currently present:
- `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods` — used by `parse_official_hourly.py`.
- `data/nts0502.ods` — DfT NTS Table NTS0502a, weekday trip start times. Used by `derive_component_profiles.py`.
- `data/nts0504.ods` — DfT NTS Table NTS0504b, average trips by day/purpose. Used by `derive_component_profiles.py`.
- `data/census-2021-apwp001.xlsx` — DZ-level workplace population. Used by `build_demographics.py` (internal nodes) and `build_census_zones.py` (external nodes).

Boundary files needed by `build_census_zones.py` (download from NISRA / OpenDataNI):
- `simulation/dz2021/DZ2021.geojson` — DZ polygon boundaries (present, gitignored).
- `simulation/sdz2021/SDZ2021.geojson` — SDZ polygon boundaries (present, gitignored).
- `simulation/dea2021/DEA2021.geojson` — DEA polygon boundaries (present, gitignored).

RoI data files for `build_wz_apportionment.py` + `build_census_zones.py` (in `data/ireland_data/`):
- `Small_Area_National_Statistical_Boundaries_2022_Ungeneralised_view_*.geojson` — 2022 SA boundaries (~410 MB).
- `Complete_set_of_Census_2022_SAPs/SAPS_2022_Small_Area_UR_171024.csv` — SA population (`T1_1AGETT`).
- `Workplace_Zones_ITM/Workplace_Zones_ITM.shp` — 2016 WZ boundaries in EPSG:2157 with workplace headcount (`T11_C1` = total workers; `T1_T` is total population, unused); used only by `build_wz_apportionment.py`.
- `cache_sa_workplace.csv` — **generated** by `build_wz_apportionment.py`; committed once computed.

---

## Model Design

### Gravity model
The assignment is **production-constrained per component** (`model.constrained_od_flows`):
`T^c_ij = K_c · p^c_i · a^c_j · f_c(d_ij) / D^c_i`, `D^c_i = Σ_k a^c_k·f_c(d_ik)`, so each origin's
trip production is fixed by its producing weight `p^c_i` and is independent of accessibility (fixes
the generation/distribution conflation). The four components (residential / commute / retail /
school), their producer/attractor weights and per-component kernels are detailed in "Four-component
flow decomposition" below. See the agent memory note `project_production_constrained_gravity`.

**Tanner deterrence kernel** (`model._tanner_kernel`; per component, own peak `P_c` + shape `BETA_c`):
`u = d/P; f(u) = u^BETA · exp(BETA·(1 − u))`.

Properties: peak f(P) = 1 (at d = P seconds), f(0) = 0, rise ~ u^BETA near origin, and an
**exponential tail** `~exp(−BETA·d/P)` for large d (characteristic decay `1/γ = P/BETA`). This is
the "keep the form (u=d/P, peak at P), swap the tail" replacement for the old rational kernel's heavy
power-law tail `~1/d^ALPHA` (Phase 1, 2026; memory `project-tanner-kernel-tld`) — the power-law tail
was ill-conditioned (the optimiser drove `ALPHA` to extremes on a flat ridge, over-concentrating
external flow). Each component has its own `P_c`, `BETA_c` (**8 shape params**; no shared BETA, no
ALPHA). Phase 2 (later) will anchor the shape to surveyed trip-length distributions instead of
free-fitting. The kernel is **model-layer** (applied to `od_dist` in `constrained_od_flows`, not
routing) — changing it needs no `build_paths` rebuild, only a re-tune.

(`model._rational_kernel` + `model.gravity_assign` retain the unconstrained power-law kernel for
back-compatibility with old pre-split param files; not used by the current model.)

Distances are least-time shortest paths (seconds). For external→internal OD pairs, the path traverses an OSRM-derived external edge (X→B, fixed weight) then the internal road network (B→J). Dijkstra selects the optimal boundary entry node for each destination.

### Stochastic route choice (probit loading)
The paths cache stores fractional link-assignment weights computed from `N_PASSES=25`
Dijkstra runs, each with log-normal edge-cost noise (CV=0.25). For each OD pair,
`link_weight[entry]` is the fraction of passes that routed through that link. Pairs
with no topological route diversity (degree-1 stubs, single-access nodes) converge to
weight=1.0 on their forced route. `od_dist` is the mean path distance across passes.
THETA is not tuned.

### External zones (big-world network)
NI is represented as a three-level hierarchy centred on Newtownards (CENTRE):

- **Core area** (DZ level): union of all DZs whose parent SDZ intersects `CORE_RADIUS` (3 km). Boundary is irregular (follows census polygon edges, not a circle).
- **SDZ external nodes**: SDZs within `SDZ_ZONE_RADIUS` (10 km) that are not in the core — one centroid node per SDZ.
- **DEA external nodes**: DEAs entirely outside `SDZ_ZONE_RADIUS` — one centroid node per DEA.

Once a DEA is broken into SDZs, all its constituent SDZs become nodes (even those beyond `SDZ_ZONE_RADIUS`). Each external node's ID **is** its census-area code (the SDZ/DEA 2021 code, a string such as `"N21000219"`) — the `id` field in `census_zones.json`; there is no separate small-integer ID. (Road node IDs are OSM integers; external node IDs are these strings.)

**Demand:** all layers are per-zone (Census 2021 DZ/SA aggregated to SDZ/DEA/ED/LEA), kept as **separate layers** (nothing is summed — each component uses its own producer/attractor, see "Four-component flow decomposition"). Producers: `commute_producers` (resident commuters) and `school_producers` (resident students), both census-derived, plus population. Attractors: `workplace` (place-of-work jobs — commute attractor), `retail_spaces` (parking spaces within the zone via `parking_demand.parking_spaces` — retail attractor), and `school_demand` (per-zone OSM school enrolment via `school_demand.assign_enrolments` — school attractor).

**Connectivity:** boundary nodes are all internal nodes with at least one road edge crossing the core polygon boundary. OSRM-derived directed edges connect each external node to its valid boundary nodes; the gravity model path distance for any external→internal pair is the OSRM leg plus the internal shortest path, computed automatically by Dijkstra on the augmented graph.

**No tuning of external zone values.** Pop/wp come directly from census data and are not adjusted by the optimizer (there are no hand-crafted city configs, ref values, or damping factors).

### Through routes
External→external OD pairs are allowed only for pairs whose OSRM route passes through at least one boundary node (i.e. genuinely transits the core area). This allowlist is auto-generated by `build_external_links.py` and stored in `data/external_links.json`. Changing the allowlist requires re-running `build_external_links.py` and `build_paths.py` (the paths cache must be rebuilt).

### Four-component flow decomposition
The gravity OD flows are split into four **production-constrained** spatial components at each
tuner evaluation (per-pair pre-K flows from `model.constrained_od_flows`, scattered onto links;
`p^c_i`/`a^c_j` = producing/attracting weight, `D^c_i = Σ_k a^c_k·f_c(d_ik)`). The three non-res
components are independent clones: a symmetric pop↔activity split, each leg per-origin-normalised,
with **no weight parameter and no self/cross term**:

- **Residential** (`flow_res`, Tanner kernel P/BETA): `T^res_ij = pop_i·pop_j·f_res/D^res,pop_i` —
  pop×pop trips. Single leg (i→j and j→i are separate OD pairs, so both directions are covered).
- **Commute** (`flow_commute`, Tanner kernel P_commute/BETA_commute): home→work producer =
  `commute_producers`, attractor = `workplace`; return work→home producer = `workplace`, attractor
  = pop — `f_com·( commprod_i·work_j/D^com,work_i + work_i·pop_j/D^com,pop_i )`.
- **Retail** (`flow_retail`, Tanner kernel P_retail/BETA_retail): home→shop producer = pop, attractor
  = `retail_spaces`; return shop→home producer = `retail_spaces`, attractor = pop —
  `f_ret·( pop_i·ret_j/D^ret,ret_i + ret_i·pop_j/D^ret,pop_i )`.
- **School** (`flow_school`, Tanner kernel P_school/BETA_school): home→school producer =
  `school_producers` (census resident students), attractor = `school_demand` (OSM school places);
  return school→home producer = `school_demand`, attractor = pop —
  `f_sch·( schoolprod_i·school_j/D^sch,sch_i + school_i·pop_j/D^sch,pop_i )`. External
  `school_demand` and `school_producers` are both populated per zone, so external school trips are
  retained intra-zonally via the self-term (below) rather than dumping into the core.

Per-leg producer weights are scaled to **vehicle-driver trips/day** by `gen_scale` so each `K_c ≈ 1`
(see "Generation pinning").

**External intra-zonal self-term (denominator-only).** Each per-origin denominator
`D^c_i = Σ_k a^c_k·f_c(d_ik)` runs over *other* zones; collapsing an external zone to one centroid
drops its `k=i` diagonal (its intra-zonal trips), so `D^c_i` is too small and the external origin's
fixed budget over-allocates to the observed core (worst for large, isolated, far DEAs). `build_intra_times.py`
OSRM-samples `M=30` uniform intra-zonal point-pairs per external zone → `data/external_intra_times.json`;
`model.constrained_od_flows` then adds `a^c_i·(1/M)·Σ_m f_c(t_im)` to each denominator (the mean kernel
over the sample, `E[f]`, not `f(mean)`). It is **denominator-only** — no link flow — and applies to
**external zones only** (internal road nodes have no zone area). Direct OSRM sampling avoids any
characteristic-distance constant, speed assumption, or zone-shape model (real per-zone times carry local
speed + network detour). Effect is **strongly kernel-tail-dependent** (sharper `BETA` under the Tanner kernel): under a sharp kernel tail an external
zone's short intra-zonal times give a large `f(t_intra)` that dominates its denominator, substantially
cutting exported external→core flow — concentrated on the near/mid SDZs that carry essentially all
core-bound flow (the far DEAs already contribute ~0 core-bound flow, so they barely move). Wired into
`build_assignment.py` and `tune_assignment.py` via `model.load_self_terms`; absent file ⇒ no self-term.
Independent of the paths cache (model-layer, not a routing input) — no `build_paths` rebuild needed.

The self-term applies to every component. External `school_demand`/`school_producers` and the
`workplace`/`retail_spaces` layers are all populated per zone, so every component's denominator carries
its intra-zonal diagonal and external zones retain school/commute/retail trips locally rather than
dumping them into the core.

Each component has its own **independent temporal shape** (`hourly_fractions.csv`, normalised so
`W_c = 1` — magnitude/split is generation's job, not the temporal profile's) and scale
(K_res, K_commute, K_retail, K_sch).
Predicted count for observation i in slot s:
`pred_i = Σ_c K_c·flow_c·(T/3600)·f_c[s]`  over  c ∈ {res, commute, retail, school}.

### Generation pinning (data-based supply)
Producer weights are carried in absolute **vehicle-driver trips/day** so each component's tuned
scale **K_c should land at ≈ 1.0** — a *verification anchor*, not a fit knob (a `K_c` away from 1
diagnoses local car-mobilisation vs the national average, to be refined later). This is a
**model-layer** change (no paths-cache rebuild). Independent of, and a prerequisite for, the
Tanner kernel (now implemented, Phase 1; see memory `project-tanner-kernel-tld`).

**Rates `ρ_c`** (`analysis/generation_rates.json`, written by `analysis/derive_generation_rates.py`
from England NTS0409a, 2023/24 avg, vehicle-driver modes = *Car or van driver + Motorcycle + Taxi
or minicab*): commute 0.178, retail 0.495, school 0.082, res 0.262 trips/person/day (sum ≈ 1.017).
Using the driver row makes it vehicles by construction.

**Purpose→component mapping — JUDGMENT ALLOCATIONS (candidate error sources, kept flagged on
purpose).** Organising principle = the attractor each component offers: workplace(jobs)→commute,
retail_spaces(**parking** = all commercial/venue)→retail, school→school, population(**homes**)→res.
NTS0409a cannot sub-split these by car-driver mode, so each below is a modelling *decision*, the
first thing to revisit if the fit is scrutinised:
- **Business → retail** (not commute): commute kept pure home↔own-workplace; business visits hit
  commercial premises (parking), not the workplace-jobs count.
- **Personal business → retail**: services/banks/medical are commercial/parking destinations.
- **Leisure split ½ retail / ½ res** (`LEISURE_RETAIL_FRAC = 0.5` in `derive_generation_rates.py`):
  venue leisure → parking(retail), visit-friends-at-home → pop↔pop(res). The 0.5 is an assumed
  fraction and the single largest such assumption (leisure is the biggest bucket).

**Per-leg producer scaling (`model.compute_generation_scales`).** res is single-leg (full `ρ_res`;
both directions are separate OD pairs); commute/retail/school are two-leg, each direction carrying
`ρ_c/2`. The per-producer rate is `r_leg = (ρ_c·share)/k_leg`, applied to the **producer term only**
inside `constrained_od_flows` (the same array's attractor use is scale-invariant — cancels in
`attr_j/D_i`, denominators untouched; `gen_scale=None` ⇒ all 1.0 ⇒ exact prior flows). The anchors
`k = Σ(producer layer)/Σ(population)` are summed **island-wide** from the node weights (the external
census nodes tile the whole island, so they recompute for any CENTRE — the transferability win):
`k_commuters`(commute_producers), `k_jobs`(workplace), `k_retail`(retail_spaces),
`k_students`(school_producers), `k_enrolment`(school_demand); `k=1` when the producer is population
(res, retail outbound). **Independent per-side:** the per-producer rate is applied to each side's
real count, so a job-/retail-rich sub-area generates more activity→home (PM-outflow) than its
residents alone — the open-region asymmetry, fed in the morning by the external zones' outbound legs.
At island scale this equals the balanced total (`Σ_leg = ρ_c·½·Pop`, since island counts = `k·Pop`);
the difference is purely the sub-regional directional balance.

**φ-share prior follow-up:** with shares now largely set by `ρ_c`, the `phi_*` priors in
`tuner_config.json` (kept as-is for now as the K-degeneracy guard) could later be re-derived from
the ρ-implied shares or relaxed — a follow-up, not part of this step.

### Direct-K convex scale solve
At each optimizer evaluation the four component scales **(K_res, K_commute, K_retail, K_sch)** are
calibrated **directly** by `solve_scales`. The temporal fractions `f_res/f_commute/f_retail/f_school`
are **pinned at the NTS profile** (`hourly_fractions.csv` `mean_fraction_*`) and never tuned — so every
observation prediction is **linear in the scales**: `pred_i = Σ_c K_c·a^c_i` with
`a^res_i = m_res_i·Th_i·f_res[s_i]` constant (commute/retail/school analogously). The inner objective —
Gaussian WLS over the official obs + Poisson **identity-link** deviance `2·Σ(n·log(n/pred)+pred−n)` over
the walking obs + a scale-share prior — is therefore **convex over K ≥ 0**, and is solved by a **damped
(Levenberg) Newton step with a backtracking line search on the full objective**: monotone by
construction, so there is **no K-collapse and no best-iterate bookkeeping**. `CALIBRATE_PROBE=1` reports
the residual global scale λ at the start params (≈1 ⇒ K at its optimum).

**Scale-share prior (degeneracy break).** `φ_commute = K_commute/ΣK`, `φ_retail = K_retail/ΣK`,
`φ_sch = K_sch/ΣK` are computed as derived ratios and penalised `~ N(phi_commute_prior,
phi_commute_std²)`, `~ N(phi_retail_prior, phi_retail_std²)`, `~ N(phi_school_prior, phi_school_std²)`
(from `tuner_config.json`). This anchors the component shares so the convex K-solve can't trade
magnitude between similarly-shaped components. It regularises the inner K-solve only and is **not**
part of the reported χ² (exactly as the old φ-prior was not). Walking obs mostly fall in slots with
`f_school≈0`, so `K_sch` is pinned almost entirely by the official school-peak hours — the joint
Poisson+Gaussian solve does not inflate the school share. (With generation pinning, the shares are
also largely set a priori by `ρ_c`; see the φ-prior follow-up note under "Generation pinning".)

**Why this is sound.** With `f` fixed the K-problem is genuinely convex (Poisson identity-link
deviance is convex in the mean), so a single small Newton solve reaches the global inner optimum.
Freezing `f` at NTS is justified by the NTS-vs-official hourly shape
match (Pearson r > 0.97 at all three sites; only a smooth ~4–6% overnight/midday bias).

Slot key: (day_type, hour), day_type = 0 (weekday), 1 (Saturday), 2 (Sunday). The pinned NTS school
profile is a sharp weekday double-peak (h08/h15), near-zero weekends.

### Observations
Observations are in count space with per-obs weights (the live counts are printed by
`tune_assignment.py` at run start):
- **Official hourly** (24 h × 3 day-types × 3 sites): from `data/official_hourly.json`;
  Gaussian error (sigma from between-weekday std, 10% floor); weight = 1/sigma².
- **Walking**: from `data/link_aadt.json`; Poisson error; weight = 1/n_eff.

### Goodness of fit
`χ²/N` = mean squared z-score over all observations. **N_eff = N** since the temporal fractions are
pinned at NTS and not fitted — no per-slot temporal df are consumed (the few global df — gravity shape
params + 4 scales — are not subtracted, following the convention that only per-slot temporal df count).
The tuner's χ²/N is pure data-fit (Gaussian + Poisson deviance); the f-prior/coupling penalty is
identically zero with `f` pinned.

`build_assignment.py` uses the four-component `compute_chi2()` (passing the commute/retail/school link-flow dicts) when the 4-component params are present in `tuned_params.json` — a **data-only** chi²/N (pure sum of squared z-scores). Since the tuner's χ²/N is now also pure data-fit (`f` pinned ⇒ no f-prior/coupling penalty), the two surfaces are directly comparable (both read the same observed links; `build_assignment` keeps the full scatter for the map). The legacy single-K `gravity_assign`/Woodbury path is used only for old pre-split param files.

**Reading "modelled flow" across reports.** The three reporting surfaces print *different
projections* of the same tuned model — they are not directly comparable line-for-line:
- `build_assignment.py` "Official count sites" block and the `"flows"` values in
  `newtownards_flows.json` → **directed daily AADT** = `Σ_c K_c·flow_c·W_c` over the four components,
  where `W_c` (`model.aadt_weights`) is the day-type-weighted (5·weekday+Sat+Sun)/7 sum of component
  `c`'s hourly fractions. With the **decoupled per-component shapes** (each normalised so `W_c ≈ 1`,
  `Σ_c W_c ≈ 4`), `K_c·flow_c` is already ≈ the component's daily AADT; the `W_c` factor is retained
  because it reads the actual slot fractions (and was load-bearing under the old partition scheme,
  where each `W_c` was a sub-1 share). (The unweighted per-component flows feed `compute_chi2`, which
  applies `f_c` itself — do not double-weight.) `newtownards_flows.json` stores the W-weighted AADT in
  `flows`/`flows_res`/`flows_commute`/`flows_retail`/`flows_school` plus an `aadt_weights` block.
  Node-based sites (508/444) sum every directed link at the node.
- `newtownards_map.html` combined layer → the same AADT but **summed over both directions** of each
  edge (`flow(u,v)+flow(v,u)`), i.e. a two-way total (~2× a single directed link).
- The tuner / `report_tune.py` fit table → **per-observation, count-space**: official rows are
  *vehicles/hour* in one (day_type, hour) slot (≈ AADT × hourly fraction), walking rows are
  reconstructed to combined AADT. Correct for goodness-of-fit; not a table of link AADTs.
- Walking "Model" column convention (display only, chi²/N unaffected): both `compute_chi2()`
  (`model.py`) and the tuner fit table show **combined directed AADT** for walking links — the sum of
  the per-component modelled AADTs `m_res+m_commute+m_retail+m_sch`.
- Fit-table street names come from the consolidated GraphML edge `name` attribute. `tune_assignment.py`
  resolves the GraphML data-key id dynamically from the `<key>` header. `report_tune.py` echoes the
  labels the tuner stored in history, so names appear only after a fresh tune run regenerates
  `tuning_history.jsonl`.

---

## Count Data

**Official hourly obs (from ODS, parsed by `parse_official_hourly.py`):**
- Site 507: A21 Bangor Road — 72 obs (24 h × 3 day-types); Gaussian sigma, 10%/15% floors
- Site 508: A48 Donaghadee Road — 72 obs
- Site 444: A20 Portaferry Road — 72 obs

Annual AADT values (507: 21,202; 508: 10,792; 444: 7,282) are retained in `model.py`
`COUNT_SITES` for `build_assignment.py` backward compatibility but are not used by the
tuner directly.

**Walking counts:** raw recorder CSVs in `data/counts/`; processed per-session AADT in
`data/counts_processed.json`, aggregated per-link (inverse-variance) in `data/link_aadt.json`.
Manual GPS-snap overrides live in `data/manual_link_overrides.json` (see `manual_assign_link.py`).
Sessions on links absorbed by the dead-end reduction (Westmount Park, Old Belfast Road) are dropped
via `EXCLUDE_LINKS` in `model.py`. The live observation count and slot breakdown are printed by
`tune_assignment.py` at run start.

---

## Tuning History

Per-run history — params, χ²/N, N obs, and notes — is logged one line per run to
`simulation/tuning_history.jsonl` (the authoritative record). Inspect it with
`simulation/restore_params.py --list` or `analysis/report_tune.py`; it is not duplicated here.

The **Tanner kernel** (Phase 1) is now in place — the rational power-law tail (which drove the
diagnosed `ALPHA` blow-up / external over-concentration) is replaced by the exponential-tail Tanner
form. A 4-component / generation-pinned / Tanner re-tune is the next step to confirm it resolves the
pathology (see memory `project-tanner-kernel-tld`); Phase 2 (TLD anchoring) is the follow-up.

A carry-forward open concern is the **school component**: it is weakly identified (walking obs mostly
fall in slots with `f_school≈0`, so it is pinned almost entirely by the official school-peak hours)
and can act as an AM/PM-peak fitter, pushing the school share up — to be re-evaluated under the Tanner
kernel (kernel tail shape / school-peak count data is the lever, not prior strength).

## Paths Cache

The paths cache (`newtownards_paths.npz`) must be rebuilt with `build_paths.py` whenever:
- The road network changes (`newtownards_consolidated.graphml`)
- External links change (`data/external_links.json` — re-run `build_external_links.py` first)
- The tuned profile (`simulation/tuned_profile.json`) or base speeds (`data/google_cache/base_speeds.json`) change
- `N_PASSES`, `PROBIT_CV`, or `PROBIT_LL_SIGMA` change

**Staleness guard (loud failure):** `build_paths.py` stamps a signature of its inputs into the
`.npz` — SHA-1 of the reduced routing graph (`newtownards_reduced.graphml`), SHA-1 of
`data/external_links.json`, SHA-1 of the tuned profile (`simulation/tuned_profile.json`) and
empirical base speeds (`data/google_cache/base_speeds.json`), and the probit noise params
`PROBIT_CV` / `PROBIT_LL_SIGMA` (`simulation/routing_config.py`). A missing profile/base-speeds
file is stamped as a `MISSING:<path>` sentinel so the comparison stays well-defined. `tune_assignment.py` and
`build_assignment.py` call `model.assert_paths_cache_fresh(cache)` right after loading the cache
and **raise `SystemExit`** if any input changed (or the cache predates the guard), naming what
changed and telling you to re-run `build_paths.py`. Helper lives in `simulation/model.py`
(`paths_cache_signature`, `assert_paths_cache_fresh`).

**Current cache format** (probit): `node_ids` covers road nodes (OSM integer IDs) + external nodes (census-area-code strings, e.g. `"N21000219"`); `link_u`/`link_v` are road-link endpoints only (external edges are not in `link_list`); `link_weight` (float32, fraction of passes using that link for each OD pair); `od_dist` (mean path distance across passes including external legs); `probit_n_passes`, `probit_cv`, `probit_ll_sigma`. `n_routed_pairs` marks the flow-carrying OD pairs occupying indices `0..n_routed_pairs-1`; the remainder are **denominator-only non-through ext→ext virtual edges** (entries in `od_src/od_dst/od_dist` but NOT in `pair_idx/link_idx` — they complete each external origin's production-constrained denominator and carry no link flow). No `pair_idx_2/3` keys — `_has_stoch = False`, THETA not tuned.

---

## Known model behaviour
- **Scale-share degeneracy:** the convex `solve_scales` could otherwise trade magnitude between
  similarly-shaped components; the `phi_commute`/`phi_retail`/`phi_school` Gaussian priors anchor the
  shares (`φ = K_c/ΣK`). With generation pinning the shares are also largely set a priori by the NTS
  rates `ρ_c`.
- **Component scales** `K_res`/`K_commute`/`K_retail`/`K_sch` are calibrated directly by the convex
  `solve_scales` at each optimizer step. With generation pinned to vehicle-driver trips/day each
  `K_c ≈ 1` (a verification anchor); a deviation diagnoses local car-mobilisation vs the national
  average. chi²/N is the reliable fit metric.
- After a structural model change (e.g. new count data or external link regeneration), a fresh tune is needed to restore fit quality.
- **External node probit loading:** all edges (road and external) receive log-normal noise each pass (CV=0.25), length-scaled by the per-edge gain `w = σ_ll/(σ_ll + CV·cost)` (`PROBIT_LL_SIGMA`, default 120 s). Route diversity for external-internal OD pairs comes from both the X→B external leg and the internal B→J portion, giving stochastic spread across similarly-weighted boundary entry points. The length-scaled gain caps a long leg's perturbation at an *absolute* sigma of `σ_ll` (~2 min) while leaving short internal edges' multiplicative noise essentially unchanged (`w≈1`), so boundary selection is driven by real time differences rather than swamped by a long external leg's noise. Default `σ_ll=120 s` is anchored to `CV × a typical in-town journey (~8 min)`; the knob lives in `simulation/routing_config.py` and is part of the paths-cache staleness signature.
- **Manual link overrides:** Use `analysis/manual_assign_link.py <session_id> <from_node> <to_node>` to assign a session to a specific directed link, bypassing GPS snap. Use when the observer stood on a parallel carriageway and the snap would land on the wrong physical road. The override is stored in `data/manual_link_overrides.json` and takes effect even if `counts_processed.json` is wiped and rebuilt. After assignment (manual or auto), `ingest_counts.py` validates each non-null count direction against the directed graph and raises `ValueError` if the edge doesn't exist.
- **Temporal profiles are pinned, not inferred:** `f_res/f_commute/f_retail/f_school` are fixed at
  the NTS `mean_fraction_*` profile (`hourly_fractions.csv`) and never tuned, so no per-slot df are
  consumed (N_eff = N). They are **independent per-component shapes** (each normalised so `W_c = 1`),
  decoupled from magnitude — the inter-component split is set by generation, not the temporal
  profiles, so they do **not** partition the aggregate profile.
- **Dead-end street absorption (ghost edges):** OSMnx `simplify_graph` treats bidirectional
  dead-end terminus nodes as degree-2 and removes them, so the dead-end edge vanishes from the
  consolidated graph; uncorrected, buildings on those streets would snap to the nearest surviving
  edge (often but not reliably the main road). `build_demographics.py` detects these absorbed termini
  by comparing raw and consolidated network nodes, reconstructs their UTM geometry from the raw
  network, and adds "ghost" edges to the STRtree; buildings snapping to a ghost edge have all their
  demand attributed to the surviving junction node (the only network entry point for that street).
  Demand-allocation only — no effect on `build_paths.py`, `model.py`, or the paths cache.
- **`tuned_params.json` structure:** the four scales `K_res`/`K_commute`/`K_retail`/`K_sch`, the 8 Tanner shape params (P, BETA, P_commute, BETA_commute, P_retail, BETA_retail, P_school, BETA_school), `"kernel": "tanner"`, and `slot_fracs_res`/`slot_fracs_commute`/`slot_fracs_retail`/`slot_fracs_school` (dicts keyed `"dt,h"`, the pinned NTS profile). `reset_gravity_params.py` regenerates this clean structure; old pre-split param files fall back to the legacy `gravity_assign` mode in `build_assignment.py`.

---

## External Zone Configuration

External zone values are now fully data-driven from Census 2021 (via `data/census_zones.json`) and OSRM routing (via `data/external_links.json`). There are no hand-crafted reference values, dampings, or city groupings to maintain.

**Gravity param refs** (`tuner_config.json` `gravity_ref` / `gravity_lambda`): anchor the L2 regularization of the shape params. The values live in `tuner_config.json` and **must not be changed without explicit approval**.

**To update external zone coverage** (e.g. after a NISRA boundary update):
1. Re-run `build_census_zones.py` (updates `data/census_zones.json`)
2. Re-run `build_demographics.py` (updates `node_weights.json`)
3. Re-run `build_external_links.py` (updates `data/external_links.json`)
4. Re-run `build_paths.py` (rebuilds paths cache with new external nodes)
5. Re-tune, then `build_assignment.py` and `build_map.py`

**External demand** — the separate attractor layers `workplace`, `retail_spaces`, and `school_demand` (and the census producers) — is measured per zone via the same estimators as internal nodes, with no external scale factors. Known limitation: OSM under-maps RoI schools (~80% of real schools mapped), so RoI external school demand runs proportionally low.

**Demand model — open items (TODO).** (1) Optional: precompute per-small-area retail/school enrolment island-wide (like `cache_sa_workplace.csv`) so the external aggregation isn't redone per CENTRE. (2) **Car-ownership mobilisation — TBC, pending NTS microdata SN 5340.** Generation currently assumes a uniform car-mobilisation level for all nodes. The intended refinement is a zone-varying car-driver trip-rate multiplier = Σ_band (zone persons per car-availability band, from NI+RoI census) × (car-**driver** trips/person/purpose for that band) — car ownership modulating the car-driver mode share of a roughly-fixed trip budget, so the 2nd-car effect is sub-linear per person. The required `(driver-trips × per-person × car-band × purpose)` cross-tab is **not a published NTS table** (NTS0205 gives the band distribution, NTS0409 per-purpose car trips/person, but not the joint) — only derivable from SN 5340 (access requested). Deliberately deferred rather than approximated from marginals (independence would discard the band×purpose interaction). Not a pipeline blocker; slots in later as a generation refinement.

---

## Google Routing-Time Calibration (offline, optional — NOT part of the main pipeline)

**Status — COMPLETED.** A calibrated per-`(road-class × speed-band)` speed-factor
profile (`simulation/tuned_profile.json`) was fit against Google and deployed (`compile_profile.py`
→ `car_roaaads.lua` → OSRM rebuilt), bringing OSRM **external-corridor** times into line with
Google (offline `predicted/Google` medians X2B 1.00, B2X 0.99, X2X 1.00). Re-tuned on the
combined v1+v2 skeleton cache and confirmed stable; the v2-combined profile is **deployed**.
**Residuals not fully resolved:** in-town
(INT) ~12% too fast (median ≈0.88) — a **turn/junction-model** gap, not base speed; specific
external corridors notably **Ballyrainey** improved but not fully matched. **Stage 2 (route
preference) explored but deferred:** class-only preference factors (`tune_preference.py`, 13
classes with urban/rural split) cannot achieve a net ranking improvement with the current data —
67 true violations vs 1033 external concordant pairs (~1:15) means any factor large enough to
resolve a violation flips 4–5× more concordant pairs. Both timing and preference calibration
share the same `_FAC`/`_PREF` two-table Lua architecture; the `_PREF` block is wired and ready
but no preference file is compiled into `car_roaaads.lua` until the conditioning problem is solved.

**Purpose / design.** OSRM was systematically *too fast* (worst in-town and on some external approach
corridors), inflating external→core flow and hurting the fit. The workflow uses **Google Routes API
as a journey-time source-of-truth** to calibrate a realistic OSRM time profile, **decoupling
impedance (travel time) from route preference (generalised cost)** — conflated in the old single
`HIGHWAY_COST_FACTOR`. The error was length-structured, confirming a **turn/junction penalty** is
needed, not just per-class speed factors (the INT residual).

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
cache in milliseconds. (Segments are labelled by exact OSM node-id → tag lookup rather than reading
back `annotation.speed`, which corrupts short urban-edge buckets.)
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

**Calibration status — COMPLETED:** empirical base speeds make the offline model a faithful proxy
for real OSRM on external corridors (verify per-leg medians ≈ 1.00–1.03); the external-focused factor
tune lands X2B/B2X/X2X medians ≈ 0.99–1.00 with physically sensible factors (motorway ~Google
free-flow, urban A/B-roads slowed), compiled to `car_roaaads.lua` and deployed. Remaining work is in
**Status** at the top of this section.

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
