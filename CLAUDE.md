# All-Ireland Gravity Traffic Model — Project Overview

A gravity-model traffic assignment pipeline for **any centre on the island of Ireland**,
calibrated against walking count data and official AADT figures. The pipeline is fully
reproducible: running the scripts in order regenerates all outputs from raw data.
**Newtownards is the current calibration centre, not the model's scope** (see the portability
rule below).

**⚠️ THE MODEL IS LOCATION-PORTABLE — read before designing any data layer.** `CENTRE`
(`simulation/zones_config.py`) is a **free parameter**: the core/study area can be relocated to
**any centre on the island of Ireland**. **Newtownards is only the current calibration centre, not
the model's scope.** Therefore **every data layer is built island-wide** (census
producers/attractors, schools, parking, road network) so a new centre needs **no per-location data
work**, and **no per-location manual step may be the mechanism** for anything (no hand-picking "the
schools near Newtownards", no verifying a *fixed handful* of core POIs, no location-specific
overrides as the primary path) — when CENTRE moves, a *different* set of areas becomes core, so any
location-specific step silently breaks the new deployment. Design every layer **and every
quality/precision step** (geocoding, snapping, estimators) for **"whichever areas fall in the
active core", uniformly island-wide**; manual overrides are only a targeted backstop.

**Failure mode to avoid (this is the one agents keep hitting):** **do not sequence, scope, defer,
or judge the relevance of any work by its impact on the Newtownards fit.** That silently
re-privileges one centre. The count fit at Newtownards is a *falsification test* of a transferable
model — not the objective — so "this barely moves the current fit" is **never** a reason to do
something coarsely, defer it, or skip a jurisdiction. **"External" nodes = areas outside the
*currently selected* core** (a per-centre partition that changes when CENTRE moves), **not** a
lesser or peripheral category: the Republic of Ireland is the *core* for some centres, NI for
others. Every jurisdiction, zone, and data layer is equally load-bearing; build and refine them
symmetrically. Before deferring/prioritizing/scoping, check: *is my criterion "the current fit"? If
so, stop and re-justify island-wide.* See agent memory `feedback-model-is-portable` and
`project-model-transferability`.

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
python3 simulation/build_schools.py          # island-wide OSM schools + per-POI enrolment → data/cache_osm_schools_island.geojson (one-off; osmctools + school_demand estimator; needs DEA boundary for NI/RoI tagging; OSM source for build_admin_schools — geocoding + third-level)
python3 simulation/build_admin_schools.py    # admin-roll school-age (NI+RoI rolls, geocoded) + OSM third-level → data/cache_admin_schools_island.geojson (the unified school-demand cache; level-tagged; resets the cache each run, applies manual_school_coords.json)
python3 simulation/geocode_school_tail.py    # Nominatim-geocode the ~2% offline tail of the admin cache (external queries; cached/resumable; run after build_admin_schools.py)
python3 simulation/build_census_zones.py     # classify NI+RoI census areas → data/census_zones.json incl. per-zone retail_spaces + per-level school_demand/school_producers (one-off; needs SDZ/DEA boundary files + cache_sa_workplace.csv + island parking/school caches)
python3 simulation/build_network.py          # build road network from local all-island .osm.pbf via osmium (core polygon + 5km bbox; needs Docker)
python3 simulation/build_demographics.py     # node weights + boundary detection + external weights (map is built separately by build_map.py)
python3 simulation/build_external_links.py   # OSRM queries → external↔boundary links + through-route allowlist (needs local OSRM)
python3 simulation/build_intra_times.py      # mass-weighted per-component intra-zonal self-term per external zone → data/external_intra_times.json (denominator self-term; needs local OSRM + road-point cache; independent of build_paths — no cache rebuild)
python3 simulation/reduce_deadends.py        # collapse residential dead-ends → newtownards_reduced.graphml + node_weights_reduced.json (consumed by build_paths/build_assignment/tune; see reduce_deadends.py row)
python3 simulation/build_paths.py            # probit stochastic paths incl. external nodes (N_PASSES=25, CV=0.25, N_WORKERS=1; build time depends on hardware)

python3 analysis/parse_official_hourly.py    # parse ODS hourly counts → data/official_hourly.json (one-off)
python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/derive_school_generation.py               # per-STUDENT school gen (escort + self-drive) from NTS microdata + DfE England age→level split → analysis/school_generation_rates.json (run before derive_generation_rates; needs data/NTS + the DfE participation CSV + node_weights)
python3 analysis/derive_generation_rates.py                # regenerate generation_rates.json — per-capita vehicle-driver trips/person/day from NTS microdata (23-cat B01 mapping); school = per-student × island students/pop; retail += pre-school escort fudge (re-run when the microdata, mapping, or school rates change)
python3 analysis/derive_component_profiles.py              # regenerate all six hourly_fractions.csv shape columns from NTS microdata (car-specific; res/commute/retail joint dow×hour; school = per-bin escort regression + age→level self-drive; needs the DfE participation CSV for the school split)

python3 analysis/tune_assignment.py                        # tune gravity params (18 double-exp willingness params, 6 production-constrained scales; external zones fixed from census)
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
| `simulation/build_wz_apportionment.py` | **(RoI data prep — one-off).** Pre-computes the WZ→SA workplace apportionment for all of RoI and writes `data/ireland_data/cache_sa_workplace.csv` (columns: `sa_code`, `workplace_pop`, `commute_car`). CSO 2016 Workplace Zone (WZ) boundaries do not align with 2022 SA boundaries; this script intersects them geometrically via `gpd.overlay(wz, sa, how="intersection")`, bypassing 2016→2022 SA boundary change codes entirely (geometry is the ground truth). Each WZ's `T11_C1` headcount ("total workers in workplace zone" — place-of-work jobs; **not** `T1_T`, which is total daytime population) is split across the intersection pieces weighted by the sum of POI weights (`EXCLUDE_AMENITY`/`POI_WEIGHTS` from `demographics_config`) falling within each piece; area-proportional fallback for zero-POI pieces. `commute_car` (the RoI car-commute attractor for `census_attractor.py`) apportions the WZ daytime-driver columns `T2_M5+T2_M6+T2_M8` (motorcycle + car driver + van; present in the shapefile) with the **same** split weights, × national work-driver share **0.9588** (nets out self-driving 3rd-level students mixed into the daytime driver count). POIs are extracted from the local PBF via the `osmctools-roaaads` Docker image (`osmfilter --keep-nodes="amenity= shop= office=" --drop-ways --drop-relations`) and cached to `data/ireland_data/cache_roi_pois.geojson`. Will reuse `osrm/edge_index/ni.o5m` if present to skip the slow PBF→o5m conversion step. Re-run only when WZ or SA boundaries change or OSM POI data is significantly stale. **Needs Docker + local PBF.** |
| `simulation/ingest_ni_census.py` | Loads NI DZ/SDZ/DEA boundaries + NISRA population + workplace into standardised GeoDataFrames consumed by `build_census_zones.py`. Public API: `load_ni_census() → (dz_gdf, sdz_gdf, dea_gdf)`. Standardised columns: `area_code`, `parent_code`, `level`, `population`, `workplace_pop`, `geometry` (in `PROJECTED_CRS`). Handles DZ→SDZ parent lookup via column or spatial join fallback; SDZ→DEA similarly. Population fetched from NISRA API (cached to `data/cache_nisra_population.csv`). Workplace from `data/census-2021-apwp001.xlsx`. |
| `simulation/ingest_roi_census.py` | Loads RoI SA/ED/LEA boundaries + CSO 2022 population + pre-computed WZ workplace into standardised GeoDataFrames consumed by `build_census_zones.py`. Public API: `load_roi_census() → (sa_gdf, ed_gdf, lea_gdf)`. Fails loud if `data/ireland_data/cache_sa_workplace.csv` is missing (run `build_wz_apportionment.py` first). ED and LEA GeoDataFrames are derived by dissolving SAs — no separate boundary file needed. Standardised columns match `ingest_ni_census.py`. (The SAPS `T1_1AGETT` *column* sums to ~2× the national population because the file carries a "State" aggregate row equal to the sum of all SAs; the per-SA `SA_PUB2022` join excludes it, so the loaded per-SA population is the correct 1×.) |
| `simulation/census_supply.py` | `load_supply() → {area_code: {commute}}` — per-small-area **commute** trip *producer* harmonised NI (DZ) + RoI (SA): residents who *drive* to work — car-driver modes only (RoI `SAP2022 T11T1` travel-to-work Car Driver + Van + Motorcycle; NI `transport_to_workplace` "Driving a car or van" + "Motorcycle, scooter or moped" + "Taxi"). WFH/not-in-employment are structurally excluded (separate travel-method categories, not selected). The car restriction matches the vehicle-driver modes in `derive_generation_rates.py`, since the model assigns car flow. RoI key = JSON-stat SA-code label (the GUID index is ignored); a national "State" aggregate row (~2M) is dropped so per-SA values are clean 1×. Source data (gitignored): `data/ireland_census/` (RoI CSO SAP JSON-stat), `data/ni_census/` (NISRA DZ CSVs). The **school producers** are now a separate per-level module (`census_school_producers.py`); the `_one`/`_ni_csv` readers here are shared with it. |
| `simulation/census_school_producers.py` | `load_school_producers() → {area_code: {primary, postprimary, tertiary}}` — per-small-area **per-level** school trip *producers* (resident students by level) harmonised NI (DZ) + RoI (SA), the producer counterpart of `school_attractor.py`'s per-level enrolment. Replaces the old lumped `census_supply` "school" producer + its childcare subtraction: pre-school childcare is excluded by construction (each level is a school-age/tertiary headcount). RoI age-band→level fractions are **data-derived** from CSO PxStat EDA42 (primary pupils by age) + EDA70 (post-primary total) against census age persons (`_F04_PRIMARY`/`_F59_PRIMARY`/`_F1014_PRIMARY`/`_F1014_POSTPRIM`/`_F1519_POSTPRIM`), scaled per-jurisdiction to the admin enrolment total. `island_enrolment_by_level()` (island totals per level) feeds `derive_generation_rates.py`. (The canonical level tuple `SCHOOL_LEVELS = ("primary","postprimary","tertiary")` lives in `model.py`.) |
| `simulation/census_attractor.py` | `load_attractor() → {area_code: car_commute_jobs}` — per-small-area car-commute *attractor* (jobs reached by car) harmonised NI (DZ) + RoI (SA), the demand counterpart of `census_supply.py`'s commute producer. Computed **once, island-wide** (jurisdiction handled internally; keyed by `area_code` so any CENTRE on the island works) and consumed by both `build_census_zones.py` (external zones — aggregated) and `build_demographics.py` (internal core — POI-distributed). **NI:** `car_jobs[DZ] = apwp001 DZ workplace total × car_share[parent SDZ]`, `car_share[SDZ] = (Driving + Motorcycle + Taxi) / Σ(all apwp035 methods incl. WFH)` (`data/ni_census/census-2021-apwp035.xlsx` sheet `SDZ`). WFH removed via the WFH-inclusive denominator — exact because `apwp001[SDZ] ≡ Σ apwp035 method columns` (verified across all 850 SDZs). DZ→SDZ from the `DZ2021_cd`/`SDZ2021_cd` attribute columns of `DZ2021.geojson` (no geometry op). **RoI:** the `commute_car` column of `cache_sa_workplace.csv` (WZ daytime drivers; fails loud if absent). |
| `simulation/build_census_zones.py` | Classifies the full island of Ireland into a three-level census hierarchy centred on `CENTRE` — works for any CENTRE on the island. Calls `load_ni_census()` (NI DZ/SDZ/DEA) and `load_roi_census()` (RoI SA/ED/LEA), concatenates the two hierarchies, then runs unified classification: small areas intersecting `CORE_RADIUS` (3 km) → core; intermediate zones in broken outer zones → SDZ/ED external nodes; non-core small areas in partially-core intermediate zones → orphan DZ/SA external nodes; outer zones outside `SDZ_ZONE_RADIUS` (10 km) → single DEA/LEA centroid nodes. Population-weighted centroids computed from constituent small areas. Outputs `data/census_zones.json`: core polygon (WGS84), external node list with `id` = census-area code (`"N21000219"` for NI, `"017001001"` for RoI SA), `level`, centroid, population, workplace_pop, `commute_attractor` (car-commute jobs, from `census_attractor.py`), `retail_spaces`, the three per-level `school_demand_<level>`, `commute_producers`, and the three per-level `school_producers_<level>`. NI codes start with `'N'`; RoI codes are pure-numeric. `retail_spaces` = sum of `parking_demand.parking_spaces` over island-parking polygons within each zone (sjoin; workplace-derived fallback for zones with no mapped parking). `school_demand_<level>` = per-zone sum of per-POI enrolment (from the unified school cache — `build_admin_schools.py`: admin-roll school-age + OSM third-level) split into primary/post-primary/tertiary by `school_attractor.py` (0 for zones with no school of that level); `school_producers_<level>` from `census_school_producers.py`. |
| `simulation/build_parking.py` | Builds the island-wide OSM parking cache → `data/cache_osm_parking_island.geojson` (gitignored), the single parking source for `build_census_zones.py` (external zones) and `build_demographics.py` (internal core). Streams parking ways from the all-island pbf via **osmctools** (reuses `ni.o5m`, then `osmfilter --keep="amenity=parking landuse=parking"`), assembles closed-way polygons (RAM-light ~0.5 GB). Saves each polygon with the tags the estimator reads (`access`, `parking`, `building`, `building:levels`, `parking:levels`, `capacity`, `fee`, `amenity`, `landuse`, `name`). **Needs Docker + the pbf/ni.o5m.** |
| `simulation/parking_demand.py` | Pure-stdlib `parking_spaces(tags, area_m2)` → estimated retail parking **spaces** for one OSM parking polygon. Recipe: exclude `access ∈ {private,no,permit}`; decks (`parking ∈ {multi-storey,underground,rooftop}` or `building=parking`) trust `capacity` (else `area×levels/30`), gate-exempt; else `capacity` only if implied `area/capacity ∈ [8,80] m²/space`, else area fallback `÷13` on-street (`street_side`/`lane`) or `÷30` otherwise. Constants in `demographics_config.py`. Destination car parks land at ~29 m²/space in both NI and RoI. Tests: `simulation/test_parking_demand.py`. |
| `simulation/build_schools.py` | Builds the island-wide OSM school cache → `data/cache_osm_schools_island.geojson` (gitignored). **No longer the model's school source** — it is the OSM input to `build_admin_schools.py` (which uses it to geocode NI admin schools by name and to supply the third-level `college`/`university` POIs). Streams `amenity=school/college/university/kindergarten` from the pbf via **osmctools** (reuses `ni.o5m`), tags each POI's jurisdiction (NI vs RoI via the DEA boundary), applies `school_demand.assign_enrolments` globally, and saves one point per kept POI with `enrolment`, `amenity`, `name`. **Needs Docker + pbf/ni.o5m + DEA boundary.** |
| `simulation/build_admin_schools.py` | Builds the **unified island school-demand cache** → `data/cache_admin_schools_island.geojson` (gitignored), the school source for `build_census_zones.py` + `build_demographics.py` (via `SCHOOL_ISLAND_CACHE`). One Point feature per school: `jurisdiction`, `level` (primary/post_primary/special/**tertiary**), `name`, `enrolment`, `geocode_method`, `matched_osm_name`, `needs_review`. **School-age** from admin rolls (NI `data/ni_data/School level - …`: primary school-age = Total − nursery − pre-school, post-primary, special from the Sex sheet; RoI `Data_on_Individual_Schools_*`: Mainstream + Special tab + post-primary, in-file coords, pure-`Boarding` dropped). **Third-level** taken from the OSM cache (`college`/`university`) but **enrolment recomputed live via `school_demand._assign_tertiary`** (curated full-time HE/FE; junk/part-time/mis-tagged-secondary dropped) — so tertiary curation applies on this run with no `build_schools.py` (Docker/pbf) rebuild. Stage-1 lumped. **NI geocoding** (rolls lack coords): match to the OSM cache by name, gated by amenity↔level + DZ containment + a same-school identity check (de-accented core tokens + role modifiers); doubtful → tail. RoI post-primary coordinate-column quirks (ITM easting / full ITM pair) recovered by projection inversion. **Manual overrides** from `data/manual_school_coords.json` (tracked) applied last. ~98% geocoded offline; the tail → `geocode_school_tail.py`. Level-tagged so the primary/post-primary/tertiary split needs no re-ingest. |
| `simulation/geocode_school_tail.py` | Geocodes the null-geometry tail of the admin cache via **Nominatim** (≤1 req/s, identifying User-Agent, cached to `data/cache_nominatim_schools.json`, resumable). Query cascade name+town → street+town+postcode → postcode/eircode centroid, first hit. NI hits validated against the stated DZ: in-DZ → trusted; near-but-outside → flagged; **>3 km → rejected** (same-name wrong match) and left for manual. Run after `build_admin_schools.py` (which resets the cache). **External queries.** |
| `simulation/school_demand.py` | `assign_enrolments(features)` → per-POI school enrolment, operating on the full island set (cross-feature clustering + institution splitting). Used by `build_schools.py` for the OSM cache and — for **third-level only** — recomputed live by `build_admin_schools.py`. **Schools:** light dedup; primary/secondary classified by `school=`/name, valued by jurisdiction-aware averages (`SCHOOL_ENROLL`). **Kindergartens (pre-school) are EXCLUDED** (the school component drops pre-school on both sides — admin-roll school-age excludes nursery, and the per-level census producers (`census_school_producers.py`) are school-age/tertiary headcounts with no pre-school by construction). **Third-level (`_assign_tertiary`) is CURATED-ONLY** (no fallback): each POI matches a curated HE institution (`INSTITUTIONS`, **full-time**: HEA 2024/25 RoI + DfE/HESA 2023/24 NI; Ulster split by campus via `CAMPUS_FT`, GB campuses excluded), a curated FE college (`FE_INSTITUTIONS` — NI DfE Table A4 per-college + Teagasc/CAFRE agri), or the RoI public-FE keep-set (national SOLAS FT total 65,851 distributed by **method (a)** total→institution→POI, `_fe_instkey` grouping). Institution total split across its matched POIs. Everything unmatched is **dropped** — OSM junk, part-time FET/adult (Youthreach/VTOS/training/`_ROI_FE_EXCL`), and second-level 'colleges' already in the admin school-age rolls. Tests: `simulation/test_school_demand.py`. |
| `simulation/school_attractor.py` | `add_level_enrolments(gdf) → gdf` — splits the unified school cache's per-POI `enrolment` into the three per-level attractor columns `enrol_primary`/`enrol_postprimary`/`enrol_tertiary` (`LEVEL_ENROL_COLS`) using each POI's `level` tag from `build_admin_schools.py` (primary/post_primary/tertiary; `special` split into primary vs post-primary by the NI primary:post-primary enrolment ratio 0.5126). The per-level school *attractor* counterpart of `census_school_producers.py`'s per-level producer. Consumed by `build_census_zones.py` (external zones — aggregated per level) and `build_demographics.py` (internal core — POI-distributed per level). |
| `simulation/build_network.py` | Builds the road graph from the **local NI `.osm.pbf`** (the same Geofabrik snapshot OSRM is built from — `demographics_config.PBF_PATH`), so road/boundary/internal node IDs share one OSM snapshot with OSRM's route node IDs. The full ~400 MB island pbf OOMs an in-process parse, so a small extract is streamed out with **osmctools** (`osmconvert` + `osmfilter`; Docker image `osmctools-roaaads`, auto-built from `simulation/osmctools.Dockerfile`, ~0.5 GB peak RAM): `osmconvert -b=<bbox> --complete-ways` (bbox = core polygon buffered by `demographics_config.BOUNDARY_BBOX_MARGIN_M` = 5 km) then `osmfilter --keep="highway=<drive set>"` (positive form of osmnx's `drive` filter), written to `simulation/_pbf_drive_extract.osm`. (osmctools is used rather than osmium-tool, whose referenced-node id-set is sized by OSM's max node id and needs several GB regardless of extract area.) `ox.graph_from_xml` reads it; `graph_from_xml` omits the `street_count` node attribute, which `consolidate_intersections` needs, so it is re-added via `ox.stats.count_streets_per_node` (without it the core under-merges). Raw graph extends 5 km beyond the core (for boundary nodes' external neighbours + `build_external_links.py` positions); the consolidated routing graph is still clipped to the core polygon, then junction-consolidated (tol 15 m) and relabelled to OSM IDs. Outputs `newtownards_network.graphml` (raw) + `newtownards_consolidated.graphml`. **Needs Docker + the pbf on disk.** |
| `simulation/build_external_links.py` | Queries a local OSRM instance (all-island extract, **biased `car_roaaads.lua` profile** — see `build_osrm_profile.py`, `http://localhost:5000`) to derive all external zone connectivity. **X→B links:** for each (external node, boundary node) ordered pair, keeps the link only if that boundary node is the first boundary node encountered in the OSRM route (i.e., it is the natural entry point into the core) (`_classify_xb_link`); duration = the route total minus its post-B travel (the penalty-inclusive time to *reach* B). Snap anomalies **fail loud**: the external origin snapping onto a core node (centroid inside/on the core), or B absent from the route (e.g. a carriageway-twin snap). **B→X links:** symmetric with X→B — keeps B→X only if no other boundary node appears in the OSRM route sequence (i.e., B is the last boundary node departed on the way to X). If another boundary node B' appears, the journey is already covered by B→B' + B'→X (`_classify_bx_link`); duration = the route total minus travel up to B's *last* occurrence (strips a start junction loop). Snap anomalies **fail loud**: the origin boundary node absent or not at sequence index 0/1, or B returning to itself after a **>500 m excursion** (`REVISIT_MAX_LOOP_M`; a nearer re-pass — e.g. B on a roundabout, which collapses to one model node — is benign and priced out). The crow-flies revisit check needs per-node geometry, fetched with a lazy second OSRM query **only** when B re-appears (rare), keeping the common path geometry-free. **Boundary→boundary exterior shortcuts:** for each ordered boundary pair (B1, B2), keeps a directed shortcut only if the OSRM route stays *outside* the core — B2 is the first boundary node reached after B1 **and** no core-interior model node is passed en route (`_classify_bb_shortcut`); duration = the route total minus its post-B2 travel **and** minus travel up to B1's last occurrence (penalty-inclusive — strips a one-way overshoot past B2 and a start junction loop; the earlier per-edge annotation sum dropped OSRM turn/signal penalties — the bug this fixes on b→b, whereas X→B/B→X already priced off the penalty-inclusive route total). A B1 return within 500 m (roundabout) is benign; a farther return fails loud (same `REVISIT_MAX_LOOP_M` / lazy-geometry mechanism as B→X). This corrects an earlier `node_seq[1]`-only "exits core first" test that mis-kept straight-through-town routes (`node_seq[1]` is usually an off-model shape point, never in `internal_node_ids`). Coordinate-snap edge cases (the origin/destination boundary node absent or misplaced in the node sequence — e.g. a dual-carriageway-twin snap) **fail loud** for review. *Limitation:* a boundary→boundary road crossing the interior on a junction-free corridor has no model junction to detect, so it can still be recorded as a shortcut — harmless (nothing is measured there) beyond a mild probit route-choice perturbation from the redundant near-duplicate of the real internal route; closing it fully would need a route-polyline-vs-core-polygon test. **Through-route allowlist:** for each ordered external-external pair, checks if any OSRM route node is a boundary node; if so adds to `allowed_through_pairs` (`_ext_ext_transits_core`, which also **fails loud** if either external centroid snapped onto a core node — centroid inside/on the core). Outputs `data/external_links.json`. ~28,000 OSRM queries; under a minute on a local instance. |
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, detects boundary nodes, adds external node weights, writes `node_weights.json` + `newtownards_demographics.geojson`. The map is built separately by `build_map.py` (run it afterwards). `--zones-only` re-reads `data/census_zones.json` and patches only the external node entries in `node_weights.json`. Shared constants (paths, OSM tag handling, map styling) live in `simulation/demographics_config.py`. **Boundary node detection:** loads core polygon from `census_zones.json` and the **pbf-sourced** raw graph (`newtownards_network.graphml` from `build_network.py`), identifies internal nodes (within core polygon), then boundary nodes = internal nodes with at least one edge going outside. Writes `boundary_node_ids`/`internal_node_ids` to `node_weights.json`. Because the raw graph comes from the same OSM snapshot as OSRM, these IDs match OSRM's route node IDs exactly (so `build_external_links.py`'s boundary/internal route-sequence matching is exact). **External node weights:** reads the external node list from `census_zones.json` and writes population + the demand layers (workplace, commute_attractor, retail_spaces, commute producers, per-level school producers, per-level school_demand) to `node_weights.json`. **Study area = core polygon (not a circle):** DZ selection and all OSM downloads (buildings/POIs/parking) are bounded by the core polygon from `census_zones.json` (extent `max_core_vertex_dist_m`, ~10.2 km), matching the road graph built by `build_network.py`. Core DZs are selected by centroid-within the polygon (recovers exactly the `n_core_dzs` core DZs) and use **full** DZ population/workplace_pop (no area-fraction clipping). OSM downloads use a circle sized to the polygon (+1 km margin); POIs and parking are then filtered to within the core polygon to avoid margin leakage (buildings are already DZ-bounded by sjoin). **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Demand layers (separate, never summed):** `node_commute_attractor` (car-commute jobs — commute attractor; per-area car_jobs from `census_attractor.load_attractor()`, island-wide and CENTRE-agnostic, distributed within each area by the same edge-snapped POI weights as `node_workplace`), `node_workplace` (all workplace jobs distributed within each DZ by POI count — kept for the retail-spaces fallback + map, **not** an attractor; currently NI-only internally) and `node_retail_spaces` (retail parking spaces via `parking_demand.parking_spaces` on the island parking cache clipped to the core, snapped to road edges — retail attractor) are written as independent layers (no `node_business_demand`); `node_commute_producers` comes from census (`census_supply`) and the three `node_school_producers_<level>` from `census_school_producers`. External nodes take each per-zone layer from `census_zones.json`. **School demand:** per-POI enrolment from the unified school cache (`build_admin_schools.py`: admin-roll school-age + OSM third-level) is split by level (`school_attractor.py`) → the three `node_school_demand_<level>` layers (internal POIs snapped to core road edges; external nodes take per-zone `school_demand_<level>` from `census_zones.json`). There is no single lumped `node_school_demand`/`node_school_producers` layer — the levels are always separate. |
| `simulation/build_map.py` | Builds the interactive folium map (`newtownards_map.html`) from artifacts written by `build_demographics.py` (`node_weights.json`, `newtownards_demographics.geojson`), the road graphs, the cached OSM POI/parking layers, and — if present — `newtownards_flows.json`. This was the old `build_demographics.py --map-only` path, now a standalone step (it always reloads POI/parking from cache). Run after `build_demographics.py`, and again after `build_assignment.py` to refresh flow layers. **Flow map layers:** combined AADT (default), plus per-component residential / commute / retail / school layers. No args (`--help` only). |
| `simulation/demographics_config.py` | Shared pure-constant config imported by `build_demographics.py`, `build_map.py` **and `build_network.py`** (file paths, OSM tag handling — `EXCLUDE_AMENITY`/`POI_WEIGHTS` — and map styling). `CENTRE` is re-exported from `zones_config.py` (not defined here). Also holds the **OSRM backend location — the single source of truth**: **`OSRM_DIR`** (the local OSRM data dir holding the `.osm.pbf`, built `.osrm` files, and `car_roaaads.lua`), which defaults to the sibling `osrm/` directory next to the repo and is overridable via the **`ROAAADS_OSRM_DIR`** environment variable for a one-time non-default layout. **`PBF_PATH`** and **`OSRM_LUA`** are derived from `OSRM_DIR`, and every OSRM consumer (`build_network`/`build_schools`/`build_parking`/`build_wz_apportionment`, the profile tooling `build_edge_index`/`compile_profile`/`build_skeleton_index`, and `build_n_of_t`) imports `OSRM_DIR`/`PBF_PATH`/`OSRM_LUA` from here rather than hardcoding a path. (Repo-relative paths are self-located: scripts derive `REPO_ROOT` from `__file__`, no hardcoded repo path anywhere.) Also holds **`BOUNDARY_BBOX_MARGIN_M`** (5 km buffer around the core polygon for `build_network.py`'s extract). The separate `NETWORK_MARGIN_M` (1 km) here sizes the OSM POI/building/parking download circle only — unrelated to the road graph. Also defines **`PROJECTED_CRS = "EPSG:2157"`** (Irish Transverse Mercator / ITM) — the single source of truth for all projected spatial operations in the pipeline. ITM covers the whole island of Ireland with uniform accuracy, avoiding UTM Zone 30N's distortion for Republic of Ireland towns west of ~6°W. All simulation and analysis scripts import this constant; `build_network.py` passes it explicitly to `ox.project_graph`. **Do not hardcode `EPSG:32630` anywhere.** Single source of truth so the split scripts don't drift. Also holds the **parking→retail-spaces estimator constants** (`PARKING_M2_PER_SPACE_OFFSTREET=30`/`_ONSTREET=13`, `PARKING_GATE_LO=8`/`PARKING_GATE_HI=80`, `PARKING_EXCLUDE_ACCESS`, `PARKING_DECK_TYPES`, `PARKING_ONSTREET_TYPES`) consumed by `parking_demand.py`, plus `PARKING_ISLAND_CACHE` and `SCHOOL_ISLAND_CACHE` (island parking/school cache paths). |
| `simulation/build_intra_times.py` | Builds the **mass-weighted, per-component** intra-zonal **self-term** per external census zone. Reconstructs each zone's member small areas (ingest-loader parent maps, as `build_census_zones.py`), then per component (res/commute/retail + 3 school levels) samples origins ∝ producer and destinations ∝ attractor **within the zone** — real POIs (parking ∝ spaces, schools ∝ per-level enrolment) for retail/school, road-snapped member-area points for res/commute (retail falls back to area-level ∝ `retail_spaces` where a zone has no mapped parking, mirroring the attractor's workplace-derived fallback) — one OSRM `/table` per batch (independent draws ⇒ S×D matrix ∝ p⊗a), histogrammed → `data/external_intra_times.json` (`{census_code: {component: {t:[bin-centre s], w:[weights Σ=1]}}}` + `_meta`). Captures **clustering** (people/jobs/schools co-locating in villages ⇒ short intra-zonal times ⇒ strong self-suppression), unlike the old uniform-in-polygon single-average which sampled empty fields; and because `p·a·f` is symmetric one histogram serves both legs (the leg-asymmetry dissolves). Reuses `build_n_of_t` (`osrm_table`/`road_point`/`build_point_cache`/`load_area_masses`/`load_poi_layers`) — self-term ≡ `n(t)` restricted to same-zone pairs; shares `data/_area_road_points.json`. `--s/--d/--batches` set the fixed-generous budget (default 50×50×8 = 20k pairs/zone-component, from a convergence probe); `--component X` re-samples just one trip type and **merges** it into the existing output (keeps the other five) — for targeted refreshes. **Run after `build_census_zones.py`, OSRM up. Independent of `build_paths.py`** — model-layer, no paths rebuild. Re-tune afterwards. |
| `simulation/build_opportunity_table.py` | **(national n(t) data prep — one-off).** Island-wide opportunity table → `data/island_opportunity_table.csv` (committed): one row per census small area (NI DZ + RoI SA, ~22.7k) with each area's producer/attractor masses (population, commute producers/attractor, per-level school producers/demand, retail parking spaces) + WGS84 centroid. The per-area aggregation `build_census_zones.py` does for external zones, run **island-wide with no core/external classification** (each small area is the unit). Reuses the identical estimators (`census_supply`/`census_school_producers`/`census_attractor`/`parking_demand`/`school_attractor`), so per-area values aggregate to the `census_zones.json` external-node values (verified). Feeds `analysis/build_n_of_t.py`. Needs the census/boundary/cache reference data (as `build_census_zones.py`); no OSRM. |
| `analysis/build_n_of_t.py` | **(national n(t) sampler).** Builds the empirical national Ireland opportunity-density-in-cost `n_Ire(t)` per purpose → `data/national_n_of_t.json`: point-cloud Monte-Carlo of `n(t)=Σ P_i·A_j·δ(c_ij−t)` (P=producer/origin, A=attractor/dest; OSRM car time on `car_roaaads.lua`), origins ∝ producer / dests ∝ attractor, batched via OSRM **`/table`** (independent draws ⇒ the B×B outer product is ∝ P_i·A_j, accumulated unweighted). Six purposes (res/commute/retail + 3 school levels); area→area dests are road-proximate points (`/nearest` snap filter), retail/school dests are the real POIs ∝ parking_spaces/enrolment. **v1 = unconstrained N(t), outbound leg only** (both flagged simplifications). **`--stratified` (default for the frozen build)** resolves the short-range head: cache K road points per area once (purpose-independent, resumable to `data/_area_road_points.json`), partition pairs by centroid-haversine into near bands (0-3/3-10/10-25 km) + a far tail, compute each band's exact P·A mass via a `cKDTree` ring sum, sample each band with its own budget (near bands generous ⇒ dense head), and reconstruct `n(t)=Σ_b M_b·ŝ_b(t)` — unbiased (banding only allocates samples; heights come from routed times × exact M_b). The naive outer-product path (`--pairs`) is retained; `--pilot --purpose P` runs one purpose + a head/geometry diagnostic plot. Needs OSRM up (island, `localhost:5000`); reads `island_opportunity_table.csv` + DZ/SA polygons. Full spec `task_empirical_n_of_t.md`; recovery-side `n_Eng` concept in memory `project-n-eng-source-geometry`. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Covers both internal road nodes and external census-area nodes. **Graph augmentation:** loads external nodes and edges from `data/external_links.json`; adds them to the routing graph before Dijkstra. External edges (X↔B, boundary shortcuts) are included in the adjacency matrix but NOT in `link_list` — they contribute to path distance but not to flow accumulation. **Probit loading:** all edges (road and external) perturbed each pass with log-normal noise `exp(eps·w)`, `eps ~ N(0, CV=0.25)`, N_PASSES=25, giving stochastic spread in boundary node selection for external→internal OD pairs with similarly-weighted entry options. **Length-scaled noise (`PROBIT_LL_SIGMA`, default 120 s ≈ 2 min):** the per-edge gain `w = σ_ll/(σ_ll + CV·cost)` ∈ (0,1] keeps the noise multiplicative for short legs (`w→1`) but saturates it to a fixed *absolute* sigma `σ_ll` for long legs (`w→σ_ll/(CV·cost)`), so a long single-edge external↔boundary leg's perturbation no longer swamps the few-minute differences between competing boundary entries. The adjusted perturbation never exceeds the pure multiplicative one, and `eps=0 ⇒ no bias`. Crossover at cost ≈ σ_ll/CV (~8 min). `PROBIT_CV` and `PROBIT_LL_SIGMA` are imported from `simulation/routing_config.py` (the gain vector is precomputed once, constant across passes). **OD pair filter:** through-routed external→external pairs (in `allowed_through_pairs`) are routed flow pairs (Dijkstra path through the core). Non-through external→external pairs (from `external_external_times`) are appended as **denominator-only** pairs — entries in `od_src/od_dst/od_dist` (distance = direct OSRM time) but NOT in `pair_idx/link_idx` and excluded from `src_groups`/probit passes, so they carry no flow; they complete each external origin's production-constrained denominator. The cache stamps `n_routed_pairs` (flow-carrying pairs occupy `0..n_routed_pairs-1`). No offscreen leg calculation. **Internal edge costs:** each road edge's routing cost is the OSRM-equivalent travel time `factor(class,band)·length·3.6/base_speed` from the Google-calibrated profile (`simulation/tuned_profile.json` + empirical `base_speeds.json`), via `simulation/edge_speed.py`. `add_edge_speeds`/`add_edge_travel_times` are still run because `deadend_collapsed` synthetic edges keep their osmnx-encoded intra-region time (factor 1.0). Internal routes are chosen on realistic time alone (no route-preference biasing); the node-based Dijkstra applies no turn penalties. Re-run if road network, external links, the tuned profile or base speeds, `N_PASSES`, `PROBIT_CV`, or `PROBIT_LL_SIGMA` change. **Performance constants:** `N_WORKERS` (default 1) controls parallel pass workers via `multiprocessing.Pool` — increase on machines with sufficient RAM (each worker uses ~100–150 MB extra); `MAX_HOPS` (default 120) caps per-pair path-trace iterations. Inner path-tracing loop is vectorised (numpy). |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, **`_modesub_kernel(d, wparams, component)`** (the production kernel `driveShare(equiv_miles(d), component)·[w·exp(−d/τs)+(1−w)·exp(−d/τl)]`, `wparams=(w, τs, τl)` — per-component driveshare curve × double-exp willingness, `component` required; imports `equiv_miles`/`driveshare` from `../analysis`), **`constrained_od_flows()`** (production-constrained per-component per-pair pre-K flows + per-origin denominators, kernel per component via each component's `(w, τs, τl)` willingness; the **three school levels** are passed as `w_school_levels`/`w_school_prod_levels` dicts keyed by `SCHOOL_LEVELS` and returned in `t_sch_by_level` — they are fully independent blocks, **each with its own per-level driveshare curve AND its own willingness** (`school_primary/postprimary/tertiary`, so the school kernel + both denominators are computed per level), with **no shared `τ_school`**, `SCHOOL_LEVELS = ("primary","postprimary","tertiary")`; optional `self_terms` (a `{component: (self_src, self_dist, self_w)}` dict) adds the **external intra-zonal self-term** `a^c_i·Σ_bin w·F_c(t_bin)` to each denominator — denominator-only, per-component, both legs share the component entry, `None` ⇒ no self-term; `return_legs=True` additionally returns the per-producer→attractor legs dict — `res`/`<comp>_out`/`<comp>_ret` — for home-anchored diagnostics like `diagnose_per_capita.py`, default off ⇒ the return signature is unchanged) and **`scatter_od_to_links()`** (the production-constrained assignment core, used by `build_assignment.py` and `tune_assignment.py`), **`load_self_terms(node_ids)`** (builds the per-component self-term arrays from `data/external_intra_times.json`'s per-zone per-component weighted time histograms; skips zones absent from `node_ids`; returns a `{component: (src,dist,w)}` dict or `None` if the file is missing), **`load_generation_rates()`** + **`compute_generation_scales(node_weights, rates)`** (the per-leg producer scales that pin generation to vehicle-driver trips/day — island per-capita anchors from the node-weight layer sums; feed `constrained_od_flows(..., gen_scale=…)`; see "Generation pinning"), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `compute_chi2()` takes the commute/retail link-flow dicts + the per-level school link-flow dicts (`link_flow_school_dicts`) + per-component `slot_fracs_*` (school per level, `slot_fracs_school_levels`); in the multi-component production-constrained mode N_eff = N (no per-slot df subtracted); legacy single-flow path otherwise. Road node IDs are OSM integers (stable); external census node IDs are census-area-code strings (e.g. `"N21000219"`) — not OSM IDs. COUNT_SITES: site 507 links 538692601↔549139252; site 508 node 136173611; site 444 node 449111329. **`WEIGHTS_FILE` and `ROUTING_GRAPH` point at the dead-end-reduced artifacts** (`node_weights_reduced.json`, `newtownards_reduced.graphml`) from `reduce_deadends.py`. EXCLUDE_LINKS: `{(181844513, 181839481)}` plus the Westmount Park and Old Belfast Road directed links (both directions) whose endpoints are absorbed by `reduce_deadends.py` and no longer exist in the reduced graph — their walking observations are discarded from calibration (regenerate this set from `deadend_broken_obs.json` if the reduction params change). |
| `simulation/build_assignment.py` | **Production-constrained** gravity assignment (via `model.constrained_od_flows` + `scatter_od_to_links`). Requires `simulation/newtownards_paths.npz`. Multi-component mode activated when the six K's + the 18 double-exp willingness params (`model.willingness_keys()`) are present in `tuned_params.json` (each school level active when its `K_<level>` > 0 and `node_school_demand_<level>` exist). Applies the per-leg `gen_scale` (`model.compute_generation_scales`) so component magnitudes are vehicle-driver trips/day. Saves `flows_res`, `flows_commute`, `flows_retail`, and the three `flows_school_primary`/`_postprimary`/`_tertiary` in `newtownards_flows.json` (`build_map.py` combines the school levels for display). Requires the multi-component params (the six K's + commute/retail/school kernels) — fails loud on an old single-K param file. External node weights come from `node_weights.json` directly (no override from tuned params). **Honours the `doubly_constrained` set** from `tuned_params.json` (Furness attraction constraint per component); it leaves `furness_max_sweeps=None` so each flagged leg **cold-converges** (exact deployed flows) rather than using the tuner's warm k-sweep — so a doubly-constrained build is substantially slower than a singly one, with the short-kernel school `IPF capped` warnings. See "Doubly-constrained option". |
| `simulation/edge_speed.py` | OSRM-equivalent internal edge-time model shared by `build_paths.py` and `reduce_deadends.py` (paths used in `model.paths_cache_signature`). `load_profile()` loads the tuned `ProfileSpec` (`simulation/tuned_profile.json`) + empirical base speeds (`data/google_cache/base_speeds.json`) — **fails loud** if the tuned profile is missing (warns if base speeds absent → analytical fallback). `edge_time_seconds(tags, length_m, spec) = factor(class,band)·length·3.6/base_speed(class,band)` reuses `profile_spec` for all bucketisation/base-speed/factor logic (same `(class×band)` buckets the deployed `car_roaaads.lua` keys on, so internal routing matches the tuned OSRM instance). Robust to `highway`/`maxspeed` stored as a list. Edge impedance only — no turn penalties (the internal Dijkstra is node-based). |
| `simulation/routing_config.py` | Holds `HIGHWAY_COST_FACTOR` + the probit noise params (`PROBIT_CV`/`PROBIT_LL_SIGMA`). `HIGHWAY_COST_FACTOR` does not drive internal routing (`build_paths.py`) or the dead-end reducer (`reduce_deadends.py`) — those use the calibrated `(class×band)` profile via `simulation/edge_speed.py`. It is used only by the legacy tooling that references it: `build_osrm_profile.py`, `build_skeleton_index.py --base-speeds`, and `skeleton_model.legacy_spec_from_highway_cost_factor`. |
| `simulation/build_osrm_profile.py` | Generates `car_roaaads.lua` — the road-class-biased OSRM car profile. Pulls the default `car.lua` from the `osrm/osrm-backend` Docker image, injects a block after the `forward_rate` assignment that divides `forward_speed`/`forward_rate` by `HIGHWAY_COST_FACTOR` (matching internal Dijkstra biasing). Re-run whenever `HIGHWAY_COST_FACTOR` changes, then re-preprocess OSRM (`osrm-extract -p car_roaaads.lua`, `osrm-partition`, `osrm-customize`). Output: `/home/matthew/Documents/CodingFun/osrm/car_roaaads.lua`. |
| `simulation/reduce_deadends.py` | Collapses "residential dead-end" regions in the consolidated routing graph to shrink node count (speeds up `build_paths.py`/tuning, enables a larger core area). A region R (entrance E ∉ R) qualifies iff: (1) R connects to the rest of the network through exactly one cut vertex E; (2) R contains no boundary node and no school-demand node (both *protected* — never absorbed — which enforces the no-boundary and zero-school rules structurally); (3) max directed journey time E→n over n∈R < `T_MAX` (default 60 routing-cost seconds); (4) total workplace+retail demand < `BIZ_CAP` (default 100; residential pop unbounded); (5) `|R| ≥ 2` (single-node spurs skipped — 1→1 saves nothing). **Algorithm:** every valid region is a protected-free connected component of H−a (H = undirected simple projection) for some articulation point a, so it enumerates all such (entrance, region) candidates, filters by constraints 2–5 + directed reachability both ways, and selects the *maximal feasible* regions (laminar family ⇒ disjoint; naturally descends into an oversized branch to find the largest collapsible sub-pockets — catches cyclic closes that leaf-pruning would miss). Each region → one super-node S (=min id, summed pop/workplace/retail/school, pop-weighted UTM centroid) joined to E by directed links E→S, S→E whose travel times are population-weighted means of the intra-region directed times. **Intra-region times use the same OSRM-equivalent `(class×band)` edge model as `build_paths.py`** (`build_cost_digraph(G, spec)` via `simulation/edge_speed.py`), so the collapse/`T_MAX` decisions match the speeds the reduced graph is later routed on. Synthetic edges use `highway="deadend_collapsed"` (factor 1.0 in `build_paths.py`) with `maxspeed`+`length` set so osmnx's `add_edge_speeds`/`add_edge_travel_times` (re-run by `build_paths.py`) reproduce that target time. **Run after `build_demographics.py` (needs pop/biz/school + boundary) and before `build_paths.py`.** Outputs (gitignored): `newtownards_reduced.graphml`, `node_weights_reduced.json`, `deadend_map.json` (provenance: super-node→absorbed nodes + times), `deadend_broken_obs.json` (observed/count links whose endpoints were eaten — **manual review before adoption**; observed-link endpoints are deliberately *not* protected). Params: `--t-max`, `--biz-cap`. **Wired into the pipeline:** `build_paths.py` (`CONS_GRAPH`), `build_assignment.py` (`CONS_GRAPH`), `tune_assignment.py` (`CONS_GRAPH`) read `newtownards_reduced.graphml`, and `model.WEIGHTS_FILE`/`ROUTING_GRAPH` point at the reduced files — so this step must run after `build_demographics.py`. The 6 absorbed walking observations (Westmount Park, Old Belfast Road) are discarded via `EXCLUDE_LINKS` in `model.py`. **Map caveat:** `build_map.py` still draws the *full* consolidated graph, so flow on collapsed interior streets is not shown on the map (demand layers are unaffected; main-road flows and the fit are unaffected). Re-mapping collapsed regions via their super-nodes is a possible follow-up. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Gravity param regularization and the scale K-prior. `gravity_lambda` + `gravity_ref` regularise the 18 double-exp willingness params (`<comp>_taus`/`_taul`/`_w` for res/commute/retail + the three independent school levels; τ's in seconds — `gravity_ref` holds each component's current `(w, τs, τl)`, no shared school kernel; seeded from the kernel fit by `sync_kernel_anchor.py`); `gravity_fixed` (an optional list of willingness keys held at their `gravity_ref` value and **excluded from the tuned vector** — for pinning individual params; unknown keys warn and are ignored, default `[]`); `K_prior_std` (per-component `res`/`commute`/`retail`/`school_primary`/`school_postprimary`/`school_tertiary`, default 0.5 — each school level has its own so they can be adjusted independently) sets the width of the **generation-anchored K-prior** `Σ_c (K_c−1)²/σ_c²` inside `solve_scales` — softly pulls each component scale toward the generation value 1 (magnitude anchor + degeneracy break in one). The anchor is fixed at 1 in code (not a config knob); only the widths σ_c are configurable. The `doubly_constrained` list (which components are attraction-constrained via Furness; default `[]` ⇒ singly everywhere — see "Doubly-constrained option") and `furness_max_sweeps` (the approximate-balancing warm-sweep budget `k`, default 12) are carried into `tuned_params.json` by `reset_gravity_params.py`. `gamma_coupling_scale` and `lambda` are unused dead keys left in the file. External demand layers (workplace, retail_spaces, per-level school_demand, commute producers, per-level school producers) are measured per zone in `census_zones.json` — no external scale factors. |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. **Imports `model.COUNT_SITES` as the single source of truth for site geometry** (`SITE_MAP`) — it stamps each site's `node`/`links` from COUNT_SITES into the output. **Re-run when the ODS file OR a COUNT_SITES site location changes** (otherwise `official_hourly.json`, which the tuner reads, drifts stale from COUNT_SITES). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. Loads manual link overrides from `data/manual_link_overrides.json`. After every new link assignment, validates each non-null count direction against the directed graph; raises `ValueError` if the edge doesn't exist. |
| `analysis/manual_assign_link.py` | CLI tool to manually assign a session to a specific directed link, bypassing GPS snap. Usage: `python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>`. Validates both nodes exist and checks count-edge consistency. Writes to `data/manual_link_overrides.json` and patches `counts_processed.json` directly. After correcting an assignment, re-run `aggregate_counts.py` then `tune_assignment.py`. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s`. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Six-component, production-constrained model:** gravity flows split into residential (`flow_res`), commute (`flow_commute`), retail (`flow_retail`), and **three school levels** (`flow_school_primary`/`_postprimary`/`_tertiary`), each production-constrained (optionally also attraction-constrained — see "Doubly-constrained option"). Tunes **18 gravity params** — each component's **double-exp willingness** `(w, τs, τl)` in the mode-substitution × willingness kernel `f(c)=driveShare(equiv_miles(c))·[w·exp(−c/τs)+(1−w)·exp(−c/τl)]` (`kernel: modesub_double`; the `driveShare` rise is shared/empirical, not fit; flat keys `<comp>_taus/_taul/_w`): res, commute, retail **and the three fully-independent school levels** (school_primary/postprimary/tertiary — each its OWN kernel, **no shared `τ_school`**; their distinct distributions come from both the per-level data and per-level shape params). External zone values are fixed from census data and are not tuned. Producer weights are scaled to vehicle-driver trips/day via `model.compute_generation_scales` (generation pinning ⇒ each `K_c ≈ 1`). **Inner calibration = direct-K convex scale solve (`solve_scales`, generic over N components):** the temporal fractions (per-level for school) are **pinned at the NTS profile** (never tuned), so with `f` fixed each prediction is linear in the six scales `(K_res, K_commute, K_retail, K_primary, K_postprimary, K_tertiary)` and the inner objective (Gaussian WLS + Poisson identity-link deviance + generation-anchored K-prior) is **convex**, solved by a damped-Newton + line-search step — **monotone, no K-collapse, no best-iterate hack**. `run_assignment` calls `model.constrained_od_flows` and scatters via the probit routing incidence. **Observed-link scatter restriction (tuner-only):** the objective reads modelled flow on only the ~230 observed links, so `run_assignment` scatters just the incidence entries landing on those links (≈32% of the ~62M), precomputed once into a compact observed-link space — bit-identical results, ~3× faster per eval (`build_assignment.py` keeps the full scatter for the map). Once the compact arrays are built the tuner frees the full-incidence cache arrays (`pair_idx`/`link_idx`/`link_weight`) to keep steady-state memory low. **Performance:** a full singly Powell run is a heavy, long-running pass (many evals) — not a quick verify. **With double-constraint active** (`doubly_constrained` non-empty in `tuned_params.json`) the cost is substantially higher: the **first eval cold-seeds the Furness `b`-cache** (all doubly legs; short-kernel school legs capped + warned), then later evals run the warm `k`-sweep at a multiple of the singly per-eval cost, so a full doubly re-tune is an overnight-scale run — see "Doubly-constrained option" in Model Design for the mechanism and the timing caveat. **Diagnostics (env-gated, no optimization, no writes, then exit):** `CALIBRATE_PROBE=1` reports the post-calibrate residual global scale λ at the start params. `SWEEP=<component>` (one of the six willingness components — `res`/`commute`/`retail`/`school_primary`/`school_postprimary`/`school_tertiary`) sweeps that component's **tail scale `τl`** over a sane grid **with all other kernels frozen at the start params**, reporting its solved `K` and the resulting χ²/N per cell — it answers only "does this component's K collapse / which way does its willingness tail lean", and its χ²/N is **conditional on the frozen others (not a joint fit, never quote as an achievable χ²/N)**. (Under the mode-substitution kernel residential is tolerated — `K_res` O(1) — at short res willingness, collapsing only at a long tail; commute is weakly identified — flat χ²/N, `K` never collapses.) |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` and `reports/slot_pulls_{id}.png`. Echoes the labels the tuner stored in history, so street names appear only after a fresh tune run regenerates `tuning_history.jsonl`. |
| `analysis/diagnose_per_capita.py` | **Fit diagnostic (read-only, no writes).** Sanity-checks the deployed model's trip generation as **true trips-per-capita for core residents**: daily car/van-driver trips whose HOME end lies in the core, per component and overall, divided by core **population** (never producers/attractors — a *true* per-capita), for comparison to travel-survey headline rates (e.g. TSNI trips/person/day by purpose). Mirrors `build_assignment.py`'s deployed setup — `model.constrained_od_flows` (honouring `doubly_constrained`) then × `K_c` × `W_c` (`model.aadt_weights`) for true daily trips, with the deployed generation scales + self-terms. **Home-anchored:** `produced` = trips with origin in core (production side), `received` = trips with destination in core (attraction side). Two-leg components (commute/retail/school — doubly-constrained) → `total = produced + received` (the two distinct legs; `produced ≈ received` since attraction is pinned). Residential is a single symmetric pop↔pop field (singly-constrained, no home/activity distinction), so summing the two would double-count the internal↔internal interaction → `total = (produced + received)/2 = CC + (CE+EC)/2` (full weight internal↔internal, half weight cross-boundary); `produced` is pinned (= ρ·K) while `received` (attraction) is **free**, so `received/produced ≠ 1` is the single-constraint signal. Core = the explicit **`internal_node_ids`** set (as `build_external_links.py`), so **fully portable** across CENTRE/radius changes (no hardcoded IDs). Uses the `return_legs=True` option of `constrained_od_flows`. Run from the repo root; prints a per-component table (`produced`/`received`/`total`/`ρ·K`/`ρ`/`K_c`). Also reports, for the **same** trip set, the trip-**length** distribution per category + pooled in **equivalent miles** (`equiv_miles(od_dist)`) — trips/day, per-capita (ties back to `total`), weighted **average** (mean) + median length, mean journey minutes, and mile-band shares — and writes a per-category TLD plot to `reports/diagnose_tld.png` (residential uses the same M-normalisation in the TLD). Writes only that PNG (no model artifacts). |
| `analysis/equiv_miles.py` | `equiv_miles(t_seconds)` → equivalent road distance in **miles** for a journey time (the seconds↔miles bridge for comparing the time-based kernel to mile-based NTS/TSNI trip-length distributions). Closed form `exp(C0+C1·ln t+C2·(ln t)²)` (zero file I/O, scalar/numpy-array, monotone) fit from the Google cache best-route `g_dur`/`g_dist` (free-flow, 989 pts, ~10% scatter floor). `--fit` re-derives the constants + writes `reports/equiv_miles.png` and prints the refreshed constants to paste back (the module constants are the single source of truth — no coefficients file is written, so nothing can drift out of sync). **Wired into the kernel** via `model._modesub_kernel` (the seconds→miles input to `driveShare`); a **local placeholder for a national average speed** — swappable through this one function. |
| `analysis/driveshare.py` | `driveshare(d_miles, component)` → vehicle-driver **mode share** at trip length `d` miles for a gravity component — the empirical 0→peak rise of the production kernel (`driveShare(equiv_miles(c))` in `model._modesub_kernel`). **Per-component** curves (the short-range walk↔drive substitution genuinely differs by purpose): closed form `PLATEAU_c·(1−exp(−(d/D0_c)^K_c))`, authoritative constants in `CURVES` — commute `(0.694,1.287,0.989)` slow rise/high plateau, retail `(0.549,0.871,1.605)` fast/low, res `(0.615,0.741,1.397)` earliest, and the **three school levels** `school_primary` `(0.847,0.913,2.155)` / `school_postprimary` `(0.427,0.990,2.307)` / `school_tertiary` `(0.320,1.697,1.998)`. `component` is **required** — no shared/legacy fallback (zero file I/O on import, scalar/numpy, monotone, `driveshare(0)=0`). Adds **zero tuned params** (driveShare stays empirical; only the willingness `(w, τs, τl)` is tuned) and de-confounds the willingness; each `PLATEAU_c` **cancels per-component** in the production constraint, so only the rise shape is load-bearing. `--fit` re-derives from the **NTS trip-level microdata** (SN 5340, via `nts_microdata`) using the shared `purpose_mapping.B01_COMPONENT` scheme: survey-weighted **binomial MLE** of the share form on trip records (each trip a Bernoulli outcome at its actual `TripDisIncSW`, weight `JJXSC×W5`), just-walk dropped exactly at record level (`TripPurpose_B01ID==17`), fit capped at 25 mi (`D0,K` stable to the cap) — prints constants + writes `reports/driveshare.png`. Non-school = **driver** share (modes {3,5,12}), 2023/24. **School = per-trip by-car share of the CHILD's own Education trip** (purpose 4; by-car = child driver OR passenger, motorcycle, taxi = modes {3,4,5,12}) — because the student is usually not the driver; level via the **same** `derive_school_generation` age→level machinery (5-10 primary, 11-15 post-primary, 16-18 DfE split, 19+ FT tertiary), tertiary pooled ex-COVID (2013-19+23-24) for sample size. A share is ride-share-invariant, so no car-sharing correction (that lives in generation's magnitude). The school run is non-monotonic (walk→car→bus at range) — the saturating form fits the rise; the bus tail is left to willingness + plateau-cancel. |
| `analysis/trip_length_dist.py` | **(TLD numerator for kernel anchoring).** Builds the six per-component **car trip-length distributions** in **miles** from the NTS microdata → `analysis/trip_length_distributions.json`. Body from the trip table (`Σ JJXSC×W5`, car-driver modes {3,5,12}); the **≥50 mi tail is LDJ-boosted** — the long-distance-journey table (`ldj`), **total-pinned to the trip-table ≥50 mi fraction** so LDJ only refines the within-tail *shape* + ~2–2.5× effective-n (contemporaneous ⇒ no year-pool drift). LDJ handling verified: `LDJPurpose_B01ID`≡`TripPurpose_B01ID` (direct `B01_COMPONENT`), LDJ car+van get a driver-fraction correction measured from the trip table, tranche-(a) [`TripID` present] reproduces the trip-table ≥50 mi mass (ratio ≈1.0), period factor ρ≈0.5. School three levels on the **child's-trip basis** (by-car {3,4,5,12}, level by age via `driveshare._load_school`; tertiary ex-COVID-pooled; no LDJ — school ≥50 mi negligible). Each component carries per-bin `share`/`density`/`eff_n`; `components[c]["distribution"]` is the recommended TLD. The **numerator of `f = TLD/n(t)`** (consumed by `fit_kernel.py`). |
| `analysis/fit_kernel.py` | **(willingness anchor — the `TLD ÷ n_Ire(t)` divide).** Divides each car TLD by the national geometry `n_Ire(t)` (`data/national_n_of_t.json`, from `build_n_of_t.py`) and divides out the fixed `driveshare` rise to recover the **willingness** `W(c)=[TLD/n]/driveshare`, then fits single- and double-exponential forms → `analysis/kernel_fit.json` + `reports/kernel_fit.png`. **Axis = OSRM seconds** (n(t)'s + the model's cost axis); the miles TLD is mapped onto it via `equiv_miles` with the density Jacobian `dd/dt`. **Consumes `national_n_of_t.json` — re-run it to re-check/anchor the kernel after any n(t) change.** **Artifact only** (no `model.py` wiring, no re-tune; the `f≠TLD` philosophy — the divide is a default *shape* to anchor, refined by local counts). Finding: willingness is two-scale — a robust fast head `τ_s`≈7–13 min + a **heavier-than-single-exp tail** (single-exp misfits, weighted-log-resid ~1.2→~0.3 with double). **The tail `τ_l` is QUALITATIVE-ONLY because `n_Ire(t)` is v1 UNCONSTRAINED (no `1/D_i`)** — the production constraint is what suppresses far trips, so dividing by an unconstrained n mis-attributes that to willingness and inflates the tail; **building the constrained (`1/D_i`) n(t) is exactly what firms it up.** Head `τ_s` is trustworthy. School per-level τ (235/338/1143 s): primary≈post-primary, tertiary confounded by its distinctive big-city (university) geometry — weak evidence on the shared-τ_school design. |
| `analysis/iterate_kernel.py` | **(constrained kernel — the production-constraint iteration that firms up `fit_kernel`'s tail).** Fixed-point iteration replacing the *unconstrained* `n_Ire(t)` with the deployed model's **constraint geometry** `Ñ`, per component, matching `model.constrained_od_flows`: **singly** components → `Ñ=Σ_{i,j}(P_i/D_i)A_j·δ` (`D_i=Σ_k A_kf`); **doubly** components → `Ñ=Σ(a_iO_i)(b_jD_j)·δ` with Furness `a,b` (`_furness_ab`, mirrors `model._furness`), then `f→Ñ[f]→fit_double([TLD/Ñ]/driveshare)→f`. **Which components are doubly is read live from `simulation/tuner_config.json` `doubly_constrained`** (single source of truth — a config change + re-run flips them, no code edit; res always singly). **Route-once-iterate-cheap:** one cached routing pass per purpose (area-level accessibility `M=128` for `a` (group-by-origin) + `b` (invert-by-dest) + the origin/dest-area-tagged stratified density for `Ñ`, resumable to `data/_kernel_iter_cache_<p>.npz`), then Phase B (Furness + reweight + fit) is seconds. Reuses `build_n_of_t.py` (sampler + point cache) + `fit_kernel.py` (willingness divide + fit). **Six independent kernels — schools per-level.** Multi-start fit avoids the `w→1` degenerate basin. → `analysis/kernel_fit_constrained.json` (per-component `constraint`, `converged`/`tail_weakly_identified` flags) + `reports/kernel_fit_constrained.png`. **Artifact only** (no `model.py` wiring, no re-tune). **Findings:** singly — res/commute/retail/school_primary converge, `τ_l` shortens (tail de-inflation, commute 4059→1964 s), school post-primary/tertiary tails weakly identified (bimodal `τ_l` → robust median). Doubly — the attraction `b` is **MC-sampled (approximate, ~5% margin, `M≈128`)** in the same spirit as the model's own approximate (`furness_max_sweeps`) balancing; it *lengthens* commute's tail (singly 1964 → doubly 3160 s) at a higher fit residual (wrms 0.21→0.59 — the approximate-b + area-level-attraction cost). Still caveated: `n_Ire` not `n_Eng`; outbound leg; finite-island truncation. Needs OSRM up. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. |
| `simulation/reset_gravity_params.py` | Reset the gravity params in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`: every `gravity_ref` shape param (the 18 `<comp>_taus/_taul/_w` double-exp willingness keys — iterates `gravity_ref`, so no rename drift), sets `kernel: modesub_double`, plus the six scales `K_res/K_commute/K_retail/K_primary/K_postprimary/K_tertiary` → 1.0. Strips dead/legacy keys (`K`, `K_biz`, the pre-split `K_sch` + `slot_fracs_school`, `W_BIZ`, `W_SCHOOL`, `P_biz`, `ALPHA_biz`, the rational-kernel tail exponents, the Tanner peak/shape params, `MU`, `SIGMA`). External params and the per-level `slot_fracs_*` are preserved. |
| `simulation/sync_kernel_anchor.py` | Seeds the tuner's willingness **anchor** from the kernel fit. Reads the per-component double-exp willingness fit (`analysis/kernel_fit.json`, from `analysis/fit_kernel.py`) and writes the 18 flat willingness keys (`<comp>_taus/_taul/_w`) into `tuner_config.json`'s `gravity_ref` — the anchor / start point + L2-pull target for `tune_assignment.py`. Only the shape `{w, τs, τl}` is carried (the fit's amplitude `A` is absorbed by K in the production constraint). Also (re)initialises `gravity_lambda` to a light uniform value on those 18 keys and strips dead single-exp `TAU_*`/`THETA` keys. **Run before `reset_gravity_params.py`** when re-anchoring to a fresh kernel fit; `--lambda` sets the reg weight. (Patch point: repoint `ANCHOR_FILE`/`_read_anchor` to the constrained-`n_Ire` iteration's params when available — everything downstream keys off `model.willingness_keys()`.) |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** Per-component temporal-**shape** profiles (168 rows = 7 days × 24 h): `mean_fraction_res`, `mean_fraction_commute`, `mean_fraction_retail`, the three per-level school shapes `mean_fraction_school_primary`/`_postprimary`/`_tertiary`, plus the aggregate `mean_fraction` (used elsewhere, e.g. `ingest_counts.py`). all six component columns are **car-specific, derived from the NTS microdata** via `analysis/derive_component_profiles.py` (res/commute/retail = joint car-driver dow×hour; school = per-bin escort regression + age→level self-drive); `mean_fraction`/`std_fraction` are inputs it preserves. Each component column is an **independent shape** normalised so its day-weighted daily sum `W_c = 1` (⇒ each column sums to 7.0 over the 168 rows) — **no** aggregate-partition constraint (magnitude/split is generation's job, see "Generation pinning"). Re-run `derive_component_profiles.py` when the NTS microdata or the purpose mapping change. |
| `analysis/derive_component_profiles.py` | Derives the **res/commute/retail** hourly-**shape** columns of `hourly_fractions.csv` **directly from the NTS trip-level microdata** (car-specific). Per component (shared `purpose_mapping.B01_COMPONENT`) it counts vehicle-driver trips (MainMode_B04ID {3,5,12}, weight `JJXSC×W5`) by **(day_type, hour)** — the joint car-driver start-time distribution — and normalises each column to Σ=7. Calendar day-of-week from `TravelWeekDay_B01ID` (Day table, joined via `DayID`; **not** `TravDay`, which is the diary-day position). Retires the old NTS0502a/0504b (all-mode) derivation and its `PROXY_0502` merged-purpose proxies, separable `V(dow)×H(hour)`, all-mode weekend fallback, Bayes-flip, and the `ρ`-weighting dependency (real trip counts weight by volume) — so the weekend shape is genuinely per-component (retail Sat-heavy, res weekend-heavy, commute weekend-light). day_type collapses weekday to Mon–Fri (the model averages them), written identically into the five weekday rows. **The profile is a pinned χ² input with no propagated uncertainty**, so every (day_type,hour) bin needs adequate stats: **adaptive per-bin year pooling** expands each bin's window through the ex-COVID tiers (2023/24 → +2018/19 → +2013–17) only until its unweighted count hits `MIN_BIN_N=100`, then stops (well-sampled daytime/weekday bins stay on 2023/24 at zero drift; only thin weekend-night bins pool wider), per-year normalised (÷ the bin's year-count) so windows are comparable (valid: weekday-peak drift ≈1 pp). A **light ±1 within-night smoothing** (hours 0–5, per day_type, night-block total preserved, daytime untouched) removes the residual noise on the intrinsically-thin res weekend-night bins → every bin ≲15% (obs σ floor). NaN start-hour/day-of-week/weight trips dropped **explicitly**. **School** is a separate path (`_school_level_grids`): a school car trip is a COUNT of when the escorting car leaves, so escort **ride-sharing matters** (unlike the driveshare *share*) — per level the car-departure timing = **escort** (individual escort trips carry no child level, so the **same household regression as the school generation rates** attributes them — run **per (day_type,hour) bin**, β_L(dt,h) is ride-share-correct and level-resolved; primary-heavy households drive at primary bell times; off-peak β<0 clamped) + **self-drive** (purpose-4 veh, age→level via `_level_shares` incl. the DfE 16-18 split), blended by the generation escort/self-drive split (`school_generation_rates.json`). Recovers real level differences (post-primary AM peak earlier than primary; tertiary spread/college-like, no sharp PM peak). **Uniform 9yr ex-COVID pool** (not adaptive — the regression runs on the whole set; thin, weekday-peak-concentrated; 2yr-vs-9yr stable for primary/post-primary, tertiary needs it), same night smoothing. `mean_fraction`/`std_fraction` preserved. This **completes the temporal migration** — all six columns microdata-derived. |
| `analysis/purpose_mapping.py` | Single source of truth for the NTS purpose→gravity-component mapping, imported by **both** `derive_generation_rates.py` and `derive_component_profiles.py` so generation and temporal use the same split. The mapping is **`B01_COMPONENT`** (23-cat `TripPurpose_B01ID`→component; rule: res iff endpoint is a home, else routed by land-use) + `B01_EXCLUDE`, carrying the judgment allocations (Business/Personal-business→retail; escorts by destination; leisure split by endpoint). |
| `analysis/derive_school_generation.py` | **Tracked output `analysis/school_generation_rates.json`.** Derives the **per-student** school generation rate (vehicle-driver trips per full-time student per day) for each level from the **NTS microdata** (via `nts_microdata`) + the **DfE England age→level split** (`data/participation-in-education-training-and-employment-age-16-to-21_2025/…allinsts…csv`, gitignored). Two parts: **escort** (adult drives child; TripPurpose_B01ID 21) via a **household regression** of escort-vehicle trips on per-household primary/secondary/tertiary/pre-school student counts (β = per-student marginal, handles ride-sharing); **self-drive** (student's own vehicle trip; purpose 4) combined per level from the per-age rate. NTS gives actual age but not level, so 16–18-yos are split secondary(sixth form)/tertiary(FE) by DfE full-time enrolment shares — **6th form = secondary, FE = tertiary**, DfE academic age → NTS actual age via the September-cutoff **half-year** offset, FT students only (matches `census_school_producers` tertiary). Firmed age-18 self-drive = pooled 2013–24 ex-COVID (`FIRM_AGE18`). Also emits the **pre-school escort per-capita magnitude** folded into retail (a documented fudge — no pre-school producers). Fudges/caveats (escort transfers well; self-drive is the least-transferable piece and dominates the small tertiary rate; England behaviour) are in the module docstring. Result: primary 0.379, post-primary 0.327, tertiary 0.049 trips/FT-student/day. Re-run when the microdata or DfE file change. |
| `analysis/derive_generation_rates.py` | **Tracked output `analysis/generation_rates.json`.** Derives per-component **per-capita** vehicle-driver trips/person/day from the **NTS microdata** (via `nts_microdata`, 2023/24, MainMode_B04ID {3,5,12}, Σ(JJXSC×W5)÷Σ W2 persons ÷7) using the 23-cat `purpose_mapping.B01_COMPONENT` mapping → commute 0.200, retail 0.581, res 0.159. **School is per-capita converted from the per-student rates** (`school_generation_rates.json`) × (island `node_school_producers_<level>`/`node_population`) so ρ_school/k_students recovers the per-student rate exactly → `school_primary` 0.0381, `school_postprimary` 0.0276, `school_tertiary` 0.0025 (tertiary far below the old enrolment-share placeholder — that overstated tertiary ~7×). Retail additionally absorbs the pre-school escort fudge. **Every rate is per-capita** (encodes an island total; the producer/attractor layer only distributes spatially). Reads `node_weights_reduced.json` for the school students/pop ratio, so it is all-Ireland-specific. Re-run when the microdata, the mapping, or the school rates change. Consumed via `model.load_generation_rates` / `model.compute_generation_scales`. |
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
`simulation/node_weights.json` keys: `node_population`, `node_workplace` (all place-of-work jobs — feeds the external retail-spaces fallback + the map; **not** the commute attractor), `node_commute_attractor` (car-commute jobs = commute attractor, internal + external; from `census_attractor.py`), `node_retail_spaces` (estimated retail parking spaces = retail attractor, internal + external), the three `node_school_demand_<level>` (per-level school places = school attractors, internal + external; from `school_attractor.py`), `node_commute_producers` (census car-driving commuters; from `census_supply.py`) and the three `node_school_producers_<level>` (census resident students by level; from `census_school_producers.py`) = trip producers, internal + external, `boundary_node_ids` (auto-detected from core polygon). There is no combined `node_business_demand` layer — the commute (`node_commute_attractor`) and retail (`node_retail_spaces`) attractors are separate. External node entries (census-area-code string IDs, e.g. `"N21000219"`) are included alongside internal OSM node IDs.
**Node ID scheme:** `build_network.py` relabels consolidated graph nodes to stable OSM IDs after junction consolidation. A consolidated junction (multiple OSM nodes merged) gets `min(osmid_original)` as its ID; a non-merged node gets `int(osmid_original)`. All road node IDs are therefore genuine OSM node IDs (in the hundreds of millions) and stable across graph regenerations. External census nodes use their **census-area-code string IDs** (SDZ/DEA 2021 codes, e.g. `"N21000219"`) — these are the `id` values in `census_zones.json`, *not* small integers (downstream code that consumes `external_links.json`/`census_zones.json` must treat external IDs as strings — `build_paths.py` does). Road node IDs are ints, external node IDs are strings; `node_to_idx` mixes both.
`simulation/newtownards_flows.json` — combined flows plus `flows_res`/`flows_commute`/`flows_retail` and the three `flows_school_primary`/`_postprimary`/`_tertiary` keys (W-weighted directed AADT) when the multi-component params are active, plus an `aadt_weights` block.
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
`data/island_opportunity_table.csv` — committed; output of `build_opportunity_table.py`. One row per island small area (NI DZ + RoI SA, ~22.7k) with producer/attractor masses + WGS84 centroid; the frozen input to `analysis/build_n_of_t.py`. Same estimators as `census_zones.json` (values aggregate to its external nodes). Re-run only if the census/parking/school estimators change.
`data/national_n_of_t.json` — committed/frozen; output of `build_n_of_t.py --stratified`. Per-purpose (res/commute/retail + 3 school levels) reconstructed `n_Ire(t)` over 30 s bins to 240 min + overflow (481 bins), unconstrained + outbound-leg (v1). Re-run only if the opportunity table, OSRM profile, or sampler change.
`data/external_links.json` — committed; output of `build_external_links.py`. Contains OSRM-derived X↔B links, boundary shortcuts, and through-route allowlist. Re-run `build_external_links.py` when boundary nodes change or OSRM data is updated.
`data/external_intra_times.json` — committed; output of `build_intra_times.py`. Per external zone, a **per-component** (res/commute/retail + 3 school levels) **mass-weighted** intra-zonal time histogram (`{t:[bin-centre s], w:[weights Σ=1]}`) for the production-suppression self-term (`model.load_self_terms` → `constrained_od_flows`; the model applies the tuned kernel to the bin centres). Committed so the model runs without re-querying OSRM. Re-run `build_intra_times.py` only when external zones change.
`data/manual_link_overrides.json` — committed so manual assignments survive a wipe of `counts_processed.json`.
`simulation/tuner_config.json` — committed as source config (gitignore exception).
`analysis/hourly_fractions.csv` — committed as source data (single authoritative version).
`analysis/generation_rates.json` — committed; per-component vehicle-driver trips/person/day (output of `derive_generation_rates.py`). Source data for generation pinning (`model.compute_generation_scales`).
`analysis/trip_length_distributions.json` — committed; output of `trip_length_dist.py`. Six per-component car trip-length distributions (miles bins: `share`/`density`/`eff_n`; non-school tail LDJ-boosted). The numerator of `f=TLD/n(t)`; consumed by `fit_kernel.py`.
`analysis/kernel_fit.json` — committed; output of `fit_kernel.py`. Per-component single- + double-exp willingness fits from the `TLD ÷ n_Ire(t)` divide; `_meta` carries the caveats (tail `τ_l` qualitative-only pending the constrained `1/D_i` n(t); the tuned-τ reference is stale/pre-mode-sub-retune). Artifact — not wired into the model.
`analysis/kernel_fit_constrained.json` — committed; output of `iterate_kernel.py`. The six double-exp willingness kernels after the `1/D_i` production-constraint iteration (constrained `τ_l`, de-inflated vs `kernel_fit.json`); per-component `converged` / `tail_weakly_identified` flags + iteration trace, `_meta` carries the method + caveats. Artifact — not wired into the model.

### Large reference data (gitignored, kept locally only)
`data/*.ods`, `data/*.xlsx`, boundary GeoJSON files — too large to commit; keep local copies.
Currently present:
- `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods` — used by `parse_official_hourly.py`.
- `data/nts0502.ods` / `data/nts0504.ods` — DfT NTS Tables NTS0502a / NTS0504b (weekday trip start times / trips by day×purpose). **No longer used** — the temporal derivation is now microdata-based; retained on disk only.
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
the generation/distribution conflation). The six components (residential / commute / retail /
school-primary / school-postprimary / school-tertiary), their producer/attractor weights and
per-component kernels are detailed in "Six-component flow decomposition" below (the three school
levels share one kernel). See the agent memory note `project_production_constrained_gravity`.

**Mode-substitution × willingness kernel** (`model._modesub_kernel(d, wparams, component)`; per component, three params `wparams = (w, τs, τl)`):
`f(c) = driveShare(equiv_miles(c)) · [w·exp(−c/τs) + (1−w)·exp(−c/τl)]`, `c` = OSRM travel time (seconds).

The kernel is **decomposed into two physically-distinct factors** instead of one free-fit shape:
- **Mode substitution** (the 0→peak rise) — `driveShare(equiv_miles(c))`: short trips are *walked,
  not driven*, so car demand is suppressed at short cost. **Empirically fixed** (not tuned) and
  **per-component** — the microdata shows the rise genuinely differs by purpose (commute rises
  slowest to the highest plateau, res earliest), so each component has its own curve
  (`analysis/driveshare.py` `CURVES`, e.g. commute `0.694·(1−exp(−(d/1.287)^0.989))`); `equiv_miles(c)`
  converts OSRM seconds→miles (`analysis/equiv_miles.py`). `f(0)=0` holds (`equiv_miles→0 ⇒
  driveShare→0`). *(All six curves — incl. the three per-level school curves — are in `CURVES` and
  wired through `model._modesub_kernel(d, wparams, component)`.)*
- **Willingness to travel** (the decay) — `w·exp(−c/τs) + (1−w)·exp(−c/τl)`: a **double
  exponential** — a fast head (weight `w`, scale `τs`) plus a heavier tail (`1−w`, `τl > τs`); the
  `TLD ÷ n_Ire` divide (`fit_kernel.py`) showed a single exponential is too light-tailed. Monotone;
  closer wins given equal opportunity. A **purpose** property → per-component, **three tuned params**
  `(w, τs, τl)` (τ's in seconds); the willingness amplitude is absorbed by K.

The peak therefore **emerges** from rise×decay rather than being a tuned "typical time" — under
production constraints (which already own trip *volume*) the kernel is pure *distribution*, so a tuned
20-min "peak" would be nonsensical. **18 tuned params** — per-component `(w, τs, τl)` for
res/commute/retail **and the three fully-independent school levels** (school_primary/postprimary/
tertiary, each its OWN kernel — **no shared `τ_school`**; flat keys `<comp>_taus/_taul/_w`). `τs`/`τl`
are the transferable / TLD-anchorable willingness scales. This replaced the Tanner kernel (whose
per-component `P_c`/`BETA_c` free-fit degenerated: `K_res` collapsing to ~0 at any sane kernel); under
the mode-substitution split `K_res` no longer collapses. The shared `driveShare` **plateau cancels** in the production
constraint (`a_j/D_i`), so it carries no magnitude — magnitude is owned entirely by generation pinning.
The kernel is **model-layer** (applied to `od_dist` in `constrained_od_flows`, not routing) — changing
it needs no `build_paths` rebuild, only a re-tune.

`equiv_miles` is currently fit to the local Google routing cache (length-dependent, ~14 mph short →
~33 mph long) as a **placeholder for a national average speed**; it is swappable via that one function.

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

**Demand:** all layers are per-zone (Census 2021 DZ/SA aggregated to SDZ/DEA/ED/LEA), kept as **separate layers** (nothing is summed — each component uses its own producer/attractor, see "Six-component flow decomposition"). Producers: `commute_producers` (car-driving resident commuters) and the three `school_producers_<level>` (resident students by level), all census-derived, plus population. Attractors: `commute_attractor` (car-commute jobs per zone via `census_attractor.py` — commute attractor; `workplace` = all jobs is retained for the retail-spaces fallback + map, not as an attractor), `retail_spaces` (parking spaces within the zone via `parking_demand.parking_spaces` — retail attractor), and the three `school_demand_<level>` (per-zone school enrolment split by level via `school_attractor.py` — school attractors).

**Connectivity:** boundary nodes are all internal nodes with at least one road edge crossing the core polygon boundary. OSRM-derived directed edges connect each external node to its valid boundary nodes; the gravity model path distance for any external→internal pair is the OSRM leg plus the internal shortest path, computed automatically by Dijkstra on the augmented graph.

**No tuning of external zone values.** Pop/wp come directly from census data and are not adjusted by the optimizer (there are no hand-crafted city configs, ref values, or damping factors).

### Through routes
External→external OD pairs are allowed only for pairs whose OSRM route passes through at least one boundary node (i.e. genuinely transits the core area). This allowlist is auto-generated by `build_external_links.py` and stored in `data/external_links.json`. Changing the allowlist requires re-running `build_external_links.py` and `build_paths.py` (the paths cache must be rebuilt).

### Six-component flow decomposition
The gravity OD flows are split into six **production-constrained** spatial components at each
tuner evaluation (per-pair pre-K flows from `model.constrained_od_flows`, scattered onto links;
`p^c_i`/`a^c_j` = producing/attracting weight, `D^c_i = Σ_k a^c_k·f_c(d_ik)`). The five non-res
components are independent clones: a symmetric pop↔activity split, each leg per-origin-normalised,
with **no weight parameter and no self/cross term**:

- **Residential** (`flow_res`, res willingness kernel): `T^res_ij = pop_i·pop_j·f_res/D^res,pop_i` —
  pop×pop trips. Single leg (i→j and j→i are separate OD pairs, so both directions are covered).
- **Commute** (`flow_commute`, commute willingness kernel): home→work producer =
  `commute_producers` (car-driving commuters), attractor = `commute_attractor` (car-commute jobs);
  return work→home producer = `commute_attractor`, attractor = `commute_producers` (the returning
  commuters' homes are distributed by resident-commuter count, **not** raw population) —
  `f_com·( commprod_i·cattr_j/D^com,cattr_i + cattr_i·commprod_j/D^com,commprod_i )`. Both sides are
  car-specific per zone (the model assigns car flow), sharpening the spatial distribution while
  generation pinning keeps the magnitude NTS-pinned. The symmetric producer↔attractor round-trip
  makes commute independent of the population layer.
- **Retail** (`flow_retail`, retail willingness kernel): home→shop producer = pop, attractor
  = `retail_spaces`; return shop→home producer = `retail_spaces`, attractor = pop —
  `f_ret·( pop_i·ret_j/D^ret,ret_i + ret_i·pop_j/D^ret,pop_i )`.
- **School — three independent levels** (`flow_school_primary`/`_postprimary`/`_tertiary`), each a
  full production-constrained component with its own producer, attractor, scale `K_<level>`,
  temporal shape **and its OWN willingness kernel** (per-level `(w, τs, τl)` — **no shared `τ_school`**;
  the distinct distributions come from both the per-level data and the per-level shape params). For each level: home→school
  producer = `school_producers_<level>` (census resident students by level, `census_school_producers.py`),
  attractor = `school_demand_<level>` (per-level enrolment, `school_attractor.py`); return school→home
  producer = `school_demand_<level>`, attractor = `school_producers_<level>` (returning students land
  where that level's resident students live, **not** by raw population) —
  `f_sch·( schoolprod^L_i·school^L_j/D^L,sch_i + school^L_i·schoolprod^L_j/D^L,schoolprod_i )`. The
  levels are fully decoupled (a change to one cannot perturb another); each level's producer↔attractor
  round-trip is independent of the population layer. External `school_demand_<level>` and `school_producers_<level>`
  are both populated per zone, so external school trips are retained intra-zonally via the self-term
  (below) rather than dumping into the core.

Per-leg producer weights are scaled to **vehicle-driver trips/day** by `gen_scale` so each `K_c ≈ 1`
(see "Generation pinning").

**External intra-zonal self-term (denominator-only, mass-weighted, per-component).** Each per-origin
denominator `D^c_i = Σ_k a^c_k·f_c(d_ik)` runs over *other* zones; collapsing an external zone to one
centroid drops its `k=i` diagonal (its intra-zonal trips), so `D^c_i` is too small and the external
origin's fixed budget over-allocates to the rest of the network (worst for large, isolated, far zones).
`build_intra_times.py` measures the diagonal as `a^c_i·S^c_i`, where **`S^c_i` is the producer×attractor
mass-weighted mean kernel over intra-zonal trips** — origins ∝ producer, destinations ∝ attractor (real
POIs for retail/school), road-snapped, within the zone → a per-component weighted time histogram in
`data/external_intra_times.json`; `model.constrained_od_flows` then adds `a^c_i·Σ_bin w·F_c(t_bin)` to
each denominator. This captures **clustering**: a sparse rural zone whose people and jobs both sit in the
same villages reads short (strong suppression) while a genuinely spread zone reads long — the old
uniform-in-polygon single-average was clustering-blind (it sampled empty fields, so sparse zones
under-suppressed and over-exported, e.g. Cusher/Clogher out-sending larger/closer towns). Because the
`p·a·f` interaction is symmetric, one `S^c` serves both legs of a component, so the previous
leg-asymmetry (0.47-vs-0.54 out/return) dissolves. It is **denominator-only** — no link flow — and
applies to **external zones only** (internal road nodes have no zone area). Effect is **kernel-dependent**
(shorter willingness `τ_c` sharpens it). Wired into `build_assignment.py` and `tune_assignment.py` via
`model.load_self_terms`; absent file ⇒ no self-term. Independent of the paths cache (model-layer) — no
`build_paths` rebuild.

The self-term applies to every component. External `school_demand_<level>`/`school_producers_<level>`
and the `workplace`/`retail_spaces` layers are all populated per zone, so every component's denominator
carries its intra-zonal diagonal and external zones retain school/commute/retail trips locally rather
than dumping them into the core.

**Doubly-constrained (Furness) option — implemented (approximate-balancing), OFF by default.**
`constrained_od_flows(..., doubly_constrained=<set>)` can additionally attraction-constrain any of
commute / retail / the three school levels (residential is held singly by design): a flagged
component's every leg is balanced (Furness/IPF) so BOTH margins hold — `Σ_j T_ij (+ self diagonal)
= gen-scaled producer_i` (production, the absolute magnitude anchor) AND `Σ_i T_ij (+ self diagonal)
∝ attractor_j` (attraction; the attractor's raw scale is normalised so `ΣD = ΣO` over the reachable
support, so only its cross-zone proportions matter). The self-term restores the intra-zonal diagonal
to **both** balancing sums (the `p·a·f` interaction is symmetric, so one histogram serves both
margins). Flow stays **linear in `K_c`** (balancing factors normalise to the raw margins, not `K`),
so the convex direct-K solve is unchanged. **Verified correct** (both margins exact to ~1e-11 vs an
independent dense reference; residential output bit-identical to the singly path).

The active set is the `doubly_constrained` key in `tuned_params.json` (default `[]` from
`tuner_config.json` ⇒ singly-constrained everywhere, zero behaviour change; both `build_assignment.py`
and `tune_assignment.py` honour it). **Why plain IPF alone won't do:** the *real* short-range kernels
make the OD system nearly decoupled, so exact IPF needs **~1000–1400 iterations per leg** (commute
≈36 s/leg) → a full converged eval is minutes and a re-tune would take weeks; warm-starting to
tolerance does **not** help (subdominant eigenvalue ≈ 1) and a looser tolerance buys only ~1.6×.
**The deployed scheme is approximate-balancing** (`furness_max_sweeps`, default 12; `furness_state`
b-cache): each leg's balancing factors are **cached across evals**, the first (cold) eval per leg
converges to seed the cache, and every later eval runs a **fixed `k` warm sweeps** ending on a
row-normalisation — so **production stays exact** and only the attraction margin is approximate
(**<1% at k≈10**, well under count noise, and self-healing as the cache tracks the slowly-drifting
`b`). **Timing:** `tune_assignment.py` uses the approximate path (the `b`-cache persists across the
Powell evals), so the **first eval is slow — it cold-seeds every doubly leg once** (the short-kernel
school legs cap at `max_iter=3000` and print a `[furness …] IPF capped …` **WARNING** — expected, not
an error), then **every later eval runs the fixed warm `k`-sweep at a multiple of the singly per-eval
cost**. Net: **a full doubly-constrained re-tune is an overnight-scale run** — far heavier than a
singly one. The tuner **preserves `doubly_constrained` / `furness_max_sweeps` in the saved
`tuned_params.json`** so the flag survives a tune. `build_assignment.py` leaves `furness_max_sweeps=None`
so the **deployed** flows converge exactly (cold-converges each doubly leg, substantially slower than
a singly build, same school warnings).
**End-to-end verified on a real doubly-constrained tune** (`diagnose_imbalance.py --sides --doubly`):
every flagged component's **internal** imbalance collapses to **0.0%** (0.6–1.9% overall, the
residual being external intra-zonal self-flow), while the singly-held `res` control stays at its
full imbalance. **Caveat:** the short-kernel school legs (`school_postprimary` τs≈90 s) balance
poorly under the fixed sweep budget (attraction residual can be tens of %); it barely matters while
those components are weakly identified / near-off, but a future tune that makes them carry real flow
would need more sweeps or a faster exact solver (Anderson / Newton-CG on the convex balancing dual)
— a possible follow-up, not needed at k≈10 today.

Each component has its own **independent temporal shape** (`hourly_fractions.csv`, normalised so
`W_c = 1` — magnitude/split is generation's job, not the temporal profile's) and scale
(K_res, K_commute, K_retail, K_primary, K_postprimary, K_tertiary).
Predicted count for observation i in slot s:
`pred_i = Σ_c K_c·flow_c·(T/3600)·f_c[s]`  over  c ∈ {res, commute, retail, school_primary, school_postprimary, school_tertiary}.

### Generation pinning (data-based supply)
Producer weights are carried in absolute **vehicle-driver trips/day** so each component's tuned
scale **K_c should land at ≈ 1.0** — a *verification anchor*, not a fit knob (a `K_c` away from 1
diagnoses local car-mobilisation vs the national average, to be refined later). This is a
**model-layer** change (no paths-cache rebuild). Independent of the kernel form (the mode-substitution
× willingness kernel; see memory `project-tanner-kernel-tld`).

**Rates `ρ_c`** (`analysis/generation_rates.json`, written by `analysis/derive_generation_rates.py`
from the **NTS microdata**, 2023/24, vehicle-driver modes MainMode_B04ID {3,5,12} = *car/van driver +
motorcycle + taxi*, via the 23-cat `purpose_mapping.B01_COMPONENT`): commute 0.200, retail 0.581
(incl. the pre-school escort fudge), res 0.159 trips/person/day. **All rates are per-capita** (each
encodes an island total; the producer/attractor layer only distributes spatially). The three school
levels — school_primary 0.0381, school_postprimary 0.0276, school_tertiary 0.0025 — are the
**per-student** escort+self-drive rates (`analysis/derive_school_generation.py`: primary 0.379,
post-primary 0.327, tertiary 0.049 per FT student/day) converted to per-capita by × island
students/pop, so `ρ_school/k_students` recovers the per-student rate exactly. Tertiary is far below
the old enrolment-share placeholder (which overstated it ~7×). Each level's ρ is independent.

**Purpose→component mapping — JUDGMENT ALLOCATIONS (candidate error sources, kept flagged on
purpose).** Organising principle = the attractor each component offers: workplace(jobs)→commute,
retail_spaces(**parking** = all commercial/venue)→retail, school→school, population(**homes**)→res.
The 23-cat microdata mapping (`purpose_mapping.B01_COMPONENT`) resolves most splits from the data
(leisure by **endpoint** — visit-friends-at-home→res, venues→retail; escorts by destination), but
two allocations remain modelling *decisions* the codes don't resolve — the first thing to revisit if
the fit is scrutinised:
- **Business / Other-work → retail** (not commute): commute kept pure home↔own-workplace; business
  visits hit commercial premises (parking), not the workplace-jobs count.
- **Personal business → retail**: services/banks/medical are commercial/parking destinations.

**Per-leg producer scaling (`model.compute_generation_scales`).** res is single-leg (full `ρ_res`;
both directions are separate OD pairs); commute/retail/school (each level) are two-leg, each direction
carrying `ρ_c/2`. The per-producer rate is `r_leg = (ρ_c·share)/k_leg`, applied to the **producer term
only** inside `constrained_od_flows` (the same array's attractor use is scale-invariant — cancels in
`attr_j/D_i`, denominators untouched; `gen_scale=None` ⇒ all 1.0 ⇒ exact prior flows). The anchors
`k = Σ(producer layer)/Σ(population)` are summed **island-wide** from the node weights (the external
census nodes tile the whole island, so they recompute for any CENTRE — the transferability win):
`k_commuters`(commute_producers), `k_jobs`(commute_attractor), `k_retail`(retail_spaces), and — **per
school level** — `k_students_<level>`(school_producers_<level>), `k_enrolment_<level>`(school_demand_<level>);
`k=1` when the producer is population (res, retail outbound). **Independent per-side:** the per-producer rate is applied to each side's
real count, so a job-/retail-rich sub-area generates more activity→home (PM-outflow) than its
residents alone — the open-region asymmetry, fed in the morning by the external zones' outbound legs.
At island scale this equals the balanced total (`Σ_leg = ρ_c·½·Pop`, since island counts = `k·Pop`);
the difference is purely the sub-regional directional balance.

**Generation anchors the K-prior directly:** because each producer weight is in vehicle-driver
trips/day, the value generation pinning predicts is `K_c ≈ 1`, so the inner solve regularises each
scale toward 1 via a **generation-anchored K-prior** (see "Direct-K convex scale solve") rather than a
share prior — the anchor *is* the generation expectation, with nothing to re-derive.

### Direct-K convex scale solve
At each optimizer evaluation the six component scales **(K_res, K_commute, K_retail, K_primary,
K_postprimary, K_tertiary)** are calibrated **directly** by `solve_scales` (generic over N components).
The temporal fractions (per-level for school) are **pinned at the NTS profile**
(`hourly_fractions.csv` `mean_fraction_*`) and never tuned — so every
observation prediction is **linear in the scales**: `pred_i = Σ_c K_c·a^c_i` with
`a^res_i = m_res_i·Th_i·f_res[s_i]` constant (commute/retail/school analogously). The inner objective —
Gaussian WLS over the official obs + Poisson **identity-link** deviance `2·Σ(n·log(n/pred)+pred−n)` over
the walking obs + a generation-anchored K-prior — is therefore **convex over K ≥ 0**, and is solved by a **damped
(Levenberg) Newton step with a backtracking line search on the full objective**: monotone by
construction, so there is **no K-collapse and no best-iterate bookkeeping**. `CALIBRATE_PROBE=1` reports
the residual global scale λ at the start params (≈1 ⇒ K at its optimum).

**Generation-anchored K-prior (magnitude anchor + degeneracy break).** Every active component is
penalised `Σ_c (K_c − 1)²/σ_c²` toward the generation value `1` (the anchor is fixed in code; the
per-component widths `σ_c` come from `tuner_config.json` `K_prior_std`, default 0.5). This does double
duty: it anchors each scale to the generation expectation (`K_c ≈ 1`) **and** pins the otherwise-flat
commute↔retail direction — so no separate share prior is needed, and res is anchored on the same
footing as the other three (not left as a free remainder). It is a clean quadratic in
K (Hessian contribution `+2/σ_c²` on the diagonal, PSD), regularises the inner K-solve only, and is
**not** part of the reported χ². Walking obs mostly fall in slots with `f_school≈0`, so the school-level
scales `K_primary`/`K_postprimary`/`K_tertiary` are pinned almost entirely by the official school-peak
hours; a soft K-prior anchors but does not force them, so a genuine school over-fit still surfaces as a
school `K_<level>` departing from 1 (a diagnostic, not a masked one).

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

`build_assignment.py` uses the multi-component `compute_chi2()` (passing the commute/retail link-flow dicts + the per-level school link-flow dicts) when the multi-component params are present in `tuned_params.json` — a **data-only** chi²/N (pure sum of squared z-scores). Since the tuner's χ²/N is now also pure data-fit (`f` pinned ⇒ no f-prior/coupling penalty), the two surfaces are directly comparable (both read the same observed links; `build_assignment` keeps the full scatter for the map).

**Reading "modelled flow" across reports.** The three reporting surfaces print *different
projections* of the same tuned model — they are not directly comparable line-for-line:
- `build_assignment.py` "Official count sites" block and the `"flows"` values in
  `newtownards_flows.json` → **directed daily AADT** = `Σ_c K_c·flow_c·W_c` over the six components,
  where `W_c` (`model.aadt_weights`) is the day-type-weighted (5·weekday+Sat+Sun)/7 sum of component
  `c`'s hourly fractions. With the **decoupled per-component shapes** (each normalised so `W_c ≈ 1`,
  `Σ_c W_c ≈ 6`), `K_c·flow_c` is already ≈ the component's daily AADT; the `W_c` factor is retained
  because it reads the actual slot fractions (and was load-bearing under the old partition scheme,
  where each `W_c` was a sub-1 share). (The unweighted per-component flows feed `compute_chi2`, which
  applies `f_c` itself — do not double-weight.) `newtownards_flows.json` stores the W-weighted AADT in
  `flows`/`flows_res`/`flows_commute`/`flows_retail`/`flows_school_primary`/`_postprimary`/`_tertiary` plus an `aadt_weights` block.
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

The current kernel is the **mode-substitution × double-exp willingness**
`f(c)=driveShare(equiv_miles(c))·[w·exp(−c/τs)+(1−w)·exp(−c/τl)]` (`kernel: modesub_double`, **18
params** — per-component `(w, τs, τl)` incl. the three independent school levels), replacing the
earlier Tanner and single-exp forms (whose free-fits degenerated — `K_res` collapsing to ~0). The
`driveShare` rise is shared/empirical (not tuned). Six components are tuned with generation pinning
(each `K_c` should land ≈ 1; a departure diagnoses local car-mobilisation or a weakly-identified
component). The optional **doubly-constrained** (Furness) attraction constraint is available per
component, OFF by default (see "Doubly-constrained option"). TLD-anchoring the willingness
(`fit_kernel.py`/`iterate_kernel.py`) and a national `equiv_miles` remain follow-ups (see memory
`project-tanner-kernel-tld`).

A carry-forward open concern is the **school component**: it is weakly identified (walking obs mostly
fall in slots with `f_school≈0`, so it is pinned almost entirely by the official school-peak hours)
and can act as an AM/PM-peak fitter, pushing the school share up (school-peak count data is the lever,
not prior strength).

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
- **Component degeneracy:** the convex `solve_scales` could otherwise trade magnitude between
  similarly-shaped components; the generation-anchored K-prior (`Σ_c (K_c−1)²/σ_c²`, `K_prior_std`
  default 0.5) pulls every scale toward the generation value 1, which both anchors magnitude and pins
  the flat commute↔retail direction.
- **Component scales** `K_res`/`K_commute`/`K_retail`/`K_primary`/`K_postprimary`/`K_tertiary` are calibrated directly by the convex
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
- **`tuned_params.json` structure:** the six scales `K_res`/`K_commute`/`K_retail`/`K_primary`/`K_postprimary`/`K_tertiary` (plus a derived `K_sch` = their sum, written for display/report only), the 18 double-exp willingness params (`<comp>_taus/_taul/_w` for res/commute/retail + the three independent school levels; τ's in seconds), `"kernel": "modesub_double"`, `slot_fracs_res`/`slot_fracs_commute`/`slot_fracs_retail`/`slot_fracs_school_primary`/`_postprimary`/`_tertiary` (dicts keyed `"dt,h"`, the pinned NTS profile), and — when double-constraint is active — `doubly_constrained` + `furness_max_sweeps`. `reset_gravity_params.py` regenerates this clean structure; `build_assignment.py` requires these multi-component params and fails loud on an old single-K param file.

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

**External demand** — the separate attractor layers `commute_attractor` (car-commute jobs), `retail_spaces`, and the three per-level `school_demand_<level>` (and the census producers `commute_producers` + per-level `school_producers_<level>`) — is measured per zone via the same estimators as internal nodes, with no external scale factors. Known limitation: OSM under-maps RoI schools (~80% of real schools mapped), so RoI external school demand runs proportionally low.

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
must be **fast** (a real OSRM re-extract is far too slow per tuning step), so a
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
