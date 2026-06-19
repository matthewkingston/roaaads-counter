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

python3 analysis/derive_component_profiles.py              # regenerate hourly_fractions.csv school column (re-run when NTS files change)

python3 analysis/tune_assignment.py                        # Stage 1: tune gravity params (9 params; ~40-50s after probit cache rebuild)
python3 analysis/tune_assignment.py --full                 # Stage 2: + external zones (31 params, ~70s after probit cache rebuild)
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
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, assigns external zone weights, builds map. `--map-only` skips demographic recomputation and rebuilds only the HTML. `--zones-only` patches boundary node weights without rebuilding. **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Business demand:** workplace population × POI count (school/college/university excluded from this layer), augmented with OSM car park polygon area (public: area/25, private: area/50 equivalent persons). **School demand:** separate `node_school_demand` layer from OSM school POIs (amenity=school/secondary_school/college/university). Enrollment (pupils) is the control total per POI: OSM `capacity` tag used where present, otherwise type-based fallback (school→300, secondary_school→900, college→2000, university→3000; DE NI 2023/24 averages). Units are pupils, so W_SCHOOL is interpretable as a trip-generation ratio relative to residential population. No DZ distribution — each school's enrollment is assigned directly to its snapped road node(s). **Flow map layers:** reads `newtownards_flows.json` and adds: (1) combined AADT layer (shown by default, blue→yellow→red, tooltip includes res/biz/school breakdown); (2) residential layer (teal, off by default); (3) business-adjacent layer (amber→red, off by default); (4) school layer (violet→purple, off by default). Component layers only appear if `flows_res`/`flows_biz`/`flows_school` keys exist. Map also includes parking (blue/red) and POI layers. **Adding a new external zone:** add an entry to `_EXT_GEO` and to `tuner_config.json` cities. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Re-run if road network changes, `HIGHWAY_COST_FACTOR` values change, or `N_PASSES`/`PROBIT_CV` change. Edge costs are travel time × a road-class multiplier (trunk/primary: ×0.67, residential/unclassified/living_street: ×1.2, others: ×1.0). Also reads `tuner_config.json` to filter external→external OD pairs. **Probit stochastic loading:** runs `N_PASSES=25` Dijkstra passes with log-normal edge-cost noise (CV=0.25), accumulates fractional link-assignment weights (`link_weight` per entry = fraction of passes using that link for that OD pair). `od_dist` is the mean path distance across passes. Build time ~30–60 min. |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (rational kernel; `return_components=True` with `w_school` provided returns `(flow_res, flow_biz, flow_school)` tuple, without `w_school` returns `(flow_res, flow_biz)`), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `gravity_assign()` accepts optional `link_weight` array (probit fractional weights). `compute_chi2()` dispatches to three-component mode when `link_flow_school_dict` provided (N_eff = N − 3·N_slots); two-component when only `link_flow_biz_dict` provided; legacy single-flow otherwise. |
| `simulation/build_assignment.py` | Gravity model assignment. Requires `simulation/newtownards_paths.npz`. Three-component mode activated when K_sch > 0 and W_SCHOOL are in `tuned_params.json` and `node_school_demand` is in `node_weights.json`. Saves `flows_res`, `flows_biz`, `flows_school` in `newtownards_flows.json`. Falls back to two-component (K_res/K_biz only) or legacy single-K for old param files. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Reference values for L2 regularization, city→node groupings, `through_route_pairs` whitelist, gravity param regularization, and temporal coupling. `lambda` regularises external zones; `gravity_lambda` + `gravity_ref` regularise P/ALPHA/BETA/W_BIZ/P_biz/ALPHA_biz; `gamma_coupling_scale` controls the per-slot aggregate coupling (γ = scale/std_f²); `phi_biz_prior` + `phi_biz_std` set the Gaussian prior on the business flow fraction (see Model Design). |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. Run once (or when the ODS file changes). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. Loads manual link overrides from `data/manual_link_overrides.json` (use `analysis/manual_assign_link.py` to add entries). After every new link assignment (manual or auto), validates that each non-null count direction maps to a real directed edge in G; raises `ValueError` if not, preventing counts from silently hitting a zero-flow phantom edge. Edges without geometry (virtual stub nodes such as Dundonald 10000) are skipped during snap candidate construction. |
| `analysis/manual_assign_link.py` | CLI tool to manually assign a session to a specific directed link, bypassing GPS snap. Usage: `python3 analysis/manual_assign_link.py <session_id> <from_node> <to_node>` where `from_node→to_node` is the observer's walking direction ("with"). Validates both nodes exist and checks count-edge consistency: a non-null count must map to a real directed edge (mirrors `ingest_counts.py` rule). Handles one-way roads correctly — the observer may walk with or against traffic; only the direction carrying non-null counts must have a real edge. Writes to `data/manual_link_overrides.json` and patches `counts_processed.json` directly. Idempotent: re-running overwrites the previous assignment; re-running `ingest_counts.py` afterwards is safe (skips already-snapped sessions) but not required. After correcting an assignment, re-run `aggregate_counts.py` then `tune_assignment.py`. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry now carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s` so the tuner can work in count space. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Three-component model:** gravity flows split into residential (pop→pop, `flow_res`), business-adjacent (pb+bb, `flow_biz`), and school (pop→school, `flow_school`) components. Stage 1 tunes 9 gravity params (W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz, W_SCHOOL, P_school, ALPHA_school); `--full` adds 14 city pop/wp + 6 dampings (31 params total). Uses 216 official hourly obs + 343 walking obs. **Alternating minimisation (5 blocks, up to 10 iters):** K-step (1D); phi_biz-step (1D, Gaussian prior from `phi_biz_prior`/`phi_biz_std` in config); phi_sch-step (1D, Gaussian prior from `phi_school_prior`/`phi_school_std`); f_res/f_biz/f_school steps (per-slot analytical via bincount) + aggregate coupling γ·(f_res+f_biz+f_school−f_agg)² per slot. School component disabled gracefully if `node_school_demand` absent from weights. **Performance estimate with probit cache:** ~40–50 s stage 1, ~70 s stage 2. `--fast` mode caps alt-min at 5 iters. |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` (summary, chi² tables, gravity params with K_res/K_biz/K_sch/phi_biz/phi_sch, external zones, slot fractions table showing f_res/f_biz/f_school with pulls) and `reports/slot_pulls_{id}.png` (two or three side-by-side heatmaps: 24 h × 3 day-types per component). History `slot_prior` entries carry 5 values: `[mean_f_agg, std_f, mean_f_res, mean_f_biz, mean_f_school]`. Old entries with 4 values handled gracefully. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. Writes all keys from `params` dict (K_res, K_biz, slot_fracs_res, slot_fracs_biz, external zone params etc.); clears stale `slot_fracs` legacy key. |
| `simulation/reset_gravity_params.py` | Reset only the gravity params (K, W_BIZ, P, ALPHA, BETA, P_biz, ALPHA_biz) in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. External zone params are preserved. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile. Includes `mean_fraction_res`, `mean_fraction_biz`, and `mean_fraction_school` columns: temporal profile priors for the three-component model. **Derived from NTS** via `analysis/derive_component_profiles.py`. Constraint: res + biz + school = agg for all 168 rows. **Summation convention:** rows sum to 7.0 (AADT normalisation). Re-run `derive_component_profiles.py` whenever NTS files change; do not edit aggregate columns by hand. |
| `analysis/derive_component_profiles.py` | Derives `mean_fraction_res`, `mean_fraction_biz`, and `mean_fraction_school` from DfT NTS data (2023–2024 rolling average). **Purpose classification:** biz = commuting + employer's business + shopping; school = education×⅕ + escort education; res = remainder. Education ÷5: the standalone child trip has far lower car generation; escort education is the adult drop-off car trip, kept at full weight. **Weekdays (NTS0502a):** per-hour shares computed for each of the three components. **Weekends (NTS0504b):** flat school_share and biz_share per Saturday/Sunday, with "Just walk" removed and education ÷5 applied; temporal shape from aggregate profile. Constraint: res + biz + school = agg exactly. Re-run whenever NTS files or purpose classification change. |

### Generated / gitignored outputs
`simulation/newtownards_paths.npz`, `simulation/node_weights.json`,
`simulation/newtownards_map.html`, `simulation/tuned_params.json` — all regenerated by the pipeline.
`simulation/node_weights.json` keys: `node_population`, `node_business_demand`, `node_school_demand`, `node_parking_equiv`, `node_effective_utm`, `boundary_node_ids`.
`simulation/newtownards_flows.json` — combined flows plus optional `flows_res`/`flows_biz`/`flows_school` keys when three-component params active.
`reports/` — generated by `report_tune.py` and `tune_assignment.py`; not tracked.

### Tracked generated outputs
`data/counts_processed.json`, `data/link_aadt.json`, `data/official_hourly.json`,
`simulation/tuning_history.jsonl` — committed so history is preserved.
`data/manual_link_overrides.json` — committed so manual assignments survive a wipe of `counts_processed.json`.
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
| 2026-06-19 | gravity | 559 | 9 | 1.3292 | three-component model (school added); phi_biz=27%, phi_sch=1.2%; school at ref |
| 2026-06-19 | full | 559 | 31 | **1.3146** | three-component full tune; phi_biz=25.6%, phi_sch=1.4%; school params at ref |

**Note on comparability:** runs from 2026-06-17 onward use the two-component model with coupling penalty terms in chi²/N; not directly comparable to earlier single-component runs. From 2026-06-19 three-component model: N_eff = 559 − 3×72 = 343 (one extra df per slot for f_school).

Current best full-tune: chi²/N = **1.3146** (559 obs, N_eff=343; three-component with probit cache, run f09a003e).
K_res=1.23e-04, K_biz=4.31e-05, K_sch=2.40e-06. phi_biz=25.6%, phi_sch=1.4%.
W_BIZ=3.82, P=117.6s, ALPHA=4.02, BETA=7.67. P_biz=83.4s, ALPHA_biz=3.66.
W_SCHOOL=1.00, P_school=600s, ALPHA_school=2.00 (at ref — school params not driven by data).
mean|z|=0.80  |z|>2: 40  |z|>3: 12.

**Confirmed working:**
- Three-component model runs end-to-end; display bug fixed (K scaling in walking display).
- School temporal prior correct: weekday h08=40.1% school, h15=34.5%; weekends ~0.5%.
- External zone tuning improved: Belfast wp now +20.7% (was +891% before ref update); Bangor wp −43.8%.

**Outstanding concerns:**
- **phi_sch=1.4%** — school component contributes very little. School params sit at gravity_ref because walking count sessions don't cover h08/h15 school-peak hours. The school kernel (P/ALPHA) is unidentifiable without school-peak observations.
- Structural outliers: `22→12 Regent Street` (z=+4.03), `23→295 Frances Street` (z=+3.94), `296→297 Nursery Road` (z=−3.70), `139→137 Portaferry Road` (z=−3.70).
- `73→70` Mill Street severe underprediction (z=−3.30; obs 26,377 vs model 2,682).
- `719→325` / `325→719` Messines Road persistent (z=−3.33/−2.53).
- Hardford Link persistent (z=−3.25/−3.19/−2.73).

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
- **`tuned_params.json` structure:** contains `K_res`, `K_biz`, `K_sch`, `W_SCHOOL`, `P_school`, `ALPHA_school`, `slot_fracs_res`, `slot_fracs_biz`, `slot_fracs_school` (dicts keyed `"dt,h"`); does **not** contain a `slot_fracs` key (legacy). Old param files without school keys fall back to two-component or legacy mode in `build_assignment.py`.

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
