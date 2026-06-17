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
python3 simulation/build_paths.py            # precompute k=3 shortest paths (~6 min)

python3 analysis/parse_official_hourly.py    # parse ODS hourly counts → data/official_hourly.json (one-off)
python3 analysis/ingest_counts.py            # process walking count CSVs → counts_processed.json
python3 analysis/aggregate_counts.py         # combine per-session AADT → link_aadt.json

python3 analysis/tune_assignment.py                        # Stage 1: tune gravity params (4 params incl. THETA, ~5 min)
python3 analysis/tune_assignment.py --full                 # Stage 2: + external zones (26 params, ~15 min)
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
| `simulation/build_demographics.py` | Downloads NISRA population, allocates to nodes, assigns external zone weights, builds map. `--map-only` skips demographic recomputation and rebuilds only the HTML. `--zones-only` patches boundary node weights without rebuilding. **Population distribution:** building centroids snapped to road edges; DZs with <3 buildings fall back to road-length weighting. **Business demand:** workplace population × POI count, augmented with OSM car park polygon area (public: area/25, private: area/50 equivalent persons). **Flow map layers:** reads `newtownards_flows.json` and adds: (1) combined AADT layer (shown by default, blue→yellow→red, tooltip includes res/biz breakdown if component flows present); (2) residential layer (teal, off by default); (3) business-adjacent layer (amber→red, off by default). Component layers only appear if `flows_res`/`flows_biz` keys exist in the flows file. Map also includes parking (blue/red) and POI layers. |
| `simulation/build_paths.py` | Precomputes all-pairs shortest paths; result cached in `newtownards_paths.npz`. Re-run if road network changes or `HIGHWAY_COST_FACTOR` values change. Edge costs are travel time × a road-class multiplier (trunk/primary: ×0.67, residential/unclassified/living_street: ×1.2, others: ×1.0) to bias routing toward major roads. Also reads `tuner_config.json` to filter which external→external OD pairs to include as through routes. Produces **k=3 alternative paths** per OD pair (k=2/k=3 via progressive edge penalisation ×10) for stochastic logit routing. Build time ~30 min (3 Dijkstra passes). |
| `simulation/model.py` | **Shared constants and functions:** `COUNT_SITES`, `EXCLUDE_LINKS`, file-path constants, `gravity_assign()` (rational kernel, `return_components=True` returns `(flow_res, flow_biz)` tuple), `site_flow()`, `compute_chi2()`, `print_chi2_table()`. `compute_chi2()` has two modes: **two-component** (when `link_flow_biz_dict` provided) uses 216 official hourly obs + count-space formula matching the tuner; **legacy** (single flow dict) uses 3 AADT obs + Woodbury correction. |
| `simulation/build_assignment.py` | Gravity model assignment. Requires `simulation/newtownards_paths.npz` (exits with an informative error if absent). Calls `gravity_assign(return_components=True)` and `compute_chi2` two-component mode when K_res/K_biz are in `tuned_params.json`. Saves `flows_res` and `flows_biz` alongside `flows` in `newtownards_flows.json` when two-component. Falls back to legacy single-K mode for old param files. |
| `simulation/edit_network.py` | Manual network edits (node deletions etc.). |
| `simulation/tuner_config.json` | **Tracked in git.** Reference values for L2 regularization, city→node groupings, `through_route_pairs` whitelist, gravity param regularization, and temporal coupling. `lambda` regularises external zones; `gravity_lambda` + `gravity_ref` regularise P/ALPHA/W_BIZ/THETA; `gamma_coupling_scale` controls the per-slot aggregate coupling (γ = scale/std_f²); `phi_prior` + `phi_std` set the Gaussian prior on the business flow fraction (see Model Design). |
| `analysis/parse_official_hourly.py` | Parses sheets 444/507/508 from the 2023 NI ODS traffic count file → `data/official_hourly.json`. Run once (or when the ODS file changes). Weekday sigma = max(between-day std, 10% relative, √count); weekend sigma = max(√count, 15% relative). The √count floor prevents unrealistically tight sigmas at overnight low-count hours. |
| `analysis/ingest_counts.py` | Reads all CSVs from `data/counts/`, snaps GPS tracks to road links, estimates per-session AADT via hourly fraction profile. Idempotent: skips already-processed sessions. **`MANUAL_LINK_OVERRIDES`** dict at the top maps session IDs to forced directed links; use when the observer was on a parallel carriageway and GPS snap would land on the wrong road (e.g. A20 Kempe Stones sessions `e644eae2`/`760b0c8e` → link 8→7). After every new link assignment (manual or auto), validates that each non-null count direction maps to a real directed edge in G; raises `ValueError` if not, preventing counts from silently hitting a zero-flow phantom edge. Edges without geometry (virtual stub nodes such as Dundonald 10000) are skipped during snap candidate construction. |
| `analysis/aggregate_counts.py` | Combines per-session AADT estimates into per-link estimates using inverse-variance weighting. Always regenerates from scratch. Each observation entry now carries `n_eff` (Jeffreys count = n + 0.5) and `duration_s` so the tuner can work in count space. Output: `data/link_aadt.json`. |
| `analysis/tune_assignment.py` | Powell's method parameter tuning. **Two-component model:** gravity flows split into residential (pop→pop, `flow_res`) and business-adjacent (pb+bb, `flow_biz`) components, each with its own temporal profile (f_s_res, f_s_biz) and scale (K_res, K_biz). Stage 1 tunes 4 gravity params (W_BIZ, P, ALPHA, THETA); `--full` adds 14 city pop/wp + 6 dampings. Uses 216 official hourly obs from `official_hourly.json` (replacing 3 AADT point estimates) plus per-session walking obs from `link_aadt.json`. **Alternating minimisation (4 blocks, 10 iters):** K-step (1D, total scale); phi-step (1D, K_biz/K ratio, Gaussian prior N(0.35, 0.15²) prevents K_biz→0); f_res-step (per-slot analytical); f_biz-step (symmetric) + aggregate coupling γ·(f_res+f_biz−2·f_agg)² per slot. **Performance (stochastic k=3):** ~300 ms/eval; ~5 min stage 1, ~15 min stage 2. |
| `analysis/report_tune.py` | Generate a structured report from a tuning history entry. Writes `reports/tune_report_{id}.txt` (summary, chi² by measurement, chi² by link, gravity params with K_res/K_biz/phi, external zones, slot fractions table showing f_res and f_biz side by side with component pulls) and `reports/slot_pulls_{id}.png` (two side-by-side heatmaps: 24 h × 3 day-types per component, colour = pull from prior). History `slot_prior` entries carry 4 values: `[mean_f_agg, std_f, mean_f_res, mean_f_biz]`. History `tuner_hyperparams` carries phi_prior, phi_std, gamma_coupling_scale, gravity_lambda, lambda at run time. |
| `simulation/restore_params.py` | Restore `tuned_params.json` from any history entry by run ID. `--list` shows all runs; partial ID prefix matching is supported. Writes all keys from `params` dict (K_res, K_biz, slot_fracs_res, slot_fracs_biz, external zone params etc.); clears stale `slot_fracs` legacy key. |
| `simulation/reset_gravity_params.py` | Reset only the gravity params (K, W_BIZ, P, ALPHA, THETA) in `tuned_params.json` to the `gravity_ref` anchors in `tuner_config.json`. External zone params are preserved. |
| `data/counts/*.csv` | Raw walking count CSVs from the recorder app. Add new files and re-run `ingest_counts.py`. |
| `analysis/hourly_fractions.csv` | **Tracked in git.** NI-average hourly fraction profile (fraction of AADT per hour×day-of-week). Now includes `mean_fraction_res` and `mean_fraction_biz` columns: synthetic priors for the two-component temporal profiles (business: AM/PM weekday peaks; residential: complement). Used as prior means in the tuner. |

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
Currently present: `data/2023-northern-ireland-traffic-count-data-in-ods-format.ods`
(used by `parse_official_hourly.py` to extract hourly counts for sites 444/507/508;
annual AADT values in `model.py` `COUNT_SITES` are no longer used by the tuner but
retained for `build_assignment.py`) and `data/census-2021-apwp001.xlsx`.

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
minimisation (10 iterations, converges in 3–5), where K_res = K·(1−phi), K_biz = K·phi.

**K-step:** 1D solve (same structure as old single-K step), using combined coefficient
`(1−phi)·c_r·f_res + phi·c_b·f_biz` per observation.

**phi-step:** 1D solve for business fraction, with Gaussian prior phi ~ N(0.35, 0.15²).
This prior prevents K_biz collapsing to 0 (the K_biz/W_BIZ degeneracy otherwise exploited
by the optimizer). phi = 0 means all-residential; phi = 1 means all-business.

**f_res-step / f_biz-step:** per-slot analytical update, anchored by component-specific
priors from `hourly_fractions.csv` columns `mean_fraction_res` / `mean_fraction_biz`
(synthetic profiles: business peaks weekday AM/PM; residential is complement).

**Aggregate coupling:** each slot also carries a penalty γ·(f_res + f_biz − 2·f_agg)²
where γ = `gamma_coupling_scale` / std_f_agg² (per-slot, from `tuner_config.json`).
This prevents the two profiles collectively drifting from the known NI aggregate profile.
`gamma_coupling_scale = 0.0` disables coupling; `1.0` gives coupling as strong as the
individual component priors.

Slot key: (day_type, hour), day_type = 0 (weekday), 1 (Saturday), 2 (Sunday).
Prior std derived from `hourly_fractions.csv` via law of total variance.

### Observations
All 374 observations are in count space with per-obs weights:
- **Official hourly** (216 obs, 24 h × 3 day-types × 3 sites): from `data/official_hourly.json`;
  Gaussian error (sigma from between-weekday std, 10% floor); weight = 1/sigma².
- **Walking** (158 obs): from `data/link_aadt.json`; Poisson error; weight = 1/n_eff.

### Goodness of fit
`χ²/N` (mean squared z-score; N=374 obs, N_eff = N − 2·N_slots = 374 − 144 = 230).
Two df lost per slot (one each for f_s_res and f_s_biz). With coupling enabled,
chi²/N includes the coupling penalty terms; pure data-fit chi²/N is lower.

`build_assignment.py` uses the two-component `compute_chi2()` when K_res/K_biz are present in tuned_params.json, giving a chi²/N directly comparable to the tuner's (minus the coupling penalty). The legacy Woodbury path is used only for old single-K param files.

---

## Count Data

**Official hourly obs (from ODS, parsed by `parse_official_hourly.py`):**
- Site 507: A21 Bangor Road — 72 obs (24 h × 3 day-types); Gaussian sigma, 10%/15% floors
- Site 508: A48 Donaghadee Road — 72 obs
- Site 444: A20 Portaferry Road — 72 obs

Annual AADT values (507: 21,202; 508: 10,792; 444: 7,282) are retained in `model.py`
`COUNT_SITES` for `build_assignment.py` backward compatibility but are no longer used
by the tuner directly.

**Walking counts:** 6 CSV files, 170 sessions, 329 per-session observations (after EXCLUDE_LINKS). Two sessions manually assigned: `e644eae2` and `760b0c8e` (A20 Kempe Stones eastbound, link 8→7; observer was on the westbound carriageway). The tuner uses per-session observations directly; per-link aggregates are retained in `link_aadt.json`.
⚠ **A re-tune is needed to incorporate the new sessions.**

**Total (after re-tune): 545 observations (216 official hourly + 329 walking) in 72 time slots.
N_eff = 545 − 2×72 = 401.**

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

**Note on comparability:** runs from 2026-06-17 onward use a different observation set
(374 obs vs 161, including temporal official hourly obs) and include coupling penalty
terms in chi²/N. They are not directly comparable to earlier single-component runs.

⚠ **PROVISIONAL STATE** — the following are outstanding before this baseline is considered validated:
1. **Full-stage (`--full`) two-component run not yet done.** External zone params are still
   at the refs from the 2026-06-16 single-component run; city populations will shift after
   a two-component `--full` re-tune.
2. **New count data (2026-06-17 ingest: 170 sessions, 165 links) not yet incorporated.**
   Current best was tuned on 158 walking obs. A re-tune with the new data is needed.
3. **sigma floor fix (√count) applied to official_hourly.json but not yet re-tuned.**
   Site 444 overnight outliers (previously z≈−5 at h04) should improve; confirmed only
   after the next tuning run.
4. **Map layers (residential/business) not yet confirmed** — depend on a `build_assignment.py`
   run with two-component params to populate `flows_res`/`flows_biz` in `newtownards_flows.json`.

Pending full validation: chi²/N = 1.4286 (374 obs, N_eff=230, 72 slots; two-component with coupling, pre-sigma-fix, pre-new-count-data).
W_BIZ=1.20, P=87s, ALPHA=3.75, THETA=0.07. phi≈13% business fraction.
Persistent structural outliers: `328→326` Comber Road (z=−3.29), `719→325` Messines Road (z=−3.27), `18→21` Hardford Link (z=−2.72).
`22→159` was a data error (snap direction bug, fixed 2026-06-15): now recorded as 159→22.
Belfast Road `20→18` persistent underprediction (z≈−2.5).

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
  ALPHA/P change). K_res and K_biz are derived from K × (1−phi) and K × phi respectively.
  K is not interpretable in isolation; chi²/N is reliable.
- After a structural model change (e.g. adding through routes or new count data), a gravity-only
  stage 1 run will show inflated chi²/N. A full `--full` re-tune is needed to restore fit quality.
- **Dundonald virtual node (added 2026-06-17):** Node 10000 is a degree-1 stub connected
  only to node 97, representing the Dundonald catchment on the A20 corridor. Alternative paths
  k=2/k=3 fall back to k=1 for all node 10000 OD pairs.
- **Manual link overrides:** `MANUAL_LINK_OVERRIDES` in `ingest_counts.py` hard-assigns specific sessions to a directed link, bypassing GPS snap. Use when the observer stood on a parallel carriageway (e.g. a dual one-way road) and the snap would land on the wrong physical road. The override is idempotent and takes effect even if `counts_processed.json` is wiped and rebuilt. After assignment (manual or auto), the script validates each non-null count direction against the directed graph and raises `ValueError` if the edge doesn't exist.
- **Snap direction bug (fixed 2026-06-15):** `ingest_counts.py` previously stored canonical
  `(min(u,v), max(u,v))` — fixed to store actual directed `(u, v)`. Only session `f56b2ce4`
  was materially affected (re-snapped from 22→159 to 159→22).
- Two temporal profiles (f_s_res, f_s_biz) are inferred per (day_type, hour) slot, each
  anchored by component-specific priors from `hourly_fractions.csv`. The aggregate coupling
  (gamma_coupling_scale / std_f²) per slot keeps their sum near 2×f_agg. With 72 slots
  and 2 df each, N_eff = N − 2×N_slots = 230.

---

## External Zone Reference Values (`tuner_config.json`)

**These values must not be changed without explicit user approval.** After a full tuning run,
updating them is something to *consider and discuss*, not an automatic step — the refs anchor
L2 regularization and changing them shifts the penalty basin for all future runs.

Last updated: 2026-06-16 full tuning run (χ²/N=1.2546, 258 obs, stochastic k=3).
Dundonald added 2026-06-17 (initial priors only, not yet tuned).

| City | Nodes | ref_pop | ref_wp | Tunable dampings |
|------|-------|---------|--------|-----------------|
| Donaghadee | 47 | 190,201 | 7,018 | — |
| Comber | 65, 617, 618, 620 | 53,571 | 2,996 | 617 (×0.38), 618 (×0.35), 620 (×0.43) |
| LowerArds | 92 | 84,500 | 5,024 | — |
| Belfast | 97, 119 | 1,034,719 | 183,661 | 119 (×0.31) |
| Dundonald | 10000 | 150,000 | 8,000 | — |
| Bangor | 98, 731 | 95,426 | 21,246 | 98 (×0.39) |
| Holywood | 99 | 3,652 | 1,203 | — |
| Millisle | 748, 749 | 2,570 | 498 | 749 (×0.47) |

---

## Worktree Convention

Background/parallel Claude work is done in `.claude/worktrees/` and
cherry-picked to `main` after review. All work is on `main`.
