"""Per-small-area school-trip *producers*, split by level: primary / post-primary / tertiary.

`load_school_producers() -> {area_code: {"primary", "postprimary", "tertiary"}}` for every NI
Data Zone + RoI Small Area, island-wide, keyed by area_code so any CENTRE on the island works
(jurisdiction handled internally — same portable pattern as census_supply / census_attractor).

This is the Phase-2 replacement for the single lumped `school` producer in census_supply.py.

Design (settled with the user — see agent memory project_school_split_phase2)
--------------------------------------------------------------------------------
The producer's ABSOLUTE total is irrelevant: in generation pinning
`trips_{c,i} = ρ_c · Σpop · (producer_i / Σproducer_c)`, the total cancels. Only two things
survive and must be right: (a) the within-jurisdiction SPATIAL distribution, and (b) the NI:RoI
relative scale. Census age bands give unreliable *magnitudes* (off enrolment by −9%..+14% for
primary/post-primary, with fence-posts at both ends: age-4 Junior Infants in AGE0-4, age-15
Transition-Year double-count), so magnitude is NOT taken from census age. Instead:

  * a per-(area, level) SPATIAL WEIGHT (shape) is built from census, then
  * each jurisdiction's per-level weight is SCALED to that jurisdiction's admin ENROLMENT
    (from the attractor cache) — putting NI and RoI producers in the same real-student units
    (the NI:RoI compatibility that matters); the island total then harmlessly cancels.

Spatial weights (SHAPES only — absolute value doesn't matter, only cross-area ratios do):
  NI  (per DZ, from NISRA census CSVs):
    primary      = compulsory_school_age_agg4 "Primary school age"
    postprimary  = compulsory_school_age_agg4 "Post-primary school age"
    tertiary     = in_full_time_education total − primary − postprimary   (>= 0)
                   (= "above-compulsory full-time students" by DZ of residence)
  RoI (per SA, SAP2022 census 5-yr age bands; band→level fractions from CSO PxStat pupils-by-age,
       EDA42/EDA70 2021-22 aligned to the 2022 census — see the constants block for provenance):
    primary      = 0.040·AGE0-4 + 1.000·AGE5-9 + 0.549·AGE10-14
    postprimary  = 0.451·AGE10-14 + 0.660·AGE15-19
    tertiary     = max(0, T10T2 "still at school or college" 15+ − 0.660·AGE15-19)
                   (15+ full-time students minus the measured 15-19 post-primary share — consistent)

Enrolment (magnitude scale + the ρ generation split) comes from the attractor cache
`cache_admin_schools_island.geojson`, summed by jurisdiction × level, with the small 'special'
level folded into primary/post-primary by the NI primary:post-primary ratio (~0.513/0.487). Using
the attractor as the single source keeps producer and attractor on the same footing (DRY).

NOTES / documented limitations
  * T10T2 = FULL-TIME students only (CSO "have you ceased full-time education?"). Full-time,
    college-based FE/PLC count; part-time and work-based (apprentices) don't. Fine here — the
    producer only needs a spatial shape and an NI:RoI scale, both robust to this.
  * Band-split fractions come from measured pupils-by-age (EDA42/EDA70, 2021-22) ÷ census persons
    (2022) — WELL-ALIGNED years, so no cohort-size drift (that pitfall bit an earlier attempt that
    mixed 2025 enrolment with 2022 persons). The measured AGE15-19 post-primary share is 0.660, i.e.
    only ~66% of 15-19s are in post-primary (senior-cycle leavers + tertiary) — notably below a naive
    100%-participation assumption; this per-SA 10-14:15-19 balance is real spatial information.
  * Age-4 primary pupils (11,931; incl. the 16.5% of Junior Infants who are age-4) are captured via
    0.040·AGE0-4. Pre-school (0-3) gets zero weight, so childcare is excluded by construction — this
    producer needs NO CHILDCARE_NATIONAL subtraction (that was a Phase-1 C02-cleanup, now obsolete).
  * AGGREGATE-ROW GOTCHA: the RoI census JSON-stat carries a national "State" SA row (== Σ of all
    SAs); dropped here via a per-cell threshold. (The attractor cache's own aggregate rows are
    already dropped upstream by build_admin_schools.)
"""
import glob
import json

from census_supply import _one, _ni_csv   # shared stdlib helpers (file glob + NISRA CSV reader)

try:
    from demographics_config import SCHOOL_ISLAND_CACHE
except Exception:                          # keep importable without the config module
    SCHOOL_ISLAND_CACHE = "data/cache_admin_schools_island.geojson"

# --- source globs -----------------------------------------------------------------------------
NI_COMPULSORY_GLOB = "data/ni_census/*compulsory_school_age_agg4*.csv"
NI_EDU_GLOB        = "data/ni_census/*in_full_time_education*.csv"
ROI_AGE_GLOB       = "data/ireland_census/SAP2022T1T1A*.json"
ROI_T10T2_GLOB     = "data/ireland_census/SAP2022T10T2SA*.json"

# RoI 5-year census bands → school levels, from the ACTUAL single-year-of-age pupil distribution
# (not a modelled ladder), aligned to the census year (2021-22 pupils vs the April-2022 census):
#   primary-by-age  : CSO PxStat EDA42 "Primary Pupils" (Age × programme), 2021-22
#   post-primary tot: CSO PxStat EDA70 "Pupils Enrolled in Second Level Schools", 2021-22 = 391,698
#   census persons  : SAP2022 T1T1A (2022, "State" aggregate row dropped)
# (API: ws.cso.ie/public/api.restful/PxStat.Data.Cube_API.ReadDataset/EDA42|EDA70/JSON-stat/2.0/en)
# Ages 10-14 are compulsory, so post-primary-aged 10-14 = persons − primary (residual); post-primary
# aged 15-19 = total post-primary − that. These reproduce the Irish Independent figure exactly
# (age-4 Junior Infants 10,482/63,583 = 16.5%) and the CSO totals (Σ primary weight ≡ EDA42 560,264,
# Σ post-primary ≡ EDA70 391,698). Used as spatial SHAPES only — magnitude comes from enrolment
# scaling in load_school_producers, so the small 2021-22-pupils vs 2022-persons offset is absorbed.
_ROI_PERSONS         = {"0_4": 295415, "5_9": 342670, "10_14": 374202, "15_19": 337628}
_ROI_PRIMARY_BY_BAND = {"0_4": 11931,  "5_9": 343006, "10_14": 205327}
_ROI_POSTPRIM_TOTAL  = 391698
_ROI_POSTPRIM_10_14  = _ROI_PERSONS["10_14"] - _ROI_PRIMARY_BY_BAND["10_14"]     # 10-14 compulsory

_F04_PRIMARY    = _ROI_PRIMARY_BY_BAND["0_4"]   / _ROI_PERSONS["0_4"]            # 0.040
_F59_PRIMARY    = min(1.0, _ROI_PRIMARY_BY_BAND["5_9"] / _ROI_PERSONS["5_9"])    # ~1.0
_F1014_PRIMARY  = _ROI_PRIMARY_BY_BAND["10_14"] / _ROI_PERSONS["10_14"]          # 0.549
_F1014_POSTPRIM = _ROI_POSTPRIM_10_14 / _ROI_PERSONS["10_14"]                    # 0.451
_F1519_POSTPRIM = (_ROI_POSTPRIM_TOTAL - _ROI_POSTPRIM_10_14) / _ROI_PERSONS["15_19"]   # 0.660

_STATE_THRESHOLD = 100_000   # any per-SA census cell above this is the national "State" aggregate
_LEVEL_KEYS = ("primary", "postprimary", "tertiary")


# ============================== enrolment (magnitude scale + ρ split) =========================
def _enrolment_by_juris_level():
    """{(jurisdiction, level): enrolment} from the attractor cache, 'special' folded into
    primary/post-primary by the NI primary:post-primary ratio. level ∈ {primary, postprimary,
    tertiary}; jurisdiction ∈ {'NI','RoI'}."""
    gj = json.load(open(_one(SCHOOL_ISLAND_CACHE)))
    lvl_map = {"primary": "primary", "post_primary": "postprimary",
               "tertiary": "tertiary", "special": "special"}
    raw = {}
    for ft in gj["features"]:
        p = ft["properties"]
        j = p.get("jurisdiction")
        lv = lvl_map.get(p.get("level"))
        if j is None or lv is None:
            continue
        raw[(j, lv)] = raw.get((j, lv), 0.0) + (p.get("enrolment") or 0.0)

    jurisdictions = {j for (j, _) in raw}
    # fold 'special' by the NI primary:post-primary ratio
    ni_p = raw.get(("NI", "primary"), 0.0)
    ni_pp = raw.get(("NI", "postprimary"), 0.0)
    r_prim = ni_p / (ni_p + ni_pp) if (ni_p + ni_pp) > 0 else 0.5
    out = {}
    for j in jurisdictions:
        sp = raw.get((j, "special"), 0.0)
        out[(j, "primary")]     = raw.get((j, "primary"), 0.0)     + sp * r_prim
        out[(j, "postprimary")] = raw.get((j, "postprimary"), 0.0) + sp * (1.0 - r_prim)
        out[(j, "tertiary")]    = raw.get((j, "tertiary"), 0.0)
    return out


def island_enrolment_by_level():
    """{level: island enrolment} — the basis for the ρ (generation) split across levels."""
    enrol = _enrolment_by_juris_level()
    return {lv: sum(v for (_, l), v in enrol.items() if l == lv) for lv in _LEVEL_KEYS}


# ============================== spatial weights (shapes) ======================================
def _ni_weights():
    """NI per-DZ spatial-weight shapes {dz: {primary, postprimary, tertiary}}."""
    primary, postprim, ft = {}, {}, {}
    for dz, label, n in _ni_csv(_one(NI_COMPULSORY_GLOB)):
        lbl = label.strip().lower()
        if "post-primary school age" in lbl:
            postprim[dz] = postprim.get(dz, 0.0) + n
        elif "primary school age" in lbl:          # after post-primary (substring guard)
            primary[dz] = primary.get(dz, 0.0) + n
    for dz, label, n in _ni_csv(_one(NI_EDU_GLOB)):
        if label.strip().lower() == "full-time student or schoolchild":
            ft[dz] = ft.get(dz, 0.0) + n
    out = {}
    for dz in set(primary) | set(postprim) | set(ft):
        p, pp = primary.get(dz, 0.0), postprim.get(dz, 0.0)
        tert = max(0.0, ft.get(dz, 0.0) - p - pp)
        out[dz] = {"primary": p, "postprimary": pp, "tertiary": tert}
    return out


def _jsonstat(path):
    """Load a JSON-stat SAP file, returning helpers for cheap SA-varying reads."""
    d = json.load(open(path))
    ids, size, val, dim = d["id"], d["size"], d["value"], d["dimension"]
    stride = [1] * len(size)
    for i in range(len(size) - 2, -1, -1):
        stride[i] = stride[i + 1] * size[i + 1]
    def codes(did):
        idx = dim[did]["category"]["index"]
        return list(idx) if isinstance(idx, list) else sorted(idx, key=lambda k: idx[k])
    def label(did):
        return dim[did]["category"]["label"]
    return ids, val, stride, codes, label


def _roi_weights():
    """RoI per-SA spatial-weight shapes {sa: {primary, postprimary, tertiary}}."""
    # age bands
    ids, val, stride, codes, label = _jsonstat(_one(ROI_AGE_GLOB))
    sa_dim, age_dim, sex_dim = "C04172V04943", "C03737V04485", "C03738V04487"
    sa_codes, sa_lab = codes(sa_dim), label(sa_dim)
    sa_stride  = stride[ids.index(sa_dim)]
    age_stride = stride[ids.index(age_dim)]
    bi = codes(sex_dim).index("B") * stride[ids.index(sex_dim)]
    ac = codes(age_dim)
    ai = {b: ac.index(b) * age_stride for b in ("AGE0-4", "AGE5-9", "AGE10-14", "AGE15-19")}
    def age_cell(sa, band):
        v = val[sa * sa_stride + ai[band] + bi]
        return v if isinstance(v, (int, float)) else 0.0
    ages = {}
    for sa, k in enumerate(sa_codes):
        a04, a59, a1014, a1519 = (age_cell(sa, "AGE0-4"), age_cell(sa, "AGE5-9"),
                                  age_cell(sa, "AGE10-14"), age_cell(sa, "AGE15-19"))
        if max(a04, a59, a1014, a1519) > _STATE_THRESHOLD:      # State aggregate row
            continue
        ages[sa_lab[k]] = (a04, a59, a1014, a1519)

    # T10T2 still-studying 15+
    ids, val, stride, codes, label = _jsonstat(_one(ROI_T10T2_GLOB))
    sa_codes, sa_lab = codes(sa_dim), label(sa_dim)
    sa_stride = stride[ids.index(sa_dim)]
    base = (codes(sex_dim).index("B") * stride[ids.index(sex_dim)]
            + codes("C03752V04501").index("SAS") * stride[ids.index("C03752V04501")])
    still = {}
    for sa, k in enumerate(sa_codes):
        v = val[sa * sa_stride + base]
        v = v if isinstance(v, (int, float)) else 0.0
        if v > _STATE_THRESHOLD:
            continue
        still[sa_lab[k]] = v

    out = {}
    for sa, (a04, a59, a1014, a1519) in ages.items():
        primary  = _F04_PRIMARY * a04 + _F59_PRIMARY * a59 + _F1014_PRIMARY * a1014
        postprim = _F1014_POSTPRIM * a1014 + _F1519_POSTPRIM * a1519
        tertiary = max(0.0, still.get(sa, 0.0) - _F1519_POSTPRIM * a1519)   # remove 15-19 secondary
        out[sa] = {"primary": primary, "postprimary": postprim, "tertiary": tertiary}
    return out


# ============================== public API ====================================================
def _scale_to_enrolment(weights, juris, enrol):
    """Scale a jurisdiction's per-level weight shapes so each level sums to its admin enrolment."""
    totals = {lv: sum(w[lv] for w in weights.values()) for lv in _LEVEL_KEYS}
    factor = {lv: (enrol[(juris, lv)] / totals[lv]) if totals[lv] > 0 else 0.0 for lv in _LEVEL_KEYS}
    return {code: {lv: w[lv] * factor[lv] for lv in _LEVEL_KEYS} for code, w in weights.items()}


def load_school_producers():
    """{area_code: {"primary","postprimary","tertiary"}} for all NI DZ + RoI SA, each jurisdiction
    scaled to its admin enrolment (NI:RoI on a real-student footing; island total cancels)."""
    enrol = _enrolment_by_juris_level()
    out = {}
    out.update(_scale_to_enrolment(_roi_weights(), "RoI", enrol))
    out.update(_scale_to_enrolment(_ni_weights(), "NI", enrol))
    return out


if __name__ == "__main__":
    enrol = _enrolment_by_juris_level()
    print("Enrolment by jurisdiction × level (special folded by NI ratio):")
    for j in ("NI", "RoI"):
        print(f"  {j:>4}: " + "  ".join(f"{lv}={enrol[(j,lv)]:,.0f}" for lv in _LEVEL_KEYS))
    isl = island_enrolment_by_level()
    tot = sum(isl.values())
    print("ρ-split (island enrolment ratios): "
          + "  ".join(f"{lv}={isl[lv]:,.0f} ({isl[lv]/tot:.3f})" for lv in _LEVEL_KEYS))

    prod = load_school_producers()
    def jsum(pred):
        return {lv: sum(v[lv] for k, v in prod.items() if pred(k)) for lv in _LEVEL_KEYS}
    ni = jsum(lambda k: str(k).startswith("N"))
    roi = jsum(lambda k: not str(k).startswith("N"))
    print(f"\nProducer sums (should equal enrolment by construction):")
    print(f"  NI : " + "  ".join(f"{lv}={ni[lv]:,.0f}" for lv in _LEVEL_KEYS))
    print(f"  RoI: " + "  ".join(f"{lv}={roi[lv]:,.0f}" for lv in _LEVEL_KEYS))
    print(f"  areas: {len(prod):,}")
    print("\nNI:RoI ratio per level (the cross-border scale that matters):")
    for lv in _LEVEL_KEYS:
        print(f"  {lv:<12}: {ni[lv]/roi[lv]:.3f}")
