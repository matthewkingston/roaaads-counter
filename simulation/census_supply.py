"""Per-small-area census trip *producers* (supply), harmonised NI (DZ) + RoI (SA).

`load_supply() -> {area_code: {"commute": float, "school": float}}` for every NI Data Zone
and RoI Small Area, consumed by build_census_zones.py (external zones) and
build_demographics.py (internal core), where they become the producing weights for the
commute and school gravity components (school wired now; commute stored for the
business→commute/retail split).

Definitions
-----------
commute_producers = resident workers who physically *drive* to work — car-driver modes only
  (the model assigns car flow), matching the vehicle-driver modes in
  analysis/derive_generation_rates.py.  WFH and not-in-employment are excluded structurally —
  they are their own travel-method categories and are simply not selected:
    RoI  SAP2022 T11T1, statistic "travel to work": Car Driver + Van + Motorcycle/scooter
         (RoI has no separate Taxi mode — it sits in "Other (incl. lorry)" and is excluded).
    NI   transport_to_workplace: "Driving a car or van" + "Motorcycle, scooter or moped"
         + "Taxi".

school_producers = student headcount (school-age + third-level — matches the school attractor,
  which is admin-roll school-age + OSM third-level; pre-school childcare excluded both sides):
    RoI  SAP2022 T11T1, statistic "travel to school, college or childcare": Total, MINUS an
         age-0-4-distributed childcare estimate (CHILDCARE_NATIONAL — see below).
    NI   in_full_time_education: "Full-time student or schoolchild" (under-compulsory-age kids are
         "No code required", so childcare is already excluded — no subtraction needed).

RoI key = the JSON-stat SA category *label* (e.g. "057103001", incl. "/NN" split composites),
which equals `SA_PUB2022` / `area_code` in ingest_roi_census. The category *index* (a CSO GUID)
is ignored.

Note: the JSON-stat carries a national "State" aggregate row (~2M, equal to the sum of all
SAs), which `_roi_supply` drops — so the per-SA values are a clean 1× and need no scaling.
"""
import glob
import json
import csv

ROI_TRAVEL_GLOB   = "data/ireland_census/SAP2022T11T1SA*.json"
ROI_AGE_GLOB      = "data/ireland_census/SAP2022T1T1A*.json"   # age-band pop (childcare proxy)
NI_TRANSPORT_GLOB = "data/ni_census/*transport_to_workplace*.csv"
NI_EDU_GLOB       = "data/ni_census/*in_full_time_education*.csv"

# RoI school producer: C02 ("travel to school, college or childcare") lumps childcare (pre-school)
# in with school/college; the model's school component must exclude it. Subtract a national
# childcare total distributed across SAs ∝ age 0-4 population (SAP2022T1T1A). Anchored to the
# CSO Census 2022 figure: ~160k children aged 0-4 travelled to childcare *or* school, less the
# ~30k age-4 junior-infants who are school (kept — they are in the admin-roll attractor) ⇒ ~130k.
# (ECCE-only would be ~106.5k but misses sub-ECCE crèche/childminder under-3s.) DOCUMENTED
# LIMITATION: the finest census age band is 5-year (0-4), so the age-3-4 pre-school cohort is only
# approximated. NI needs no equivalent — in_full_time_education excludes under-compulsory-age kids.
CHILDCARE_NATIONAL = 130_000

def _one(glob_pat):
    hits = glob.glob(glob_pat)
    if not hits:
        raise FileNotFoundError(f"No census file matching {glob_pat}")
    return hits[0]


def _roi_pop_0_4():
    """{SA label: population aged 0-4 (both sexes)} from SAP2022T1T1A, plus the national total.
    Drops the national 'State' aggregate row (>100k)."""
    d = json.load(open(_one(ROI_AGE_GLOB)))
    size, val, dim = d["size"], d["value"], d["dimension"]
    def cat(did):
        c = dim[did]["category"]; idx = c["index"]
        return (list(idx) if isinstance(idx, list) else sorted(idx, key=lambda k: idx[k])), c["label"]
    sa_keys, sa_lab = cat("C04172V04943")
    age_keys, _     = cat("C03737V04485")
    sex_keys, _     = cat("C03738V04487")
    nSA, nAGE, nSEX = size[2], size[3], size[4]   # leading STATISTIC/TLIST dims are size 1
    ai = age_keys.index("AGE0-4"); bi = sex_keys.index("B")
    out, tot = {}, 0.0
    for sa, k in enumerate(sa_keys):
        v = val[(sa * nAGE + ai) * nSEX + bi]
        v = v if isinstance(v, (int, float)) else 0.0
        if v > 100_000:        # State aggregate
            continue
        out[sa_lab[k]] = v; tot += v
    return out, tot


def _roi_supply():
    """RoI per-SA producers from SAP2022 T11T1. → {area_code: (commute, school)}.
    school = C02 (travel to school/college/childcare) minus an age-0-4-distributed childcare
    estimate (CHILDCARE_NATIONAL)."""
    d = json.load(open(_one(ROI_TRAVEL_GLOB)))
    size, val, dim, idd = d["size"], d["value"], d["dimension"], d["id"]

    def cat(did):
        c = dim[did]["category"]; idx = c["index"]
        keys = list(idx) if isinstance(idx, list) else sorted(idx, key=lambda k: idx[k])
        return keys, c["label"]

    stat_keys, stat_lab = cat("STATISTIC")
    sa_keys,   sa_lab   = cat("C04172V04943")
    mean_keys, mean_lab = cat("C03767V04516")
    nSA, nM = size[2], size[3]

    def find(keys, lab, pred):
        for i, k in enumerate(keys):
            if pred(lab[k].lower()):
                return i
        raise KeyError("category not found")

    si_work = find(stat_keys, stat_lab, lambda s: "travel to work" in s and "school" not in s)
    si_sch  = find(stat_keys, stat_lab, lambda s: "school" in s)
    mi_tot  = find(mean_keys, mean_lab, lambda s: s == "total")
    # Car-driver modes (the model assigns car flow). No separate Taxi in RoI.
    mi_cd   = find(mean_keys, mean_lab, lambda s: s == "car driver")
    mi_van  = find(mean_keys, mean_lab, lambda s: s == "van")
    mi_moto = find(mean_keys, mean_lab, lambda s: "motorcycle" in s)

    def cell(si, sa, mi):
        v = val[(si * nSA + sa) * nM + mi]
        return v if isinstance(v, (int, float)) else 0.0

    # Childcare subtraction: distribute CHILDCARE_NATIONAL across SAs ∝ age 0-4 population.
    pop04, pop04_tot = _roi_pop_0_4()
    f_childcare = (CHILDCARE_NATIONAL / pop04_tot) if pop04_tot > 0 else 0.0

    out = {}
    for sa, k in enumerate(sa_keys):
        commute = (cell(si_work, sa, mi_cd) + cell(si_work, sa, mi_van)
                   + cell(si_work, sa, mi_moto))
        school_raw = cell(si_sch, sa, mi_tot)
        # The SA dimension includes one national "State" aggregate row (~2M) equal to the
        # sum of all real SAs; a real SA never exceeds a few thousand. Drop it so the
        # output is clean per-SA (and naive totals aren't 2×).
        if commute > 100_000 or school_raw > 100_000:
            continue
        childcare = f_childcare * pop04.get(sa_lab[k], 0.0)
        school = max(0.0, school_raw - childcare)        # exclude pre-school childcare
        out[sa_lab[k]] = (max(commute, 0.0), school)
    return out


def _ni_csv(path):
    """Yield (dz_code, label, count) from a NISRA people-by-X CSV."""
    with open(path, newline="") as f:
        r = csv.reader(f)
        next(r)  # header
        for row in r:
            if len(row) < 5:
                continue
            try:
                yield row[0], row[3], float(row[4])
            except ValueError:
                continue


def _ni_supply():
    """NI per-DZ producers from transport_to_workplace + in_full_time_education (1×)."""
    commute = {}
    # Car-driver modes only (WFH/no-code are separate labels, structurally excluded).
    _CAR = ("driving a car or van", "motorcycle, scooter or moped", "taxi")
    for dz, label, n in _ni_csv(_one(NI_TRANSPORT_GLOB)):
        if label.lower() not in _CAR:
            continue
        commute[dz] = commute.get(dz, 0.0) + n   # car driver + motorcycle + taxi
    school = {}
    for dz, label, n in _ni_csv(_one(NI_EDU_GLOB)):
        if label.lower() == "full-time student or schoolchild":
            school[dz] = school.get(dz, 0.0) + n
    out = {}
    for dz in set(commute) | set(school):
        out[dz] = (commute.get(dz, 0.0), school.get(dz, 0.0))
    return out


def load_supply():
    """{area_code: {"commute": float, "school": float}} for all NI DZ + RoI SA."""
    out = {}
    for code, (c, s) in {**_roi_supply(), **_ni_supply()}.items():
        out[code] = {"commute": c, "school": s}
    return out


if __name__ == "__main__":
    _p04, _p04tot = _roi_pop_0_4()
    print(f"RoI age 0-4: {_p04tot:,.0f}; childcare target {CHILDCARE_NATIONAL:,} "
          f"(f={CHILDCARE_NATIONAL/_p04tot:.3f})")
    roi = _roi_supply()
    ni = _ni_supply()
    rc = sum(c for c, s in roi.values()); rs = sum(s for c, s in roi.values())
    nc = sum(c for c, s in ni.values()); ns = sum(s for c, s in ni.values())
    print(f"RoI SAs: {len(roi):,}  commute={rc:,.0f}  school={rs:,.0f}")
    print(f"NI DZs:  {len(ni):,}  commute={nc:,.0f}  school={ns:,.0f}")
    print(f"TOTAL    commute={rc+nc:,.0f}  school={rs+ns:,.0f}")
