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

school_producers = student headcount (all ages, incl. third-level — matches the
  university-inclusive school attractor):
    RoI  SAP2022 T11T1, statistic "travel to school, college or childcare": Total.
    NI   in_full_time_education: "Full-time student or schoolchild".

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
NI_TRANSPORT_GLOB = "data/ni_census/*transport_to_workplace*.csv"
NI_EDU_GLOB       = "data/ni_census/*in_full_time_education*.csv"

def _one(glob_pat):
    hits = glob.glob(glob_pat)
    if not hits:
        raise FileNotFoundError(f"No census file matching {glob_pat}")
    return hits[0]


def _roi_supply():
    """RoI per-SA producers from SAP2022 T11T1 (÷2). → {area_code: (commute, school)}."""
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

    out = {}
    for sa, k in enumerate(sa_keys):
        commute = (cell(si_work, sa, mi_cd) + cell(si_work, sa, mi_van)
                   + cell(si_work, sa, mi_moto))
        school  = cell(si_sch, sa, mi_tot)
        # The SA dimension includes one national "State" aggregate row (~2M) equal to the
        # sum of all real SAs; a real SA never exceeds a few thousand. Drop it so the
        # output is clean per-SA (and naive totals aren't 2×).
        if commute > 100_000 or school > 100_000:
            continue
        out[sa_lab[k]] = (max(commute, 0.0), max(school, 0.0))
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
    roi = _roi_supply()
    ni = _ni_supply()
    rc = sum(c for c, s in roi.values()); rs = sum(s for c, s in roi.values())
    nc = sum(c for c, s in ni.values()); ns = sum(s for c, s in ni.values())
    print(f"RoI SAs: {len(roi):,}  commute={rc:,.0f}  school={rs:,.0f}")
    print(f"NI DZs:  {len(ni):,}  commute={nc:,.0f}  school={ns:,.0f}")
    print(f"TOTAL    commute={rc+nc:,.0f}  school={rs+ns:,.0f}")
