"""Per-small-area census trip *producers* (supply), harmonised NI (DZ) + RoI (SA).

`load_supply() -> {area_code: {"commute": float}}` for every NI Data Zone and RoI Small Area,
consumed by build_census_zones.py (external zones) and build_demographics.py (internal core), where
it becomes the producing weight for the commute gravity component.

(School producers moved to census_school_producers.py in Phase-2 — per-level
primary/post-primary/tertiary, enrolment-anchored — replacing the old lumped C02-minus-childcare
producer that used to live here. This module now provides the commute producer plus the shared
stdlib helpers `_one` / `_ni_csv`, which census_school_producers imports.)

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

RoI key = the JSON-stat SA category *label* (e.g. "057103001", incl. "/NN" split composites),
which equals `SA_PUB2022` / `area_code` in ingest_roi_census. The category *index* (a CSO GUID)
is ignored. The JSON-stat carries a national "State" aggregate row (~2M, = the sum of all SAs),
which `_roi_supply` drops so the per-SA values are a clean 1×.
"""
import glob
import json
import csv

ROI_TRAVEL_GLOB   = "data/ireland_census/SAP2022T11T1SA*.json"
NI_TRANSPORT_GLOB = "data/ni_census/*transport_to_workplace*.csv"


def _one(glob_pat):
    hits = glob.glob(glob_pat)
    if not hits:
        raise FileNotFoundError(f"No census file matching {glob_pat}")
    return hits[0]


def _roi_supply():
    """RoI per-SA car-commute producers from SAP2022 T11T1 → {area_code: commute}.
    commute = "travel to work" Car Driver + Van + Motorcycle (car-driver modes; the model assigns
    car flow — RoI has no separate Taxi mode)."""
    d = json.load(open(_one(ROI_TRAVEL_GLOB)))
    size, val, dim = d["size"], d["value"], d["dimension"]

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
        # Drop the national "State" aggregate SA row (~2M; a real SA is a few thousand).
        if commute > 100_000:
            continue
        out[sa_lab[k]] = max(commute, 0.0)
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
    """NI per-DZ car-commute producers from transport_to_workplace → {dz: commute}.
    Car-driver modes only (WFH/no-code are separate labels, structurally excluded)."""
    commute = {}
    _CAR = ("driving a car or van", "motorcycle, scooter or moped", "taxi")
    for dz, label, n in _ni_csv(_one(NI_TRANSPORT_GLOB)):
        if label.lower() in _CAR:
            commute[dz] = commute.get(dz, 0.0) + n   # car driver + motorcycle + taxi
    return commute


def load_supply():
    """{area_code: {"commute": float}} for all NI DZ + RoI SA."""
    out = {}
    for code, c in {**_roi_supply(), **_ni_supply()}.items():
        out[code] = {"commute": c}
    return out


if __name__ == "__main__":
    roi = _roi_supply()
    ni = _ni_supply()
    rc, nc = sum(roi.values()), sum(ni.values())
    print(f"RoI SAs: {len(roi):,}  commute={rc:,.0f}")
    print(f"NI DZs:  {len(ni):,}  commute={nc:,.0f}")
    print(f"TOTAL    commute={rc + nc:,.0f}")
