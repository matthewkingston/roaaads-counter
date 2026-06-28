"""Assign school enrolment (pupils/students) to OSM school POIs — the school-demand
attractor, shared by build_schools.py (which writes the island cache) and, via that
cache, by build_demographics.py (internal core) and build_census_zones.py (external).

Mirrors parking_demand.py but the logic is cross-feature (clustering / institution
splitting), so the entry point `assign_enrolments(features)` takes the FULL island
feature set and returns per-feature enrolment for the kept features. Doing it globally
keeps each third-level institution's total split consistently across all its POIs.

Scheme (see CLAUDE.md / plan):
  • school / kindergarten  → light targeted dedup (drop same-name node+way dupes and
    unnamed sub-buildings within a named campus cluster; keep distinct co-located
    schools), each kept feature → flat fallback (school 300, kindergarten 40).
    [primary/secondary split is a planned follow-up that refines the school value.]
  • college / university   → match name to a curated institution with a SOURCED total
    enrolment (HEA / HESA / institutional); split that total equally across the
    institution's matched POIs (no area-dedup — a dense campus's many POIs then carry
    more weight). Unmatched → flat fallback (university 300, college 700). Obvious
    non-teaching POIs (research stations, accommodation, single ancillary buildings,
    language centres, etc.) are dropped.
"""
import re
import collections
import numpy as np
from scipy.spatial import cKDTree

# Per-type fallback enrolment when there is no curated institution match.
# (school handled separately by jurisdiction-aware primary/secondary classification.)
FALLBACK = {"college": 700, "university": 300, "kindergarten": 40}

# Jurisdiction-aware primary/secondary average enrolment (SOURCED averages = total
# pupils / number of schools; see CLAUDE.md). NI schools are notably larger than RoI.
#   NI  : DE NI 2023/24 — primary 169,463/~800 ≈ 210 ; post-primary 156,399/~190 ≈ 820
#   RoI : DoE 2023/24   — primary 546,758/~3,240 ≈ 170 ; post-primary 416,620/~723 ≈ 575
SCHOOL_ENROLL = {
    ("NI",  "primary"):   210, ("NI",  "secondary"): 820,
    ("RoI", "primary"):   170, ("RoI", "secondary"): 575,
}
SEN_ENROLL = 80   # special-education schools (jurisdiction-agnostic, small)

CLUSTER_M = 150.0   # co-location radius for school/kindergarten dedup

# Primary vs secondary name patterns (NI + RoI). school= tag is preferred when present.
_PRIM_RE = re.compile(
    r"national school|\bN\.?S\b|\bB\.?N\.?S\b|\bG\.?N\.?S\b|primary|\bP\.?S\b|"
    r"gaelscoil|bunscoil|\bscoil\b|infant|nursery|pre-?school", re.I)
_SEC_RE = re.compile(
    r"secondary|\bcollege\b|gramm?ar|community (school|college)|comprehensive|"
    r"high school|post.?primary|col(á|a)iste|pobalscoil|vocational|academy", re.I)


def classify_school(school_tag, name):
    """Return 'primary' | 'secondary' | 'sen' | 'unknown' from the school= tag (preferred)
    then name patterns."""
    st = (school_tag or "").lower()
    if "secondary" in st:
        return "secondary"
    if "primary" in st:
        return "primary"
    if "special" in st or "sen" == st:
        return "sen"
    n = name or ""
    p, s = bool(_PRIM_RE.search(n)), bool(_SEC_RE.search(n))
    if s and not p:
        return "secondary"
    if p and not s:
        return "primary"
    return "unknown"   # → primary value (the majority class)


def school_enrolment(feature):
    """Jurisdiction-aware enrolment for one school/kindergarten feature."""
    if feature["amenity"] == "kindergarten":
        return FALLBACK["kindergarten"]
    cls = classify_school(feature.get("school"), feature.get("name"))
    if cls == "sen":
        return SEN_ENROLL
    key = "secondary" if cls == "secondary" else "primary"   # unknown → primary
    return SCHOOL_ENROLL[(feature.get("juris", "RoI"), key)]

# Curated third-level institutions: (key, total enrolment, name regex). SOURCED figures
# (HEA 2024/25 RoI totals, HESA 2024/25 NI, institutional reports — see CLAUDE.md). Order
# matters: first regex match wins, so keep patterns specific and non-overlapping.
INSTITUTIONS = [
    ("UCD",         34715, r"University College Dublin|\bUCD\b"),
    ("TUDublin",    28500, r"TU Dublin|Technological University,? (of )?Dublin|\bDIT\b|Dublin Institute of Technology"),
    ("Ulster",      32915, r"Ulster University|University of Ulster|Magee College|Jordanstown"),
    ("QUB",         25080, r"Queen.?s University|\bQUB\b"),
    ("UCC",         24195, r"University College Cork|\bUCC\b"),
    ("Trinity",     20490, r"Trinity College|\bTCD\b|The Lir"),
    ("ATU",         20418, r"Atlantic Technological|\bATU\b|Galway[- ]Mayo|\bGMIT\b|Letterkenny Institute"),
    ("DCU",         20377, r"Dublin City University|\bDCU\b"),
    ("Galway",      19335, r"University of Galway|NUI,? Galway|\bNUIG\b|National University of Ireland,? Galway|Cairnes"),
    ("SETU",        18338, r"South East Technological|\bSETU\b|Waterford Institute of Technology|\bWIT\b|Institute of Technology,? Carlow"),
    ("MTU",         18000, r"Munster Technological|\bMTU\b|Cork Institute of Technology|\bCIT\b|Crawford College|National Maritime College|Institute of Technology,? Tralee"),
    ("Limerick",    18000, r"University of Limerick"),
    ("Maynooth",    16000, r"Maynooth|Froebel"),
    ("TUS",         14000, r"Technological University of the Shannon|\bTUS\b|Athlone Institute|\bAIT\b|Limerick Institute of Technology|\bLIT\b|Limerick School of Art"),
    ("DKIT",         5400, r"Dundalk Institute|\bDk?IT\b"),
    ("MIC",          5000, r"Mary Immaculate"),
    ("RCSI",         4094, r"Royal College of Surgeons|\bRCSI\b"),
    ("IADT",         2500, r"Institute of Art,? Design|\bIADT\b|D(ú|u)n Laoghaire Institute"),
    ("Stranmillis",  1430, r"Stranmillis"),
    ("StMarysBel",   1060, r"St\.? Mary.?s University College"),
]
_INST = [(k, n, re.compile(p, re.I)) for k, n, p in INSTITUTIONS]
_INST_TOTAL = {k: n for k, n, _ in INSTITUTIONS}

# Obvious non-teaching POIs to drop (only applied to UNMATCHED third-level features —
# curated institutions are never dropped).
_NONTEACHING = re.compile(
    r"research (station|centre|center|facility)|Mace Head|Moorepark|MaREI|Food Research|"
    r"Reservoir House|Riddel Hall|Merriman House|\b(Á|A)ras\b|Enterprise Centre|Questum|"
    r"Accomm?odation|Halls? of Residence|Student Village|Law Society|Garda College|"
    r"Sports (Field|Ground|Centre|Broombridge)|Nurse Education|Medical Academy",
    re.I,
)


def _norm(name):
    return re.sub(r"[^a-z0-9]", "", name.lower()) if name else None


def _match_institution(name):
    if not name:
        return None
    for key, _total, rx in _INST:
        if rx.search(name):
            return key
    return None


def _union_find_clusters(xy, radius):
    if len(xy) == 0:
        return []
    tree = cKDTree(xy)
    parent = list(range(len(xy)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in tree.query_pairs(radius):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    groups = collections.defaultdict(list)
    for i in range(len(xy)):
        groups[find(i)].append(i)
    return list(groups.values())


def _dedup_schools(feats):
    """Light targeted dedup for school/kindergarten. Returns [(feat, enrolment)]."""
    if not feats:
        return []
    xy = np.array([[f["x"], f["y"]] for f in feats])
    out = []
    for cluster in _union_find_clusters(xy, CLUSTER_M):
        named = {}      # normalised name -> feature index (keep one per distinct name)
        unnamed = []
        for i in cluster:
            nm = _norm(feats[i]["name"])
            if nm:
                named.setdefault(nm, i)
            else:
                unnamed.append(i)
        keep = list(named.values())
        if not named and unnamed:           # all-unnamed cluster → keep one
            keep = [unnamed[0]]
        # (unnamed dropped when a named feature shares the cluster: sub-buildings/mobiles)
        for i in keep:
            out.append((feats[i], school_enrolment(feats[i])))
    return out


def _assign_tertiary(feats):
    """Curated-institution split + fallback + non-teaching drop. Returns [(feat, enrol)]."""
    groups = collections.defaultdict(list)
    out = []
    for f in feats:
        inst = _match_institution(f["name"])
        if inst:
            groups[inst].append(f)
        elif f["name"] and _NONTEACHING.search(f["name"]):
            continue                        # drop obvious non-teaching
        else:
            out.append((f, FALLBACK[f["amenity"]]))
    for inst, members in groups.items():
        share = _INST_TOTAL[inst] / len(members)
        for f in members:
            out.append((f, share))
    return out


def assign_enrolments(features):
    """features: list of dicts with keys amenity, name, x, y (x/y projected metres).
    Returns list of (feature, enrolment) for the KEPT features."""
    sch = [f for f in features if f["amenity"] in ("school", "kindergarten")]
    ter = [f for f in features if f["amenity"] in ("college", "university")]
    return _dedup_schools(sch) + _assign_tertiary(ter)
