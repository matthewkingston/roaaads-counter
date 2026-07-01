"""Assign school enrolment (pupils/students) to OSM school POIs — the school-demand
attractor, shared by build_schools.py (which writes the island cache) and, via that
cache, by build_demographics.py (internal core) and build_census_zones.py (external).

Mirrors parking_demand.py but the logic is cross-feature (clustering / institution
splitting), so the entry point `assign_enrolments(features)` takes the FULL island
feature set and returns per-feature enrolment for the kept features. Doing it globally
keeps each third-level institution's total split consistently across all its POIs.

Scheme (see CLAUDE.md / plan):
  • school                  → light targeted dedup (drop same-name node+way dupes and
    unnamed sub-buildings within a named campus cluster; keep distinct co-located
    schools), each kept feature valued by jurisdiction-aware primary/secondary average.
    **Kindergartens (pre-school) are excluded** — the school component excludes pre-school
    on both sides (admin-roll school-age drops nursery; the producer subtracts childcare).
  • college / university   → **curated-only** third-level (see `_assign_tertiary`): match to a
    curated HE institution (SOURCED full-time enrolment — HEA 2024/25 RoI / DfE-HESA 2023/24 NI),
    a curated FE college (NI DfE Table A4 per-college + Teagasc/CAFRE agri), or the RoI public-FE
    keep-set (national SOLAS total distributed by method (a)). Each institution's total is split
    across its matched POIs (HE Ulster split by campus). Anything unmatched is DROPPED — OSM junk,
    part-time FET/adult, and second-level 'colleges' already counted in the admin-roll school-age.
"""
import re
import collections
import numpy as np
from scipy.spatial import cKDTree

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
    """Jurisdiction-aware enrolment for one school feature (kindergartens are excluded upstream)."""
    cls = classify_school(feature.get("school"), feature.get("name"))
    if cls == "sen":
        return SEN_ENROLL
    key = "secondary" if cls == "secondary" else "primary"   # unknown → primary
    return SCHOOL_ENROLL[(feature.get("juris", "RoI"), key)]

# Curated third-level (HE) institutions: (key, full-time enrolment, name regex). Order matters:
# first regex match wins, so keep patterns specific and non-overlapping.
#   RoI: full-time 2024/25 (HEA Students Summary by Institute × Mode); national total 215,580.
#   NI:  full-time 2023/24 (DfE/HESA "Enrolments 2324 tables" Table 8b). Ulster is island-only —
#        GB campuses (Birmingham/London) excluded, split by campus via CAMPUS_FT from Table 8c
#        (Belfast 10,300 + Coleraine 2,925 + Magee 4,915 = 18,140).
INSTITUTIONS = [
    # RoI: full-time enrolment 2024/25 (HEA Students Summary by Institute × Mode of Study).
    ("UCD",         27255, r"University College Dublin|\bUCD\b"),
    ("TUDublin",    19090, r"TU Dublin|Technological University,? (of )?Dublin|\bDIT\b|Dublin Institute of Technology"),
    ("UCC",         20055, r"University College Cork|\bUCC\b"),
    ("Trinity",     19225, r"Trinity College|\bTCD\b|The Lir"),
    ("Galway",      17175, r"University of Galway|NUI,? Galway|\bNUIG\b|National University of Ireland,? Galway|Cairnes"),
    ("DCU",         17205, r"Dublin City University|\bDCU\b"),
    ("Limerick",    16030, r"University of Limerick"),
    ("Maynooth",    14025, r"Maynooth|Froebel"),
    ("ATU",         13885, r"Atlantic Technological|\bATU\b|Galway[- ]Mayo|\bGMIT\b|Letterkenny Institute"),
    ("MTU",         11165, r"Munster Technological|\bMTU\b|Cork Institute of Technology|\bCIT\b|Crawford College|National Maritime College|Institute of Technology,? Tralee"),
    ("SETU",        10745, r"South East Technological|\bSETU\b|Waterford Institute of Technology|\bWIT\b|Institute of Technology,? Carlow"),
    ("TUS",          9620, r"Technological University of the Shannon|\bTUS\b|Athlone Institute|\bAIT\b|Limerick Institute of Technology|\bLIT\b|Limerick School of Art"),
    ("DKIT",         4160, r"Dundalk Institute|\bDk?IT\b"),
    ("MIC",          4030, r"Mary Immaculate"),
    ("NCI",          4040, r"National College of Ireland"),
    ("RCSI",         3475, r"Royal College of Surgeons|\bRCSI\b"),
    ("IADT",         2015, r"Institute of Art,? Design|\bIADT\b|D(ú|u)n Laoghaire Institute"),
    ("NCAD",         1220, r"National College of Art"),
    ("StAngela",     1170, r"St\.? Angela.?s College"),
    # Private/specialist HEIs not in the HEA return. Marino/RIAM full-time sourced; Griffith/DBS
    # publish no FT/PT split → total × HE full-time ratio 0.773 (likely an upper bound; PT-heavy).
    ("DBS",          6950, r"Dublin Business School|\bDBS\b"),
    ("Griffith",     6200, r"Griffith College"),
    ("MarinoIoE",     750, r"Marino Institute of Education"),
    ("RIAM",          200, r"Royal Irish Academy of Music|\bRIAM\b"),
    # NI: full-time enrolment 2023/24 (DfE/HESA). Ulster is island-only (GB campuses excluded)
    # and split by campus via CAMPUS_FT.
    ("Ulster",      18140, r"Ulster University|University of Ulster|Magee College|Jordanstown"),
    ("QUB",         21345, r"Queen.?s University|\bQUB\b"),
    ("Stranmillis",   875, r"Stranmillis"),
    ("StMarysBel",    880, r"St\.? Mary.?s University College"),
]
_INST = [(k, n, re.compile(p, re.I)) for k, n, p in INSTITUTIONS]
_INST_TOTAL = {k: n for k, n, _ in INSTITUTIONS}

# Per-campus full-time enrolment for institutions with a known campus split (else the institution
# total is spread evenly across its POIs). Each institution → [(campus regex, FT enrolment)]; a POI
# matching a campus regex takes that campus's enrolment (÷ POIs matching it); an institution POI
# matching NO campus regex gets 0 (stale/closed campus, e.g. Ulster Jordanstown → moved to Belfast).
CAMPUS_FT = {
    "Ulster": [                                    # DfE/HESA Table 8c, 2023/24, island campuses only
        (re.compile(r"Belfast", re.I), 10300),
        (re.compile(r"Coleraine", re.I), 2925),
        (re.compile(r"Magee", re.I), 4915),
    ],
}

# ── Further Education ──────────────────────────────────────────────────────────
# FE colleges with a per-college full-time total (split across their matched campus POIs, like
# HE). NI: DfE Table A4 FT 2023/24. Agri colleges (Teagasc RoI + CAFRE NI) are funded separately
# (not in the SOLAS/DfE FE totals) → included here (tiny component). Kildalton ~1,200 and CAFRE
# 523 (FT HE 2023/24) sourced; other Teagasc colleges ~400 est. NB CAFRE's 523 are on UU/QUB-
# validated courses so HESA may already count them under Ulster/QUB — a ~523 double-count accepted
# to place them at the CAFRE campuses.
FE_INSTITUTIONS = [
    ("BelfastMet",   4840, r"Belfast Metropolitan|Titanic Quarter|Castlereagh Camp|Castlereagh College|Millfield Campus|Springvale|\be3\b"),
    ("NorthernReg",  3325, r"Northern Regional College"),
    ("NorthWestReg", 2610, r"North[- ]?West Regional College"),
    ("SouthEastReg", 4465, r"South Eastern Regional College|\bSERC\b"),
    ("SouthernReg",  4080, r"Southern Regional College"),
    ("SouthWest",    3620, r"South West College"),
    ("CAFRE",         523, r"CAFRE|Greenmount|Loughry"),            # NI agri FT HE 2023/24 (sourced); Enniskillen→SWC
    ("Kildalton",    1200, r"Kildalton"),                           # RoI Teagasc (sourced)
    ("ClonakiltyAg",  400, r"Clonakilty Agric"),                    # RoI Teagasc (est)
    ("Ballyhaise",    400, r"Ballyhaise"),                          # RoI Teagasc (est)
    ("MountbellewAg", 350, r"Mountbellew"),                         # RoI agri (est)
    ("SalesianAg",    200, r"Salesian Agric"),                      # RoI agri (est)
]
_FE = [(k, n, re.compile(p, re.I)) for k, n, p in FE_INSTITUTIONS]
_FE_TOTAL = {k: n for k, n, _ in FE_INSTITUTIONS}

# RoI public FE (Colleges of FE / FET Colleges / named institutes) — only a NATIONAL full-time
# total is published (SOLAS 2023), so it is distributed by method (a) total→institution→POI (see
# CLAUDE.md; POI-count is a mapping artifact, so per-POI splitting is avoided). Keep pattern first,
# then exclude part-time/adult (FET/Youthreach/VTOS/BTEI/training/adult) which are not the daily
# full-time cohort. Private colleges and mis-tagged second-level fall through → dropped.
ROI_FE_TOTAL = 65851
_ROI_FE_KEEP = re.compile(
    r"College of Further Education|Institute of Further Education|Further Education Institute|"
    r"\bCollege of FET\b|\bFET College\b|Cavan Institute|Monaghan Institute|Galway Technical Institute|"
    r"Limerick City College|Liberties College|Pearse College|Crumlin College|Killester|Plunket College|"
    r"Rathmines College|Ballsbridge College|Whitehall College|Inchicore College|Ringsend College|"
    r"Ballyfermot College|Bray Institute|Col(á|a)iste (Dh(ú|u)laigh|(Í|I)de|Stiof(á|a)in)", re.I)
_ROI_FE_EXCL = re.compile(
    r"Youthreach|Training Centre|Training Centres|\bVTOS\b|\bBTEI\b|Adult (Education|Learning)|"
    r"Education Centre|Learning Centre|Community College", re.I)

def _fe_instkey(nm):
    """Normalise an RoI FE POI name to its institution (drop campus suffixes) for method (a)."""
    s = re.sub(r"\s*[-,(].*$", "", nm)                       # drop after '-', ',', '('
    s = re.sub(r"\b(Central |Main )?Campus\b.*$", "", s, flags=re.I)
    return re.sub(r"[^a-z0-9]", "", s.lower())

def _match_fe(name):
    for key, _t, rx in _FE:
        if rx.search(name):
            return key
    return None

def _is_roi_public_fe(name):
    return bool(_ROI_FE_KEEP.search(name)) and not _ROI_FE_EXCL.search(name)




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
    """Third-level (college/university) enrolment. Curated-only — anything not matched to a
    curated HE institution, a curated FE college, or the RoI public-FE keep-set is DROPPED
    (removes the old flat fallback: OSM junk, part-time FET/adult, and mis-tagged second-level
    'colleges' already in the admin school-age rolls). Returns [(feat, enrol)].
      HE / FE-named: per-institution total split across matched POIs (HE Ulster via CAMPUS_FT).
      RoI public FE: national total distributed by method (a) total→institution→POI."""
    he = collections.defaultdict(list)
    fe = collections.defaultdict(list)
    roi_fe = []
    out = []
    for f in feats:
        nm = f["name"] or ""
        inst = _match_institution(nm)
        if inst:
            he[inst].append(f); continue
        fk = _match_fe(nm)
        if fk:
            fe[fk].append(f); continue
        if _is_roi_public_fe(nm):
            roi_fe.append(f); continue
        # else: DROP (junk / part-time FET-adult / mis-tagged second-level / private / unmatched)

    for inst, members in he.items():
        campuses = CAMPUS_FT.get(inst)
        if campuses:                              # per-campus (Ulster); POI matching none → 0
            for f in members:
                nm = f["name"] or ""
                val = 0.0
                for rx, enrol in campuses:
                    if rx.search(nm):
                        val = enrol / sum(1 for g in members if rx.search(g["name"] or ""))
                        break
                out.append((f, val))
        else:
            share = _INST_TOTAL[inst] / len(members)
            for f in members:
                out.append((f, share))

    for fk, members in fe.items():                # FE with per-college total
        share = _FE_TOTAL[fk] / len(members)
        for f in members:
            out.append((f, share))

    if roi_fe:                                    # RoI public FE: method (a)
        inst_groups = collections.defaultdict(list)
        for f in roi_fe:
            inst_groups[_fe_instkey(f["name"])].append(f)
        per_inst = ROI_FE_TOTAL / len(inst_groups)
        for members in inst_groups.values():
            per_poi = per_inst / len(members)
            for f in members:
                out.append((f, per_poi))
    return out


def assign_enrolments(features):
    """features: list of dicts with keys amenity, name, x, y (x/y projected metres).
    Returns list of (feature, enrolment) for the KEPT features."""
    sch = [f for f in features if f["amenity"] == "school"]   # kindergartens (pre-school) excluded
    ter = [f for f in features if f["amenity"] in ("college", "university")]
    return _dedup_schools(sch) + _assign_tertiary(ter)
