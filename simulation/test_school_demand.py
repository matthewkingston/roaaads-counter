"""Unit tests for the school-demand estimator (simulation/school_demand.py).
Run from simulation/:  python3 simulation/test_school_demand.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from school_demand import (assign_enrolments, _INST_TOTAL, _FE_TOTAL, ROI_FE_TOTAL,
                            SCHOOL_ENROLL, SEN_ENROLL)

fails = []


def feat(am, name, x, y, juris="RoI", school=None):
    return dict(amenity=am, name=name, x=float(x), y=float(y), juris=juris, school=school)


def total(kept):
    return sum(e for f, e in kept)


def check(name, got, want, tol=0.5):
    ok = abs(got - want) <= tol
    print(f"  [{'ok' if ok else 'FAIL'}] {name}: got {got:.1f}, want {want:.1f}")
    if not ok:
        fails.append(name)


RP, RS = SCHOOL_ENROLL[("RoI", "primary")], SCHOOL_ENROLL[("RoI", "secondary")]
NP, NS = SCHOOL_ENROLL[("NI", "primary")], SCHOOL_ENROLL[("NI", "secondary")]

# 1. Two distinct-named primary schools at same site → both kept (2× RoI primary)
check("distinct co-located primaries",
      total(assign_enrolments([feat("school", "Junior NS", 0, 0), feat("school", "Senior NS", 5, 5)])),
      2 * RP)

# 2. Named school + unnamed sub-building (mobile) → 1 primary
check("named + unnamed sub-building",
      total(assign_enrolments([feat("school", "St X PS", 0, 0), feat("school", None, 3, 3)])), RP)

# 3. node+way same name → 1
check("same-name duplicate",
      total(assign_enrolments([feat("school", "Dup PS", 0, 0), feat("school", "Dup PS", 2, 2)])), RP)

# 4. Jurisdiction-aware: NI secondary (by name) vs RoI secondary (by school= tag)
check("NI secondary by name",
      total(assign_enrolments([feat("school", "Foyle College", 0, 0, juris="NI")])), NS)
check("RoI secondary by school= tag",
      total(assign_enrolments([feat("school", "Some School", 0, 0, juris="RoI", school="secondary")])), RS)
check("NI primary by school= tag",
      total(assign_enrolments([feat("school", "Some PS", 0, 0, juris="NI", school="primary")])), NP)

# 5. SEN → small value; unknown-name → primary
check("SEN school", total(assign_enrolments([feat("school", "X Special School", 0, 0, school="special_education_needs")])), SEN_ENROLL)
check("unknown name → primary", total(assign_enrolments([feat("school", "Surf School", 0, 0, juris="RoI")])), RP)

# 6. University split: Trinity POIs → total split across them
r = assign_enrolments([feat("university", "Trinity College Dublin", 0, 0),
                       feat("university", "Goldsmith Hall (Trinity College Dublin)", 500, 0),
                       feat("university", "Trinity Biomedical Sciences Institute (Trinity College, Dublin)", 1000, 0)])
check("Trinity split sums to total", total(r), _INST_TOTAL["Trinity"])
check("Trinity per-POI share", r[0][1], _INST_TOTAL["Trinity"] / 3)

# 7. Curated-only: unmatched third-level is DROPPED (no fallback)
check("unmatched university dropped", total(assign_enrolments([feat("university", "Some Random Uni", 0, 0)])), 0)
check("unmatched college dropped", total(assign_enrolments([feat("college", "Random Private Academy", 0, 0)])), 0)

# 8. Non-teaching / junk dropped (curated-only)
check("research station dropped", total(assign_enrolments([feat("university", "Moorepark Food Research Centre", 0, 0)])), 0)
check("accommodation dropped", total(assign_enrolments([feat("university", "Rock Mills Student Accommodation", 0, 0)])), 0)

# 9. kindergarten (pre-school) EXCLUDED
check("kindergarten excluded", total(assign_enrolments([feat("kindergarten", "Little Tots", 0, 0)])), 0)

# 10. college-tagged curated HE institution still matches
check("MIC as college matches curated", total(assign_enrolments([feat("college", "Mary Immaculate College", 0, 0)])), _INST_TOTAL["MIC"])

# 11. Ulster campus split (Belfast campus → its campus FT)
check("Ulster Belfast campus FT",
      total(assign_enrolments([feat("university", "Ulster University, Belfast Campus", 0, 0, juris="NI")])), 10300)

# 12. NI FE college → per-college total (DfE Table A4)
check("Belfast Met FE total",
      total(assign_enrolments([feat("college", "Belfast Metropolitan College", 0, 0, juris="NI")])), _FE_TOTAL["BelfastMet"])

# 13. RoI public FE, single POI → full national total (1 institution × 1 POI, method a)
check("RoI public FE method-a",
      total(assign_enrolments([feat("college", "Ballyfermot College of Further Education", 0, 0)])), ROI_FE_TOTAL)

# 14. mis-tagged second-level 'college' dropped (avoids double-count with admin rolls)
check("community college (2nd-level) dropped",
      total(assign_enrolments([feat("college", "Clonturk Community College", 0, 0)])), 0)

print()
if fails:
    print(f"FAILED {len(fails)}: {fails}"); sys.exit(1)
print("All school_demand tests passed.")
