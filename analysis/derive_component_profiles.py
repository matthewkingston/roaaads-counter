"""
Derive per-component temporal SHAPE profiles for the four-component gravity model.

Under generation pinning (analysis/derive_generation_rates.py) each component's
absolute magnitude AND the inter-component split are set by data (ρ_c × producers,
K_c ≈ 1).  The temporal profiles f_c therefore carry ONLY each component's relative
time-of-day-and-week SHAPE, derived independently — they no longer partition a
shared aggregate.

For each component c the weekly profile is
    f_c(dow, h) = V_c(dow) · H_c(daytype(dow), h)
  - H_c = the within-day hourly shape (Σ_h H = 1):
       weekday:  ρ-weighted blend of the component's purposes' weekday hourly
                 distributions (NTS0502a × the aggregate hourly volume);
       weekend:  the aggregate weekend hourly shape (no per-purpose weekend hourly
                 data exists), shared across components.
  - V_c(dow) = the component's relative daily volume Mon–Sun, a ρ-weighted blend of
       its purposes' day-of-week distributions (NTS0504b), normalised so
       Σ_{7 days} V_c = 7  ⇒  the day-weighted daily sum W_c = 1.

W_c = 1 makes K_c·m_c the component's daily AADT directly (so K_c ≈ 1 with the
generation-pinned m_c) and matches the existing "rows sum to 7" convention; the
magnitude/split lives entirely in generation.  The old res+commute+retail+school =
agg partition constraint is intentionally DROPPED.

Mapping: analysis/purpose_mapping.py (shared with derive_generation_rates.py).
ρ_p (per-purpose car-driver rates) come from analysis/generation_rates.json.

Sources (DfT NTS, England, 2023/24 rolling avg):
  NTS0502a (data/nts0502.ods) — weekday trip start-time × purpose; rows sum to 100,
       i.e. P(purpose | hour).
  NTS0504b (data/nts0504.ods) — trips/person/year by day-of-week × purpose (Mon–Sun).
NTS0502a merges some purposes, so each purpose's weekday hourly SHAPE uses a
best-available proxy column (PROXY_0502 below) — an approximation in the shape only
(magnitude is exact from generation).

Usage:  python3 analysis/derive_component_profiles.py
Overwrites mean_fraction_{res,commute,retail,school_primary,school_postprimary,school_tertiary}
in analysis/hourly_fractions.csv (primary/post-primary = escort-education shape; tertiary = commute
shape — three explicit copied columns).
Re-run when the NTS files, the purpose mapping, or generation_rates.json change.
"""

import csv, json, os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from purpose_mapping import COMPONENT_PURPOSES, COMPONENTS, CANONICAL_PURPOSES

NTS0502_FILE = "data/nts0502.ods"
NTS0504_FILE = "data/nts0504.ods"
FRACS_FILE   = "analysis/hourly_fractions.csv"
GEN_RATES    = "analysis/generation_rates.json"
NTS_YEAR     = "2023 to 2024"

# canonical purpose -> NTS0502a column for its weekday hourly SHAPE proxy.  The table
# merges purposes, so several share a proxy column (shape only; magnitude is from
# generation): personal_business + other_escort -> col 7 "Other work, other escort and
# personal business"; leisure -> col 8 "Visiting friends, entertainment and sport"
# (dominant leisure); other -> col 9 "Holiday, day trip and other".  Education (col 4,
# the non-vehicle child trip) is NOT used — school uses the Escort-education shape (col 5).
PROXY_0502 = {
    "commuting": 2, "business": 3, "education_escort": 5, "shopping": 6,
    "personal_business": 7, "other_escort": 7, "leisure": 8, "other": 9,
}

# canonical purpose -> NTS0504b column(s) for its day-of-week volume.
DOW_0504 = {
    "commuting": [2], "business": [3], "education_escort": [5], "shopping": [6],
    "other_escort": [7], "personal_business": [8],
    "leisure": [9, 10, 11, 12],   # visit-home + visit-elsewhere + sport/ent + holiday
    "other": [14],
}
_DAY_ORDER = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
              "Saturday", "Sunday"]


def _load_0502_pct():
    """{col_idx: [P(purpose|hour) for h in 0..23]} from NTS0502a (rows sum to 100)."""
    raw = pd.read_excel(NTS0502_FILE, sheet_name="NTS0502a_start_time_by_purpose",
                        header=None, engine="odf")
    d = raw.iloc[6:].copy(); d.columns = range(d.shape[1])
    yr = d[d[0] == NTS_YEAR].copy()
    yr = yr[yr[1] != "All day"].copy()
    if len(yr) != 24:
        sys.exit(f"ERROR: expected 24 hourly rows in NTS0502a '{NTS_YEAR}', got {len(yr)}")
    yr["hour"] = yr[1].str[:2].astype(int)
    yr = yr.sort_values("hour").reset_index(drop=True)
    cols = sorted(set(PROXY_0502.values()))
    for c in cols:
        yr[c] = pd.to_numeric(yr[c], errors="coerce").fillna(0.0)
    return {c: yr[c].tolist() for c in cols}


def _load_0504_dayfrac():
    """{canonical_purpose: [frac(dow) for dow 0..6]} — purpose's Mon–Sun volume
    distribution from NTS0504b, summing to 1."""
    raw = pd.read_excel(NTS0504_FILE, sheet_name="NTS0504b_day_purpose",
                        header=None, engine="odf")
    d = raw.iloc[6:].copy(); d.columns = range(d.shape[1])
    yr = d[d[0] == NTS_YEAR].copy()
    allcols = sorted({c for cols in DOW_0504.values() for c in cols})
    for c in allcols:
        yr[c] = pd.to_numeric(yr[c], errors="coerce").fillna(0.0)
    rowbyday = {row[1]: row for _, row in yr.iterrows()}
    missing = [day for day in _DAY_ORDER if day not in rowbyday]
    if missing:
        sys.exit(f"ERROR: NTS0504b '{NTS_YEAR}' missing days {missing}")
    out = {}
    for p, cols in DOW_0504.items():
        trips = [sum(float(rowbyday[day][c]) for c in cols) for day in _DAY_ORDER]
        s = sum(trips)
        out[p] = [t / s for t in trips] if s > 0 else [1.0 / 7] * 7
    return out


def main():
    rho = json.load(open(GEN_RATES))["purpose_rates"]   # canonical -> ρ_p

    # ── Aggregate hourly profile (input column) ──────────────────────────────
    rows = list(csv.DictReader(open(FRACS_FILE, newline="")))
    fieldnames = rows[0].keys() if rows else []
    mf = {(int(r["day_of_week"]), int(r["hour"].split(":")[0])): float(r["mean_fraction"])
          for r in rows}
    # weekday aggregate hourly weight P(hour) (avg Mon–Fri, normalised) → converts
    # NTS0502a's P(purpose|hour) into P(hour|purpose).
    mfwd = [sum(mf[(d, h)] for d in range(5)) / 5 for h in range(24)]
    s = sum(mfwd); mfwd = [x / s for x in mfwd]
    # weekend hourly shapes (shared across components — no per-purpose weekend data).
    def _norm_day(dow):
        v = [mf[(dow, h)] for h in range(24)]; t = sum(v)
        return [x / t for x in v] if t > 0 else [1.0 / 24] * 24
    H_weekend = {5: _norm_day(5), 6: _norm_day(6)}

    # ── Per-purpose weekday hourly distribution g_p(h) = P(hour|purpose) ──────
    col_pct = _load_0502_pct()
    g = {}
    for p in CANONICAL_PURPOSES:
        raw = [col_pct[PROXY_0502[p]][h] * mfwd[h] for h in range(24)]
        t = sum(raw)
        g[p] = [x / t for x in raw] if t > 0 else [1.0 / 24] * 24

    # ── Component weekday hourly shape Hwd_c(h): ρ-weighted blend, Σ_h = 1 ─────
    Hwd = {}
    for comp, terms in COMPONENT_PURPOSES.items():
        blend = [0.0] * 24
        for p, w in terms:
            wt = w * rho[p]
            for h in range(24):
                blend[h] += wt * g[p][h]
        t = sum(blend)
        Hwd[comp] = [x / t for x in blend] if t > 0 else [1.0 / 24] * 24

    # ── Component day-of-week volume V_c(dow): ρ-weighted, Σ_7 = 7 (⇒ W_c=1) ──
    dayfrac = _load_0504_dayfrac()
    V = {}
    for comp, terms in COMPONENT_PURPOSES.items():
        vol = [0.0] * 7
        for p, w in terms:
            wt = w * rho[p]
            for dow in range(7):
                vol[dow] += wt * dayfrac[p][dow]
        t = sum(vol)
        V[comp] = [x / t * 7 for x in vol] if t > 0 else [1.0] * 7

    # ── Assemble f_c(dow, h) = V_c(dow) · H_c(daytype, h) ────────────────────
    fc = {comp: {} for comp in COMPONENTS}
    for comp in COMPONENTS:
        for dow in range(7):
            H = Hwd[comp] if dow < 5 else H_weekend[dow]
            for h in range(24):
                fc[comp][(dow, h)] = V[comp][dow] * H[h]

    # Per-level school shapes (settled design): primary + post-primary share the escort-education
    # school-run shape; tertiary uses the commute shape (self-driven — spread AM / later PM).
    # Written as three explicit copied columns so each level's temporal is self-contained in the CSV
    # (no relying on documentation to know primary≡post-primary≡escort and tertiary≡commute).
    fc["school_primary"]     = dict(fc["school"])
    fc["school_postprimary"] = dict(fc["school"])
    fc["school_tertiary"]    = dict(fc["commute"])

    # ── Write component columns (drop legacy biz + lumped-school columns) ──────
    out_components = ["res", "commute", "retail",
                      "school_primary", "school_postprimary", "school_tertiary"]
    for r in rows:
        dow = int(r["day_of_week"]); h = int(r["hour"].split(":")[0])
        for comp in out_components:
            r[f"mean_fraction_{comp}"] = f"{fc[comp][(dow, h)]:.10f}"
    comp_cols = [f"mean_fraction_{c}" for c in out_components]
    base_cols = [c for c in fieldnames
                 if c not in comp_cols and c not in ("mean_fraction_biz", "mean_fraction_school")]
    out_cols = base_cols + comp_cols
    with open(FRACS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote component shape columns → {FRACS_FILE}")

    # ── Diagnostics + verification ───────────────────────────────────────────
    print("\nWeekday hourly shape (peak hours, % of the day at each hour):")
    print(f"  {'h':>3} {'commute':>8} {'retail':>7} {'school':>7} {'res':>7}")
    for h in range(6, 20):
        print(f"  {h:02d}: {100*Hwd['commute'][h]:7.1f}% {100*Hwd['retail'][h]:6.1f}% "
              f"{100*Hwd['school'][h]:6.1f}% {100*Hwd['res'][h]:6.1f}%")
    print("\nDay-of-week volume V_c (Mon→Sun, Σ=7):")
    for comp in ("commute", "retail", "school", "res"):
        print(f"  {comp:8s} " + " ".join(f"{v:4.2f}" for v in V[comp]))

    print("\nVerification (W_c uses the Mon–Fri average, as the model collapses dow→day_type):")
    ok = True
    for comp in COMPONENTS:
        col_sum = sum(fc[comp].values())                       # over 168 rows
        wd_avg = sum(sum(fc[comp][(d, h)] for h in range(24))
                     for d in range(5)) / 5.0                  # avg weekday daily sum
        wc = (5 * wd_avg
              + sum(fc[comp][(5, h)] for h in range(24))
              + sum(fc[comp][(6, h)] for h in range(24))) / 7.0
        flag = "" if abs(col_sum - 7) < 1e-6 and abs(wc - 1) < 1e-6 else "  <-- FAIL"
        ok = ok and not flag
        print(f"  {comp:8s} Σ_168 = {col_sum:.4f} (→7)   W_c = {wc:.4f} (→1){flag}")
    neg = [(c, k) for c in COMPONENTS for k, v in fc[c].items() if v < 0]
    print("  all non-negative" if not neg else f"  NEGATIVE at {neg[:5]}")
    print("  ✓ all checks pass" if ok and not neg else "  ✗ CHECK FAILED")


if __name__ == "__main__":
    main()
