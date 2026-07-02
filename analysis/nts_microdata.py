"""Generic ingestion + exploration layer for the NTS trip-level microdata (UKDS SN 5340).

This module is deliberately **domain-agnostic**: it knows how to *read* the NTS
End-User-Licence STATA files and *decode* their coded variables, and nothing about
this project's gravity components, trip-purpose groupings, vehicle-mode sets, or any
Tier-1/Tier-2 modelling choice.  Those belong in the derivation scripts that import
this module — keep them out of here.

Layout it expects (override the root with the ``NTS_DIR`` env var):

    data/NTS/
      stata/stata13/<table>_eul_2002-2024.dta      the data tables (STATA 13)
      mrdoc/excel/*response_levels*.xlsx            code -> label lookup
      mrdoc/excel/*lookup_table*eul*.xlsx           variable descriptions + year availability

The tables are large (the ``trip`` table is ~600 MB and is NOT sorted by year), so
every reader streams in chunks and never materialises more than it must.  Value
labels are **not** embedded in the .dta files (``StataReader.value_labels()`` is
empty); the ``*response_levels*`` workbook is the authoritative code book and is what
:func:`value_labels` reads.

**Special / unclassified codes — the module never decides for you.**  There are no
total/aggregate rows here (record-level data, one real trip per row; no "All"/"Total"
code level), but coded variables do carry non-substantive levels — missing/NA and the
``-8`` unclassified / ``-10`` 'DEAD' sentinels — exposed by :func:`special_codes`.
Whether to include or exclude them is question-, definition- and use-case-dependent,
so nothing here filters them implicitly: loaders keep every row, :func:`weighted_counts`
shows every code and adds a ``special`` flag column, and exclusion is an explicit,
loud opt-in via :func:`drop_special`.

Quick tour (see ``python3 analysis/nts_microdata.py --help`` for the CLI):

    import nts_microdata as nts
    nts.list_tables()                       # available tables + sizes
    nts.list_variables("trip")              # variables in a table (+ descriptions)
    nts.value_labels("TripPurpose_B01ID")   # {code: label}
    nts.years_available("TripDisIncSW")     # survey years the variable exists in
    df = nts.load("trip", columns=[...], years=[2023, 2024])   # streamed, filtered
    for chunk in nts.iter_chunks("trip", columns=[...]): ...    # fold-your-own
    nts.weighted_counts(df, by="TripPurpose_B01ID")            # Σ W5 per group
"""

from __future__ import annotations

import argparse
import functools
import glob
import os
import sys
from typing import Iterable, Iterator, Sequence

import pandas as pd


# --------------------------------------------------------------------------- paths

def _default_root() -> str:
    # repo_root/data/NTS  (this file lives in repo_root/analysis/)
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(os.path.dirname(here), "data", "NTS")


NTS_DIR = os.environ.get("NTS_DIR", _default_root())
STATA_DIR = os.path.join(NTS_DIR, "stata", "stata13")
_EXCEL_DIR = os.path.join(NTS_DIR, "mrdoc", "excel")

# The standard NTS trip weight (weighted travel sample, incl. household weight).  It
# is only a *default* here — every weighted helper takes an overridable ``weight=``
# argument; this module hard-codes no modelling weight choice.
DEFAULT_WEIGHT = "W5"


def _resolve_root() -> None:
    if not os.path.isdir(STATA_DIR):
        sys.exit(f"ERROR: NTS STATA dir not found: {STATA_DIR}\n"
                 f"       set NTS_DIR to the folder containing stata/stata13/ + mrdoc/")


def _one_glob(pattern: str, exclude: str | None = None) -> str:
    hits = sorted(glob.glob(os.path.join(_EXCEL_DIR, pattern)))
    if exclude:
        hits = [h for h in hits if exclude not in os.path.basename(h)]
    if not hits:
        sys.exit(f"ERROR: no lookup workbook matching {pattern!r} in {_EXCEL_DIR}")
    return hits[0]


# ------------------------------------------------------------------------- tables

def table_path(table: str) -> str:
    """Absolute path to a table's .dta, given a short stem (e.g. ``"trip"``)."""
    direct = os.path.join(STATA_DIR, table) if table.endswith(".dta") else None
    if direct and os.path.isfile(direct):
        return direct
    hits = sorted(glob.glob(os.path.join(STATA_DIR, f"{table}_*.dta")))
    if not hits:
        hits = sorted(glob.glob(os.path.join(STATA_DIR, f"{table}*.dta")))
    if not hits:
        sys.exit(f"ERROR: no .dta table matching {table!r} in {STATA_DIR}")
    return hits[0]


def list_tables() -> pd.DataFrame:
    """All .dta tables as a DataFrame (stem, size_mb, path), largest first."""
    _resolve_root()
    rows = []
    for p in glob.glob(os.path.join(STATA_DIR, "*.dta")):
        stem = os.path.basename(p).split("_eul_")[0].split("_2002")[0]
        rows.append({"table": stem, "size_mb": round(os.path.getsize(p) / 2**20, 1),
                     "file": os.path.basename(p)})
    return (pd.DataFrame(rows).sort_values("size_mb", ascending=False)
            .reset_index(drop=True))


@functools.lru_cache(maxsize=None)
def _reader(table: str) -> pd.io.stata.StataReader:
    return pd.io.stata.StataReader(table_path(table))


def table_columns(table: str) -> list[str]:
    """Column names present in a table's .dta (header read only, no data)."""
    return list(_reader(table).variable_labels().keys())


# ------------------------------------------------------------ variable descriptions

@functools.lru_cache(maxsize=1)
def _var_lookup() -> pd.DataFrame:
    _resolve_root()
    path = _one_glob("*lookup_table*.xlsx", exclude="response_levels")
    df = pd.read_excel(path, sheet_name="Main Table Variables")
    # Preserve the integer per-year column headers (2002..2024); strip only the
    # string labels, so years_available() can find them.
    df.columns = [c if isinstance(c, int) else str(c).strip() for c in df.columns]
    return df


def variable_info(table: str | None = None) -> pd.DataFrame:
    """Variable dictionary (Table, Variable, Data Type, Description, per-year cols).

    Filter to one table's variables with ``table=`` (case-insensitive match on the
    lookup's Table column, e.g. ``"trip"`` -> ``"Trip"``).
    """
    df = _var_lookup()
    if table is not None:
        df = df[df["Table"].astype(str).str.lower() == table.lower()]
    return df.reset_index(drop=True)


def _year_columns() -> list[int]:
    return [c for c in _var_lookup().columns if isinstance(c, int)]


def years_available(variable: str) -> list[int]:
    """Survey years for which ``variable`` is populated (from the lookup flags)."""
    df = _var_lookup()
    row = df[df["Variable"] == variable]
    if row.empty:
        return []
    yr = _year_columns()
    r = row.iloc[0]
    return [y for y in yr if pd.to_numeric(r[y], errors="coerce") == 1]


def list_variables(table: str) -> pd.DataFrame:
    """Variables actually in a table's .dta, merged with lookup descriptions.

    Descriptions are drawn from this table's rows in the lookup (a Variable name can
    recur across tables, so filtering by table avoids row-multiplying the merge); if
    the table isn't named in the lookup, falls back to a de-duplicated global merge.
    """
    cols = pd.DataFrame({"Variable": table_columns(table)})
    info = variable_info(table)[["Variable", "Data Type", "Description"]]
    if info.empty:
        info = (variable_info()[["Variable", "Data Type", "Description"]]
                .drop_duplicates("Variable"))
    else:
        info = info.drop_duplicates("Variable")
    return cols.merge(info, on="Variable", how="left")


# ---------------------------------------------------------------- code -> label maps

@functools.lru_cache(maxsize=1)
def _resp_lookup() -> pd.DataFrame:
    _resolve_root()
    path = _one_glob("*response_levels*.xlsx")
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]
    return df


@functools.lru_cache(maxsize=None)
def value_labels(variable: str) -> dict[int, str]:
    """{code: label} for a coded variable, from the response-levels code book.

    Includes NTS sentinel/negative codes (e.g. -8 'not answered', -10 'DEAD') as they
    appear in the book — the caller decides how to treat them.  Empty dict for a
    variable with no coded levels (a plain numeric like ``TripDisIncSW``).
    """
    df = _resp_lookup()
    sub = df[df["Variable"] == variable]
    out: dict[int, str] = {}
    for _, r in sub.iterrows():
        try:
            code = int(r["ID"])
        except (ValueError, TypeError):
            continue
        desc = r.get("Desc")
        out[code] = "" if pd.isna(desc) else str(desc).strip()
    return out


@functools.lru_cache(maxsize=None)
def special_codes(variable: str) -> dict[int, str]:
    """Non-substantive codes for a variable, from the code book's Part-2 block.

    ``Part == 2`` in the response-levels workbook flags every level that is NOT a real
    response: the negative sentinels (``-8`` unclassified, ``-10`` 'DEAD') and the
    positive missing/NA levels (e.g. a mode's ``'NA (public)'``).  This dataset has
    **no** total/aggregate levels — every row of the microdata is one real trip, and
    the coded variables carry no "All …"/"Total" sum level — but were a variable ever
    to carry one it would live in this same block, so this is the general guard: a
    caller summing ``ΣW5`` over a coded field should decide, explicitly, whether to
    exclude ``special_codes(col)`` first (see :func:`drop_special`).
    """
    df = _resp_lookup()
    part = pd.to_numeric(df["Part"], errors="coerce")
    sub = df[(df["Variable"] == variable) & (part == 2)]
    out: dict[int, str] = {}
    for _, r in sub.iterrows():
        try:
            code = int(r["ID"])
        except (ValueError, TypeError):
            continue
        desc = r.get("Desc")
        out[code] = "" if pd.isna(desc) else str(desc).strip()
    return out


def drop_special(df: pd.DataFrame, columns: str | Sequence[str],
                 weight: str | None = DEFAULT_WEIGHT, verbose: bool = True):
    """Drop rows carrying a non-substantive code in any of ``columns`` — loudly.

    Returns ``(filtered_df, report)`` where ``report[col] = (rows_removed, Σweight)``.
    This is deliberately explicit and off by default everywhere else: nothing in this
    module silently discards sentinel/NA (or would-be total) rows.  Inspect
    :func:`special_codes` for a column before trusting a sum over it, then call this to
    exclude them on purpose.
    """
    cols = [columns] if isinstance(columns, str) else list(columns)
    keep = pd.Series(True, index=df.index)
    report: dict[str, tuple] = {}
    for col in cols:
        sc = set(special_codes(col))
        if not sc or col not in df.columns:
            continue
        bad = df[col].isin(sc)
        if bad.any():
            w = float(df.loc[bad, weight].sum()) if (weight and weight in df) else None
            report[col] = (int(bad.sum()), w)
        keep &= ~bad
    if verbose and report:
        for col, (nrow, w) in report.items():
            wtxt = f", Σ{weight}={w:,.0f}" if w is not None else ""
            print(f"  drop_special: {col} removed {nrow:,} rows{wtxt} "
                  f"(codes {sorted(special_codes(col))})")
    return df[keep].copy(), report


def decode(df: pd.DataFrame, columns: Sequence[str] | None = None,
           suffix: str = "_label") -> pd.DataFrame:
    """Add ``<col><suffix>`` label columns for every coded column that has a code book.

    Non-destructive (keeps the numeric codes).  With ``columns=None`` it tries every
    column in ``df``; pass a list to decode only some.  Columns with no code book are
    skipped silently.
    """
    df = df.copy()
    for col in (columns if columns is not None else df.columns):
        if col not in df.columns:
            continue
        labels = value_labels(col)
        if labels:
            df[f"{col}{suffix}"] = df[col].map(labels)
    return df


# ----------------------------------------------------------------------- data loading

def iter_chunks(table: str, columns: Sequence[str] | None = None,
                years: Iterable[int] | None = None, chunksize: int = 500_000,
                year_col: str = "SurveyYear",
                convert_categoricals: bool = False) -> Iterator[pd.DataFrame]:
    """Stream a table in row chunks, optionally column-subset and year-filtered.

    Memory-safe entry point for whole-file passes over the large tables.  ``columns``
    is passed to ``read_stata`` so only those columns are materialised; if ``years``
    is given, ``year_col`` is force-included for the filter and each chunk is filtered
    before it is yielded.
    """
    _resolve_root()
    read_cols = list(columns) if columns is not None else None
    if years is not None and read_cols is not None and year_col not in read_cols:
        read_cols = read_cols + [year_col]
    yrs = set(int(y) for y in years) if years is not None else None

    it = pd.read_stata(table_path(table), columns=read_cols, chunksize=chunksize,
                       convert_categoricals=convert_categoricals)
    for chunk in it:
        if yrs is not None:
            chunk = chunk[chunk[year_col].isin(yrs)]
        if not len(chunk):
            continue
        yield chunk


def load(table: str, columns: Sequence[str] | None = None,
         years: Iterable[int] | None = None, chunksize: int = 500_000,
         year_col: str = "SurveyYear", convert_categoricals: bool = False,
         max_rows: int | None = None) -> pd.DataFrame:
    """Materialise a (column-subset, year-filtered) table into one DataFrame.

    Convenience wrapper over :func:`iter_chunks` for when the filtered result fits in
    memory (almost always true once you pick a few columns and recent years).  Guard
    against accident: loading the full ``trip`` table with no ``columns`` and no
    ``years`` refuses unless ``max_rows`` is set.
    """
    if columns is None and years is None and max_rows is None:
        sys.exit("refusing to load a full table with no columns/years filter — pass "
                 "columns=, years=, or max_rows= (streaming: use iter_chunks).")
    parts, n = [], 0
    for chunk in iter_chunks(table, columns, years, chunksize, year_col,
                             convert_categoricals):
        if max_rows is not None and n + len(chunk) > max_rows:
            parts.append(chunk.iloc[: max_rows - n])
            break
        parts.append(chunk)
        n += len(chunk)
    return (pd.concat(parts, ignore_index=True) if parts
            else pd.DataFrame(columns=list(columns) if columns else None))


# -------------------------------------------------------------- weighted aggregation

def weighted_counts(df: pd.DataFrame, by: str | Sequence[str],
                    weight: str | None = DEFAULT_WEIGHT,
                    decode_labels: bool = True) -> pd.DataFrame:
    """Σ(weight) per group of ``by`` (+ unweighted n).  ``weight=None`` -> plain n.

    Generic groupby helper — knows nothing about which variable ``by`` is.  With
    ``decode_labels`` it appends a label column for any coded ``by`` field.
    """
    by_list = [by] if isinstance(by, str) else list(by)
    g = df.groupby(by_list, dropna=False)
    out = g.size().rename("n").to_frame()
    if weight is not None and weight in df.columns:
        out[weight] = g[weight].sum()
    out = out.reset_index().sort_values(out.columns[-1], ascending=False)
    if decode_labels:
        for col in by_list:
            labels = value_labels(col)
            if labels:
                out[f"{col}_label"] = out[col].map(labels)
    # Flag (never drop) non-substantive codes so include/exclude stays the caller's
    # explicit decision — a group is 'special' if any of its by-values is a Part-2 code.
    special = pd.Series(False, index=out.index)
    for col in by_list:
        sc = set(special_codes(col))
        if sc:
            special |= out[col].isin(sc)
    if special.any():
        out["special"] = special.values
    return out.reset_index(drop=True)


def weighted_crosstab(df: pd.DataFrame, index: str, columns: str,
                      weight: str | None = DEFAULT_WEIGHT) -> pd.DataFrame:
    """Σ(weight) pivot of ``index`` × ``columns`` (unweighted counts if ``weight=None``)."""
    if weight is not None and weight in df.columns:
        return pd.pivot_table(df, index=index, columns=columns, values=weight,
                              aggfunc="sum", fill_value=0.0)
    return pd.crosstab(df[index], df[columns])


# ---------------------------------------------------------------------------- CLI

def _print_df(df: pd.DataFrame, n: int | None = None) -> None:
    with pd.option_context("display.max_rows", n or 200, "display.width", 160,
                           "display.max_colwidth", 50):
        print(df if n is None else df.head(n))


def main(argv: Sequence[str] | None = None) -> None:
    ap = argparse.ArgumentParser(
        description="Explore the NTS SN 5340 microdata (generic; no modelling logic).")
    ap.add_argument("--tables", action="store_true", help="list all tables + sizes")
    ap.add_argument("--vars", metavar="TABLE", help="list a table's variables + descriptions")
    ap.add_argument("--levels", metavar="VARIABLE", help="print code->label map for a variable")
    ap.add_argument("--years", metavar="VARIABLE", help="print survey years a variable exists in")
    ap.add_argument("--head", metavar="TABLE", help="show first rows of a table")
    ap.add_argument("--columns", help="comma-separated columns for --head")
    ap.add_argument("--filter-years", help="comma-separated survey years for --head")
    ap.add_argument("--n", type=int, default=15, help="rows for --head (default 15)")
    ap.add_argument("--counts", metavar="TABLE", help="weighted Σ(W5) counts by --by")
    ap.add_argument("--by", help="grouping variable(s) for --counts (comma-separated)")
    args = ap.parse_args(argv)

    cols = args.columns.split(",") if args.columns else None
    fyears = [int(y) for y in args.filter_years.split(",")] if args.filter_years else None

    if args.tables:
        _print_df(list_tables())
    if args.vars:
        _print_df(list_variables(args.vars))
    if args.levels:
        labels = value_labels(args.levels)
        if not labels:
            print(f"(no coded levels for {args.levels} — likely a plain numeric)")
        for code, lab in sorted(labels.items()):
            print(f"  {code:>4}: {lab}")
    if args.years:
        print(f"{args.years}: {years_available(args.years)}")
    if args.head:
        df = load(args.head, columns=cols, years=fyears, max_rows=args.n)
        _print_df(decode(df) if cols else df, args.n)
    if args.counts:
        if not args.by:
            sys.exit("--counts needs --by VARIABLE")
        by = args.by.split(",")
        need = list(dict.fromkeys(by + [DEFAULT_WEIGHT] + (["SurveyYear"] if fyears else [])))
        df = load(args.counts, columns=need, years=fyears)
        _print_df(weighted_counts(df, by if len(by) > 1 else by[0]))

    if not any([args.tables, args.vars, args.levels, args.years, args.head, args.counts]):
        ap.print_help()


if __name__ == "__main__":
    main()
