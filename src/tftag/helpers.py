from __future__ import annotations

import os
import re
import pandas as pd


def parse_genes_arg(genes: str | None) -> list[str] | None:
    """
    Parse a `genes` argument, accepting:
      - None / empty -> None (meaning: all genes in the GTF DB)
      - path to a file -> one gene ID per non-empty line
      - comma-separated list -> split on commas
      - single ID -> returned as a 1-element list

    Notes
    -----
    This function does not validate gene IDs against the GTF DB; that happens downstream.
    """
    if genes is None:
        return None

    genes = str(genes).strip()
    if not genes:
        return None

    if os.path.exists(genes):
        with open(genes) as fh:
            return [ln.strip() for ln in fh if ln.strip()]

    # comma-separated list (also works for a single gene ID)
    return [g.strip() for g in genes.split(",") if g.strip()]


def _infer_max_mismatch_from_columns(df: pd.DataFrame) -> int | None:
    """
    Infer the largest mismatch bin present in columns named n_mm<k>.
    Returns None if no such columns exist.
    """
    mm = []
    for c in df.columns:
        m = re.fullmatch(r"n_mm(\d+)", c)
        if m:
            mm.append(int(m.group(1)))
    return max(mm) if mm else None


def filter_by_offtarget_mismatch(
    df: pd.DataFrame,
    min_mismatch: int | None,
    *,
    max_mismatches: int | None = None,
) -> pd.DataFrame:
    """
    Filter guides by an off-target mismatch threshold using Cas-OFFinder summary bins.

    Definitions
    ----------
    n_mmK = number of genomic matches with exactly K mismatches to the query sequence.
            (n_mm0 includes the intended on-target hit; for most guides n_mm0 >= 1)

    Semantics
    ---------
    - min_mismatch=None:
        return df unchanged (no filtering).

    - min_mismatch = 0:
        keep guides with exactly one perfect match:
          n_mm0 == 1

    - min_mismatch = 2:
        keep guides with exactly one perfect match AND no 1-mismatch off-targets:
          n_mm0 == 1 and n_mm1 == 0

    - min_mismatch = N (>0):
        keep guides with:
          n_mm0 == 1 and n_mm1..n_mm(N-1) == 0

    - min_mismatch = -1 (strict uniqueness):
        keep guides with:
          n_mm0 == 1 and n_mm1..n_mm(max_mismatches) == 0

    Notes
    -----
    - This function assumes off-target enumeration ran and mismatch bins are present.
    - If max_mismatches is not provided, we infer it from existing n_mm* columns.
    """
    if df.empty:
        return df

    if min_mismatch is None:
        return df

    min_mismatch = int(min_mismatch)

    # Determine how many bins we can/should enforce.
    if max_mismatches is None:
        max_mismatches = _infer_max_mismatch_from_columns(df)
    if max_mismatches is None:
        raise KeyError("No n_mm* columns found (did you run Cas-OFFinder and merge specificity?)")

    # Always require exactly one perfect match.
    if "n_mm0" not in df.columns:
        raise KeyError("Missing column 'n_mm0' (did you run Cas-OFFinder and merge specificity?)")
    keep = (df["n_mm0"] == 1)

    # Strict mode: disallow *all* off-target bins we have.
    if min_mismatch == -1:
        cols = [f"n_mm{k}" for k in range(1, max_mismatches + 1)]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise KeyError(f"Missing mismatch-bin columns: {missing}")
        keep &= (df[cols] == 0).all(axis=1)
        return df.loc[keep].copy()

    # min_mismatch < 0 but not -1: treat as no filter beyond n_mm0==1
    if min_mismatch < 0:
        return df.loc[keep].copy()

    # Enforce that bins strictly below the threshold are zero.
    # Example: min_mismatch=2 enforces n_mm1==0.
    upper = min(min_mismatch - 1, max_mismatches)
    if upper >= 1:
        cols = [f"n_mm{k}" for k in range(1, upper + 1)]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise KeyError(f"Missing mismatch-bin columns: {missing}")
        keep &= (df[cols] == 0).all(axis=1)

    return df.loc[keep].copy()