from __future__ import annotations

import os
import re

import pandas as pd


def parse_genes_arg(genes: str | None) -> list[str] | None:
    """
    Parse a gene-selection argument.

    Accepted forms
    --------------
    - None / empty string:
        return None, meaning "use all genes in the GTF database"

    - path to an existing file:
        return one non-empty gene ID per line

    - comma-separated string:
        return split gene IDs

    - single gene ID:
        return a one-element list

    This function only parses user input. It does not validate that gene IDs
    exist in the annotation database.
    """
    if genes is None:
        return None

    genes = str(genes).strip()
    if not genes:
        return None

    if os.path.exists(genes):
        with open(genes, encoding="utf-8") as fh:
            return [line.strip() for line in fh if line.strip()]

    return [gene.strip() for gene in genes.split(",") if gene.strip()]


def _infer_max_mismatch_from_columns(df: pd.DataFrame) -> int | None:
    """
    Infer the largest mismatch bin present in columns named n_mm<k>.
    """
    bins: list[int] = []

    for col in df.columns:
        match = re.fullmatch(r"n_mm(\d+)", col)
        if match:
            bins.append(int(match.group(1)))

    return max(bins) if bins else None


def filter_by_offtarget_mismatch(
    df: pd.DataFrame,
    min_mismatch: int | None,
    *,
    max_mismatches: int | None = None,
) -> pd.DataFrame:
    """
    Filter guides by Cas-OFFinder mismatch-bin summary columns.

    Definitions
    -----------
    n_mmK:
        Number of genomic matches with exactly K mismatches to the query.
        n_mm0 includes the intended on-target hit. For a unique guide, n_mm0
        should usually be exactly 1.

    Semantics
    ---------
    - min_mismatch is None or 0:
        return df unchanged; no off-target filtering

    - min_mismatch = 2:
        keep guides with n_mm0 == 1 and n_mm1 == 0

    - min_mismatch = 3:
        keep guides with n_mm0 == 1 and n_mm1 == 0 and n_mm2 == 0

    - min_mismatch = N (>0):
        keep guides with n_mm0 == 1 and no off-targets in bins below N

    - min_mismatch = -1:
        strict uniqueness within searched mismatch space:
        n_mm0 == 1 and all n_mm1..n_mm(max_mismatches) == 0

    Notes
    -----
    This function assumes Cas-OFFinder has run and n_mm* columns are present.
    """
    if df.empty:
        return df

    if min_mismatch is None:
        return df

    min_mismatch = int(min_mismatch)

    if min_mismatch == 0:
        return df

    if max_mismatches is None:
        max_mismatches = _infer_max_mismatch_from_columns(df)

    if max_mismatches is None:
        raise KeyError("No n_mm* columns found. Did you run and merge Cas-OFFinder results?")

    if "n_mm0" not in df.columns:
        raise KeyError("Missing column 'n_mm0'. Did you run and merge Cas-OFFinder results?")

    keep = pd.to_numeric(df["n_mm0"], errors="coerce").fillna(0).eq(1)

    if min_mismatch == -1:
        cols = [f"n_mm{k}" for k in range(1, max_mismatches + 1)]
    elif min_mismatch > 0:
        upper = min(min_mismatch - 1, max_mismatches)
        cols = [f"n_mm{k}" for k in range(1, upper + 1)]
    else:
        raise ValueError("min_mismatch must be -1, 0, a positive integer, or None")

    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise KeyError(f"Missing mismatch-bin columns: {missing}")

    if cols:
        keep &= df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).eq(0).all(axis=1)

    return df.loc[keep].copy()

def offtarget_keep_mask(
    df: pd.DataFrame,
    min_mismatch: int | None,
    *,
    max_mismatches: int | None = None,
) -> pd.Series:
    kept = filter_by_offtarget_mismatch(
        df,
        min_mismatch,
        max_mismatches=max_mismatches,
    )
    return df.index.isin(kept.index)