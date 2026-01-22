"""
Selection utilities.

Purpose
-------
Reduce multiple candidate guides per (gene_id, tag) to a single guide, using
a user-selected ranking strategy:

  - mode="all"     : return all rows unchanged
  - mode="closest" : prefer smallest cut_distance; tie-break by higher rs3_score;
                     then fewer off-target hits (if present)
  - mode="rs3"     : prefer highest rs3_score; tie-break by smaller cut_distance;
                     then fewer off-target hits (if present)

Notes
-----
- This function is intentionally robust to pipelines that did not compute
  off-targets (NaNs in n_hits / n_mm*).
- We enforce deterministic selection by adding a final stable tie-breaker column.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name that exists in df.columns, else None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def select_one_per_tag(df: pd.DataFrame, mode: str = "all") -> pd.DataFrame:
    """
    Reduce multiple guides per (gene_id, tag) to one guide per (gene_id, tag).

    Parameters
    ----------
    df:
      Candidate guide table.
    mode:
      "all" | "closest" | "rs3"

    Returns
    -------
    DataFrame
      A filtered view/copy containing at most one row per (gene_id, tag).
    """
    if df.empty or mode == "all":
        return df

    if mode not in ("closest", "rs3"):
        raise ValueError("mode must be one of: all, closest, rs3")

    # --- Grouping columns (enforce your requirement: 1 of N1, 1 of N2, etc.) ---
    if "tag" not in df.columns:
        raise KeyError("select_one_per_tag requires a 'tag' column")

    if "gene_id" in df.columns:
        group_cols = ["gene_id", "tag"]
    else:
        # Conservative fallbacks if users pass legacy schema
        gid = _first_present(df, ["FB-id", "FB_id", "gene", "name", "gene_symbol"])
        if gid is None:
            raise KeyError("select_one_per_tag could not find a gene identifier column (expected 'gene_id')")
        group_cols = [gid, "tag"]

    # Work on a copy to avoid mutating caller's df
    df2 = df.copy()

    # --- Ensure ranking columns exist ---
    if "cut_distance" not in df2.columns:
        df2["cut_distance"] = np.nan
    if "rs3_score" not in df2.columns:
        df2["rs3_score"] = np.nan

    sort_cols: list[str] = []
    ascending: list[bool] = []

    if mode == "closest":
        sort_cols += ["cut_distance", "rs3_score"]
        ascending += [True, False]
    else:  # mode == "rs3"
        sort_cols += ["rs3_score", "cut_distance"]
        ascending += [False, True]

    # --- Off-target tie-breakers (prefer fewer hits) ---
    # Prefer n_hits if it exists and has at least one non-NA value.
    if "n_hits" in df2.columns and df2["n_hits"].notna().any():
        sort_cols.append("n_hits")
        ascending.append(True)
    else:
        # Fall back to whatever n_mm* columns exist (in increasing mismatch order)
        mm_cols = [c for c in df2.columns if c.startswith("n_mm")]
        # Sort numerically by the suffix, so n_mm0, n_mm1, ...
        def _mm_key(c: str) -> int:
            try:
                return int(c.replace("n_mm", ""))
            except Exception:
                return 10**9

        mm_cols = sorted(mm_cols, key=_mm_key)

        # If present, prefer fewer off-targets at small mismatch counts.
        # (We generally do NOT want to include n_mm0 here as a "fewer is better"
        # criterion, because n_mm0 should be exactly 1 for designable guides, and
        # lower could reflect missing enumeration.)
        for c in mm_cols:
            if c == "n_mm0":
                continue
            sort_cols.append(c)
            ascending.append(True)

    # --- Deterministic final tie-breaker ---
    # Use a stable identifier if available; otherwise fall back to genomic coords.
    stable_id = _first_present(df2, ["spacer", "grna_seq_23"])
    if stable_id is not None:
        sort_cols.append(stable_id)
        ascending.append(True)
    else:
        # Genomic fallback (only add those that exist)
        for c in ["chromosome", "protospacer_start", "pam_start", "cut_pos"]:
            if c in df2.columns:
                sort_cols.append(c)
                ascending.append(True)

    # Sort and pick first within each group
    df2 = df2.sort_values(sort_cols, ascending=ascending, na_position="last", kind="mergesort")
    # mergesort is stable: preserves order for ties, which + our stable_id improves determinism
    return df2.drop_duplicates(group_cols, keep="first").reset_index(drop=True)