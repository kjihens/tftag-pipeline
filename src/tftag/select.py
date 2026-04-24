"""
Guide selection utilities.

Purpose
-------
Reduce multiple candidate guides per (gene_id, tag) to a single guide using
transparent, deterministic ranking.

Selection modes
---------------
- mode="all"
    Return all rows unchanged.

- mode="closest"
    Prefer smallest cut_distance, then higher rs3_score, then fewer off-targets.

- mode="rs3"
    Prefer highest rs3_score, then smaller cut_distance, then fewer off-targets.

- mode="score"
    Use a composite guide-selection score with explicit last-resort tiers.

Composite score logic
---------------------
The score mode combines:
- lower cut distance
- higher RS3 efficiency score
- fewer off-targets, weighted by mismatch class
- preference for guides that do not require homology-arm editing

Important last-resort rules
---------------------------
Some guides should not be selected unless no better option exists.

We therefore use selection_tier before numerical score:

  tier 0: normal usable guide
  tier 1: n_mm0 > 1
          multiple perfect genomic matches; avoid unless no tier 0 guide exists
  tier 2: rejected=True
          usually means guide would require editing a non-coding homology arm;
          use only as a true last resort

Within each tier, higher selection_score is better.

This means a guide with very high RS3 cannot accidentally outrank a safer guide
from a better tier.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate column present in df.columns."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _mm_cols(df: pd.DataFrame) -> list[str]:
    """Return n_mm* columns sorted by numeric mismatch suffix."""
    cols = [c for c in df.columns if c.startswith("n_mm")]

    def key(c: str) -> int:
        try:
            return int(c.replace("n_mm", ""))
        except Exception:
            return 10**9

    return sorted(cols, key=key)


def _group_columns(df: pd.DataFrame) -> list[str]:
    """
    Return columns used to define one biological target.

    Normally this is:
      gene_id + tag

    where tag is N1, N2, C1, C2, etc.
    """
    if "tag" not in df.columns:
        raise KeyError("select_one_per_tag requires a 'tag' column")

    if "gene_id" in df.columns:
        return ["gene_id", "tag"]

    fallback = _first_present(df, ["FB-id", "FB_id", "gene", "name", "gene_symbol"])
    if fallback is None:
        raise KeyError("select_one_per_tag could not find a gene identifier column")
    return [fallback, "tag"]


def add_guide_selection_score(
    df: pd.DataFrame,
    *,
    max_cut_distance: int = 30,
    w_cut: float = 40.0,
    w_rs3: float = 25.0,
    w_offtarget: float = 25.0,
    w_no_edit: float = 10.0,
    coding_edit_penalty: float = 5.0,
) -> pd.DataFrame:
    """
    Add guide-selection fields.

    Adds
    ----
    selection_tier:
      0 = normal guide
      1 = multiple perfect genomic matches, n_mm0 > 1
      2 = rejected guide, usually non-coding HA edit required

    selection_score:
      Composite score. Higher is better within a tier.

    selection_warning:
      Human-readable warning for guides that are not ideal but may be selected
      if no better guide exists for that terminus.

    Notes
    -----
    This function does not remove any rows. It annotates them for downstream
    selection.
    """
    out = df.copy()

    # Defensive defaults for optional columns.
    if "cut_distance" not in out.columns:
        out["cut_distance"] = np.nan
    if "rs3_score" not in out.columns:
        out["rs3_score"] = np.nan

    if "n_mm0" not in out.columns:
        out["n_mm0"] = 1
    for k in range(1, 5):
        if f"n_mm{k}" not in out.columns:
            out[f"n_mm{k}"] = 0

    if "requires_edit_arm" not in out.columns:
        out["requires_edit_arm"] = "none"
    if "rejected" not in out.columns:
        out["rejected"] = False
    if "reject_reason" not in out.columns:
        out["reject_reason"] = "none"

    # ------------------------------------------------------------
    # 1. Last-resort tiering
    # ------------------------------------------------------------
    out["selection_tier"] = 0

    n_mm0 = pd.to_numeric(out["n_mm0"], errors="coerce").fillna(1)
    multi_perfect = n_mm0 > 1

    rejected = out["rejected"].fillna(False).astype(bool)

    out.loc[multi_perfect, "selection_tier"] = np.maximum(
        out.loc[multi_perfect, "selection_tier"],
        1,
    )

    out.loc[rejected, "selection_tier"] = np.maximum(
        out.loc[rejected, "selection_tier"],
        2,
    )

    # ------------------------------------------------------------
    # 2. Cut-distance component
    # ------------------------------------------------------------
    cut = pd.to_numeric(out["cut_distance"], errors="coerce").fillna(max_cut_distance)
    cut_norm = 1.0 - (cut.clip(0, max_cut_distance) / max_cut_distance)

    # ------------------------------------------------------------
    # 3. RS3 efficiency component
    # ------------------------------------------------------------
    rs3 = pd.to_numeric(out["rs3_score"], errors="coerce")

    if rs3.notna().any():
        if rs3.min(skipna=True) >= 0 and rs3.max(skipna=True) <= 1:
            rs3_norm = rs3.fillna(rs3.median())
        else:
            # Robust fallback if RS3 scores are not naturally 0..1.
            rs3_norm = rs3.rank(pct=True).fillna(0.5)
    else:
        rs3_norm = pd.Series(0.5, index=out.index)

    # ------------------------------------------------------------
    # 4. Off-target component
    # ------------------------------------------------------------
    # More severe mismatch classes are penalised more strongly.
    # n_mm0 is handled as a tier rule, not as a normal score component.
    n_mm1 = pd.to_numeric(out["n_mm1"], errors="coerce").fillna(0)
    n_mm2 = pd.to_numeric(out["n_mm2"], errors="coerce").fillna(0)
    n_mm3 = pd.to_numeric(out["n_mm3"], errors="coerce").fillna(0)
    n_mm4 = pd.to_numeric(out["n_mm4"], errors="coerce").fillna(0)

    off_penalty_raw = (
        16.0 * n_mm1
        + 8.0 * n_mm2
        + 4.0 * n_mm3
        + 2.0 * n_mm4
    )

    off_score = 1.0 / (1.0 + off_penalty_raw)

    # ------------------------------------------------------------
    # 5. Homology-arm edit component
    # ------------------------------------------------------------
    requires_edit = out["requires_edit_arm"].fillna("none").astype(str)

    no_edit_score = (requires_edit == "none").astype(float)

    # Mild penalty for coding-arm edits. Non-coding-arm edit cases should already
    # have rejected=True and therefore tier 2.
    coding_edit_pen = ((requires_edit != "none") & (~rejected)).astype(float)

    # ------------------------------------------------------------
    # 6. Composite score
    # ------------------------------------------------------------
    out["selection_score"] = (
        w_cut * cut_norm
        + w_rs3 * rs3_norm
        + w_offtarget * off_score
        + w_no_edit * no_edit_score
        - coding_edit_penalty * coding_edit_pen
    )

    # ------------------------------------------------------------
    # 7. Human-readable warnings
    # ------------------------------------------------------------
    warnings: list[str] = []
    for _, row in out.iterrows():
        w: list[str] = []

        try:
            if float(row.get("n_mm0", 1)) > 1:
                w.append("multiple perfect genomic matches; use only if no unique guide is available")
        except Exception:
            pass

        if bool(row.get("rejected", False)):
            reason = row.get("reject_reason", "requires non-coding HA edit")
            if reason in (None, "", "none"):
                reason = "requires non-coding HA edit"
            w.append(f"last-resort guide: {reason}")

        warnings.append("; ".join(w) if w else "none")

    out["selection_warning"] = warnings

    return out


def _add_deterministic_tiebreakers(
    df: pd.DataFrame,
    sort_cols: list[str],
    ascending: list[bool],
) -> tuple[list[str], list[bool]]:
    """
    Add stable final tie-breakers so repeated runs select the same guide.
    """
    stable_id = _first_present(df, ["grna_seq_23", "spacer"])
    if stable_id is not None:
        sort_cols.append(stable_id)
        ascending.append(True)
        return sort_cols, ascending

    for c in ["chromosome", "grna_23_start", "protospacer_start", "pam_start", "cut_pos"]:
        if c in df.columns:
            sort_cols.append(c)
            ascending.append(True)

    return sort_cols, ascending


def select_one_per_tag(df: pd.DataFrame, mode: str = "all") -> pd.DataFrame:
    """
    Reduce multiple guides per (gene_id, tag) to one guide per terminus.

    Parameters
    ----------
    df:
      Candidate guide table.

    mode:
      "all" | "closest" | "rs3" | "score"

    Returns
    -------
    pd.DataFrame
      A filtered copy containing at most one row per (gene_id, tag), unless
      mode="all".
    """
    if df.empty or mode == "all":
        return df

    if mode not in ("closest", "rs3", "score"):
        raise ValueError("mode must be one of: all, closest, rs3, score")

    group_cols = _group_columns(df)
    df2 = df.copy()

    # Ensure common ranking columns exist.
    if "cut_distance" not in df2.columns:
        df2["cut_distance"] = np.nan
    if "rs3_score" not in df2.columns:
        df2["rs3_score"] = np.nan

    # ------------------------------------------------------------
    # Composite scoring mode
    # ------------------------------------------------------------
    if mode == "score":
        df2 = add_guide_selection_score(df2)

        sort_cols = [
            "selection_tier",
            "selection_score",
            "cut_distance",
            "rs3_score",
        ]
        ascending = [
            True,   # lower tier is better
            False,  # higher score is better
            True,   # closer cut is better
            False,  # higher RS3 is better
        ]

        sort_cols, ascending = _add_deterministic_tiebreakers(df2, sort_cols, ascending)

        df2 = df2.sort_values(
            sort_cols,
            ascending=ascending,
            na_position="last",
            kind="mergesort",
        )

        return df2.drop_duplicates(group_cols, keep="first").reset_index(drop=True)

    # ------------------------------------------------------------
    # Legacy/simple modes
    # ------------------------------------------------------------
    sort_cols: list[str] = []
    ascending: list[bool] = []

    if mode == "closest":
        sort_cols += ["cut_distance", "rs3_score"]
        ascending += [True, False]

    elif mode == "rs3":
        sort_cols += ["rs3_score", "cut_distance"]
        ascending += [False, True]

    # Prefer fewer off-targets as tie-breakers if present.
    if "n_hits" in df2.columns and df2["n_hits"].notna().any():
        sort_cols.append("n_hits")
        ascending.append(True)
    else:
        for c in _mm_cols(df2):
            if c == "n_mm0":
                continue
            sort_cols.append(c)
            ascending.append(True)

    sort_cols, ascending = _add_deterministic_tiebreakers(df2, sort_cols, ascending)

    df2 = df2.sort_values(
        sort_cols,
        ascending=ascending,
        na_position="last",
        kind="mergesort",
    )

    return df2.drop_duplicates(group_cols, keep="first").reset_index(drop=True)