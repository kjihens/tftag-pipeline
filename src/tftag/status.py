"""
Guide status helpers.

These functions keep all tested guides in the output and record why a guide
was rejected instead of silently dropping rows.
"""

from __future__ import annotations

import pandas as pd

from .utils import merge_warn


RETAINED = "retained"
REJECTED = "rejected"
NO_PAM_FOUND = "no_pam_found"


def initialise_guide_status(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add guide-status columns if missing.

    guide_status:
      retained | rejected | no_pam_found

    guide_reject_reason:
      semicolon-separated rejection reasons, or "none"
    """
    out = df.copy()

    if "guide_status" not in out.columns:
        out["guide_status"] = RETAINED

    if "guide_reject_reason" not in out.columns:
        out["guide_reject_reason"] = "none"

    if "guide_found" not in out.columns:
        out["guide_found"] = True

    return out


def mark_rejected(
    df: pd.DataFrame,
    mask,
    reason: str,
) -> pd.DataFrame:
    """
    Mark rows as rejected and append a rejection reason.

    Existing rejection reasons are preserved.
    """
    out = df.copy()
    mask = pd.Series(mask, index=out.index).fillna(False).astype(bool)

    out.loc[mask, "guide_status"] = REJECTED
    out.loc[mask, "guide_found"] = True
    out.loc[mask, "guide_reject_reason"] = out.loc[mask, "guide_reject_reason"].apply(
        lambda x: merge_warn(x, reason)
    )

    return out


def retained_mask(df: pd.DataFrame) -> pd.Series:
    """Return boolean mask for currently retained guide rows."""
    if "guide_status" not in df.columns:
        return pd.Series(True, index=df.index)

    return df["guide_status"].fillna(RETAINED).eq(RETAINED)