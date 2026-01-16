from __future__ import annotations
import os
import pandas as pd


def parse_genes_arg(genes: str | None) -> list[str] | None:
    """
    Parse `genes` argument which may be:
      - None
      - path to a file (one gene per line)
      - comma-separated list
      - single gene id
    """
    if genes is None:
        return None
    genes = str(genes).strip()
    if not genes:
        return None
    if os.path.exists(genes):
        with open(genes) as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    return [g.strip() for g in genes.split(",") if g.strip()]

def filter_by_offtarget_mismatch(df: pd.DataFrame, min_mismatch: int | None) -> pd.DataFrame:
    """
    Filter guides by off-target mismatch threshold.

    Rules:
      - Always require n_mm0 == 1 (exact on-target match present once).
      - For i in [1, min_mismatch-1], require n_mmi == 0.
      - min_mismatch <= 0 disables filtering EXCEPT:
      - min_mismatch == -1 enforces strict uniqueness:
            n_mm0 == 1 and all other n_mm* == 0

    Examples:
      min_mismatch=None -> no filtering
      min_mismatch=0    -> no filtering
      min_mismatch=1    -> require n_mm0 == 1
      min_mismatch=2    -> require n_mm0 == 1 and n_mm1 == 0
      min_mismatch=3    -> require n_mm0 == 1, n_mm1 == 0, n_mm2 == 0
      min_mismatch=-1   -> require ONLY the on-target hit (no off-targets at any mismatch)
    """

    if min_mismatch is None:
        return df

    min_mismatch = int(min_mismatch)

    # Always require exactly one perfect match
    if "n_mm0" not in df.columns:
        raise KeyError("Missing column 'n_mm0' (did you run Cas-OFFinder?)")

    keep = (df["n_mm0"] == 1)

    # Strict mode: no off-targets of any kind
    if min_mismatch == -1:
        for col in df.columns:
            if col.startswith("n_mm") and col != "n_mm0":
                keep &= (df[col] == 0)
        return df.loc[keep].copy()

    # No filtering
    if min_mismatch <= 0:
        return df

    # Disallow mismatches below threshold
    for i in range(1, min_mismatch):
        col = f"n_mm{i}"
        if col not in df.columns:
            raise KeyError(f"Missing column '{col}'")
        keep &= (df[col] == 0)

    return df.loc[keep].copy()


