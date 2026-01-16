import numpy as np
import pandas as pd

def select_one_per_tag(df: pd.DataFrame, mode: str) -> pd.DataFrame:
    """
    Reduce multiple guides per (gene_id, tag) to one guide per (gene_id, tag),
    according to mode: 'all'|'closest'|'rs3'.
    """
    if df.empty or mode == "all":
        return df

    # Determine grouping columns robustly
    if "gene_id" in df.columns:
        group_cols = ["gene_id", "tag"]
    elif "FB-id" in df.columns:
        group_cols = ["FB-id", "tag"]
    else:
        group_cols = ["name", "tag"]

    # Ensure required columns exist
    if "cut_distance" not in df.columns:
        df = df.copy()
        df["cut_distance"] = np.nan
    if "rs3_score" not in df.columns:
        df = df.copy()
        df["rs3_score"] = np.nan

    sort_cols = []
    ascending = []

    if mode == "closest":
        sort_cols = ["cut_distance", "rs3_score"]
        ascending = [True, False]
    elif mode == "rs3":
        sort_cols = ["rs3_score", "cut_distance"]
        ascending = [False, True]
    else:
        raise ValueError(f"Unknown mode: {mode}")

    # Off-target tie-breakers if present
    if "n_hits" in df.columns:
        sort_cols.append("n_hits")
        ascending.append(True)
    elif "n_mm1" in df.columns:
        # prefer fewer close mismatches if you don't have n_hits
        sort_cols.extend(["n_mm0", "n_mm1", "n_mm2"])
        ascending.extend([True, True, True])

    df2 = df.sort_values(sort_cols, ascending=ascending, na_position="last")
    return df2.drop_duplicates(group_cols, keep="first")