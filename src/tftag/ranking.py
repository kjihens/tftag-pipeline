
"""
Optional guide selection policies.
"""
from __future__ import annotations
import pandas as pd

def retain_one_per_tag(df: pd.DataFrame,
                       prefer_max: bool = True) -> pd.DataFrame:
    x = df.copy()
    key = ["name", "feature", "tag"]
    x["rk_cut"] = x.groupby(key)["cut_distance"].rank(method="min", ascending=True)
    x["rk_rs3"] = x.groupby(key)["rs3_score"].rank(method="min", ascending=False)
    ccol = "cclmoff_max" if prefer_max and "cclmoff_max" in x.columns else ("cclmoff_sum" if "cclmoff_sum" in x.columns else None)
    x["rk_ccl"] = x.groupby(key)[ccol].rank(method="min", ascending=True) if ccol else 0
    x["score"] = x["rk_cut"] - x["rk_rs3"] + x["rk_ccl"]
    idx = x.groupby(key)["score"].idxmin()
    return x.loc[idx].sort_values(key).reset_index(drop=True)
