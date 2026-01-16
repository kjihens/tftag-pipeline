from __future__ import annotations

import os
import pandas as pd
import numpy as np

from . import offtarget


def enumerate_offtargets_cas_offinder(
    candidates: pd.DataFrame,
    genome_fasta_path: str,
    *,
    out_dir: str = "out",
    run_id: str | None = None,
    cas_offinder_bin: str = "cas-offinder",
    device_spec: str = "C",
    mismatches: int = 4,
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run Cas-OFFinder and return:
      (hits_df, specificity_summary_df)

    - hits_df: one row per off-target hit
    - specificity_summary_df: one row per spacer, with n_hits and n_mm0..n_mmN
    """
    if candidates.empty:
        return pd.DataFrame(), pd.DataFrame(columns=["spacer", "n_hits"])

    os.makedirs(out_dir, exist_ok=True)
    prefix = run_id or "run"

    inp_path = os.path.join(out_dir, f"{prefix}_cas_offinder_input.txt")
    out_path = os.path.join(out_dir, f"{prefix}_cas_offinder_hits.txt")

    inp = offtarget.write_cas_offinder_input(
        candidates,
        genome_fasta_path,
        pam_pattern=pam_pattern,
        outfile=inp_path,
        mismatches=mismatches,
    )
    out = offtarget.run_cas_offinder(
        inp,
        cas_offinder_bin=cas_offinder_bin,
        device_spec=device_spec,
        output_file=out_path,
    )

    hits = offtarget.parse_cas_offinder_output(out)
    spec = offtarget.summarize_specificity(hits)

    if not spec.empty and not spec["spacer"].is_unique:
        raise RuntimeError("Specificity summary is not unique per spacer; check parser/aggregator.")

    return hits, spec


def merge_specificity(
    candidates: pd.DataFrame,
    spec: pd.DataFrame,
    *,
    mismatches: int = 4,
) -> pd.DataFrame:
    """
    Left-merge specificity summary into candidates and fill missing with zeros.
    Ensures columns: n_hits, n_mm0..n_mmN exist and are int.
    """
    df = candidates.copy()

    if spec is not None and not spec.empty:
        df = df.merge(spec, on="spacer", how="left")

    # ensure columns exist
    if "n_hits" not in df.columns:
        df["n_hits"] = 0
    for k in range(0, mismatches + 1):
        col = f"n_mm{k}"
        if col not in df.columns:
            df[col] = 0

    cols = ["n_hits"] + [f"n_mm{k}" for k in range(0, mismatches + 1)]
    df[cols] = df[cols].fillna(0).astype(int)
    return df

def filter_no_offtargets(df: pd.DataFrame, mismatches: int) -> pd.DataFrame:
    """
    Keep only guides with:
      - exactly one mm0 hit (assumed on-target)
      - zero hits for mm1..mm<mismatches>
    """
    if df.empty:
        return df

    mm_cols = [f"n_mm{k}" for k in range(0, mismatches + 1)]
    missing = [c for c in mm_cols if c not in df.columns]
    if missing:
        return df.iloc[0:0].copy()

    mask = (df["n_mm0"] == 1)
    for k in range(1, mismatches + 1):
        mask &= (df[f"n_mm{k}"] == 0)

    return df.loc[mask].copy()