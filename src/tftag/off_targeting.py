"""
Pipeline-level off-target utilities.

Responsibilities
---------------
- Run Cas-OFFinder once per pipeline run (not per guide row).
- Parse hits and summarise mismatch-bin counts per *spacer*.
- Merge the summary into the candidate table with consistent n_mm* columns.
- Provide common filters (e.g., "no off-targets detected").
"""
from __future__ import annotations

import os
import pandas as pd

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

    Returns
    -------
    hits_df:
      One row per off-target hit.

    spec_df:
      One row per 20-nt spacer with columns:
        spacer, n_hits, n_mm0..n_mm{mismatches}

    Notes
    -----
    - We summarise by *spacer* (20-nt) even if Cas-OFFinder queries were 23-mers.
    - n_mm0 is expected to be >= 1 for most guides (includes the on-target hit).
    """
    # Return empty but schema-consistent frames for downstream logic.
    if candidates.empty:
        cols = ["spacer", "n_hits"] + [f"n_mm{k}" for k in range(mismatches + 1)]
        return pd.DataFrame(), pd.DataFrame(columns=cols)

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

    # IMPORTANT: mismatch bins must match the mismatches value used for the Cas-OFFinder run.
    spec = offtarget.summarize_specificity(hits, max_mismatches=mismatches)

    if not spec.empty and not spec["spacer"].is_unique:
        raise RuntimeError("Specificity summary is not unique per spacer; check parser/aggregator.")

    return hits, spec


def merge_specificity(
    candidates: pd.DataFrame,
    spec: pd.DataFrame | None,
    *,
    mismatches: int = 4,
) -> pd.DataFrame:
    """
    Left-merge specificity summary into candidates and fill missing with zeros.

    Guarantees the presence (and integer dtype) of:
      n_hits, n_mm0..n_mm{mismatches}

    Assumption
    ----------
    This function should be called only when off-target enumeration was run.
    If you call it with spec=None to mean "not run", you will get zeros, which
    are not semantically correct.
    """
    df = candidates.copy()

    if spec is not None and not spec.empty:
        df = df.merge(spec, on="spacer", how="left")

    # Ensure all expected columns exist (important for downstream filters/selection).
    if "n_hits" not in df.columns:
        df["n_hits"] = 0
    for k in range(mismatches + 1):
        col = f"n_mm{k}"
        if col not in df.columns:
            df[col] = 0

    cols = ["n_hits"] + [f"n_mm{k}" for k in range(mismatches + 1)]
    df[cols] = df[cols].fillna(0).astype(int)
    return df
