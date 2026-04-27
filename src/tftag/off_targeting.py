"""
Pipeline-level off-target utilities.

Responsibilities
----------------
- Run Cas-OFFinder once per batch of unique spacers (not per guide row).
- Parse hits and summarise mismatch-bin counts per *spacer*.
- Merge the summary into the candidate table with consistent n_mm* columns.
- Provide common filters (e.g., "no off-targets detected").
"""
from __future__ import annotations

import math
import os
import pandas as pd
from tqdm.auto import tqdm

from . import offtarget


def enumerate_offtargets_cas_offinder(
    candidates: pd.DataFrame,
    genome_fasta_path: str,
    *,
    outdir: str = "out",
    run_id: str | None = None,
    cas_offinder_bin: str = "cas-offinder",
    device_spec: str = "C",
    mismatches: int = 4,
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
    batch_size: int = 500,
    show_progress: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run Cas-OFFinder in batches and return:
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
    - Batching is used to provide a progress bar and avoid one giant monolithic run.
    """
    # Return empty but schema-consistent frames for downstream logic.
    if candidates.empty:
        cols = ["spacer", "n_hits"] + [f"n_mm{k}" for k in range(mismatches + 1)]
        return pd.DataFrame(), pd.DataFrame(columns=cols)

    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")

    os.makedirs(outdir, exist_ok=True)
    prefix = run_id or "run"

    # Cas-OFFinder only needs unique spacers as queries.
    keep_cols = [c for c in ["spacer", "pam_seq", "grna_seq_23"] if c in candidates.columns]
    if "spacer" not in keep_cols:
        raise KeyError("Candidates must contain 'spacer' for off-target enumeration.")

    unique_guides = candidates[keep_cols].drop_duplicates(subset=["spacer"]).reset_index(drop=True)

    n_guides = len(unique_guides)
    n_batches = math.ceil(n_guides / batch_size)

    all_hits = []

    batch_iter = range(n_batches)
    if show_progress:
        batch_iter = tqdm(batch_iter, total=n_batches, desc="Cas-OFFinder", leave=False)

    for batch_idx in batch_iter:
        start = batch_idx * batch_size
        end = min((batch_idx + 1) * batch_size, n_guides)
        chunk = unique_guides.iloc[start:end].copy()

        inp_path = os.path.join(outdir, f"{prefix}_cas_offinder_input_batch{batch_idx:04d}.txt")
        out_path = os.path.join(outdir, f"{prefix}_cas_offinder_hits_batch{batch_idx:04d}.txt")

        inp = offtarget.write_cas_offinder_input(
            chunk,
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
        if not hits.empty:
            all_hits.append(hits)

    if all_hits:
        hits = pd.concat(all_hits, ignore_index=True)
    else:
        hits = pd.DataFrame()

    # IMPORTANT: mismatch bins must match the mismatches value used for the Cas-OFFinder run.
    spec = offtarget.summarize_specificity(
    hits,
    candidates=candidates,
    max_mismatches=mismatches,
)
    
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

    # Ensure all expected global columns exist.
    if "n_hits" not in df.columns:
        df["n_hits"] = 0

    expected_cols = ["n_hits"]

    for k in range(mismatches + 1):
        expected_cols.append(f"n_mm{k}")
        if f"n_mm{k}" not in df.columns:
            df[f"n_mm{k}"] = 0

        for suffix in ["same_chr", "other_chr"]:
            col = f"n_mm{k}_{suffix}"
            expected_cols.append(col)
            if col not in df.columns:
                df[col] = 0

    df[expected_cols] = df[expected_cols].fillna(0).astype(int)
    return df