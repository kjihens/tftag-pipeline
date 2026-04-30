"""
Pipeline-level off-target utilities.

Responsibilities
----------------
This module sits between the main TFTag pipeline and the lower-level
Cas-OFFinder wrapper functions in `offtarget.py`.

It is responsible for:
- running Cas-OFFinder in batches of unique guide sequences
- parsing and combining hit files across batches
- summarising off-target hits per spacer
- preserving chromosome-aware off-target counts
- merging off-target summaries back into the full candidate-guide table

Design notes
------------
Cas-OFFinder is run on unique spacers, not every candidate row. This avoids
re-running identical off-target searches for guides that may occur in multiple
transcript/terminus contexts.

Specificity is summarised per 20 nt spacer. The summary is then merged back
onto all candidate rows that use that spacer.
"""

from __future__ import annotations

import math
import os

import pandas as pd
from tqdm.auto import tqdm

from . import offtarget


def _empty_specificity_summary(mismatches: int, *, chromosome_aware: bool = True) -> pd.DataFrame:
    """
    Return an empty specificity summary with the expected output schema.

    Keeping the schema stable avoids downstream special cases when no hits are
    returned or when no guides are present.
    """
    cols = ["spacer", "n_hits"]

    for k in range(mismatches + 1):
        cols.append(f"n_mm{k}")

        if chromosome_aware:
            cols.append(f"n_mm{k}_same_chr")
            cols.append(f"n_mm{k}_other_chr")

    return pd.DataFrame(columns=cols)


def _unique_guide_queries(candidates: pd.DataFrame) -> pd.DataFrame:
    """
    Build the unique query table passed to Cas-OFFinder.

    Required
    --------
    spacer:
      20 nt guide spacer used as the primary merge key.

    Optional but useful
    -------------------
    pam_seq / grna_seq_23:
      Used by `offtarget.write_cas_offinder_input()` depending on its expected
      input schema.

    chromosome:
      Retained so chromosome-aware specificity summaries can classify hits as
      same-chromosome or other-chromosome relative to the intended target.
    """
    if "spacer" not in candidates.columns:
        raise KeyError("Candidates must contain 'spacer' for off-target enumeration.")

    keep_cols = [
        c
        for c in ["spacer", "pam_seq", "grna_seq_23", "chromosome"]
        if c in candidates.columns
    ]

    return (
        candidates[keep_cols]
        .drop_duplicates(subset=["spacer"])
        .reset_index(drop=True)
    )


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
    Run Cas-OFFinder in batches and return hits plus per-spacer summaries.

    Returns
    -------
    hits_df:
        One row per Cas-OFFinder hit.

    spec_df:
        One row per spacer with global mismatch counts:

            spacer, n_hits, n_mm0, n_mm1, ...

        and, when target chromosomes are available:

            n_mm0_same_chr, n_mm0_other_chr, ...

    Notes
    -----
    - `n_mm0` is expected to be at least 1 for most guides because the intended
      on-target site is usually recovered as a perfect match.
    - `n_mm0 > 1` indicates multiple perfect genomic matches and is handled
      downstream as a last-resort selection tier.
    - Batching is primarily used to provide progress reporting and avoid a
      single monolithic Cas-OFFinder invocation.
    """
    if mismatches < 0:
        raise ValueError("mismatches must be >= 0")

    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")

    if candidates.empty:
        return pd.DataFrame(), _empty_specificity_summary(mismatches)

    os.makedirs(outdir, exist_ok=True)
    prefix = run_id or "run"

    unique_guides = _unique_guide_queries(candidates)
    n_guides = len(unique_guides)
    n_batches = math.ceil(n_guides / batch_size)

    all_hits: list[pd.DataFrame] = []

    batch_iter = range(n_batches)
    if show_progress:
        batch_iter = tqdm(
            batch_iter,
            total=n_batches,
            desc="Cas-OFFinder",
            leave=False,
        )

    for batch_idx in batch_iter:
        start = batch_idx * batch_size
        end = min((batch_idx + 1) * batch_size, n_guides)
        chunk = unique_guides.iloc[start:end].copy()

        inp_path = os.path.join(
            outdir,
            f"{prefix}_cas_offinder_input_batch{batch_idx:04d}.txt",
        )
        out_path = os.path.join(
            outdir,
            f"{prefix}_cas_offinder_hits_batch{batch_idx:04d}.txt",
        )

        input_file = offtarget.write_cas_offinder_input(
            chunk,
            genome_fasta_path,
            pam_pattern=pam_pattern,
            outfile=inp_path,
            mismatches=mismatches,
        )

        output_file = offtarget.run_cas_offinder(
            input_file,
            cas_offinder_bin=cas_offinder_bin,
            device_spec=device_spec,
            output_file=out_path,
        )

        hits = offtarget.parse_cas_offinder_output(output_file)

        if not hits.empty:
            all_hits.append(hits)

    hits_df = (
        pd.concat(all_hits, ignore_index=True)
        if all_hits
        else pd.DataFrame()
    )

    spec_df = offtarget.summarize_specificity(
        hits_df,
        candidates=candidates,
        max_mismatches=mismatches,
    )

    if spec_df.empty:
        spec_df = _empty_specificity_summary(mismatches)

    if not spec_df.empty and not spec_df["spacer"].is_unique:
        raise RuntimeError(
            "Specificity summary is not unique per spacer; check parser/aggregator."
        )

    return hits_df, spec_df


def merge_specificity(
    candidates: pd.DataFrame,
    spec: pd.DataFrame | None,
    *,
    mismatches: int = 4,
) -> pd.DataFrame:
    """
    Merge per-spacer off-target summaries into the candidate guide table.

    Guarantees these global columns:
        n_hits, n_mm0, n_mm1, ..., n_mm{mismatches}

    Also guarantees chromosome-aware columns:
        n_mm0_same_chr, n_mm0_other_chr, ...

    Important
    ---------
    This function should only be used after off-target enumeration has actually
    been run. If off-target enumeration is deliberately skipped, the pipeline
    should populate these fields with NaN rather than zero, because zero would
    incorrectly imply that no off-targets were found.
    """
    if mismatches < 0:
        raise ValueError("mismatches must be >= 0")

    if "spacer" not in candidates.columns:
        raise KeyError("merge_specificity requires candidates to contain 'spacer'.")

    df = candidates.copy()

    if spec is not None and not spec.empty:
        df = df.merge(spec, on="spacer", how="left")

    expected_cols = ["n_hits"]

    if "n_hits" not in df.columns:
        df["n_hits"] = 0

    for k in range(mismatches + 1):
        global_col = f"n_mm{k}"
        same_col = f"n_mm{k}_same_chr"
        other_col = f"n_mm{k}_other_chr"

        expected_cols.extend([global_col, same_col, other_col])

        for col in [global_col, same_col, other_col]:
            if col not in df.columns:
                df[col] = 0

    df[expected_cols] = df[expected_cols].fillna(0).astype(int)

    return df