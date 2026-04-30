"""
Homology-arm utilities for TFTag.

Responsibilities
----------------
- Check whether a guide locus can support complete homology arms and validation
  primer windows.
- Build 240 bp HAL/HAR homology arms around start/stop codons.
- Return arm sequences in gene orientation.

Coordinate conventions
----------------------
- Genomic coordinates are 1-based inclusive.
- HAL/HAR are defined in gene orientation:
    HAL = upstream/left arm in gene orientation
    HAR = downstream/right arm in gene orientation
"""

from __future__ import annotations

from typing import Tuple

import pandas as pd
from tqdm.auto import tqdm

from .utils import clamp_range, get_sequence
from .utils import merge_warn


def arm_ranges(
    feature: str,
    strand: str,
    codon_start: int,
    codon_end: int,
    *,
    arm_len: int = 240,
) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """
    Return genomic ranges for HAL and HAR.

    Returns
    -------
    (HAL, HAR)
        Each is a 1-based inclusive genomic interval `(start, end)`.
    """
    if feature not in ("start_codon", "stop_codon"):
        raise ValueError("feature must be 'start_codon' or 'stop_codon'")

    if strand == "+":
        hal = (codon_start - arm_len, codon_start - 1)
        har = (codon_end + 1, codon_end + arm_len)
    elif strand == "-":
        hal = (codon_end + 1, codon_end + arm_len)
        har = (codon_start - arm_len, codon_start - 1)
    else:
        raise ValueError("strand must be '+' or '-'")

    return hal, har


def prefilter_designable(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    *,
    hal_len: int = 240,
    har_len: int = 240,
    primer_upstream_window: tuple[int, int] = (100, 300),
    primer_downstream_window: tuple[int, int] = (100, 300),
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Mark guide rows as designable if arms and primer windows fit within contig.

    Adds/updates:
    - designable
    - skip_reason
    """
    df = gRNA_df.copy()

    if "skip_reason" not in df.columns:
        df["skip_reason"] = "none"

    df["designable"] = True

    required = ["chromosome", "feature", "gene_strand", "codon_start", "codon_end"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"prefilter_designable missing required columns: {missing}")

    def in_bounds(start: int, end: int, chrom_len: int) -> bool:
        return start >= 1 and end <= chrom_len and end >= start

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Prefiltering designable loci", leave=False)

    for row in iterator:
        idx = row.Index

        chrom = row.chromosome
        feature = row.feature
        strand = row.gene_strand
        codon_start = int(row.codon_start)
        codon_end = int(row.codon_end)
        chrom_len = len(fasta_dict[chrom].seq)

        (hal_s, hal_e), (har_s, har_e) = arm_ranges(
            feature,
            strand,
            codon_start,
            codon_end,
        )

        up_lo, up_hi = primer_upstream_window
        dn_lo, dn_hi = primer_downstream_window

        if strand == "+":
            up_start, up_end = hal_s - up_hi, hal_s - up_lo
            dn_start, dn_end = har_e + dn_lo, har_e + dn_hi
        else:
            up_start, up_end = hal_e + up_lo, hal_e + up_hi
            dn_start, dn_end = har_s - dn_hi, har_s - dn_lo

        reasons: list[str] = []

        if not in_bounds(hal_s, hal_e, chrom_len):
            reasons.append("HAL out of contig")
        if not in_bounds(har_s, har_e, chrom_len):
            reasons.append("HAR out of contig")
        if (hal_e - hal_s + 1) < hal_len:
            reasons.append(f"HAL < {hal_len} bp")
        if (har_e - har_s + 1) < har_len:
            reasons.append(f"HAR < {har_len} bp")
        if not in_bounds(up_start, up_end, chrom_len):
            reasons.append("Upstream primer window out of contig")
        if not in_bounds(dn_start, dn_end, chrom_len):
            reasons.append("Downstream primer window out of contig")

        if reasons:
            df.at[idx, "designable"] = False
            df.at[idx, "skip_reason"] = merge_warn(df.at[idx, "skip_reason"], reasons)

    return df


def add_homology_arms(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    *,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Add HAL/HAR genomic coordinates and gene-oriented arm sequences.

    Adds:
    - HALs, HALe, HARs, HARe
    - HAL_seq_gene, HAR_seq_gene
    """
    df = gRNA_df.copy()

    if "warnings" not in df.columns:
        df["warnings"] = "none"

    if "designable" not in df.columns:
        df["designable"] = False

    required = ["chromosome", "feature", "gene_strand", "codon_start", "codon_end"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"add_homology_arms missing required columns: {missing}")

    designable_mask = df["designable"].fillna(False).astype(bool)
    if not designable_mask.any():
        return df

    work = df.loc[designable_mask].copy()

    iterator = work.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(work), desc="Computing homology arms", leave=False)

    contig_len_cache: dict[str, int] = {}

    for row in iterator:
        idx = row.Index

        chrom = row.chromosome
        feature = row.feature
        strand = row.gene_strand
        codon_start = int(row.codon_start)
        codon_end = int(row.codon_end)

        (hal_s, hal_e), (har_s, har_e) = arm_ranges(
            feature,
            strand,
            codon_start,
            codon_end,
        )

        chrom_len = contig_len_cache.setdefault(chrom, len(fasta_dict[chrom].seq))

        # Defensive clamp; designability prefilter should already ensure this is unnecessary.
        hal_s, hal_e = clamp_range(hal_s, hal_e, chrom_len)
        har_s, har_e = clamp_range(har_s, har_e, chrom_len)

        work.at[idx, "HALs"] = hal_s
        work.at[idx, "HALe"] = hal_e
        work.at[idx, "HARs"] = har_s
        work.at[idx, "HARe"] = har_e

        hal_seq = get_sequence(fasta_dict, chrom, hal_s, hal_e, strand)
        har_seq = get_sequence(fasta_dict, chrom, har_s, har_e, strand)

        work.at[idx, "HAL_seq_gene"] = hal_seq
        work.at[idx, "HAR_seq_gene"] = har_seq

        warnings: list[str] = []

        if len(hal_seq) < 240:
            warnings.append(f"HAL truncated to {len(hal_seq)} bp")
        if len(har_seq) < 240:
            warnings.append(f"HAR truncated to {len(har_seq)} bp")

        if warnings:
            work.at[idx, "warnings"] = merge_warn(work.at[idx, "warnings"], warnings)

    out = df.copy()
    for col in work.columns:
        out.loc[work.index, col] = work[col]

    return out