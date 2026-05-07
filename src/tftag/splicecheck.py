"""
Splice-site and intron-overlap annotation for TFTag guides.

Purpose
-------
Annotate whether a candidate gRNA overlaps intronic sequence or splice-sensitive
regions near exon/intron boundaries.

This is useful because a guide can be close to a start/stop codon while part of
its protospacer or PAM overlaps:

- intron
- exon/intron junction
- splice donor core
- splice acceptor core
- broader splice donor/acceptor region

This module does not reject guides. It annotates them for scoring, warnings, and
manual review.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable

import pandas as pd
from tqdm.auto import tqdm

from .utils import merge_warn


def _feature_transcript_ids(feat) -> set[str]:
    """Extract transcript-like IDs from a gffutils feature."""
    ids: set[str] = set()

    for key in ("transcript_id", "Parent", "parent"):
        for value in feat.attributes.get(key, []):
            if value:
                ids.add(str(value))

    return ids


def _interval_overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Return overlap length between two 1-based inclusive intervals."""
    start = max(int(a_start), int(b_start))
    end = min(int(a_end), int(b_end))
    return max(0, end - start + 1)


def _interval_overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Return True if two 1-based inclusive intervals overlap."""
    return _interval_overlap_len(a_start, a_end, b_start, b_end) > 0


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or directly adjacent 1-based inclusive intervals."""
    sorted_intervals = sorted((int(s), int(e)) for s, e in intervals if int(e) >= int(s))

    if not sorted_intervals:
        return []

    merged: list[tuple[int, int]] = []
    cur_s, cur_e = sorted_intervals[0]

    for s, e in sorted_intervals[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e

    merged.append((cur_s, cur_e))
    return merged


def _clip_interval(start: int, end: int, chrom_len: int | None = None) -> tuple[int, int] | None:
    """
    Clip an interval to positive genomic coordinates.

    If chrom_len is supplied, also clips to chromosome length.
    """
    s = max(1, int(start))
    e = int(end)

    if chrom_len is not None:
        e = min(e, int(chrom_len))

    if e < s:
        return None

    return s, e


def _exons_by_transcript_for_gene(gene_id: str, db) -> dict[str, list[tuple[int, int]]]:
    """
    Return exon intervals grouped by transcript ID for one gene.

    Uses transcript_id / Parent attributes rather than relying entirely on
    gffutils parent-child traversal, which is more robust for GTF files.
    """
    try:
        gene = db[gene_id]
    except Exception:
        return {}

    tx_to_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)

    for feat in db.children(gene, featuretype="exon", order_by="start"):
        tx_ids = _feature_transcript_ids(feat)

        for tx_id in tx_ids:
            tx_to_exons[tx_id].append((int(feat.start), int(feat.end)))

    return {
        tx_id: _merge_intervals(intervals)
        for tx_id, intervals in tx_to_exons.items()
        if intervals
    }


def _infer_introns_from_exons(
    exons: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """
    Infer introns from sorted exon intervals.

    For exon1=(100,200), exon2=(300,400), intron=(201,299).
    """
    exons = _merge_intervals(exons)
    introns: list[tuple[int, int]] = []

    for (_, prev_end), (next_start, _) in zip(exons, exons[1:]):
        intron_s = prev_end + 1
        intron_e = next_start - 1

        if intron_e >= intron_s:
            introns.append((intron_s, intron_e))

    return introns


def _splice_windows_for_transcript(
    *,
    exons: list[tuple[int, int]],
    strand: str,
    donor_core_intronic: int = 2,
    acceptor_core_intronic: int = 2,
    donor_region_intronic: int = 6,
    acceptor_region_intronic: int = 20,
) -> dict[str, list[tuple[int, int]]]:
    """
    Build splice-site windows for a transcript.

    Definitions
    -----------
    donor core:
        Exon boundary plus first 2 intronic bases.

    acceptor core:
        Last 2 intronic bases plus exon boundary.

    donor region:
        Broader donor-sensitive region, including first 6 intronic bases.

    acceptor region:
        Broader acceptor-sensitive region, including last 20 intronic bases.

    Coordinate notes
    ----------------
    Genomic coordinates are always returned on the '+' coordinate axis.
    Gene/transcript strand determines which side of an intron is donor vs acceptor.
    """
    exons = _merge_intervals(exons)
    introns = _infer_introns_from_exons(exons)

    windows = {
        "splice_donor_core": [],
        "splice_acceptor_core": [],
        "splice_donor_region": [],
        "splice_acceptor_region": [],
        "exon_intron_junction": [],
        "intron": introns,
    }

    for intron_s, intron_e in introns:
        if strand == "+":
            donor_boundary = intron_s - 1
            acceptor_boundary = intron_e + 1

            donor_core = (donor_boundary, intron_s + donor_core_intronic - 1)
            acceptor_core = (intron_e - acceptor_core_intronic + 1, acceptor_boundary)

            donor_region = (donor_boundary, intron_s + donor_region_intronic - 1)
            acceptor_region = (intron_e - acceptor_region_intronic + 1, acceptor_boundary)

        elif strand == "-":
            # On the minus strand, transcription runs from high to low coordinates.
            # Therefore donor is at the high-coordinate end of the intron and
            # acceptor is at the low-coordinate end.
            donor_boundary = intron_e + 1
            acceptor_boundary = intron_s - 1

            donor_core = (intron_e - donor_core_intronic + 1, donor_boundary)
            acceptor_core = (acceptor_boundary, intron_s + acceptor_core_intronic - 1)

            donor_region = (intron_e - donor_region_intronic + 1, donor_boundary)
            acceptor_region = (acceptor_boundary, intron_s + acceptor_region_intronic - 1)

        else:
            raise ValueError("strand must be '+' or '-'")

        windows["splice_donor_core"].append(donor_core)
        windows["splice_acceptor_core"].append(acceptor_core)
        windows["splice_donor_region"].append(donor_region)
        windows["splice_acceptor_region"].append(acceptor_region)

        # Compact generic exon/intron-junction windows. These are deliberately
        # small and strand-independent.
        windows["exon_intron_junction"].append((intron_s - 1, intron_s + 1))
        windows["exon_intron_junction"].append((intron_e - 1, intron_e + 1))

    return {
        key: _merge_intervals(value)
        for key, value in windows.items()
    }


def _build_gene_splice_model(
    gene_id: str,
    db,
    *,
    strand: str,
    donor_core_intronic: int,
    acceptor_core_intronic: int,
    donor_region_intronic: int,
    acceptor_region_intronic: int,
) -> dict[str, list[tuple[int, int]]]:
    """
    Build a unioned splice/intron model for all transcripts of one gene.

    The unioned model is conservative: if any transcript has an intron/splice
    boundary overlapping the guide, the guide is annotated.
    """
    tx_to_exons = _exons_by_transcript_for_gene(gene_id, db)

    union: dict[str, list[tuple[int, int]]] = {
        "intron": [],
        "splice_donor_core": [],
        "splice_acceptor_core": [],
        "splice_donor_region": [],
        "splice_acceptor_region": [],
        "exon_intron_junction": [],
    }

    for exons in tx_to_exons.values():
        tx_windows = _splice_windows_for_transcript(
            exons=exons,
            strand=strand,
            donor_core_intronic=donor_core_intronic,
            acceptor_core_intronic=acceptor_core_intronic,
            donor_region_intronic=donor_region_intronic,
            acceptor_region_intronic=acceptor_region_intronic,
        )

        for key, intervals in tx_windows.items():
            union[key].extend(intervals)

    return {
        key: _merge_intervals(intervals)
        for key, intervals in union.items()
    }


def _overlap_any(
    query_start: int,
    query_end: int,
    intervals: list[tuple[int, int]],
) -> tuple[bool, int]:
    """Return whether query overlaps any interval and the total overlap length."""
    total = 0

    for s, e in intervals:
        total += _interval_overlap_len(query_start, query_end, s, e)

    return total > 0, total


def annotate_splice_overlap(
    guides_df: pd.DataFrame,
    db,
    *,
    donor_core_intronic: int = 2,
    acceptor_core_intronic: int = 2,
    donor_region_intronic: int = 6,
    acceptor_region_intronic: int = 20,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Annotate gRNA rows for intron and splice-region overlap.

    Required columns
    ----------------
    - gene_id
    - gene_strand
    - chromosome
    - grna_23_start
    - grna_23_end

    Added columns
    -------------
    - grna_overlaps_intron
    - grna_intronic_bases
    - grna_overlaps_exon_intron_junction
    - grna_overlaps_splice_donor_core
    - grna_overlaps_splice_acceptor_core
    - grna_overlaps_splice_donor_region
    - grna_overlaps_splice_acceptor_region
    - grna_overlaps_splice_core
    - grna_overlaps_splice_region
    - splice_warning
    """
    df = guides_df.copy()

    if df.empty:
        return df

    required = [
        "gene_id",
        "gene_strand",
        "chromosome",
        "grna_23_start",
        "grna_23_end",
    ]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"annotate_splice_overlap missing required columns: {missing}")

    for col, default in [
        ("grna_overlaps_intron", False),
        ("grna_intronic_bases", 0),
        ("grna_overlaps_exon_intron_junction", False),
        ("grna_overlaps_splice_donor_core", False),
        ("grna_overlaps_splice_acceptor_core", False),
        ("grna_overlaps_splice_donor_region", False),
        ("grna_overlaps_splice_acceptor_region", False),
        ("grna_overlaps_splice_core", False),
        ("grna_overlaps_splice_region", False),
        ("splice_warning", "none"),
    ]:
        if col not in df.columns:
            df[col] = default

    if "warnings" not in df.columns:
        df["warnings"] = "none"

    model_cache: dict[tuple[str, str], dict[str, list[tuple[int, int]]]] = {}

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Annotating splice overlap", leave=False)

    for row in iterator:
        idx = row.Index

        gene_id = str(row.gene_id)
        strand = str(row.gene_strand)
        key = (gene_id, strand)

        if key not in model_cache:
            model_cache[key] = _build_gene_splice_model(
                gene_id,
                db,
                strand=strand,
                donor_core_intronic=donor_core_intronic,
                acceptor_core_intronic=acceptor_core_intronic,
                donor_region_intronic=donor_region_intronic,
                acceptor_region_intronic=acceptor_region_intronic,
            )

        model = model_cache[key]

        q_start = int(min(row.grna_23_start, row.grna_23_end))
        q_end = int(max(row.grna_23_start, row.grna_23_end))

        overlaps_intron, intronic_bases = _overlap_any(q_start, q_end, model["intron"])
        overlaps_junction, _ = _overlap_any(q_start, q_end, model["exon_intron_junction"])
        overlaps_donor_core, _ = _overlap_any(q_start, q_end, model["splice_donor_core"])
        overlaps_acceptor_core, _ = _overlap_any(q_start, q_end, model["splice_acceptor_core"])
        overlaps_donor_region, _ = _overlap_any(q_start, q_end, model["splice_donor_region"])
        overlaps_acceptor_region, _ = _overlap_any(q_start, q_end, model["splice_acceptor_region"])

        overlaps_core = overlaps_donor_core or overlaps_acceptor_core
        overlaps_region = overlaps_donor_region or overlaps_acceptor_region

        df.at[idx, "grna_overlaps_intron"] = overlaps_intron
        df.at[idx, "grna_intronic_bases"] = int(intronic_bases)
        df.at[idx, "grna_overlaps_exon_intron_junction"] = overlaps_junction
        df.at[idx, "grna_overlaps_splice_donor_core"] = overlaps_donor_core
        df.at[idx, "grna_overlaps_splice_acceptor_core"] = overlaps_acceptor_core
        df.at[idx, "grna_overlaps_splice_donor_region"] = overlaps_donor_region
        df.at[idx, "grna_overlaps_splice_acceptor_region"] = overlaps_acceptor_region
        df.at[idx, "grna_overlaps_splice_core"] = overlaps_core
        df.at[idx, "grna_overlaps_splice_region"] = overlaps_region

        warning_items: list[str] = []

        if overlaps_core:
            warning_items.append("gRNA overlaps splice donor/acceptor core")
        elif overlaps_region:
            warning_items.append("gRNA overlaps broader splice donor/acceptor region")

        if overlaps_intron:
            warning_items.append(f"gRNA overlaps intronic sequence ({intronic_bases} bp)")

        if overlaps_junction:
            warning_items.append("gRNA overlaps exon-intron junction")

        if warning_items:
            splice_warning = "; ".join(warning_items)
            df.at[idx, "splice_warning"] = splice_warning
            df.at[idx, "warnings"] = merge_warn(df.at[idx, "warnings"], splice_warning)

    return df