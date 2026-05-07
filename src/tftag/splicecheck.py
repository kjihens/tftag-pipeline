"""
Splice-site and edit-risk annotation for TFTag.

Purpose
-------
This module annotates whether actual donor blocking edits overlap intronic
sequence or splice-sensitive regions.

Important design principle
--------------------------
We do NOT penalise a guide merely because the gRNA 23-mer overlaps an intron
from another splice variant. The relevant question is:

    Does the actual donor edit mutate intronic sequence or a splice site
    in the transcript(s) covered by this terminus?

Therefore the main public function is:

    annotate_edit_splice_risk()

It should be called after apply_silent_edits().
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


def _merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent 1-based inclusive intervals."""
    sorted_intervals = sorted(
        (int(start), int(end))
        for start, end in intervals
        if int(end) >= int(start)
    )

    if not sorted_intervals:
        return []

    merged: list[tuple[int, int]] = []
    cur_start, cur_end = sorted_intervals[0]

    for start, end in sorted_intervals[1:]:
        if start <= cur_end + 1:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end

    merged.append((cur_start, cur_end))
    return merged


def _infer_introns_from_exons(
    exons: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """
    Infer intron intervals from transcript exon intervals.

    Example
    -------
    exon1 = 100..200
    exon2 = 300..400

    inferred intron = 201..299
    """
    exons = _merge_intervals(exons)
    introns: list[tuple[int, int]] = []

    for (_, previous_exon_end), (next_exon_start, _) in zip(exons, exons[1:]):
        intron_start = previous_exon_end + 1
        intron_end = next_exon_start - 1

        if intron_end >= intron_start:
            introns.append((intron_start, intron_end))

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
    Build intron and splice-window intervals for one transcript.

    Definitions
    -----------
    splice_donor_core:
        Exon boundary plus first 2 intronic bases.

    splice_acceptor_core:
        Last 2 intronic bases plus exon boundary.

    splice_donor_region:
        Broader donor-sensitive region, including first 6 intronic bases.

    splice_acceptor_region:
        Broader acceptor-sensitive region, including last 20 intronic bases.

    Notes
    -----
    Coordinates are always returned on the genomic '+' coordinate axis.
    The transcript strand determines which side of each intron is donor vs
    acceptor.
    """
    exons = _merge_intervals(exons)
    introns = _infer_introns_from_exons(exons)

    windows: dict[str, list[tuple[int, int]]] = {
        "intron": introns,
        "splice_donor_core": [],
        "splice_acceptor_core": [],
        "splice_donor_region": [],
        "splice_acceptor_region": [],
        "exon_intron_junction": [],
    }

    for intron_start, intron_end in introns:
        if strand == "+":
            donor_boundary = intron_start - 1
            acceptor_boundary = intron_end + 1

            donor_core = (
                donor_boundary,
                intron_start + donor_core_intronic - 1,
            )
            acceptor_core = (
                intron_end - acceptor_core_intronic + 1,
                acceptor_boundary,
            )

            donor_region = (
                donor_boundary,
                intron_start + donor_region_intronic - 1,
            )
            acceptor_region = (
                intron_end - acceptor_region_intronic + 1,
                acceptor_boundary,
            )

        elif strand == "-":
            # On the minus strand, transcription runs high -> low.
            # The donor is therefore at the high-coordinate end of the intron,
            # and the acceptor is at the low-coordinate end.
            donor_boundary = intron_end + 1
            acceptor_boundary = intron_start - 1

            donor_core = (
                intron_end - donor_core_intronic + 1,
                donor_boundary,
            )
            acceptor_core = (
                acceptor_boundary,
                intron_start + acceptor_core_intronic - 1,
            )

            donor_region = (
                intron_end - donor_region_intronic + 1,
                donor_boundary,
            )
            acceptor_region = (
                acceptor_boundary,
                intron_start + acceptor_region_intronic - 1,
            )

        else:
            raise ValueError("strand must be '+' or '-'")

        windows["splice_donor_core"].append(donor_core)
        windows["splice_acceptor_core"].append(acceptor_core)
        windows["splice_donor_region"].append(donor_region)
        windows["splice_acceptor_region"].append(acceptor_region)

        # Generic small junction windows. These are mainly diagnostic.
        windows["exon_intron_junction"].append((intron_start - 1, intron_start + 1))
        windows["exon_intron_junction"].append((intron_end - 1, intron_end + 1))

    return {
        key: _merge_intervals(intervals)
        for key, intervals in windows.items()
    }


def _parse_semicolon_positions(value) -> list[int]:
    """
    Parse semicolon-separated genomic positions.

    Accepts:
      - "none"
      - ""
      - NaN
      - "12345"
      - "12345;12350"
    """
    if value is None or pd.isna(value):
        return []

    text = str(value).strip()
    if not text or text == "none":
        return []

    positions: list[int] = []

    for item in text.split(";"):
        item = item.strip()
        if not item or item == "none":
            continue

        try:
            positions.append(int(item))
        except ValueError:
            continue

    return positions


def _parse_terminus_transcripts(value) -> list[str]:
    """Parse semicolon-separated transcript labels from terminus_transcripts."""
    if value is None or pd.isna(value):
        return []

    text = str(value).strip()
    if not text or text == "none":
        return []

    return [item.strip() for item in text.split(";") if item.strip()]


def _transcript_label_to_id_map(gene_id: str, db) -> dict[str, str]:
    """
    Map transcript labels back to transcript IDs for one gene.

    FlyBase GTF example
    -------------------
    transcript_id     = FBtr0078054
    transcript_symbol = Polr1B-RA

    terminus_transcripts usually stores transcript_symbol values, so this map
    allows us to recover the corresponding transcript_id.
    """
    try:
        gene = db[gene_id]
    except Exception:
        return {}

    mapping: dict[str, str] = {}

    transcript_types = [
        "transcript",
        "mRNA",
        "ncRNA",
        "lnc_RNA",
        "tRNA",
        "rRNA",
        "snRNA",
        "snoRNA",
    ]

    for tx in db.children(gene, featuretype=transcript_types, order_by="start"):
        tx_id_values = tx.attributes.get("transcript_id", [])
        tx_symbol_values = tx.attributes.get("transcript_symbol", [])
        tx_name_values = tx.attributes.get("transcript_name", [])

        possible_ids = {str(tx.id)}
        for value in tx_id_values:
            if value:
                possible_ids.add(str(value))

        labels = set(possible_ids)

        for value in tx_symbol_values:
            if value:
                labels.add(str(value))

        for value in tx_name_values:
            if value:
                labels.add(str(value))

        canonical_id = str(tx_id_values[0]) if tx_id_values else str(tx.id)

        for label in labels:
            mapping[label] = canonical_id

    return mapping


def _exons_for_transcript(
    gene_id: str,
    transcript_id: str,
    db,
) -> list[tuple[int, int]]:
    """
    Return exon intervals for a specific transcript ID.

    This uses feature attributes rather than only relying on gffutils hierarchy.
    That is more robust for GTF files where transcript membership is encoded via
    transcript_id attributes.
    """
    try:
        gene = db[gene_id]
    except Exception:
        return []

    exons: list[tuple[int, int]] = []

    for feat in db.children(gene, featuretype="exon", order_by="start"):
        tx_ids = _feature_transcript_ids(feat)

        if transcript_id in tx_ids:
            exons.append((int(feat.start), int(feat.end)))

    return _merge_intervals(exons)


def _build_relevant_transcript_splice_model(
    gene_id: str,
    gene_strand: str,
    transcript_labels: list[str],
    db,
    *,
    donor_core_intronic: int,
    acceptor_core_intronic: int,
    donor_region_intronic: int,
    acceptor_region_intronic: int,
) -> dict[str, list[tuple[int, int]]]:
    """
    Build intron/splice windows for transcripts relevant to this terminus.

    This avoids false positives caused by other isoforms of the same gene.
    For example, an exon in transcript RA may be intronic in transcript RB.
    If the terminus only tags RA, RB-specific introns should not be used to
    penalise the design.
    """
    label_to_id = _transcript_label_to_id_map(gene_id, db)

    transcript_ids: list[str] = []
    for label in transcript_labels:
        transcript_ids.append(label_to_id.get(label, label))

    union: dict[str, list[tuple[int, int]]] = {
        "intron": [],
        "splice_donor_core": [],
        "splice_acceptor_core": [],
        "splice_donor_region": [],
        "splice_acceptor_region": [],
        "exon_intron_junction": [],
    }

    for transcript_id in transcript_ids:
        exons = _exons_for_transcript(gene_id, transcript_id, db)
        if not exons:
            continue

        tx_windows = _splice_windows_for_transcript(
            exons=exons,
            strand=gene_strand,
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


def _point_overlaps_any(
    pos: int,
    intervals: list[tuple[int, int]],
) -> bool:
    """Return True if genomic position lies inside any interval."""
    return any(start <= int(pos) <= end for start, end in intervals)


def annotate_edit_splice_risk(
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
    Annotate whether actual donor edit positions overlap introns/splice sites.

    This function is edit-position based. It does not look at the entire gRNA.
    It only evaluates bases that were actually mutated in the donor arm.

    Required columns
    ----------------
    - gene_id
    - gene_strand
    - terminus_transcripts
    - HAL_mutation_positions
    - HAR_mutation_positions

    Added columns
    -------------
    - edit_overlaps_intron
    - edit_overlaps_splice_donor_core
    - edit_overlaps_splice_acceptor_core
    - edit_overlaps_splice_core
    - edit_overlaps_splice_donor_region
    - edit_overlaps_splice_acceptor_region
    - edit_overlaps_splice_region
    - edit_splice_risk
    - edit_splice_warning

    Risk interpretation
    -------------------
    splice_core_edit:
        An edit hits the donor/acceptor core. This is the strongest warning and
        may be worth rejecting.

    splice_region_edit:
        An edit hits a broader donor/acceptor region but not the core.

    intronic_edit:
        An edit lies inside an intron of a relevant transcript.

    none:
        No edit-splice risk detected.
    """
    df = guides_df.copy()

    if df.empty:
        return df

    required = [
        "gene_id",
        "gene_strand",
        "terminus_transcripts",
        "HAL_mutation_positions",
        "HAR_mutation_positions",
    ]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"annotate_edit_splice_risk missing required columns: {missing}")

    for column, default in [
        ("edit_overlaps_intron", False),
        ("edit_overlaps_splice_donor_core", False),
        ("edit_overlaps_splice_acceptor_core", False),
        ("edit_overlaps_splice_core", False),
        ("edit_overlaps_splice_donor_region", False),
        ("edit_overlaps_splice_acceptor_region", False),
        ("edit_overlaps_splice_region", False),
        ("edit_splice_risk", "none"),
        ("edit_splice_warning", "none"),
    ]:
        if column not in df.columns:
            df[column] = default

    if "warnings" not in df.columns:
        df["warnings"] = "none"

    model_cache: dict[
        tuple[str, str, tuple[str, ...]],
        dict[str, list[tuple[int, int]]],
    ] = {}

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(
            iterator,
            total=len(df),
            desc="Annotating edit splice risk",
            leave=False,
        )

    for row in iterator:
        idx = row.Index

        mutation_positions = (
            _parse_semicolon_positions(getattr(row, "HAL_mutation_positions"))
            + _parse_semicolon_positions(getattr(row, "HAR_mutation_positions"))
        )

        if not mutation_positions:
            continue

        gene_id = str(row.gene_id)
        gene_strand = str(row.gene_strand)
        transcript_labels = _parse_terminus_transcripts(row.terminus_transcripts)

        if not transcript_labels:
            warning = "cannot assess edit splice risk: no terminus_transcripts annotation"
            df.at[idx, "edit_splice_warning"] = warning
            df.at[idx, "warnings"] = merge_warn(df.at[idx, "warnings"], warning)
            continue

        cache_key = (gene_id, gene_strand, tuple(sorted(transcript_labels)))

        if cache_key not in model_cache:
            model_cache[cache_key] = _build_relevant_transcript_splice_model(
                gene_id,
                gene_strand,
                transcript_labels,
                db,
                donor_core_intronic=donor_core_intronic,
                acceptor_core_intronic=acceptor_core_intronic,
                donor_region_intronic=donor_region_intronic,
                acceptor_region_intronic=acceptor_region_intronic,
            )

        model = model_cache[cache_key]

        intron_hit = any(
            _point_overlaps_any(position, model["intron"])
            for position in mutation_positions
        )
        donor_core_hit = any(
            _point_overlaps_any(position, model["splice_donor_core"])
            for position in mutation_positions
        )
        acceptor_core_hit = any(
            _point_overlaps_any(position, model["splice_acceptor_core"])
            for position in mutation_positions
        )
        donor_region_hit = any(
            _point_overlaps_any(position, model["splice_donor_region"])
            for position in mutation_positions
        )
        acceptor_region_hit = any(
            _point_overlaps_any(position, model["splice_acceptor_region"])
            for position in mutation_positions
        )

        core_hit = donor_core_hit or acceptor_core_hit
        region_hit = donor_region_hit or acceptor_region_hit

        df.at[idx, "edit_overlaps_intron"] = intron_hit
        df.at[idx, "edit_overlaps_splice_donor_core"] = donor_core_hit
        df.at[idx, "edit_overlaps_splice_acceptor_core"] = acceptor_core_hit
        df.at[idx, "edit_overlaps_splice_core"] = core_hit
        df.at[idx, "edit_overlaps_splice_donor_region"] = donor_region_hit
        df.at[idx, "edit_overlaps_splice_acceptor_region"] = acceptor_region_hit
        df.at[idx, "edit_overlaps_splice_region"] = region_hit

        if core_hit:
            risk = "splice_core_edit"
        elif region_hit:
            risk = "splice_region_edit"
        elif intron_hit:
            risk = "intronic_edit"
        else:
            risk = "none"

        df.at[idx, "edit_splice_risk"] = risk

        if risk != "none":
            warning = f"blocking edit overlaps {risk.replace('_', ' ')}"
            df.at[idx, "edit_splice_warning"] = warning
            df.at[idx, "warnings"] = merge_warn(df.at[idx, "warnings"], warning)

    return df