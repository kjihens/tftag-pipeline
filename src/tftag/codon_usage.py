"""
Codon-usage utilities for TFTag.

Purpose
-------
Build, cache, load, and query codon-usage tables from a GTF/GFF annotation and
genome FASTA.

The table can then be used by editing.py to choose synonymous blocking mutations
that preserve codon-usage preference as closely as possible.

Default fallback
----------------
If codon-usage optimisation is not requested, editing.py can instead choose
synonymous codons that best preserve GC content.
"""

from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from typing import Any

import pandas as pd
from Bio.Seq import Seq
from tqdm.auto import tqdm

from .utils import get_sequence


STANDARD_CODON_TABLE = 1


def translate_codon(codon: str) -> str:
    """Translate one codon using the standard genetic code."""
    codon = str(codon).upper()
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return "X"
    return str(Seq(codon).translate(table=STANDARD_CODON_TABLE))


def codon_gc_fraction(codon: str) -> float:
    """Return GC fraction for a 3-nt codon."""
    codon = str(codon).upper()
    if len(codon) != 3:
        return 0.0
    return sum(base in "GC" for base in codon) / 3.0


def synonymous_codons(codon: str) -> list[str]:
    """
    Return all synonymous codons for the amino acid encoded by `codon`.

    Stop codons are included as synonymous to other stop codons, although editing
    code should normally avoid stop-codon edits unless explicitly intended.
    """
    codon = str(codon).upper()
    aa = translate_codon(codon)

    if aa == "X":
        return []

    all_codons = [
        a + b + c
        for a in "ACGT"
        for b in "ACGT"
        for c in "ACGT"
    ]

    return sorted(c for c in all_codons if translate_codon(c) == aa)


def _feature_transcript_ids(feat) -> set[str]:
    """Extract transcript-like IDs from a gffutils feature."""
    ids: set[str] = set()

    for key in ("transcript_id", "Parent", "parent"):
        for value in feat.attributes.get(key, []):
            if value:
                ids.add(str(value))

    return ids


def _transcript_like_children(gene, db):
    """Return transcript-like children for a gene."""
    tx_types = [
        "transcript",
        "mRNA",
        "ncRNA",
        "lnc_RNA",
        "tRNA",
        "rRNA",
        "snRNA",
        "snoRNA",
    ]

    txs = list(db.children(gene, featuretype=tx_types, order_by="start"))
    if txs:
        return txs

    # Fallback for annotations where transcript containers are not explicit.
    out = []
    for child in db.children(gene, level=1, order_by="start"):
        has_cds = any(True for _ in db.children(child, featuretype="CDS"))
        if has_cds:
            out.append(child)

    return out


def _cds_features_for_transcript(gene, tx, db) -> list[Any]:
    """
    Return CDS features belonging to a transcript.

    Uses transcript_id / Parent attributes because GTF files often encode
    transcript membership through attributes rather than strict hierarchy.
    """
    tx_ids = _feature_transcript_ids(tx)
    tx_ids.add(str(tx.id))

    cds_features = []

    for feat in db.children(gene, featuretype="CDS", order_by="start"):
        feat_tx_ids = _feature_transcript_ids(feat)

        if tx_ids & feat_tx_ids:
            cds_features.append(feat)

    return cds_features


def _cds_sequence_for_transcript(gene, tx, db, fasta_dict) -> str:
    """
    Reconstruct transcript CDS sequence in coding orientation.

    CDS segments are sorted in transcript order:
      + strand: increasing genomic coordinate
      - strand: decreasing genomic coordinate
    """
    cds_features = _cds_features_for_transcript(gene, tx, db)

    if not cds_features:
        return ""

    strand = cds_features[0].strand

    if strand == "+":
        cds_features = sorted(cds_features, key=lambda f: int(f.start))
    elif strand == "-":
        cds_features = sorted(cds_features, key=lambda f: int(f.start), reverse=True)
    else:
        return ""

    pieces: list[str] = []

    for feat in cds_features:
        pieces.append(
            get_sequence(
                fasta_dict,
                feat.seqid,
                int(feat.start),
                int(feat.end),
                strand,
            )
        )

    seq = "".join(pieces).upper()

    # Trim trailing incomplete codon if annotation length is not divisible by 3.
    usable_len = (len(seq) // 3) * 3
    return seq[:usable_len]


def build_codon_usage_table(
    db,
    fasta_dict,
    *,
    genes: list[str] | None = None,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Build a codon-usage table from annotated CDS features.

    Parameters
    ----------
    db:
        gffutils FeatureDB.

    fasta_dict:
        FASTA dictionary from genome_io.load_fasta_dict().

    genes:
        Optional list of gene IDs. If None, all genes in the DB are used.

    Returns
    -------
    pd.DataFrame
        Columns:
        - amino_acid
        - codon
        - count
        - amino_acid_total
        - frequency
        - gc_fraction
    """
    if genes is None:
        genes = [gene.id for gene in db.features_of_type("gene")]

    counts: Counter[str] = Counter()

    iterator = genes
    if show_progress:
        iterator = tqdm(genes, total=len(genes), desc="Building codon usage", leave=False)

    for gene_id in iterator:
        try:
            gene = db[gene_id]
        except Exception:
            continue

        txs = _transcript_like_children(gene, db)

        for tx in txs:
            cds_seq = _cds_sequence_for_transcript(gene, tx, db, fasta_dict)

            if not cds_seq:
                continue

            for i in range(0, len(cds_seq), 3):
                codon = cds_seq[i : i + 3].upper()

                if len(codon) != 3:
                    continue

                if any(base not in "ACGT" for base in codon):
                    continue

                aa = translate_codon(codon)

                # Skip invalid/ambiguous codons, but keep stop codons if present.
                if aa == "X":
                    continue

                counts[codon] += 1

    aa_totals: dict[str, int] = defaultdict(int)

    for codon, count in counts.items():
        aa_totals[translate_codon(codon)] += int(count)

    rows: list[dict[str, Any]] = []

    for codon in sorted(counts):
        aa = translate_codon(codon)
        total = aa_totals[aa]

        rows.append(
            {
                "amino_acid": aa,
                "codon": codon,
                "count": int(counts[codon]),
                "amino_acid_total": int(total),
                "frequency": float(counts[codon] / total) if total else 0.0,
                "gc_fraction": codon_gc_fraction(codon),
            }
        )

    return pd.DataFrame(rows)


def load_codon_usage_table(path: str) -> pd.DataFrame:
    """Load a cached codon-usage table."""
    return pd.read_csv(path, sep="\t")


def write_codon_usage_table(
    table: pd.DataFrame,
    path: str,
    *,
    metadata: dict[str, Any] | None = None,
) -> None:
    """
    Write codon-usage table plus optional JSON metadata.

    Metadata is written to:
      <path>.json
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    table.to_csv(path, sep="\t", index=False)

    if metadata is not None:
        with open(path + ".json", "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2)


def get_or_build_codon_usage_table(
    *,
    db,
    fasta_dict,
    cache_path: str,
    genes: list[str] | None = None,
    force_recompute: bool = False,
    metadata: dict[str, Any] | None = None,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Load a cached codon-usage table, or build and cache one if needed.
    """
    if os.path.exists(cache_path) and not force_recompute:
        return load_codon_usage_table(cache_path)

    table = build_codon_usage_table(
        db,
        fasta_dict,
        genes=genes,
        show_progress=show_progress,
    )

    write_codon_usage_table(
        table,
        cache_path,
        metadata=metadata,
    )

    return table


def codon_usage_lookup(table: pd.DataFrame | None) -> dict[str, float]:
    """
    Convert codon-usage table into codon -> frequency mapping.

    Returns an empty dict if table is None or empty.
    """
    if table is None or table.empty:
        return {}

    required = {"codon", "frequency"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Codon-usage table missing required columns: {sorted(missing)}")

    return {
        str(row.codon).upper(): float(row.frequency)
        for row in table.itertuples(index=False)
    }


def choose_synonymous_codon(
    original_codon: str,
    required_base_index: int,
    required_new_base: str | None = None,
    *,
    mode: str = "gc",
    usage_lookup: dict[str, float] | None = None,
) -> str | None:
    """
    Choose a synonymous codon under a position constraint.

    Parameters
    ----------
    original_codon:
        Original 3-nt codon.

    required_base_index:
        0-based codon position that must change.

    required_new_base:
        If supplied, candidate codons must have this base at required_base_index.
        If None, candidates only need to differ from the original base at that index.

    mode:
        "gc" or "usage".

    usage_lookup:
        codon -> amino-acid-normalised frequency.

    Returns
    -------
    str | None
        Best synonymous codon, or None if no valid synonymous change exists.
    """
    original_codon = str(original_codon).upper()
    required_new_base = required_new_base.upper() if required_new_base else None

    if mode not in ("gc", "usage"):
        raise ValueError("mode must be one of: gc, usage")

    if len(original_codon) != 3:
        return None

    if required_base_index not in (0, 1, 2):
        return None

    if any(base not in "ACGT" for base in original_codon):
        return None

    candidates = []

    for candidate in synonymous_codons(original_codon):
        if candidate == original_codon:
            continue

        if required_new_base is not None:
            if candidate[required_base_index] != required_new_base:
                continue
        else:
            if candidate[required_base_index] == original_codon[required_base_index]:
                continue

        candidates.append(candidate)

    if not candidates:
        return None

    if mode == "usage" and usage_lookup:
        original_usage = usage_lookup.get(original_codon, 0.0)

        return min(
            candidates,
            key=lambda codon: (
                abs(usage_lookup.get(codon, 0.0) - original_usage),
                abs(codon_gc_fraction(codon) - codon_gc_fraction(original_codon)),
                codon,
            ),
        )

    # Default: preserve GC content as closely as possible.
    original_gc = codon_gc_fraction(original_codon)

    return min(
        candidates,
        key=lambda codon: (
            abs(codon_gc_fraction(codon) - original_gc),
            codon,
        ),
    )