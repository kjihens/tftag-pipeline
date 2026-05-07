"""
strain-aware variant checks for TFTag.

Purpose
-------
Given a reference-designed gRNA (23-mer: 20 bp protospacer + NGG PAM),
reconstruct the corresponding sequence in a *specific injection strain*
(using its VCF), and assess whether the guide is still valid.

Core rules
----------
1. If either G in the PAM (NGG) is mutated → guide is rejected.
2. Mutation in PAM "N" is ignored.
3. Mutations in the protospacer (positions 1–20) are reported.
4. Optionally, require exact identity between reference and strain.

Key assumption
--------------
The guides_df contains *reference-derived* coordinates and sequences,
and we reconstruct the strain sequence by applying VCF variants on top
of the reference genome.

Coordinate conventions
----------------------
- All genomic coordinates are 1-based inclusive.
- pysam uses 0-based half-open intervals → we convert carefully.
"""

from __future__ import annotations

import os
from typing import Dict, Any, List

import pandas as pd
import pysam
from tqdm.auto import tqdm

from .utils import revcomp


# ---------------------------------------------------------------------
# VCF / reference compatibility check
# ---------------------------------------------------------------------

def check_vcf_reference_compatibility(
    vcf_path: str,
    fasta_path: str,
    max_records: int | None = 10000,
) -> dict:
    """
    Sanity check: does the VCF match the reference FASTA?

    This verifies:
    - contig naming overlap
    - contig length consistency
    - REF alleles match the FASTA

    This is *critical* before doing any variant-based reconstruction.
    """
    vcf = pysam.VariantFile(vcf_path)
    fasta = pysam.FastaFile(fasta_path)

    vcf_contigs = dict(vcf.header.contigs)
    fasta_refs = set(fasta.references)
    fasta_lengths = dict(zip(fasta.references, fasta.lengths))

    contigs_only_in_vcf = sorted(set(vcf_contigs) - fasta_refs)
    contigs_only_in_fasta = sorted(fasta_refs - set(vcf_contigs))

    # detect length mismatches for shared contigs
    length_mismatches = []
    for contig in sorted(set(vcf_contigs).intersection(fasta_refs)):
        vcf_len = vcf.header.contigs[contig].length
        fa_len = fasta_lengths[contig]
        if vcf_len != fa_len:
            length_mismatches.append((contig, vcf_len, fa_len))

    n_records_checked = 0
    n_ref_mismatches = 0
    mismatch_examples = []

    # check REF allele correctness
    for rec in vcf:
        if rec.contig not in fasta_refs:
            continue

        ref_from_fasta = fasta.fetch(
            rec.contig,
            rec.pos - 1,
            rec.pos - 1 + len(rec.ref)
        ).upper()

        if ref_from_fasta != rec.ref.upper():
            n_ref_mismatches += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append(
                    (rec.contig, rec.pos, rec.ref, ref_from_fasta)
                )

        n_records_checked += 1
        if max_records and n_records_checked >= max_records:
            break

    return {
        "vcf_path": vcf_path,
        "contigs_only_in_vcf": contigs_only_in_vcf,
        "contigs_only_in_fasta": contigs_only_in_fasta,
        "length_mismatches": length_mismatches,
        "n_records_checked": n_records_checked,
        "n_ref_mismatches": n_ref_mismatches,
        "mismatch_examples": mismatch_examples,
    }


def ensure_fetchable_vcf(vcf_path: str) -> None:
    """
    Enforce that VCF supports interval queries (required for performance).

    Requirements:
    - bgzipped (.vcf.gz)
    - tabix (.tbi) or CSI index

    Without this, per-guide variant lookup would be prohibitively slow.
    """
    if not vcf_path.endswith(".vcf.gz"):
        raise ValueError(
            f"strain VCF {vcf_path!r} must be bgzipped (.vcf.gz). "
            f"Convert with: bgzip -c input.vcf > input.vcf.gz"
        )

    tbi_path = vcf_path + ".tbi"
    csi_path = vcf_path + ".csi"

    if not (os.path.exists(tbi_path) or os.path.exists(csi_path)):
        raise ValueError(
            f"strain VCF {vcf_path!r} is missing index. "
            f"Run: tabix -p vcf {vcf_path}"
        )

    # quick open test
    pysam.VariantFile(vcf_path)


# ---------------------------------------------------------------------
# Core sequence reconstruction
# ---------------------------------------------------------------------

def _apply_variants_to_interval(
    vcf: pysam.VariantFile,
    fasta: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
) -> tuple[str, List[dict[str, Any]]]:
    """
    Reconstruct sequence for [start, end] by applying VCF variants.

    Strategy:
    - Walk through overlapping VCF records
    - Stitch together:
        reference segments + variant alleles
    - Return:
        - reconstructed sequence (reference orientation)
        - list of variant metadata

    Simplifications:
    - Heterozygous (0/1) treated as ALT (conservative)
    - Only first ALT allele used
    """

    ref_seq = fasta.fetch(chrom, start - 1, end).upper()
    variants: List[dict[str, Any]] = []

    records = []
    for rec in vcf.fetch(chrom, start - 1, end):
        rec_start = rec.pos
        rec_end = rec.pos + len(rec.ref) - 1

        if rec_end < start or rec_start > end:
            continue

        records.append(rec)

    if not records:
        return ref_seq, variants

    pieces: list[str] = []
    cursor = start

    for rec in sorted(records, key=lambda r: r.pos):
        rec_start = rec.pos
        rec_end = rec.pos + len(rec.ref) - 1

        if rec_start < cursor:
            variants.append({
                "pos": rec.pos,
                "ref": rec.ref.upper(),
                "alt": ",".join(str(a).upper() for a in (rec.alts or [])),
                "gt": "overlap_skipped",
            })
            continue

        if cursor < rec_start:
            pieces.append(fasta.fetch(chrom, cursor - 1, rec_start - 1).upper())

        sample = next(iter(rec.samples.values())) if rec.samples else None
        gt = sample.get("GT") if sample else None

        alt = None
        gt_label = "NA"

        if gt == (0, 0):
            gt_label = "0/0"
        elif gt is not None and any(a not in (0, None) for a in gt):
            alt_indices = [a for a in gt if a not in (0, None)]
            alt_idx = alt_indices[0] if alt_indices else None
            alt = (
                rec.alts[alt_idx - 1]
                if alt_idx and rec.alts and alt_idx <= len(rec.alts)
                else None
            )
            gt_label = "/".join("." if x is None else str(x) for x in gt)
        else:
            gt_label = str(gt)

        pieces.append(str(alt).upper() if alt is not None else rec.ref.upper())

        variants.append({
            "pos": rec.pos,
            "ref": rec.ref.upper(),
            "alt": ",".join(str(a).upper() for a in (rec.alts or [])),
            "gt": gt_label,
        })

        cursor = rec_end + 1

    if cursor <= end:
        pieces.append(fasta.fetch(chrom, cursor - 1, end).upper())

    return "".join(pieces), variants

def _orient_to_guide(seq: str, strand: str) -> str:
    """
    Convert reference-oriented sequence into guide orientation.

    Important:
    - guides are defined in their own 5'→3' orientation
    - minus-strand guides require reverse-complement
    """
    if strand == "+":
        return seq
    if strand == "-":
        return revcomp(seq)
    raise ValueError(f"Invalid grna_strand: {strand!r}")


def _guide_interval_23mer(row) -> tuple[int, int]:
    """
    Compute genomic span of the full 23-mer.

    For + strand:
        protospacer ----> PAM
    For - strand:
        PAM <---- protospacer

    Therefore:
    + : start = protospacer_start, end = pam_end
    - : start = pam_start, end = protospacer_end
    """
    strand = row["grna_strand"]
    if strand == "+":
        return int(row["protospacer_start"]), int(row["pam_end"])
    if strand == "-":
        return int(row["pam_start"]), int(row["protospacer_end"])
    raise ValueError(f"Invalid grna_strand: {strand!r}")


# ---------------------------------------------------------------------
# PAM + mutation logic
# ---------------------------------------------------------------------

def _pam_gg_mutated(ref23: str, strain23: str) -> bool:
    """
    Critical rule: if PAM GG is mutated → guide will not cut.

    Positions (1-based):
        21 = N
        22 = G
        23 = G
    """
    if len(ref23) != 23 or len(strain23) != 23:
        return True
    return (ref23[21] != strain23[21]) or (ref23[22] != strain23[22])


def _nonpam_mutated(ref23: str, strain23: str) -> bool:
    """
    Detect mutations in protospacer only (positions 1–20).
    """
    if len(ref23) != 23 or len(strain23) != 23:
        return True
    return ref23[:20] != strain23[:20]


def _diff_positions(ref23: str, strain23: str) -> str:
    """
    Return human-readable mutation summary:
    e.g. "5:A>G;17:C>T"
    """
    if len(ref23) != len(strain23):
        return f"length:{len(ref23)}>{len(strain23)}"

    return ";".join(
        f"{i}:{r}>{s}"
        for i, (r, s) in enumerate(zip(ref23, strain23), start=1)
        if r != s
    )


# ---------------------------------------------------------------------
# Main annotation
# ---------------------------------------------------------------------

def _empty_strain_record(status: str, strain_name: str | None = None) -> dict:
    return {
        "strain_name": strain_name or "",
        "strain_seq_23": "",
        "strain_seq_matches_ref": False,
        "strain_pam_gg_mutated": True,
        "strain_nonpam_mutated": False,
        "strain_variant_positions": "none",
        "strain_variant_descriptions": "none",
        "strain_check_status": status,
    }

def annotate_strain_variants(
    guides_df: pd.DataFrame,
    strain_vcfs: Dict[str, str],
    chrom_to_strain: Dict[str, str],
    reference_fasta_path: str,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Main entry point.

    For each guide:
    1. Determine which strain applies (chromosome → strain mapping)
    2. Reconstruct the 23-mer sequence in that strain
    3. Compare against reference
    4. Annotate mutation properties
    """
    df = guides_df.copy()

    fasta = pysam.FastaFile(reference_fasta_path)
    vcf_by_strain = {
        name: pysam.VariantFile(path)
        for name, path in strain_vcfs.items()
    }

    try:

        records = []

        iterator = df.iterrows()
        if show_progress:
            iterator = tqdm(iterator, total=len(df), desc="strain annotation", leave=False)

        for _, row in iterator:
            chrom = row["chromosome"]
            strand = row["grna_strand"]
            ref23 = str(row["grna_seq_23"]).upper()

            strain_name = chrom_to_strain.get(chrom)

            if strain_name is None:
                records.append(_empty_strain_record("no_strain_for_chromosome"))
                continue

            vcf = vcf_by_strain.get(strain_name)
            if vcf is None:
                records.append({"strain_check_status": "strain_vcf_not_loaded"})
                continue

            start, end = _guide_interval_23mer(row)

            strain_ref, vars_meta = _apply_variants_to_interval(
                vcf, fasta, chrom, start, end
            )

            strain23 = _orient_to_guide(strain_ref, strand)

            records.append({
                "strain_name": strain_name,
                "strain_seq_23": strain23,
                "strain_seq_matches_ref": strain23 == ref23,
                "strain_pam_gg_mutated": _pam_gg_mutated(ref23, strain23),
                "strain_nonpam_mutated": _nonpam_mutated(ref23, strain23),
                "strain_variant_positions": _diff_positions(ref23, strain23),
                "strain_variant_descriptions": ";".join(
                    f"{v['pos']}:{v['ref']}>{v['alt']}({v['gt']})"
                    for v in vars_meta
                ),
                "strain_check_status": "ok" if len(strain23) == len(ref23) else "length_changed",
            })

        ann = pd.DataFrame(records)
        return pd.concat([df.reset_index(drop=True), ann.reset_index(drop=True)], axis=1)
    
    finally:
        fasta.close()
        for vcf in vcf_by_strain.values():
            vcf.close() 