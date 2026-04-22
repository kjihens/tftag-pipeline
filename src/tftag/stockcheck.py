"""
Stock VCF / reference FASTA compatibility checks.
"""
from __future__ import annotations

import pysam


def check_vcf_reference_compatibility(
    vcf_path: str,
    fasta_path: str,
    max_records: int | None = 10000,
) -> dict:
    """
    Check whether a stock VCF is compatible with the reference FASTA.

    Returns a summary dict with:
      - contigs_only_in_vcf
      - contigs_only_in_fasta
      - length_mismatches
      - n_records_checked
      - n_ref_mismatches
      - mismatch_examples
    """
    vcf = pysam.VariantFile(vcf_path)
    fasta = pysam.FastaFile(fasta_path)

    vcf_contigs = dict(vcf.header.contigs)
    fasta_refs = set(fasta.references)
    fasta_lengths = dict(zip(fasta.references, fasta.lengths))

    contigs_only_in_vcf = sorted(set(vcf_contigs) - fasta_refs)
    contigs_only_in_fasta = sorted(fasta_refs - set(vcf_contigs))

    length_mismatches = []
    for contig in sorted(set(vcf_contigs).intersection(fasta_refs)):
        vcf_len = vcf.header.contigs[contig].length
        fa_len = fasta_lengths[contig]
        if vcf_len != fa_len:
            length_mismatches.append((contig, vcf_len, fa_len))

    n_records_checked = 0
    n_ref_mismatches = 0
    mismatch_examples = []

    for rec in vcf:
        if rec.contig not in fasta_refs:
            continue

        ref_from_fasta = fasta.fetch(
            rec.contig,
            rec.pos - 1,
            rec.pos - 1 + len(rec.ref)
        ).upper()

        ref_from_vcf = rec.ref.upper()

        n_records_checked += 1
        if ref_from_fasta != ref_from_vcf:
            n_ref_mismatches += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append(
                    (rec.contig, rec.pos, ref_from_vcf, ref_from_fasta)
                )

        if max_records is not None and n_records_checked >= max_records:
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