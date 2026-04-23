"""
Stock-aware variant checks for TFTag.

Design principles:
-----------------
- Work directly on the 23-mer target sequence
- Reject guides if PAM GG is mutated
- Ignore mutation in PAM "N"
- Report all other mutations
- Optionally filter to exact matches only
"""

from __future__ import annotations

from typing import Dict, Any, List
import os
import pysam
import pandas as pd
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
    Require a bgzipped + tabix-indexed VCF for interval queries.
    """
    if not vcf_path.endswith(".vcf.gz"):
        raise ValueError(
            f"Stock VCF {vcf_path!r} must be bgzipped (.vcf.gz) for interval queries. "
            f"Convert with: bgzip -c input.vcf > input.vcf.gz"
        )

    tbi_path = vcf_path + ".tbi"
    csi_path = vcf_path + ".csi"
    if not (os.path.exists(tbi_path) or os.path.exists(csi_path)):
        raise ValueError(
            f"Stock VCF {vcf_path!r} is missing a tabix/csi index. "
            f"Create one with: tabix -p vcf {vcf_path}"
        )

    try:
        vcf = pysam.VariantFile(vcf_path)
        # just touching header/contigs is enough to verify it opens
        _ = list(vcf.header.contigs.keys())[:1]
    except Exception as e:
        raise ValueError(
            f"Could not open indexed stock VCF {vcf_path!r}: {type(e).__name__}: {e}"
        ) from e
    

# ---------------------------------------------------------------------
# Core sequence reconstruction
# ---------------------------------------------------------------------

def _apply_variants_to_interval(
    vcf: pysam.VariantFile,
    fasta: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
) -> tuple[str, List[dict]]:
    """
    Reconstruct stock sequence across interval [start, end] (1-based inclusive).

    Simplified genotype handling:
    - 1/1 -> apply ALT
    - 0/1 -> apply ALT (conservative)
    - 0/0 -> keep REF
    """

    ref_seq = fasta.fetch(chrom, start - 1, end).upper()
    variants = []

    records = []
    for rec in vcf.fetch(chrom, start - 1, end):
        rec_start = rec.pos
        rec_end = rec.pos + len(rec.ref) - 1
        if rec_end < start or rec_start > end:
            continue
        records.append(rec)

    if not records:
        return ref_seq, variants

    pieces = []
    cursor = start

    for rec in sorted(records, key=lambda r: r.pos):

        rec_start = rec.pos
        rec_end = rec.pos + len(rec.ref) - 1

        # add reference segment before variant
        if cursor < rec_start:
            pieces.append(fasta.fetch(chrom, cursor - 1, rec_start - 1).upper())

        # genotype handling
        sample = next(iter(rec.samples.values())) if rec.samples else None
        gt = sample.get("GT") if sample else None

        alt = None
        if gt is None:
            pass
        elif gt == (0, 0):
            pass
        elif 1 in gt:
            alt = rec.alts[0] if rec.alts else None

        if alt is None:
            pieces.append(rec.ref.upper())
        else:
            pieces.append(str(alt).upper())

        variants.append({
            "pos": rec.pos,
            "ref": rec.ref.upper(),
            "alt": ",".join(str(a) for a in (rec.alts or [])),
        })

        cursor = rec_end + 1

    if cursor <= end:
        pieces.append(fasta.fetch(chrom, cursor - 1, end).upper())

    return "".join(pieces), variants


def _orient_to_guide(seq: str, strand: str) -> str:
    return seq if strand == "+" else revcomp(seq)


# ---------------------------------------------------------------------
# PAM + mutation logic
# ---------------------------------------------------------------------

def _pam_gg_mutated(ref23: str, stock23: str) -> bool:
    """
    Reject if either G in NGG is mutated.
    Positions:
      21 = N
      22 = G
      23 = G
    """
    return (ref23[21] != stock23[21]) or (ref23[22] != stock23[22])


def _nonpam_mutated(ref23: str, stock23: str) -> bool:
    """
    Only consider protospacer (positions 1–20).
    Ignore PAM entirely.
    """
    return ref23[:20] != stock23[:20]


def _diff_positions(ref23: str, stock23: str) -> str:
    out = []
    for i, (r, s) in enumerate(zip(ref23, stock23), start=1):
        if r != s:
            out.append(f"{i}:{r}>{s}")
    return ";".join(out)


# ---------------------------------------------------------------------
# Main annotation
# ---------------------------------------------------------------------

def annotate_stock_variants(
    guides_df: pd.DataFrame,
    stock_vcfs: Dict[str, str],
    chrom_to_stock: Dict[str, str],
    reference_fasta_path: str,
    show_progress: bool = True,
) -> pd.DataFrame:

    df = guides_df.copy()

    fasta = pysam.FastaFile(reference_fasta_path)
    vcf_by_stock = {
        name: pysam.VariantFile(path)
        for name, path in stock_vcfs.items()
    }

    records = []

    iterator = df.iterrows()
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Stock annotation", leave=False)

    for _, row in iterator:

        chrom = row["chromosome"]
        strand = row["gRNA_strand"]
        ref23 = str(row["gRNA_seq"]).upper()

        stock_name = chrom_to_stock.get(chrom)

        if stock_name is None:
            records.append({
                "stock_name": None,
                "stock_seq_23": "",
                "stock_seq_matches_ref": False,
                "stock_pam_gg_mutated": False,
                "stock_nonpam_mutated": False,
                "stock_variant_positions": "",
                "stock_variant_descriptions": "",
                "stock_check_status": "no_stock_for_chromosome",
            })
            continue

        vcf = vcf_by_stock.get(stock_name)

        start = int(row["gRNA_start"])
        end = int(row["gRNA_end"])

        stock_ref, vars_meta = _apply_variants_to_interval(
            vcf,
            fasta,
            chrom,
            start,
            end,
        )

        stock23 = _orient_to_guide(stock_ref, strand)

        exact = stock23 == ref23
        pam_bad = _pam_gg_mutated(ref23, stock23)
        nonpam = _nonpam_mutated(ref23, stock23)

        records.append({
            "stock_name": stock_name,
            "stock_seq_23": stock23,
            "stock_seq_matches_ref": exact,
            "stock_pam_gg_mutated": pam_bad,
            "stock_nonpam_mutated": nonpam,
            "stock_variant_positions": _diff_positions(ref23, stock23),
            "stock_variant_descriptions": ";".join(
                f"{v['pos']}:{v['ref']}>{v['alt']}" for v in vars_meta
            ),
            "stock_check_status": "ok",
        })

    ann = pd.DataFrame(records)
    return pd.concat([df.reset_index(drop=True), ann.reset_index(drop=True)], axis=1)