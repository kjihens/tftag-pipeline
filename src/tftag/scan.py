"""
Scan for NGG PAMs and emit candidate guides around codons.
"""
from __future__ import annotations
import re
from typing import Dict, List, Tuple
import pandas as pd
from tqdm.auto import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
from .utils import get_sequence, clamp_range

# + strand: 20 nt spacer followed by "NGG" (we capture spacer + PAM triplet)
PAM_REGEX_PLUS  = re.compile(r"(?=(.{" + str(20) + r"})(.GG))")
# - strand: "CCN" followed by 20 nt on +; reverse-complement gives spacer+NGG
PAM_REGEX_MINUS = re.compile(r"(?=(CC.)(.{" + str(20) + r"}))")


def _define_gRNA_space(
    fasta_dict: Dict[str, SeqIO.SeqRecord],
    chrom: str,
    codon_start: int,
    codon_end: int,
    window_up: int = 30,
    window_down: int = 30,
) -> Tuple[str, int, int]:
    # Window around codon in genomic coordinates (1-based inclusive) returned in '+' orientation.
    s = min(codon_start, codon_end) - window_up
    e = max(codon_start, codon_end) + window_down
    clen = len(fasta_dict[chrom].seq)
    s, e = clamp_range(s, e, clen)
    seq_plus = get_sequence(fasta_dict, chrom, s, e, "+")
    return seq_plus, s, e


def _cut_distance(codon_s: int, codon_e: int, cut_pos: int) -> int:
    if (codon_s - 1) <= cut_pos <= (codon_e + 1):
        return 0
    return min(abs(cut_pos - codon_s), abs(cut_pos - codon_e))


def _find_PAMs(seq_plus: str, s0: int, row) -> List[dict]:
    """
    Find + and - strand gRNAs within seq_plus (which is '+' oriented).
    Emit dict rows with explicit protospacer & PAM coords (1-based inclusive).
    """
    out: List[dict] = []

    chrom = row.chromosome
    gene_id = row.gene_id
    gene_symbol = row.gene_symbol
    feature = row.feature
    terminus = row.terminus
    tag = row.tag
    gene_strand = row.gene_strand
    codon_start = int(row.codon_start)
    codon_end = int(row.codon_end)

    codon_s = min(codon_start, codon_end)
    codon_e = max(codon_start, codon_end)

    # + strand: spacer(20) + PAM(3)
    for m in PAM_REGEX_PLUS.finditer(seq_plus):
        spacer_start_win = m.start(1)  # 0-based in window
        spacer = m.group(1).upper()
        pam_seq = m.group(2).upper()   # actual 3-mer
        grna_seq_23 = spacer + pam_seq

        prot_start = s0 + spacer_start_win
        prot_end   = prot_start + 19
        pam_start  = prot_end + 1
        pam_end    = prot_end + 3

        # SpCas9 cut: 3 bp upstream of PAM on target strand.
        # In 1-based coordinates, cut position often represented as the base immediately before the cut.
        # We store cut_pos as the last base before the cut on the genomic '+' strand for '+' guides.
        cut_pos = pam_start - 4

        out.append(dict(
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            terminus=terminus,
            feature=feature,
            tag=tag,
            chromosome=chrom,
            gene_strand=gene_strand,
            codon_start=codon_start,
            codon_end=codon_end,
            grna_strand="+",
            spacer=spacer,
            pam_seq=pam_seq,
            protospacer_start=prot_start,
            protospacer_end=prot_end,
            pam_start=pam_start,
            pam_end=pam_end,
            grna_seq_23=grna_seq_23,
            cut_pos=cut_pos,
            cut_distance=_cut_distance(codon_s, codon_e, cut_pos),
        ))

    # - strand: detect CCN + 20 on '+'; reverse-complement to gRNA orientation
    for m in PAM_REGEX_MINUS.finditer(seq_plus):
        cc_start_win = m.start(1)  # 0-based in window, points to first C of CCN
        seg_plus = seq_plus[cc_start_win : cc_start_win + 23].upper()  # 23 nt on '+' (PAM-like then spacer)
        grna_seq_23 = str(Seq(seg_plus).reverse_complement()).upper()
        spacer = grna_seq_23[:20]
        pam_seq = grna_seq_23[20:23]

        # On '+' strand, seg_plus spans [g0, g0+22]
        g0 = s0 + cc_start_win
        g1 = g0 + 22

        # For '-' guide, PAM is at low coords on '+' (CCN region), but in guide orientation PAM is at end.
        # Genomic coords:
        pam_start = g0
        pam_end   = g0 + 2
        prot_start = g0 + 3
        prot_end   = g1

        # For '-' guides, cut is 3 bp upstream of PAM on the '-' strand,
        # which corresponds to 3 bp downstream of PAM on '+'.
        cut_pos = pam_end + 3

        out.append(dict(
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            terminus=terminus,
            feature=feature,
            tag=tag,
            chromosome=chrom,
            gene_strand=gene_strand,
            codon_start=codon_start,
            codon_end=codon_end,
            grna_strand="-",
            spacer=spacer,
            pam_seq=pam_seq,
            protospacer_start=prot_start,
            protospacer_end=prot_end,
            pam_start=pam_start,
            pam_end=pam_end,
            grna_seq_23=grna_seq_23,
            cut_pos=cut_pos,
            cut_distance=_cut_distance(codon_s, codon_e, cut_pos),
        ))

    return out


def scan_for_guides(attribute_df: pd.DataFrame, fasta_dict, window_up: int = 30, window_down: int = 30) -> pd.DataFrame:
    rows: List[dict] = []

    iterator = attribute_df.itertuples(index=False)
    iterator = tqdm(iterator, total=len(attribute_df), desc="Scanning PAMs", leave=False)

    for row in iterator:
        seq_plus, s0, _ = _define_gRNA_space(fasta_dict, row.chromosome, int(row.codon_start), int(row.codon_end),
                                            window_up=window_up, window_down=window_down)
        rows.extend(_find_PAMs(seq_plus, s0, row))

    df = pd.DataFrame(rows)
    if not df.empty:
        df["cut_distance"] = pd.to_numeric(df["cut_distance"], errors="coerce")
    return df