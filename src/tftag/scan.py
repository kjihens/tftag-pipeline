"""
Scan for NGG PAMs and emit candidate guides around codons.

Outputs are **explicit and strand-aware**:
- protospacer_start / protospacer_end : genomic 1-based inclusive (always on '+' coordinate axis)
- pam_start / pam_end                 : genomic 1-based inclusive (always on '+' coordinate axis)
- grna_seq_23                         : 23-mer in *guide orientation* (5'->3' of the gRNA; PAM at end)
- spacer                              : first 20 nt of grna_seq_23
- pam_seq                             : last 3 nt of grna_seq_23 (typically NGG)
- grna_strand                         : '+' or '-' (strand of the target site on the genome)

Note on coordinate conventions
-----------------------------
We use 1-based inclusive coordinates everywhere to match the rest of the pipeline.
"""
from __future__ import annotations

import re
from typing import Dict, List, Tuple, Mapping

import pandas as pd
from tqdm.auto import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utils import get_sequence, clamp_range


SPACER_LEN = 20
PAM_LEN = 3
GRNA_LEN = SPACER_LEN + PAM_LEN  # 23


# + strand: 20 nt spacer followed by "NGG" (captures spacer + PAM triplet)
PAM_REGEX_PLUS = re.compile(rf"(?=(.{{{SPACER_LEN}}})(.GG))")

# - strand: "CCN" followed by 20 nt on '+'; reverse-complement gives spacer+NGG in guide orientation
PAM_REGEX_MINUS = re.compile(rf"(?=(CC.)(.{{{SPACER_LEN}}}))")


def _define_grna_space(
    fasta_dict: Mapping[str, SeqRecord],
    chrom: str,
    codon_start: int,
    codon_end: int,
    *,
    window_up: int = 30,
    window_down: int = 30,
) -> Tuple[str, int, int]:
    """
    Extract a window around the codon on the genomic '+' orientation.

    Returns
    -------
    seq_plus : str
        Window sequence on '+' strand.
    s : int
        Genomic start (1-based inclusive) of window.
    e : int
        Genomic end (1-based inclusive) of window.
    """
    s = min(codon_start, codon_end) - int(window_up)
    e = max(codon_start, codon_end) + int(window_down)

    clen = len(fasta_dict[chrom].seq)
    s, e = clamp_range(s, e, clen)

    seq_plus = get_sequence(fasta_dict, chrom, s, e, "+")
    return seq_plus, s, e


def _cut_distance(codon_s: int, codon_e: int, cut_pos: int) -> int:
    """Distance from cut position to codon interval; 0 if overlapping the codon neighbourhood."""
    if (codon_s - 1) <= cut_pos <= (codon_e + 1):
        return 0
    return min(abs(cut_pos - codon_s), abs(cut_pos - codon_e))


def _find_pams(seq_plus: str, window_start: int, row) -> List[dict]:
    """
    Find + and - strand SpCas9 sites within `seq_plus` (which is '+' oriented).
    Emit dict rows with explicit protospacer & PAM coords (1-based inclusive).

    Parameters
    ----------
    seq_plus:
      Sequence window on '+' strand.
    window_start:
      Genomic 1-based inclusive start coordinate of seq_plus.
    row:
      Attribute row (namedtuple from itertuples) with gene metadata and codon coords.
    """
    out: List[dict] = []

    chrom = row.chromosome
    codon_start = int(row.codon_start)
    codon_end = int(row.codon_end)
    codon_s = min(codon_start, codon_end)
    codon_e = max(codon_start, codon_end)

    # ---------- + strand guides ----------
    # Pattern: [20 spacer][NGG]
    for m in PAM_REGEX_PLUS.finditer(seq_plus):
        spacer_start_win = m.start(1)  # 0-based in seq_plus

        spacer = m.group(1).upper()
        pam_seq = m.group(2).upper()          # actual 3-mer (NGG)
        grna_seq_23 = spacer + pam_seq        # in guide orientation already

        # Convert window offsets to genomic coordinates (1-based inclusive)
        prot_start = window_start + spacer_start_win
        prot_end = prot_start + (SPACER_LEN - 1)
        pam_start = prot_end + 1
        pam_end = prot_end + PAM_LEN

        # SpCas9 cut position: 3 bp upstream of PAM on target strand.
        # For '+' guides, cut occurs between (pam_start-4) and (pam_start-3).
        # We store cut_pos as the last base before the cut on '+' coordinates.
        cut_pos = pam_start - 4

        out.append(
            dict(
                gene_id=row.gene_id,
                gene_symbol=row.gene_symbol,
                terminus=row.terminus,
                feature=row.feature,
                tag=row.tag,
                chromosome=chrom,
                gene_strand=row.gene_strand,
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
            )
        )

    # ---------- - strand guides ----------
    # Detect CCN + 20 on '+'; reverse-complement the 23-mer to get guide orientation (spacer+NGG).
    for m in PAM_REGEX_MINUS.finditer(seq_plus):
        cc_start_win = m.start(1)  # 0-based in seq_plus; first base of CCN

        # Ensure we can take a full 23-mer from the window
        if cc_start_win + GRNA_LEN > len(seq_plus):
            continue

        seg_plus = seq_plus[cc_start_win : cc_start_win + GRNA_LEN].upper()
        grna_seq_23 = str(Seq(seg_plus).reverse_complement()).upper()
        spacer = grna_seq_23[:SPACER_LEN]
        pam_seq = grna_seq_23[SPACER_LEN:GRNA_LEN]

        # Genomic coordinates of the 23-mer segment on '+' are:
        g0 = window_start + cc_start_win
        g1 = g0 + (GRNA_LEN - 1)

        # For '-' guide sites:
        # On '+' coordinates, the PAM-like CCN is at low coords (g0..g0+2),
        # and the protospacer is at higher coords (g0+3..g1).
        pam_start = g0
        pam_end = g0 + (PAM_LEN - 1)
        prot_start = g0 + PAM_LEN
        prot_end = g1

        # For '-' guides, the cut is 3 bp upstream of PAM on the '-' strand,
        # which corresponds to 3 bp downstream of PAM on '+' coordinates.
        # Cut occurs between (pam_end+3) and (pam_end+4); store last base before cut.
        cut_pos = pam_end + 3

        out.append(
            dict(
                gene_id=row.gene_id,
                gene_symbol=row.gene_symbol,
                terminus=row.terminus,
                feature=row.feature,
                tag=row.tag,
                chromosome=chrom,
                gene_strand=row.gene_strand,
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
            )
        )

    return out


def scan_for_guides(
    attribute_df: pd.DataFrame,
    fasta_dict: Mapping[str, SeqRecord],
    *,
    window_up: int = 30,
    window_down: int = 30,
    show_progress: bool = True,
) -> pd.DataFrame:
    """Scan all attribute rows and return candidate guides as a DataFrame."""
    rows: List[dict] = []

    iterator = attribute_df.itertuples(index=False)
    if show_progress:
        iterator = tqdm(iterator, total=len(attribute_df), desc="Scanning PAMs", leave=False)

    for row in iterator:
        seq_plus, s0, _ = _define_grna_space(
            fasta_dict,
            row.chromosome,
            int(row.codon_start),
            int(row.codon_end),
            window_up=window_up,
            window_down=window_down,
        )
        rows.extend(_find_pams(seq_plus, s0, row))

    df = pd.DataFrame(rows)
    if not df.empty:
        df["cut_distance"] = pd.to_numeric(df["cut_distance"], errors="coerce")
    return df