"""
Tiny shared utilities (1-based inclusive coordinates).
"""
from __future__ import annotations
from typing import Tuple, Mapping
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def revcomp(seq: str) -> str:
    """Reverse-complement DNA sequence (BioPython handles IUPAC codes)."""
    return str(Seq(seq).reverse_complement()).upper()


def clamp_range(start: int, end: int, chrom_len: int) -> Tuple[int, int]:
    """Clamp a 1-based inclusive interval to [1, chrom_len]."""
    s = max(1, int(start))
    e = min(int(end), int(chrom_len))
    if e < s:
        e = s
    return s, e


def get_sequence(
    fasta_dict: Mapping[str, SeqRecord],
    chrom: str,
    start: int,
    end: int,
    strand: str,
    *,
    clamp: bool = False,
) -> str:
    """
    Extract sequence using 1-based inclusive coordinates.
    Returns 5'->3' sequence on requested strand.
    """
    if chrom not in fasta_dict:
        raise KeyError(f"Chromosome '{chrom}' not found in FASTA dictionary.")

    chrom_len = len(fasta_dict[chrom].seq)
    s, e = int(start), int(end)

    if clamp:
        s, e = clamp_range(s, e, chrom_len)
    else:
        if s < 1 or e > chrom_len or e < s:
            raise ValueError(f"Out-of-bounds interval {chrom}:{s}-{e} (len={chrom_len}).")

    seq = str(fasta_dict[chrom].seq[s - 1 : e]).upper()
    if strand == "-":
        seq = revcomp(seq)
    elif strand != "+":
        raise ValueError("strand must be '+' or '-'")
    return seq