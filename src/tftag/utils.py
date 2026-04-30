"""
Small shared utilities for TFTag.

Coordinate convention
---------------------
Unless stated otherwise, genomic intervals are 1-based inclusive.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
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
    Extract a sequence using 1-based inclusive coordinates.

    Returns the sequence 5'→3' on the requested strand.
    """
    if chrom not in fasta_dict:
        raise KeyError(f"Chromosome {chrom!r} not found in FASTA dictionary.")

    chrom_len = len(fasta_dict[chrom].seq)
    s, e = int(start), int(end)

    if clamp:
        s, e = clamp_range(s, e, chrom_len)
    elif s < 1 or e > chrom_len or e < s:
        raise ValueError(f"Out-of-bounds interval {chrom}:{s}-{e} (len={chrom_len}).")

    seq = str(fasta_dict[chrom].seq[s - 1 : e]).upper()

    if strand == "+":
        return seq
    if strand == "-":
        return revcomp(seq)

    raise ValueError("strand must be '+' or '-'")


def merge_warn(prev, new) -> str:
    """
    Merge warning strings into a semicolon-separated unique list.

    Parameters
    ----------
    prev:
        Existing warning value. Usually "none", empty, or a semicolon-separated string.

    new:
        New warning string or iterable of warning strings.

    Returns
    -------
    str
        Semicolon-separated warning string, or "none".
    """
    if isinstance(new, str):
        new_list = [new] if new and new != "none" else []
    elif isinstance(new, Iterable):
        new_list = [str(w) for w in new if w and str(w) != "none"]
    else:
        new_list = []

    prev_list = (
        []
        if prev in (None, "", "none")
        else [w.strip() for w in str(prev).split(";") if w.strip() and w.strip() != "none"]
    )

    merged = sorted(set(prev_list + new_list))
    return "; ".join(merged) if merged else "none"