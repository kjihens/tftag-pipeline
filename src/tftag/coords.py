"""
Coordinate helpers (1-based inclusive) for protospacer/PAM/gRNA spans.
All functions accept a dict-like row (pd.Series, dict) and use .get().
"""
from __future__ import annotations
from typing import Optional, Tuple

SPACER_LEN = 20
PAM_LEN = 3


def pam_coords(row) -> Tuple[Optional[int], Optional[int]]:
    """
    Return PAM genomic coords (1-based inclusive), if available in row.
    Prefer explicit pam_start/pam_end columns; otherwise infer from gRNA span and strand.
    """
    ps = row.get("pam_start")
    pe = row.get("pam_end")
    if ps is not None and pe is not None and not _isna(ps) and not _isna(pe):
        return int(ps), int(pe)

    # Fallback inference if only gRNA_start/gRNA_end present
    gs = row.get("gRNA_start") or row.get("grna_start")
    ge = row.get("gRNA_end") or row.get("grna_end")
    strand = row.get("gRNA_strand") or row.get("grna_strand")
    if _isna(gs) or _isna(ge) or strand is None:
        return None, None

    gs = int(gs); ge = int(ge)
    length = abs(ge - gs) + 1

    if strand == "+":
        if length == 23:
            return ge - 2, ge
        elif length == 20:
            return ge + 1, ge + 3
    elif strand == "-":
        if length == 23:
            return gs, gs + 2
        elif length == 20:
            return gs - 3, gs - 1

    # Unknown
    return None, None


def protospacer_coords(row) -> Tuple[Optional[int], Optional[int]]:
    """
    Return protospacer genomic coords (1-based inclusive), if available in row.
    Prefer explicit protospacer_start/protospacer_end; otherwise infer from gRNA span and strand.
    """
    ps = row.get("protospacer_start")
    pe = row.get("protospacer_end")
    if ps is not None and pe is not None and not _isna(ps) and not _isna(pe):
        return int(ps), int(pe)

    gs = row.get("gRNA_start") or row.get("grna_start")
    ge = row.get("gRNA_end") or row.get("grna_end")
    strand = row.get("gRNA_strand") or row.get("grna_strand")
    if _isna(gs) or _isna(ge) or strand is None:
        return None, None

    gs = int(gs); ge = int(ge)
    length = abs(ge - gs) + 1

    if strand == "+":
        if length == 23:
            return gs, ge - 3
        elif length == 20:
            return gs, ge
    elif strand == "-":
        if length == 23:
            return gs + 3, ge
        elif length == 20:
            return gs, ge

    return None, None


def _isna(x) -> bool:
    try:
        import pandas as pd
        return pd.isna(x)
    except Exception:
        return x is None