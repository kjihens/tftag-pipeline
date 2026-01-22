"""
Coordinate helpers (1-based inclusive) for protospacer/PAM/gRNA spans.

All functions accept a dict-like row (pd.Series, dict) and use .get().

Design goal
-----------
Prefer explicit coordinates when present:
  - pam_start / pam_end
  - protospacer_start / protospacer_end

Fall back to inference only when necessary, using:
  - gRNA_start/gRNA_end (or grna_start/grna_end)
  - gRNA_strand/grna_strand

Notes
-----
- SPACER_LEN = 20, PAM_LEN = 3
- gRNA spans may represent either:
    - 20 nt protospacer only, or
    - 23 nt protospacer+PAM
"""
from __future__ import annotations

from typing import Optional, Tuple, Any

SPACER_LEN = 20
PAM_LEN = 3


def _isna(x: Any) -> bool:
    """pandas-aware NA check with a safe fallback."""
    try:
        import pandas as pd
        return pd.isna(x)
    except Exception:
        return x is None


def _first_present(row, *keys: str) -> Any:
    """
    Return the first value in row[key] that is present (not missing/NA).
    This is safer than `row.get(a) or row.get(b)` because 0/"" are valid values.
    """
    for k in keys:
        if k in row:
            v = row.get(k)
            if not _isna(v):
                return v
    return None


def _get_strand(row) -> Optional[str]:
    """Return '+' or '-', else None."""
    s = _first_present(row, "grna_strand", "gRNA_strand")
    if s is None:
        return None
    s = str(s).strip()
    return s if s in ("+", "-") else None


def pam_coords(row) -> Tuple[Optional[int], Optional[int]]:
    """
    Return PAM genomic coordinates (1-based inclusive).

    Preference order:
      1) explicit pam_start/pam_end
      2) infer from gRNA span + strand, if span is 20 or 23
    """
    ps = _first_present(row, "pam_start")
    pe = _first_present(row, "pam_end")
    if ps is not None and pe is not None:
        return int(ps), int(pe)

    gs = _first_present(row, "grna_start", "gRNA_start")
    ge = _first_present(row, "grna_end", "gRNA_end")
    strand = _get_strand(row)
    if gs is None or ge is None or strand is None:
        return None, None

    gs = int(gs)
    ge = int(ge)
    length = abs(ge - gs) + 1

    if strand == "+":
        if length == SPACER_LEN + PAM_LEN:        # 23
            return ge - (PAM_LEN - 1), ge
        if length == SPACER_LEN:                  # 20
            return ge + 1, ge + PAM_LEN
    else:  # strand == "-"
        if length == SPACER_LEN + PAM_LEN:        # 23
            return gs, gs + (PAM_LEN - 1)
        if length == SPACER_LEN:                  # 20
            return gs - PAM_LEN, gs - 1

    return None, None


def protospacer_coords(row) -> Tuple[Optional[int], Optional[int]]:
    """
    Return protospacer genomic coordinates (1-based inclusive).

    Preference order:
      1) explicit protospacer_start/protospacer_end
      2) infer from gRNA span + strand, if span is 20 or 23
    """
    ps = _first_present(row, "protospacer_start")
    pe = _first_present(row, "protospacer_end")
    if ps is not None and pe is not None:
        return int(ps), int(pe)

    gs = _first_present(row, "grna_start", "gRNA_start")
    ge = _first_present(row, "grna_end", "gRNA_end")
    strand = _get_strand(row)
    if gs is None or ge is None or strand is None:
        return None, None

    gs = int(gs)
    ge = int(ge)
    length = abs(ge - gs) + 1

    if strand == "+":
        if length == SPACER_LEN + PAM_LEN:        # 23
            return gs, ge - PAM_LEN
        if length == SPACER_LEN:                  # 20
            return gs, ge
    else:  # strand == "-"
        if length == SPACER_LEN + PAM_LEN:        # 23
            return gs + PAM_LEN, ge
        if length == SPACER_LEN:                  # 20
            return gs, ge

    return None, None