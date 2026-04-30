"""
Coordinate helpers for guide, protospacer, and PAM intervals.

Coordinate convention
---------------------
All genomic coordinates are 1-based inclusive.

Design goal
-----------
Prefer explicit coordinates when available:

- pam_start / pam_end
- protospacer_start / protospacer_end

Fallback inference is only used for older/intermediate schemas containing:

- grna_start / grna_end or gRNA_start / gRNA_end
- grna_strand or gRNA_strand

Current TFTag scan output should normally contain explicit coordinates, so this
module mainly provides compatibility and centralised coordinate interpretation.
"""

from __future__ import annotations

from typing import Any

import pandas as pd


SPACER_LEN = 20
PAM_LEN = 3
GRNA_LEN = SPACER_LEN + PAM_LEN


def _isna(value: Any) -> bool:
    """Return True for None/NA/NaN-like values."""
    try:
        return bool(pd.isna(value))
    except Exception:
        return value is None


def _first_present(row, *keys: str) -> Any:
    """
    Return the first non-missing value among candidate row keys.

    This avoids using `or`, which would incorrectly skip valid falsy values.
    """
    for key in keys:
        if key in row:
            value = row.get(key)
            if not _isna(value):
                return value

    return None


def _get_strand(row) -> str | None:
    """Return guide strand '+' or '-', else None."""
    strand = _first_present(row, "grna_strand", "gRNA_strand")

    if strand is None:
        return None

    strand = str(strand).strip()
    return strand if strand in ("+", "-") else None


def pam_coords(row) -> tuple[int | None, int | None]:
    """
    Return PAM genomic coordinates.

    Preference order:
    1. explicit pam_start / pam_end
    2. infer from gRNA span + strand if span length is 20 or 23
    """
    pam_start = _first_present(row, "pam_start")
    pam_end = _first_present(row, "pam_end")

    if pam_start is not None and pam_end is not None:
        return int(pam_start), int(pam_end)

    grna_start = _first_present(row, "grna_start", "gRNA_start")
    grna_end = _first_present(row, "grna_end", "gRNA_end")
    strand = _get_strand(row)

    if grna_start is None or grna_end is None or strand is None:
        return None, None

    grna_start = int(grna_start)
    grna_end = int(grna_end)
    length = abs(grna_end - grna_start) + 1

    if strand == "+":
        if length == GRNA_LEN:
            return grna_end - (PAM_LEN - 1), grna_end
        if length == SPACER_LEN:
            return grna_end + 1, grna_end + PAM_LEN

    if strand == "-":
        if length == GRNA_LEN:
            return grna_start, grna_start + (PAM_LEN - 1)
        if length == SPACER_LEN:
            return grna_start - PAM_LEN, grna_start - 1

    return None, None


def protospacer_coords(row) -> tuple[int | None, int | None]:
    """
    Return protospacer genomic coordinates.

    Preference order:
    1. explicit protospacer_start / protospacer_end
    2. infer from gRNA span + strand if span length is 20 or 23
    """
    prot_start = _first_present(row, "protospacer_start")
    prot_end = _first_present(row, "protospacer_end")

    if prot_start is not None and prot_end is not None:
        return int(prot_start), int(prot_end)

    grna_start = _first_present(row, "grna_start", "gRNA_start")
    grna_end = _first_present(row, "grna_end", "gRNA_end")
    strand = _get_strand(row)

    if grna_start is None or grna_end is None or strand is None:
        return None, None

    grna_start = int(grna_start)
    grna_end = int(grna_end)
    length = abs(grna_end - grna_start) + 1

    if strand == "+":
        if length == GRNA_LEN:
            return grna_start, grna_end - PAM_LEN
        if length == SPACER_LEN:
            return grna_start, grna_end

    if strand == "-":
        if length == GRNA_LEN:
            return grna_start + PAM_LEN, grna_end
        if length == SPACER_LEN:
            return grna_start, grna_end

    return None, None