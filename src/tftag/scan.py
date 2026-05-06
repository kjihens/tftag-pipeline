"""
Scan for NGG PAMs and emit candidate guides around codons.

Outputs are explicit and strand-aware:
- protospacer_start / protospacer_end : genomic 1-based inclusive
- pam_start / pam_end                 : genomic 1-based inclusive
- grna_seq_23                         : 23-mer in guide orientation
- spacer                              : first 20 nt of grna_seq_23
- pam_seq                             : last 3 nt of grna_seq_23
- grna_strand                         : '+' or '-'

Important design choice
-----------------------
Each scanned guide inherits the full metadata of the terminus row it came from.
This ensures annotate.py columns such as terminus_transcripts,
terminus_isoforms, and potential_readthrough propagate into the candidate table.

Window logic
------------
The biological search criterion is based on cut distance from the codon, not on
whether the full 23-mer fits inside the requested +/- window.

Therefore we extract an expanded sequence window but only emit guides whose cut
site is within max_cut_distance of the codon interval.

Coordinate convention
---------------------
All genomic coordinates are 1-based inclusive.
"""

from __future__ import annotations

import re
from typing import List, Mapping

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm.auto import tqdm

from .utils import clamp_range, get_sequence


SPACER_LEN = 20
PAM_LEN = 3
GRNA_LEN = SPACER_LEN + PAM_LEN


PAM_REGEX_PLUS = re.compile(rf"(?=(.{{{SPACER_LEN}}})(.GG))")
PAM_REGEX_MINUS = re.compile(rf"(?=(CC.)(.{{{SPACER_LEN}}}))")


def _define_grna_space(
    fasta_dict: Mapping[str, SeqRecord],
    chrom: str,
    codon_start: int,
    codon_end: int,
    *,
    window_up: int = 30,
    window_down: int = 30,
    flank_for_full_grna: int = GRNA_LEN,
) -> tuple[str, int, int]:
    """
    Extract an expanded genomic window around the codon.

    The extraction window is deliberately larger than the requested biological
    search window. This prevents missing valid guides whose cut site is close to
    the codon but whose full 20 bp spacer extends beyond the original +/- window.
    """
    s = min(codon_start, codon_end) - int(window_up) - int(flank_for_full_grna)
    e = max(codon_start, codon_end) + int(window_down) + int(flank_for_full_grna)

    chrom_len = len(fasta_dict[chrom].seq)
    s, e = clamp_range(s, e, chrom_len)

    seq_plus = get_sequence(fasta_dict, chrom, s, e, "+")
    return seq_plus, s, e


def _cut_distance(codon_s: int, codon_e: int, cut_pos: int) -> int:
    """
    Return distance from SpCas9 cut position to codon interval.

    A cut immediately adjacent to or inside the codon is treated as distance 0.
    """
    if (codon_s - 1) <= cut_pos <= (codon_e + 1):
        return 0

    return min(abs(cut_pos - codon_s), abs(cut_pos - codon_e))


def _row_metadata(row) -> dict:
    """
    Convert an attribute-row namedtuple into metadata copied to every guide row.
    """
    return row._asdict().copy()


def _find_pams(
    seq_plus: str,
    window_start: int,
    row,
    *,
    max_cut_distance: int,
) -> List[dict]:
    """
    Find + and - strand SpCas9 sites in an expanded '+'-oriented sequence window.

    Only guides with cut_distance <= max_cut_distance are emitted.
    """
    out: List[dict] = []

    chrom = row.chromosome
    codon_start = int(row.codon_start)
    codon_end = int(row.codon_end)
    codon_s = min(codon_start, codon_end)
    codon_e = max(codon_start, codon_end)

    base = _row_metadata(row)

    # ---------- + strand guides ----------
    for match in PAM_REGEX_PLUS.finditer(seq_plus):
        spacer_start_win = match.start(1)

        spacer = match.group(1).upper()
        pam_seq = match.group(2).upper()
        grna_seq_23 = spacer + pam_seq

        prot_start = window_start + spacer_start_win
        prot_end = prot_start + (SPACER_LEN - 1)
        pam_start = prot_end + 1
        pam_end = prot_end + PAM_LEN

        # For + guides, the cut occurs 3 bp upstream of PAM.
        # Stored as the last genomic base before the cut.
        cut_pos = pam_start - 4
        cut_distance = _cut_distance(codon_s, codon_e, cut_pos)

        if cut_distance > max_cut_distance:
            continue

        guide = base.copy()
        guide.update(
            {
                "grna_strand": "+",
                "spacer": spacer,
                "pam_seq": pam_seq,
                "protospacer_start": prot_start,
                "protospacer_end": prot_end,
                "pam_start": pam_start,
                "pam_end": pam_end,
                "grna_seq_23": grna_seq_23,
                "cut_pos": cut_pos,
                "cut_distance": cut_distance,
            }
        )
        out.append(guide)

    # ---------- - strand guides ----------
    for match in PAM_REGEX_MINUS.finditer(seq_plus):
        cc_start_win = match.start(1)

        if cc_start_win + GRNA_LEN > len(seq_plus):
            continue

        seg_plus = seq_plus[cc_start_win : cc_start_win + GRNA_LEN].upper()
        grna_seq_23 = str(Seq(seg_plus).reverse_complement()).upper()

        spacer = grna_seq_23[:SPACER_LEN]
        pam_seq = grna_seq_23[SPACER_LEN:GRNA_LEN]

        g0 = window_start + cc_start_win
        g1 = g0 + (GRNA_LEN - 1)

        pam_start = g0
        pam_end = g0 + (PAM_LEN - 1)
        prot_start = g0 + PAM_LEN
        prot_end = g1

        # For - guides, the cut is downstream of PAM on '+' coordinates.
        # Stored as the last genomic base before the cut.
        cut_pos = pam_end + 3
        cut_distance = _cut_distance(codon_s, codon_e, cut_pos)

        if cut_distance > max_cut_distance:
            continue

        guide = base.copy()
        guide.update(
            {
                "grna_strand": "-",
                "spacer": spacer,
                "pam_seq": pam_seq,
                "protospacer_start": prot_start,
                "protospacer_end": prot_end,
                "pam_start": pam_start,
                "pam_end": pam_end,
                "grna_seq_23": grna_seq_23,
                "cut_pos": cut_pos,
                "cut_distance": cut_distance,
            }
        )
        out.append(guide)

    return out


def scan_for_guides(
    attribute_df: pd.DataFrame,
    fasta_dict: Mapping[str, SeqRecord],
    *,
    window_up: int = 30,
    window_down: int = 30,
    max_cut_distance: int | None = None,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Scan all terminus rows and return candidate guides.

    Parameters
    ----------
    window_up, window_down:
        Used to define the expanded extraction window.

    max_cut_distance:
        Maximum allowed distance between Cas9 cut site and codon interval.
        If None, defaults to max(window_up, window_down).
    """
    if max_cut_distance is None:
        max_cut_distance = max(int(window_up), int(window_down))

    if max_cut_distance < 0:
        raise ValueError("max_cut_distance must be >= 0")

    rows: List[dict] = []

    iterator = attribute_df.itertuples(index=False)
    if show_progress:
        iterator = tqdm(
            iterator,
            total=len(attribute_df),
            desc="Scanning PAMs",
            leave=False,
        )

    for row in iterator:
        seq_plus, window_start, _ = _define_grna_space(
            fasta_dict,
            row.chromosome,
            int(row.codon_start),
            int(row.codon_end),
            window_up=window_up,
            window_down=window_down,
        )

        rows.extend(
            _find_pams(
                seq_plus,
                window_start,
                row,
                max_cut_distance=max_cut_distance,
            )
        )

    df = pd.DataFrame(rows)

    if not df.empty:
        df["cut_distance"] = pd.to_numeric(df["cut_distance"], errors="coerce")

    return df


def add_grna_23_coordinates(candidates: pd.DataFrame) -> pd.DataFrame:
    """
    Add genomic start/end coordinates for the full 23-mer target sequence.

    For '+' guides:
        protospacer_start ... pam_end

    For '-' guides:
        pam_start ... protospacer_end
    """
    out = candidates.copy()

    if out.empty:
        out["grna_23_start"] = pd.Series(dtype="Int64")
        out["grna_23_end"] = pd.Series(dtype="Int64")
        return out

    required = [
        "grna_strand",
        "protospacer_start",
        "protospacer_end",
        "pam_start",
        "pam_end",
    ]
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f"add_grna_23_coordinates missing required columns: {missing}")

    plus = out["grna_strand"] == "+"

    out["grna_23_start"] = out["pam_start"]
    out["grna_23_end"] = out["protospacer_end"]

    out.loc[plus, "grna_23_start"] = out.loc[plus, "protospacer_start"]
    out.loc[plus, "grna_23_end"] = out.loc[plus, "pam_end"]

    return out