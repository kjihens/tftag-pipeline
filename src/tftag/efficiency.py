"""
On-target efficiency scoring (Rule Set 3).

This module:
  1) builds the Rule Set 3 required 30-mer sequence for each guide:
       [4 nt upstream][20 nt protospacer][3 nt PAM][3 nt downstream]
     in *target orientation* (i.e. the gRNA strand orientation).
  2) scores those 30-mers using the `rs3` Python package.

Coordinate conventions
----------------------
- protospacer_coords() and pam_coords() are assumed to return genomic coordinates
  (1-based inclusive) regardless of strand.
- get_sequence() returns the 5'->3' sequence on the requested strand; for '-' it
  reverse-complements the genomic slice.

Therefore:
- For '+' guides, we fetch the genomic slice on '+'.
- For '-' guides, we fetch the genomic slice on '-' so that the returned sequence
  is already in target orientation.
"""
from __future__ import annotations

from typing import Any, Tuple, Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from .utils import get_sequence
from .coords import protospacer_coords, pam_coords


def _get_grna_strand(row: dict[str, Any]) -> Optional[str]:
    """
    Support historical column names:
      - 'grna_strand' (preferred)
      - 'gRNA_strand' (legacy)
    """
    s = row.get("grna_strand", None)
    if s is None:
        s = row.get("gRNA_strand", None)
    if s in ("+", "-"):
        return s
    return None


def _build_rs3_30mer(row: dict[str, Any], fasta_dict) -> Tuple[Optional[str], list[str]]:
    """
    Build RS3 30-mer in target orientation:
      [4 nt upstream][20 nt protospacer][PAM(3)][3 nt downstream]

    Uses explicit protospacer/pam coords (genomic, 1-based inclusive).
    Returns (seq30 or None, warnings list).
    """
    chrom = row.get("chromosome")
    if chrom is None:
        return None, ["RS3 skipped: missing chromosome"]

    grna_strand = _get_grna_strand(row)
    if grna_strand is None:
        return None, ["RS3 skipped: missing/invalid gRNA strand"]

    prot_s, prot_e = protospacer_coords(row)
    pam_s, pam_e = pam_coords(row)
    if prot_s is None or pam_s is None:
        return None, ["RS3 skipped: missing protospacer/pam coords"]

    # Construct the genomic interval that corresponds to:
    #   4 bp upstream of protospacer ... 3 bp downstream of PAM
    #
    # '+' guide on genome: protospacer then PAM in increasing genomic coords.
    # '-' guide on genome: PAM occurs at lower genomic coords than protospacer,
    #   but we fetch on '-' strand so the returned sequence is in target orientation.
    if grna_strand == "+":
        start_30 = prot_s - 4
        end_30 = pam_e + 3
    else:  # grna_strand == "-"
        start_30 = pam_s - 3
        end_30 = prot_e + 4

    clen = len(fasta_dict[chrom].seq)
    if start_30 < 1 or end_30 > clen:
        return None, [f"RS3 skipped: 30-mer out of contig ({chrom}:{start_30}-{end_30}, len={clen})"]

    seq30 = get_sequence(fasta_dict, chrom, start_30, end_30, grna_strand)
    if len(seq30) != 30:
        return None, [f"RS3 skipped: 30-mer length {len(seq30)} (expected 30)"]

    warns: list[str] = []

    # Sanity check: in RS3 definition, PAM bases occupy positions 25–27 (1-based),
    # i.e. indices [24:27] in 0-based Python slicing.
    pam_triplet = seq30[24:27]
    if len(pam_triplet) == 3 and pam_triplet[1:3] != "GG":
        warns.append(f"RS3 30-mer PAM not NGG at expected positions (found {pam_triplet})")

    return seq30, warns


def _merge_warn(prev: Any, new: list[str]) -> str:
    """Merge warning strings into a single semicolon-delimited string."""
    prev_list = [] if prev in (None, "", "none") else [w.strip() for w in str(prev).split(";") if w.strip()]
    new_list = [w for w in new if w]
    merged = sorted(set(prev_list + new_list))
    return "; ".join(merged) if merged else "none"


def score_rs3(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    *,
    tracrRNA: str = "Hsu2013",
    show_progress: bool = True,
    batch_size: int = 1024,
) -> pd.DataFrame:
    """
    Add RS3 scoring columns:
      - rs3_30mer
      - rs3_tracrRNA
      - rs3_score
      - warnings (appended)

    Notes
    -----
    - If a 'designable' column exists, non-designable rows are skipped for RS3 30-mer build.
    - If the rs3 package is unavailable, rs3_score is left NaN and a warning is added.
    """
    df = gRNA_df.copy()
    if "warnings" not in df.columns:
        df["warnings"] = "none"

    n = len(df)
    seq30_list: list[Optional[str]] = [None] * n
    row_warns: list[list[str]] = [[] for _ in range(n)]

    # Iterate positionally to avoid any index-label surprises.
    records = df.to_dict(orient="records")
    iterator = enumerate(records)
    if show_progress:
        iterator = tqdm(iterator, total=n, desc="RS3: build 30-mers", leave=False)

    has_designable = "designable" in df.columns

    for pos, row in iterator:
        if has_designable and not bool(df["designable"].iat[pos]):
            # Keep rs3_30mer None; do not add extra warnings (prefilter already explains).
            continue

        s, ws = _build_rs3_30mer(row, fasta_dict)
        seq30_list[pos] = s
        if ws:
            row_warns[pos].extend(ws)

    df["rs3_30mer"] = seq30_list
    df["rs3_tracrRNA"] = tracrRNA

    # Score valid 30-mers with rs3 in batches.
    scores: list[float] = [np.nan] * n
    try:
        from rs3.seq import predict_seq  # type: ignore

        valid_pos = [i for i, s in enumerate(seq30_list) if s is not None]
        if valid_pos:
            pbar = tqdm(total=len(valid_pos), desc="RS3: scoring", leave=False) if show_progress else None

            for start in range(0, len(valid_pos), int(batch_size)):
                chunk = valid_pos[start : start + int(batch_size)]
                seq_chunk = [seq30_list[i] for i in chunk]  # type: ignore[list-item]
                preds = predict_seq(seq_chunk, sequence_tracr=tracrRNA)

                for j, i in enumerate(chunk):
                    scores[i] = float(preds[j])

                if pbar:
                    pbar.update(len(chunk))

            if pbar:
                pbar.close()

    except Exception as e:
        # Keep this short but informative for deployment debugging.
        msg = f"RS3_unavailable ({type(e).__name__}: {str(e)[:200]})"
        for i in range(n):
            row_warns[i].append(msg)

    df["rs3_score"] = scores

    # Merge warnings back into the DataFrame.
    df["warnings"] = [_merge_warn(df["warnings"].iat[i], row_warns[i]) for i in range(n)]
    return df