"""
On-target efficiency scoring (Rule Set 3).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from tqdm import tqdm
from .utils import get_sequence, clamp_range
from .coords import protospacer_coords, pam_coords


def _build_rs3_30mer(row, fasta_dict):
    """
    Build 30-mer in target orientation:
      [4 nt upstream][20 nt spacer][PAM(3)][3 nt downstream]
    Uses explicit protospacer/pam coords (genomic, 1-based inclusive).
    """
    warns = []
    chrom = row["chromosome"]
    grna_strand = row["grna_strand"]

    prot_s, prot_e = protospacer_coords(row)
    pam_s, pam_e = pam_coords(row)
    if prot_s is None or pam_s is None:
        return None, ["RS3 skipped: missing protospacer/pam coords"]

    if grna_strand == "+":
        # protospacer then PAM on the + strand
        start_30 = prot_s - 4
        end_30   = pam_e + 3

    elif grna_strand == "-":
        # PAM then protospacer in genomic coords; fetch reverse-complement
        start_30 = pam_s - 3
        end_30   = prot_e + 4
        
    else:
        return None, [f"RS3 skipped: invalid grna_strand '{grna_strand}'"]

    clen = len(fasta_dict[chrom].seq)
    if start_30 < 1 or end_30 > clen:
        return None, [f"RS3 skipped: 30-mer out of contig ({chrom}:{start_30}-{end_30}, len={clen})"]

    seq30 = get_sequence(fasta_dict, chrom, start_30, end_30, grna_strand)
    if len(seq30) != 30:
        return None, [f"RS3 skipped: 30-mer length {len(seq30)} (expected 30)"]

    # PAM should occupy positions 25–27 (1-based) => indices 24:27 (0-based slice)
    pam_triplet = seq30[24:27]
    if pam_triplet[1:3] != "GG":
        warns.append(f"RS3 30-mer PAM not NGG at expected positions (found {pam_triplet})")

    return seq30, warns


def score_rs3(gRNA_df: pd.DataFrame, fasta_dict, tracrRNA: str = "Hsu2013", show_progress: bool = True, batch_size: int = 1024) -> pd.DataFrame:
    df = gRNA_df.copy()
    if "warnings" not in df.columns:
        df["warnings"] = "none"

    seq30_list = [None] * len(df)
    row_warns = [[] for _ in range(len(df))]

    # iterate positionally to avoid index-label bugs
    iterator = enumerate(df.to_dict(orient="records"))
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="RS3: build 30-mers", leave=False)

    for pos, row in iterator:
        if "designable" in df.columns and not bool(df.iloc[pos]["designable"]):
            seq30_list[pos] = None
            continue
        s, ws = _build_rs3_30mer(row, fasta_dict)
        seq30_list[pos] = s
        row_warns[pos].extend(ws)

    df["rs3_30mer"] = seq30_list
    df["rs3_tracrRNA"] = tracrRNA

    # Score
    scores = [np.nan] * len(df)
    try:
        from rs3.seq import predict_seq
        valid_pos = [i for i, s in enumerate(seq30_list) if s is not None]
        if valid_pos:
            pbar = tqdm(total=len(valid_pos), desc="RS3: scoring", leave=False) if show_progress else None
            for start in range(0, len(valid_pos), batch_size):
                chunk = valid_pos[start : start + batch_size]
                seq_chunk = [seq30_list[i] for i in chunk]
                preds = predict_seq(seq_chunk, sequence_tracr=tracrRNA)
                for j, i in enumerate(chunk):
                    scores[i] = float(preds[j])
                if pbar:
                    pbar.update(len(chunk))
            if pbar:
                pbar.close()
    except Exception as e:
        msg = f"RS3_unavailable ({type(e).__name__})"
        for i in range(len(df)):
            row_warns[i].append(msg)

    df["rs3_score"] = scores

    # merge warnings
    def _merge_warn(prev, new):
        prev_list = [] if prev in (None, "", "none") else [w.strip() for w in str(prev).split(";") if w.strip()]
        new_list = [w for w in new if w]
        return "; ".join(sorted(set(prev_list + new_list))) if (prev_list or new_list) else "none"

    df["warnings"] = [_merge_warn(df.iloc[i]["warnings"], row_warns[i]) for i in range(len(df))]
    return df