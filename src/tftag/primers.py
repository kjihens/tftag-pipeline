"""
Validation primer design utilities for TFTag.

Design
------
For each retained guide row, design:
- PCR1 forward primer upstream of HAL in gene orientation
- PCR2 reverse primer downstream of HAR in gene orientation

Primer3 is run through progressively relaxed parameter rounds to maximise the
chance of finding usable validation primers.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import primer3
from tqdm.auto import tqdm

from .homology import arm_ranges
from .utils import get_sequence, merge_warn


def _count_trailing_gc(primer: str) -> int:
    """Count consecutive G/C bases from the 3' end of a primer."""
    count = 0

    for base in reversed(primer.upper()):
        if base in ("G", "C"):
            count += 1
        else:
            break

    return count


def _design_primers_api(template: str, p3_args: dict) -> dict:
    """
    Call primer3-py using whichever API name is available.

    primer3-py has used both:
    - bindings.design_primers
    - bindings.designPrimers
    """
    bindings = primer3.bindings
    seq_args = {
        "SEQUENCE_ID": "template",
        "SEQUENCE_TEMPLATE": template,
    }

    if hasattr(bindings, "design_primers"):
        return bindings.design_primers(seq_args, p3_args)

    return bindings.designPrimers(seq_args, p3_args)


def _pick_single_primer(
    seq_template: str,
    *,
    which: str,
    rounds: list[dict],
    max_end_gc: int = 3,
    show_reason: bool = False,
):
    """
    Pick one primer from a sequence template.

    Multiple parameter rounds are tried from strict to permissive. The first
    acceptable primer from the first successful round is returned.
    """
    which = which.lower()
    if which not in ("left", "right"):
        raise ValueError("which must be 'left' or 'right'")

    last_res = None

    for round_idx, params in enumerate(rounds, start=1):
        p3_args = {
            "PRIMER_OPT_SIZE": params["size_opt"],
            "PRIMER_MIN_SIZE": params["size_min"],
            "PRIMER_MAX_SIZE": params["size_max"],
            "PRIMER_OPT_TM": params["tm_opt"],
            "PRIMER_MIN_TM": params["tm_min"],
            "PRIMER_MAX_TM": params["tm_max"],
            "PRIMER_MIN_GC": params["gc_min"],
            "PRIMER_MAX_GC": params["gc_max"],
            "PRIMER_MAX_HAIRPIN_TH": params["hairpin_max"],
            "PRIMER_MAX_POLY_X": params["poly_x_max"],
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_NUM_RETURN": 20,
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_PICK_LEFT_PRIMER": 1 if which == "left" else 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1 if which == "right" else 0,
            "PRIMER_GC_CLAMP": params["clamp_min"],
        }

        res = _design_primers_api(seq_template, p3_args)
        last_res = res

        key_prefix = "PRIMER_LEFT" if which == "left" else "PRIMER_RIGHT"
        n_returned = int(res.get(f"{key_prefix}_NUM_RETURNED", 0) or 0)

        candidates: list[tuple[str, float, float, int]] = []

        for i in range(n_returned):
            seq = res[f"{key_prefix}_{i}_SEQUENCE"]
            tm = float(res.get(f"{key_prefix}_{i}_TM", float("nan")))
            gc = float(res.get(f"{key_prefix}_{i}_GC_PERCENT", float("nan")))
            length = len(seq)

            trailing_gc = _count_trailing_gc(seq)
            if trailing_gc < params["clamp_min"] or trailing_gc > max_end_gc:
                continue

            candidates.append((seq, tm, gc, length))

        if candidates:
            seq, tm, gc, length = candidates[0]
            return seq, tm, gc, length, f"round{round_idx} ({which})"

    reason = "no primer found" if not show_reason else str(last_res)
    return None, None, None, None, reason


def design_validation_primers(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    *,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Design validation primers for retained guide rows.

    Adds:
    - PCR1_forward_primer_seq
    - PCR1_forward_primer_Tm
    - PCR1_forward_primer_GC
    - PCR1_forward_primer_len
    - PCR2_reverse_primer_seq
    - PCR2_reverse_primer_Tm
    - PCR2_reverse_primer_GC
    - PCR2_reverse_primer_len
    """
    rounds = [
        dict(
            tm_opt=60.0,
            tm_min=58.0,
            tm_max=62.0,
            gc_min=40.0,
            gc_max=60.0,
            size_opt=20,
            size_min=18,
            size_max=22,
            clamp_min=1,
            hairpin_max=24.0,
            poly_x_max=4,
        ),
        dict(
            tm_opt=60.0,
            tm_min=57.0,
            tm_max=63.0,
            gc_min=35.0,
            gc_max=65.0,
            size_opt=20,
            size_min=18,
            size_max=23,
            clamp_min=1,
            hairpin_max=48.0,
            poly_x_max=5,
        ),
        dict(
            tm_opt=60.0,
            tm_min=55.0,
            tm_max=65.0,
            gc_min=30.0,
            gc_max=70.0,
            size_opt=20,
            size_min=17,
            size_max=28,
            clamp_min=1,
            hairpin_max=48.0,
            poly_x_max=6,
        ),
    ]

    df = gRNA_df.copy()

    if "warnings" not in df.columns:
        df["warnings"] = "none"

    if "designable" not in df.columns:
        df["designable"] = False

    required = ["chromosome", "feature", "gene_strand", "codon_start", "codon_end"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"design_validation_primers missing required columns: {missing}")

    df["PCR1_forward_primer_seq"] = ""
    df["PCR1_forward_primer_Tm"] = np.nan
    df["PCR1_forward_primer_GC"] = np.nan
    df["PCR1_forward_primer_len"] = np.nan

    df["PCR2_reverse_primer_seq"] = ""
    df["PCR2_reverse_primer_Tm"] = np.nan
    df["PCR2_reverse_primer_GC"] = np.nan
    df["PCR2_reverse_primer_len"] = np.nan

    designable_mask = df["designable"].astype("boolean").fillna(False)
    bad_mask = ~designable_mask

    df.loc[bad_mask, "warnings"] = df.loc[bad_mask, "warnings"].apply(
        lambda w: merge_warn(w, "skipped: prefilter_designable")
    )

    iterator = df.iterrows()
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Designing validation primers", leave=False)

    for idx, row in iterator:
        if not designable_mask.loc[idx]:
            continue

        chrom = row["chromosome"]
        feature = row["feature"]
        gene_strand = row["gene_strand"]
        codon_start = int(row["codon_start"])
        codon_end = int(row["codon_end"])

        chrom_len = len(fasta_dict[chrom].seq)

        warnings = (
            []
            if df.at[idx, "warnings"] in ("", "none")
            else [w.strip() for w in str(df.at[idx, "warnings"]).split(";") if w.strip()]
        )

        (hal_s, hal_e), (har_s, har_e) = arm_ranges(
            feature,
            gene_strand,
            codon_start,
            codon_end,
        )

        if gene_strand == "+":
            up_start, up_end = hal_s - 300, hal_s - 100
            dn_start, dn_end = har_e + 100, har_e + 300
        else:
            up_start, up_end = hal_e + 100, hal_e + 300
            dn_start, dn_end = har_s - 300, har_s - 100

        out_of_bounds = (
            up_start < 1
            or up_end > chrom_len
            or dn_start < 1
            or dn_end > chrom_len
            or up_end < up_start
            or dn_end < dn_start
        )

        if out_of_bounds:
            warnings.append("Primer window out of contig despite designable=True")
            df.at[idx, "warnings"] = "; ".join(sorted(set(warnings)))
            continue

        upstream_seq = get_sequence(fasta_dict, chrom, up_start, up_end, gene_strand)
        downstream_seq = get_sequence(fasta_dict, chrom, dn_start, dn_end, gene_strand)

        f_seq, f_tm, f_gc, f_len, _ = _pick_single_primer(
            upstream_seq,
            which="left",
            rounds=rounds,
            max_end_gc=3,
        )

        if f_seq is None:
            warnings.append("No PCR1 forward primer found after relaxation")
        else:
            df.at[idx, "PCR1_forward_primer_seq"] = f_seq
            df.at[idx, "PCR1_forward_primer_Tm"] = f_tm
            df.at[idx, "PCR1_forward_primer_GC"] = f_gc
            df.at[idx, "PCR1_forward_primer_len"] = f_len

        r_seq, r_tm, r_gc, r_len, _ = _pick_single_primer(
            downstream_seq,
            which="right",
            rounds=rounds,
            max_end_gc=3,
        )

        if r_seq is None:
            warnings.append("No PCR2 reverse primer found after relaxation")
        else:
            df.at[idx, "PCR2_reverse_primer_seq"] = r_seq
            df.at[idx, "PCR2_reverse_primer_Tm"] = r_tm
            df.at[idx, "PCR2_reverse_primer_GC"] = r_gc
            df.at[idx, "PCR2_reverse_primer_len"] = r_len

        df.at[idx, "warnings"] = "; ".join(sorted(set(warnings))) if warnings else "none"

    return df