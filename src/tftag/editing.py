"""
Silent-edit utilities for TFTag donor design.

Responsibilities
----------------
- Decide whether an HDR donor arm can recreate the gRNA target site.
- Mark which arm, if any, requires mutation.
- Apply synonymous edits where possible to prevent re-cutting.

Editing rule
------------
An arm requires mutation only if it contains:
- the full PAM
- at least N PAM-proximal protospacer bases

If the donor does not include the PAM, it cannot recreate a functional SpCas9
target site, so no edit is required.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from tqdm.auto import tqdm

from .utils import merge_warn


def _genomic_to_coding_index(
    pos: int,
    arm_start: int,
    arm_end: int,
    gene_strand: str,
) -> Optional[int]:
    """Map genomic coordinate to 0-based index within a gene-oriented arm."""
    if pos < arm_start or pos > arm_end:
        return None
    if gene_strand == "+":
        return pos - arm_start
    if gene_strand == "-":
        return arm_end - pos
    raise ValueError("gene_strand must be '+' or '-'")


def _coding_index_to_genomic(
    idx: int,
    arm_start: int,
    arm_end: int,
    gene_strand: str,
) -> int:
    """Inverse of `_genomic_to_coding_index()`."""
    if gene_strand == "+":
        return arm_start + idx
    if gene_strand == "-":
        return arm_end - idx
    raise ValueError("gene_strand must be '+' or '-'")


def _indices_within_arm(
    arm_start: int,
    arm_end: int,
    gene_strand: str,
    positions: List[int],
) -> List[int]:
    """Return gene-oriented arm indices for genomic positions inside the arm."""
    out: list[int] = []
    for pos in positions:
        idx = _genomic_to_coding_index(pos, arm_start, arm_end, gene_strand)
        if idx is not None:
            out.append(idx)
    return out


def _translate_codon(codon: str) -> str:
    """Translate a codon using the standard genetic code."""
    return str(Seq(codon).translate(table=1))


def _codon_start_for_idx(idx: int, *, anchor: int) -> int:
    """
    Return codon start index for a sequence index using a known frame anchor.

    `anchor` is an index known to be codon phase 0 in the gene-oriented arm.
    """
    return idx - ((idx - anchor) % 3)


def _try_synonymous_change(
    coding_seq: str,
    idx: int,
    *,
    anchor: int,
) -> Optional[Tuple[str, Dict]]:
    """
    Try one synonymous nucleotide substitution at `idx`.

    Returns
    -------
    (mutated_sequence, mutation_record) or None
    """
    bases = ["A", "C", "G", "T"]

    if idx < 0 or idx >= len(coding_seq) or coding_seq[idx] not in bases:
        return None

    codon_start = _codon_start_for_idx(idx, anchor=anchor)
    codon_end = codon_start + 3

    if codon_start < 0 or codon_end > len(coding_seq):
        return None

    codon = list(coding_seq[codon_start:codon_end])
    aa_orig = _translate_codon("".join(codon))

    for base in bases:
        if base == coding_seq[idx]:
            continue

        codon_mut = codon.copy()
        codon_mut[idx - codon_start] = base

        if _translate_codon("".join(codon_mut)) == aa_orig:
            seq = list(coding_seq)
            seq[idx] = base

            return (
                "".join(seq),
                {
                    "idx": idx,
                    "from": coding_seq[idx],
                    "to": base,
                    "codon_from": "".join(codon),
                    "codon_to": "".join(codon_mut),
                    "codon_start": codon_start,
                    "anchor": anchor,
                },
            )

    return None


def _two_synonymous_changes_near(
    coding_seq: str,
    candidate_indices: List[int],
    *,
    anchor: int,
) -> Optional[Tuple[str, List[Dict]]]:
    """
    Try to apply two synonymous changes using candidate indices in priority order.
    """
    seq = coding_seq
    muts: list[dict] = []

    for idx in candidate_indices:
        attempt = _try_synonymous_change(seq, idx, anchor=anchor)
        if attempt is None:
            continue

        seq, rec = attempt
        muts.append(rec)

        if len(muts) == 2:
            return seq, muts

    return None


def _fmt_mut_locations(
    chrom: str,
    gene_strand: str,
    arm_start: int,
    arm_end: int,
    muts: List[Dict],
) -> str:
    """Format genomic edit locations as e.g. A-3L:25668456-G."""
    if not muts:
        return "none"

    def comp(base: str) -> str:
        return base.translate(str.maketrans("ACGT", "TGCA"))

    items = []

    for mut in muts:
        pos = _coding_index_to_genomic(mut["idx"], arm_start, arm_end, gene_strand)

        base_from = mut["from"] if gene_strand == "+" else comp(mut["from"])
        base_to = mut["to"] if gene_strand == "+" else comp(mut["to"])

        items.append(f"{base_from}-{chrom}:{pos}-{base_to}")

    return ", ".join(items)


def _frame_anchor_for_arm(*, feature: str, target_arm: str, arm_len: int) -> int:
    """
    Define codon-frame anchor within a gene-oriented arm sequence.

    start_codon:
        HAR starts exactly after the start codon, so anchor = 0.

    stop_codon:
        HAL ends exactly before the stop codon, so the last codon start is
        arm_len - 3.
    """
    if feature == "start_codon" and target_arm == "HAR":
        return 0

    if feature == "stop_codon" and target_arm == "HAL":
        return arm_len - 3

    raise ValueError(
        f"Unsupported frame anchor request: feature={feature}, target_arm={target_arm}"
    )


def choose_arm_for_mutation(
    gRNA_df: pd.DataFrame,
    *,
    protospacer_overlap_len: int = 13,
    coding_only: bool = True,
    show_progress: bool = True,
    debug: bool = False,
    debug_n: int = 20,
) -> pd.DataFrame:
    """
    Determine whether HAL or HAR needs mutation to prevent Cas9 re-cutting.

    Adds:
    - requires_edit_arm
    - rejected
    - reject_reason
    - HAL_pam_in_arm / HAR_pam_in_arm
    - HAL_pam_proximal_overlap / HAR_pam_proximal_overlap
    """
    if protospacer_overlap_len < 1:
        raise ValueError("protospacer_overlap_len must be >= 1")

    df = gRNA_df.copy()

    if "warnings" not in df.columns:
        df["warnings"] = "none"

    for col, default in [
        ("requires_edit_arm", "none"),
        ("rejected", False),
        ("reject_reason", "none"),
        ("HAL_pam_in_arm", False),
        ("HAR_pam_in_arm", False),
        ("HAL_pam_proximal_overlap", 0),
        ("HAR_pam_proximal_overlap", 0),
    ]:
        if col not in df.columns:
            df[col] = default

    if "designable" not in df.columns:
        df["designable"] = True

    designable_mask = df["designable"].astype("boolean").fillna(False)
    if not designable_mask.any():
        return df

    work = df.loc[designable_mask].copy()

    def interval_contains(a_s: int, a_e: int, b_s: int, b_e: int) -> bool:
        return a_s <= b_s and b_e <= a_e

    def interval_overlaps(a_s: int, a_e: int, b_s: int, b_e: int) -> bool:
        return not (a_e < b_s or b_e < a_s)

    def count_pam_proximal_spacer_bases_in_arm(
        arm_s: int,
        arm_e: int,
        prot_s: int,
        prot_e: int,
        grna_strand: str,
    ) -> int:
        """Count contiguous PAM-proximal protospacer bases present in an arm."""
        if grna_strand == "+":
            positions = range(prot_e, prot_s - 1, -1)
        elif grna_strand == "-":
            positions = range(prot_s, prot_e + 1)
        else:
            raise ValueError(f"Invalid grna_strand: {grna_strand!r}")

        count = 0
        for pos in positions:
            if arm_s <= pos <= arm_e:
                count += 1
            else:
                break

        return count

    iterator = work.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(work), desc="Choosing edit arm", leave=False)

    for n, row in enumerate(iterator):
        idx = row.Index

        if pd.isna(getattr(row, "HALs", np.nan)) or pd.isna(getattr(row, "HARs", np.nan)):
            work.at[idx, "requires_edit_arm"] = "none"
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                "Missing homology arm coordinates; cannot choose edit arm",
            )
            continue

        hal_s, hal_e = int(row.HALs), int(row.HALe)
        har_s, har_e = int(row.HARs), int(row.HARe)

        rowd = row._asdict()

        try:
            pam_s, pam_e = int(rowd["pam_start"]), int(rowd["pam_end"])
            prot_s, prot_e = int(rowd["protospacer_start"]), int(rowd["protospacer_end"])
        except Exception:
            work.at[idx, "requires_edit_arm"] = "none"
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                "Missing PAM/protospacer coordinates; cannot choose edit arm",
            )
            continue

        grna_strand = rowd.get("grna_strand") or rowd.get("gRNA_strand")
        if grna_strand not in ("+", "-"):
            work.at[idx, "requires_edit_arm"] = "none"
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                "Missing/invalid grna_strand; cannot choose edit arm",
            )
            continue

        codon_s = int(row.codon_start)
        codon_e = int(row.codon_end)

        pam_in_hal = interval_contains(hal_s, hal_e, pam_s, pam_e)
        pam_in_har = interval_contains(har_s, har_e, pam_s, pam_e)

        hal_overlap = (
            count_pam_proximal_spacer_bases_in_arm(
                hal_s, hal_e, prot_s, prot_e, grna_strand
            )
            if pam_in_hal
            else 0
        )
        har_overlap = (
            count_pam_proximal_spacer_bases_in_arm(
                har_s, har_e, prot_s, prot_e, grna_strand
            )
            if pam_in_har
            else 0
        )

        work.at[idx, "HAL_pam_in_arm"] = pam_in_hal
        work.at[idx, "HAR_pam_in_arm"] = pam_in_har
        work.at[idx, "HAL_pam_proximal_overlap"] = hal_overlap
        work.at[idx, "HAR_pam_proximal_overlap"] = har_overlap

        hal_requires = pam_in_hal and hal_overlap >= protospacer_overlap_len
        har_requires = pam_in_har and har_overlap >= protospacer_overlap_len

        required_coding_arm = "HAR" if row.feature == "start_codon" else "HAL"

        if hal_requires and har_requires:
            target_arm = "none"
            msg = f"Internal error: both arms satisfy PAM+{protospacer_overlap_len}bp rule"
            work.at[idx, "rejected"] = True
            work.at[idx, "reject_reason"] = msg
            work.at[idx, "warnings"] = merge_warn(work.at[idx, "warnings"], msg)
        elif hal_requires:
            target_arm = "HAL"
        elif har_requires:
            target_arm = "HAR"
        else:
            target_arm = "none"

        work.at[idx, "requires_edit_arm"] = target_arm

        if coding_only and target_arm != "none" and target_arm != required_coding_arm:
            work.at[idx, "rejected"] = True
            work.at[idx, "reject_reason"] = (
                f"requires edits in non-coding arm {target_arm} for {row.feature}; rejected"
            )

        if debug and n < debug_n:
            print(
                f"[DEBUG choose_arm_for_mutation] row={n} idx={idx}\n"
                f"  gene_id={row.gene_id} gene_symbol={getattr(row, 'gene_symbol', 'NA')}\n"
                f"  feature={row.feature} terminus={getattr(row, 'terminus', 'NA')}\n"
                f"  chromosome={row.chromosome} gene_strand={row.gene_strand} grna_strand={grna_strand}\n"
                f"  codon=({codon_s},{codon_e})\n"
                f"  HAL=({hal_s},{hal_e}) len={hal_e - hal_s + 1}\n"
                f"  HAR=({har_s},{har_e}) len={har_e - har_s + 1}\n"
                f"  PAM=({pam_s},{pam_e}) len={pam_e - pam_s + 1}\n"
                f"  PROT=({prot_s},{prot_e}) len={prot_e - prot_s + 1}\n"
                f"  pam_overlaps_codon={interval_overlaps(pam_s, pam_e, codon_s, codon_e)}\n"
                f"  prot_overlaps_codon={interval_overlaps(prot_s, prot_e, codon_s, codon_e)}\n"
                f"  pam_in_HAL={pam_in_hal} pam_in_HAR={pam_in_har}\n"
                f"  HAL_pam_proximal_overlap={hal_overlap}\n"
                f"  HAR_pam_proximal_overlap={har_overlap}\n"
                f"  HAL_requires={hal_requires} HAR_requires={har_requires}\n"
                f"  required_coding_arm={required_coding_arm}\n"
                f"  chosen_requires_edit_arm={target_arm}\n"
                f"  rejected={work.at[idx, 'rejected']}\n"
                f"  reject_reason={work.at[idx, 'reject_reason']}\n"
            )

    if debug:
        print("[DEBUG choose_arm_for_mutation] summary:")
        print(work["requires_edit_arm"].value_counts(dropna=False))
        print(work["rejected"].value_counts(dropna=False))

    out = df.copy()
    out.update(work)
    return out


def apply_silent_edits(
    gRNA_df: pd.DataFrame,
    *,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Apply synonymous edits to the required coding arm where possible.

    Adds:
    - HAL_seq_mut
    - HAR_seq_mut
    - HAL_mutation_location
    - HAR_mutation_location
    - edit_arm
    """
    df = gRNA_df.copy()

    if "warnings" not in df.columns:
        df["warnings"] = "none"

    for col, default in [
        ("HAL_seq_mut", ""),
        ("HAR_seq_mut", ""),
        ("HAL_mutation_location", "none"),
        ("HAR_mutation_location", "none"),
        ("edit_arm", "none"),
    ]:
        if col not in df.columns:
            df[col] = default

    if "designable" not in df.columns:
        df["designable"] = True

    designable_mask = df["designable"].astype("boolean").fillna(False)
    if not designable_mask.any():
        return df

    work = df.loc[designable_mask].copy()

    iterator = work.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(work), desc="Applying silent edits", leave=False)

    for row in iterator:
        idx = row.Index
        rowd = row._asdict()

        chrom = rowd["chromosome"]
        gene_strand = rowd["gene_strand"]
        feature = rowd["feature"]

        hal_gene = rowd.get("HAL_seq_gene")
        har_gene = rowd.get("HAR_seq_gene")

        if isinstance(hal_gene, str) and hal_gene:
            work.at[idx, "HAL_seq_mut"] = hal_gene
        if isinstance(har_gene, str) and har_gene:
            work.at[idx, "HAR_seq_mut"] = har_gene

        if bool(rowd.get("rejected", False)):
            continue

        target_arm = rowd.get("requires_edit_arm", "none")
        if target_arm in (None, "", "none"):
            continue

        if target_arm == "HAR":
            arm_start, arm_end = int(rowd["HARs"]), int(rowd["HARe"])
            arm_seq = rowd.get("HAR_seq_gene", "")
        else:
            arm_start, arm_end = int(rowd["HALs"]), int(rowd["HALe"])
            arm_seq = rowd.get("HAL_seq_gene", "")

        if not isinstance(arm_seq, str) or not arm_seq:
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                f"Missing {target_arm}_seq_gene; cannot edit",
            )
            continue

        grna_strand = rowd.get("grna_strand") or rowd.get("gRNA_strand")
        if grna_strand not in ("+", "-"):
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                "Missing/invalid grna_strand; cannot apply edits",
            )
            continue

        try:
            pam_s, pam_e = int(rowd["pam_start"]), int(rowd["pam_end"])
            prot_s, prot_e = int(rowd["protospacer_start"]), int(rowd["protospacer_end"])
        except Exception:
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                "Missing PAM/protospacer coordinates; cannot apply edits",
            )
            continue

        anchor = _frame_anchor_for_arm(
            feature=feature,
            target_arm=target_arm,
            arm_len=len(arm_seq),
        )

        hal_muts: list[dict] = []
        har_muts: list[dict] = []

        def append_mut(record: dict) -> None:
            if target_arm == "HAR":
                har_muts.append(record)
            else:
                hal_muts.append(record)

        did_edit = False

        # ---------- (1) Try single synonymous mutation in PAM GG only ----------
        # PAM in guide orientation occupies positions:
        #   21 = N
        #   22 = G
        #   23 = G
        #
        # Mutating position 21 does not break SpCas9 NGG recognition, so we only
        # attempt synonymous edits at guide-PAM indices 1 and 2.
        
        pam_seq_guide = str(rowd.get("pam_seq", "")).upper()

        if len(pam_seq_guide) == 3:
            pam_gg_guide_indices = [
                guide_idx
                for guide_idx, base in enumerate(pam_seq_guide)
                if guide_idx in (1, 2) and base == "G"
            ]

            for guide_idx in pam_gg_guide_indices:
                if grna_strand == "+":
                    gpos = pam_s + guide_idx
                else:
                    gpos = pam_e - guide_idx

                arm_idx = _genomic_to_coding_index(
                    gpos,
                    arm_start,
                    arm_end,
                    gene_strand,
                )

                if arm_idx is None:
                    continue

                attempt = _try_synonymous_change(
                    arm_seq,
                    arm_idx,
                    anchor=anchor,
                )

                if attempt is not None:
                    arm_seq, rec = attempt
                    append_mut(rec)
                    did_edit = True
                    break

        if not did_edit:
            prot_positions = list(range(prot_s, prot_e + 1))
            prot_idxs = _indices_within_arm(
                arm_start,
                arm_end,
                gene_strand,
                prot_positions,
            )

            if prot_idxs:
                pam_edge_pos = prot_e if grna_strand == "+" else prot_s
                pam_edge_idx = _genomic_to_coding_index(
                    pam_edge_pos,
                    arm_start,
                    arm_end,
                    gene_strand,
                )

                if pam_edge_idx is not None:
                    prot_idxs = sorted(prot_idxs, key=lambda i: (abs(i - pam_edge_idx), i))

                for prot_idx in prot_idxs:
                    attempt = _try_synonymous_change(arm_seq, prot_idx, anchor=anchor)
                    if attempt is not None:
                        arm_seq, rec = attempt
                        append_mut(rec)
                        did_edit = True
                        break

        if not did_edit:
            prot_positions = list(range(prot_s, prot_e + 1))
            prot_idxs = _indices_within_arm(
                arm_start,
                arm_end,
                gene_strand,
                prot_positions,
            )

            if prot_idxs:
                pam_edge_pos = prot_e if grna_strand == "+" else prot_s
                pam_edge_idx = _genomic_to_coding_index(
                    pam_edge_pos,
                    arm_start,
                    arm_end,
                    gene_strand,
                )

                if pam_edge_idx is not None:
                    prot_idxs = sorted(prot_idxs, key=lambda i: (abs(i - pam_edge_idx), i))

                attempt = _two_synonymous_changes_near(arm_seq, prot_idxs, anchor=anchor)

                if attempt is not None:
                    arm_seq, muts = attempt
                    for rec in muts:
                        append_mut(rec)
                    did_edit = True

        work.at[idx, "edit_arm"] = target_arm

        if target_arm == "HAR":
            work.at[idx, "HAR_seq_mut"] = arm_seq
        else:
            work.at[idx, "HAL_seq_mut"] = arm_seq

        work.at[idx, "HAL_mutation_location"] = _fmt_mut_locations(
            chrom,
            gene_strand,
            int(rowd["HALs"]),
            int(rowd["HALe"]),
            hal_muts,
        )
        work.at[idx, "HAR_mutation_location"] = _fmt_mut_locations(
            chrom,
            gene_strand,
            int(rowd["HARs"]),
            int(rowd["HARe"]),
            har_muts,
        )

        if not did_edit:
            work.at[idx, "warnings"] = merge_warn(
                work.at[idx, "warnings"],
                f"No synonymous mutation could be applied in target {target_arm} arm",
            )

    out = df.copy()
    out.update(work)
    return out