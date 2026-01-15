"""
Design utilities for TFTag:
- Build 240 bp homology arms (HAL/HAR) around start/stop codons.
- Apply *silent* edits to prevent re-cutting by the on-target gRNA:
  1) Prefer a single synonymous change that breaks the PAM (rank 1).
  2) Otherwise, one synonymous change within the protospacer (rank 2).
  3) Otherwise, two synonymous changes within the protospacer near the PAM edge (rank 3).

Conventions
-----------
- Genomic coordinates are 1-based inclusive.
- Homology arm sequences (HAL_seq_gene/HAR_seq_gene) are returned in *gene orientation*
  (i.e. reverse-complemented for genes on the '-' strand). This keeps codon indexing
  consistent across strands.

Expected input schema (post-scan)
---------------------------------
Required columns for most functions:
  - chromosome
  - feature: 'start_codon' | 'stop_codon'
  - gene_strand: '+' | '-'
  - codon_start, codon_end  (1-based inclusive)
  - designable (bool)  [added by prefilter_designable]

For editing functions, you also need:
  - grna_strand: '+' | '-'
  - protospacer_start, protospacer_end (1-based inclusive)
  - pam_start, pam_end (1-based inclusive)

This module relies on coords.py for PAM/protospacer coordinate access:
  from .coords import pam_coords, protospacer_coords
"""
from __future__ import annotations

from typing import Tuple, Optional, List, Dict

import numpy as np
import pandas as pd
import primer3
from Bio.Seq import Seq
from tqdm import tqdm

from .utils import get_sequence, clamp_range
from .coords import pam_coords, protospacer_coords


# -----------------------------
# Arm coordinate helpers
# -----------------------------
def _arm_ranges(feature: str, strand: str, codon_start: int, codon_end: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """
    Return genomic ranges for HAL and HAR (240 bp each) around a start/stop codon.
    Ranges are in genomic coordinates (1-based inclusive).

    HAL is the "left" arm in *gene orientation*; HAR is the "right" arm.
    """
    if feature not in ("start_codon", "stop_codon"):
        raise ValueError("feature must be 'start_codon' or 'stop_codon'")
    if strand == "+":
        HAL = (codon_start - 240, codon_start - 1)
        HAR = (codon_end + 1, codon_end + 240)
    elif strand == "-":
        HAL = (codon_end + 1, codon_end + 240)
        HAR = (codon_start - 240, codon_start - 1)
    else:
        raise ValueError("strand must be '+' or '-'")
    return HAL, HAR


def _genomic_to_coding_index(pos: int, arm_start: int, arm_end: int, gene_strand: str) -> Optional[int]:
    """Map a genomic position to 0-based index within a gene-oriented arm; return None if outside."""
    if pos < arm_start or pos > arm_end:
        return None
    if gene_strand == "+":
        return pos - arm_start
    else:
        return arm_end - pos


def _coding_index_to_genomic(idx: int, arm_start: int, arm_end: int, gene_strand: str) -> int:
    """Inverse of _genomic_to_coding_index."""
    if gene_strand == "+":
        return arm_start + idx
    else:
        return arm_end - idx


def _indices_within_arm(arm_start: int, arm_end: int, gene_strand: str, positions: List[int]) -> List[int]:
    """Return gene-oriented indices of genomic positions that fall within the arm interval."""
    out: List[int] = []
    for pos in positions:
        idx = _genomic_to_coding_index(pos, arm_start, arm_end, gene_strand)
        if idx is not None:
            out.append(idx)
    return out


# -----------------------------
# Warnings helper
# -----------------------------
def _merge_warn(prev, new) -> str:
    """
    Merge warnings into a semicolon-separated unique list.
    `new` may be a list[str] or a single str.
    """
    if isinstance(new, str):
        new_list = [new] if new else []
    else:
        new_list = [w for w in new if w]

    prev_list = [] if prev in (None, "", "none") else [w.strip() for w in str(prev).split(";") if w.strip()]
    merged = sorted(set(prev_list + new_list))
    return "; ".join(merged) if merged else "none"


# -----------------------------
# Synonymous edit helpers
# -----------------------------
def _translate_codon(codon: str) -> str:
    return str(Seq(codon).translate(table=1))


def _try_synonymous_change(coding_seq: str, idx: int) -> Optional[Tuple[str, Dict]]:
    """Attempt a single-nucleotide synonymous change at coding index `idx`."""
    bases = ["A", "C", "G", "T"]
    if idx < 0 or idx >= len(coding_seq) or coding_seq[idx] not in bases:
        return None

    codon_start = (idx // 3) * 3
    codon_end = codon_start + 3
    if codon_end > len(coding_seq):
        return None

    codon = list(coding_seq[codon_start:codon_end])
    aa_orig = _translate_codon("".join(codon))

    for b in bases:
        if b == coding_seq[idx]:
            continue
        codon_mut = codon.copy()
        codon_mut[idx - codon_start] = b
        if _translate_codon("".join(codon_mut)) == aa_orig:
            cs = list(coding_seq)
            cs[idx] = b
            return (
                "".join(cs),
                {
                    "idx": idx,
                    "from": coding_seq[idx],
                    "to": b,
                    "codon_from": "".join(codon),
                    "codon_to": "".join(codon_mut),
                    "codon_start": codon_start,
                },
            )
    return None


def _two_synonymous_changes_near(coding_seq: str, candidate_indices: List[int]) -> Optional[Tuple[str, List[Dict]]]:
    """
    Try to apply two synonymous changes using the provided candidate indices (in order).
    Returns (mutated_seq, [mut1, mut2]) or None.
    """
    seq = coding_seq
    muts: List[Dict] = []
    for idx in candidate_indices:
        attempt = _try_synonymous_change(seq, idx)
        if attempt is None:
            continue
        seq, rec = attempt
        muts.append(rec)
        if len(muts) == 2:
            return seq, muts
    return None


def _fmt_mut_locations(chrom: str, gene_strand: str, arm_start: int, arm_end: int, muts: List[Dict]) -> str:
    """Return comma-separated genomic edits like A-3L:25668456-G (uses genomic base letters)."""
    if not muts:
        return "none"

    def comp(b: str) -> str:
        table = str.maketrans("ACGT", "TGCA")
        return b.translate(table)

    items = []
    for m in muts:
        pos = _coding_index_to_genomic(m["idx"], arm_start, arm_end, gene_strand)
        # m['from']/['to'] are coding-oriented; convert to genomic letters on minus
        b_from = m["from"] if gene_strand == "+" else comp(m["from"])
        b_to = m["to"] if gene_strand == "+" else comp(m["to"])
        items.append(f"{b_from}-{chrom}:{pos}-{b_to}")
    return ", ".join(items)


# -----------------------------
# Prefilter
# -----------------------------
def prefilter_designable(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    *,
    hal_len: int = 240,
    har_len: int = 240,
    primer_upstream_window: tuple[int, int] = (100, 300),
    primer_downstream_window: tuple[int, int] = (100, 300),
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Mark rows as designable only if BOTH homology arms and BOTH validation-primer windows
    lie fully within the chromosome (no truncation). Avoids later clamping.

    Requires:
      chromosome, feature, gene_strand, codon_start, codon_end

    Adds:
      designable: bool
      skip_reason: "none" or semicolon-separated reasons
    """
    df = gRNA_df.copy()
    if "skip_reason" not in df.columns:
        df["skip_reason"] = "none"
    df["designable"] = True

    required = ["chromosome", "feature", "gene_strand", "codon_start", "codon_end"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"prefilter_designable missing required columns: {missing}")

    iterator = df.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Prefiltering designable loci", leave=False)

    def in_bounds(s: int, e: int, clen: int) -> bool:
        return (s >= 1) and (e <= clen) and (e >= s)

    for idx, row in iterator:
        chrom = row.chromosome
        feature = row.feature
        strand = row.gene_strand
        cs = int(row.codon_start)
        ce = int(row.codon_end)
        clen = len(fasta_dict[chrom].seq)

        (HALs, HALe), (HARs, HARe) = _arm_ranges(feature, strand, cs, ce)

        up_lo, up_hi = primer_upstream_window
        dn_lo, dn_hi = primer_downstream_window

        if strand == "+":
            up_start, up_end = HALs - up_hi, HALs - up_lo
            dn_start, dn_end = HARe + dn_lo, HARe + dn_hi
        else:
            up_start, up_end = HALe + up_lo, HALe + up_hi
            dn_start, dn_end = HARs - dn_hi, HARs - dn_lo

        reasons = []
        if not in_bounds(HALs, HALe, clen):
            reasons.append("HAL out of contig")
        if not in_bounds(HARs, HARe, clen):
            reasons.append("HAR out of contig")
        if (HALe - HALs + 1) < hal_len:
            reasons.append(f"HAL < {hal_len} bp")
        if (HARe - HARs + 1) < har_len:
            reasons.append(f"HAR < {har_len} bp")
        if not in_bounds(up_start, up_end, clen):
            reasons.append("Upstream primer window out of contig")
        if not in_bounds(dn_start, dn_end, clen):
            reasons.append("Downstream primer window out of contig")

        if reasons:
            df.at[idx, "designable"] = False
            prev = [] if df.at[idx, "skip_reason"] in ("", "none") else [
                x.strip() for x in str(df.at[idx, "skip_reason"]).split(";") if x.strip()
            ]
            df.at[idx, "skip_reason"] = "; ".join(sorted(set(prev + reasons)))

    return df


# -----------------------------
# Homology arms
# -----------------------------
def add_homology_arms(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Compute HAL/HAR genomic ranges (clamped) and gene-oriented sequences.

    Requires:
      chromosome, feature, gene_strand, codon_start, codon_end, designable

    Adds:
      HALs, HALe, HARs, HARe,
      HAL_seq_gene, HAR_seq_gene,
      warnings
    """
    df = gRNA_df.copy()
    if "warnings" not in df.columns:
        df["warnings"] = "none"

    required = ["chromosome", "feature", "gene_strand", "codon_start", "codon_end", "designable"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"add_homology_arms missing required columns: {missing}")

    if not df["designable"].any():
        return df

    work = df[df["designable"]].copy()
    iterator = work.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(work), desc="Computing homology arms", leave=False)

    contig_len_cache: dict[str, int] = {}

    for idx, row in iterator:
        chrom = row.chromosome
        feature = row.feature
        strand = row.gene_strand
        cs = int(row.codon_start)
        ce = int(row.codon_end)

        (HALs, HALe), (HARs, HARe) = _arm_ranges(feature, strand, cs, ce)
        clen = contig_len_cache.setdefault(chrom, len(fasta_dict[chrom].seq))

        # clamp as a safety net; designable prefilter should already guarantee in-bounds
        HALs, HALe = clamp_range(HALs, HALe, clen)
        HARs, HARe = clamp_range(HARs, HARe, clen)

        work.at[idx, "HALs"] = HALs
        work.at[idx, "HALe"] = HALe
        work.at[idx, "HARs"] = HARs
        work.at[idx, "HARe"] = HARe

        HAL_seq_gene = get_sequence(fasta_dict, chrom, HALs, HALe, strand)
        HAR_seq_gene = get_sequence(fasta_dict, chrom, HARs, HARe, strand)
        work.at[idx, "HAL_seq_gene"] = HAL_seq_gene
        work.at[idx, "HAR_seq_gene"] = HAR_seq_gene

        warns = []
        if len(HAL_seq_gene) < 240:
            warns.append(f"HAL truncated to {len(HAL_seq_gene)} bp (chromosome end)")
        if len(HAR_seq_gene) < 240:
            warns.append(f"HAR truncated to {len(HAR_seq_gene)} bp (chromosome end)")
        if warns:
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], warns)

    out = df.copy()
    out.update(work)
    return out


# -----------------------------
# Choose edit arm
# -----------------------------
def choose_arm_for_mutation(
    gRNA_df: pd.DataFrame,
    protospacer_overlap_len: int = 13,
    coding_only: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Determine which arm (HAL/HAR) needs modification based on containment of:
      PAM + protospacer_overlap_len PAM-proximal bases of protospacer.

    Requires:
      designable, HALs/HALe/HARs/HARe present (from add_homology_arms),
      plus PAM/protospacer coords and grna_strand.

    Adds/updates:
      requires_edit_arm: 'HAL'|'HAR'|'none'
      rejected: bool
      reject_reason: text
      warnings (if missing coords)
    """
    df = gRNA_df.copy()
    if "warnings" not in df.columns:
        df["warnings"] = "none"

    for col, default in [("requires_edit_arm", "none"), ("rejected", False), ("reject_reason", "none")]:
        if col not in df.columns:
            df[col] = default

    if "designable" in df.columns and not df["designable"].any():
        return df

    work = df[df.get("designable", True)].copy()
    iterator = work.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(work), desc="Choosing edit arm", leave=False)

    for idx, row in iterator:
        if pd.isna(getattr(row, "HALs", np.nan)) or pd.isna(getattr(row, "HARs", np.nan)):
            work.at[idx, "requires_edit_arm"] = "none"
            continue

        HALs, HALe = int(row.HALs), int(row.HALe)
        HARs, HARe = int(row.HARs), int(row.HARe)

        rowd = row._asdict()
        pam_s, pam_e = pam_coords(rowd)
        prot_s, prot_e = protospacer_coords(rowd)
        if pam_s is None or prot_s is None:
            work.at[idx, "requires_edit_arm"] = "none"
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], ["Missing PAM/protospacer coords; cannot choose edit arm"])
            continue

        grna_strand = rowd.get("grna_strand") or rowd.get("gRNA_strand")
        if grna_strand not in ("+", "-"):
            work.at[idx, "requires_edit_arm"] = "none"
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], ["Missing/invalid grna_strand; cannot choose edit arm"])
            continue

        if grna_strand == "+":
            window_s = max(prot_e - (protospacer_overlap_len - 1), prot_s)
            window_e = pam_e
        else:
            window_s = pam_s
            window_e = min(prot_s + (protospacer_overlap_len - 1), prot_e)

        def contains(a_s, a_e, w_s, w_e):
            return a_s <= w_s and w_e <= a_e

        in_HAL = contains(HALs, HALe, window_s, window_e)
        in_HAR = contains(HARs, HARe, window_s, window_e)

        # HAL and HAR should be disjoint; if both true something is wrong
        if in_HAL and in_HAR:
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], ["Internal error: edit window contained in both arms"])
            target_arm = "none"
        else:
            target_arm = "HAL" if in_HAL else ("HAR" if in_HAR else "none")

        work.at[idx, "requires_edit_arm"] = target_arm

        if coding_only and target_arm != "none":
            required_arm = "HAR" if row.feature == "start_codon" else "HAL"
            if target_arm != required_arm:
                work.at[idx, "rejected"] = True
                work.at[idx, "reject_reason"] = f"requires edits in non-coding arm {target_arm} for {row.feature}; rejected"

    out = df.copy()
    out.update(work)
    return out


# -----------------------------
# Apply silent edits
# -----------------------------
def apply_silent_edits(
    gRNA_df: pd.DataFrame,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Apply PAM/protospacer silent edits to the required coding arm when needed & allowed.

    Requires (for rows you want to edit):
      - HALs/HALe/HARs/HARe
      - HAL_seq_gene/HAR_seq_gene
      - gene_strand
      - grna_strand
      - protospacer/pam coords
      - requires_edit_arm (from choose_arm_for_mutation)
      - rejected (optional)

    Adds/updates:
      - HAL_seq_mut, HAR_seq_mut
      - HAL mutation location, HAR mutation location
      - edit_arm
      - warnings (appends if no synonymous edits possible)
    """
    df = gRNA_df.copy()
    if "warnings" not in df.columns:
        df["warnings"] = "none"

    # Ensure columns exist
    for col, default in [
        ("HAL_seq_mut", ""),
        ("HAR_seq_mut", ""),
        ("HAL mutation location", "none"),
        ("HAR mutation location", "none"),
        ("edit_arm", "none"),
    ]:
        if col not in df.columns:
            df[col] = default

    if "designable" in df.columns and not df["designable"].any():
        return df

    work = df[df.get("designable", True)].copy()
    iterator = work.itertuples(index=True)
    if show_progress:
        iterator = tqdm(iterator, total=len(work), desc="Applying silent edits", leave=False)

    for idx, row in iterator:
        rowd = row._asdict()
        chrom = rowd["chromosome"]
        gene_strand = rowd["gene_strand"]

        # Default: copy sequences through
        HAL_gene = rowd.get("HAL_seq_gene")
        HAR_gene = rowd.get("HAR_seq_gene")

        if HAL_gene is not None:
            work.at[idx, "HAL_seq_mut"] = HAL_gene
        if HAR_gene is not None:
            work.at[idx, "HAR_seq_mut"] = HAR_gene

        # If missing prerequisites or rejected, keep originals
        if rowd.get("rejected", False):
            continue

        target_arm = rowd.get("requires_edit_arm", "none")
        if target_arm in (None, "", "none"):
            continue

        # Arm sequence and bounds in gene orientation
        if target_arm == "HAR":
            arm_start, arm_end = int(rowd["HARs"]), int(rowd["HARe"])
            arm_seq = rowd["HAR_seq_gene"]
        else:
            arm_start, arm_end = int(rowd["HALs"]), int(rowd["HALe"])
            arm_seq = rowd["HAL_seq_gene"]

        if not isinstance(arm_seq, str) or not arm_seq:
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], [f"Missing {target_arm}_seq_gene; cannot edit"])
            continue

        grna_strand = rowd.get("grna_strand") or rowd.get("gRNA_strand")
        if grna_strand not in ("+", "-"):
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], ["Missing/invalid grna_strand; cannot apply edits"])
            continue

        pam_s, pam_e = pam_coords(rowd)
        prot_s, prot_e = protospacer_coords(rowd)
        if pam_s is None or prot_s is None:
            work.at[idx, "warnings"] = _merge_warn(work.at[idx, "warnings"], ["Missing PAM/protospacer coords; cannot apply edits"])
            continue

        hal_muts: List[Dict] = []
        har_muts: List[Dict] = []
        did_single_pam = False
        did_single_prot = False
        did_double = False

        # (1) Single synonymous edit in PAM (prefer breaking a G in NGG if possible).
        # Note: We do NOT explicitly validate coding frame here; we simply ensure synonymous at AA level.
        pam_positions = list(range(int(pam_s), int(pam_e) + 1))
        pam_idxs = _indices_within_arm(arm_start, arm_end, gene_strand, pam_positions)

        # If PAM fully within arm: length 3
        if len(pam_idxs) == 3:
            # Heuristic: mutate PAM-proximal "GG" in the guide orientation.
            # When grna_strand == gene_strand, gene-oriented sequence matches guide target strand orientation,
            # so GG correspond to indices [1,2] within the PAM triplet in gene orientation.
            # Otherwise, attempt [0,1] first.
            if grna_strand == gene_strand:
                candidate_pam_indices = [pam_idxs[1], pam_idxs[2]]
            else:
                candidate_pam_indices = [pam_idxs[0], pam_idxs[1]]

            for pi in candidate_pam_indices:
                attempt = _try_synonymous_change(arm_seq, pi)
                if attempt is not None:
                    arm_seq, rec = attempt
                    (har_muts if target_arm == "HAR" else hal_muts).append(rec)
                    did_single_pam = True
                    break

        # (2) Single synonymous edit in protospacer (try nearest to PAM edge first if possible)
        if not did_single_pam:
            prot_positions = list(range(int(prot_s), int(prot_e) + 1))
            prot_idxs = _indices_within_arm(arm_start, arm_end, gene_strand, prot_positions)

            if prot_idxs:
                # Try indices closest to PAM edge (in gene-oriented indexing)
                pam_edge_pos = int(pam_s) if grna_strand == gene_strand else int(pam_e)
                pam_edge_idx = _genomic_to_coding_index(pam_edge_pos, arm_start, arm_end, gene_strand)
                if pam_edge_idx is not None:
                    prot_idxs = sorted(prot_idxs, key=lambda i: (abs(i - pam_edge_idx), i))

                for pi in prot_idxs:
                    attempt = _try_synonymous_change(arm_seq, pi)
                    if attempt is not None:
                        arm_seq, rec = attempt
                        (har_muts if target_arm == "HAR" else hal_muts).append(rec)
                        did_single_prot = True
                        break

        # (3) Two synonymous edits in protospacer (near PAM if possible)
        if not (did_single_pam or did_single_prot):
            prot_positions = list(range(int(prot_s), int(prot_e) + 1))
            prot_idxs = _indices_within_arm(arm_start, arm_end, gene_strand, prot_positions)
            if prot_idxs:
                pam_edge_pos = int(pam_s) if grna_strand == gene_strand else int(pam_e)
                pam_edge_idx = _genomic_to_coding_index(pam_edge_pos, arm_start, arm_end, gene_strand)
                if pam_edge_idx is not None:
                    prot_idxs = sorted(prot_idxs, key=lambda i: (abs(i - pam_edge_idx), i))

                attempt2 = _two_synonymous_changes_near(arm_seq, prot_idxs)
                if attempt2 is not None:
                    arm_seq, muts = attempt2
                    (har_muts if target_arm == "HAR" else hal_muts).extend(muts)
                    did_double = True

        # Commit mutated sequences
        work.at[idx, "edit_arm"] = target_arm
        if target_arm == "HAR":
            work.at[idx, "HAR_seq_mut"] = arm_seq
            work.at[idx, "HAL_seq_mut"] = rowd.get("HAL_seq_gene", work.at[idx, "HAL_seq_mut"])
        else:
            work.at[idx, "HAL_seq_mut"] = arm_seq
            work.at[idx, "HAR_seq_mut"] = rowd.get("HAR_seq_gene", work.at[idx, "HAR_seq_mut"])

        # Record locations
        work.at[idx, "HAL mutation location"] = _fmt_mut_locations(
            chrom, gene_strand, int(rowd["HALs"]), int(rowd["HALe"]), hal_muts
        )
        work.at[idx, "HAR mutation location"] = _fmt_mut_locations(
            chrom, gene_strand, int(rowd["HARs"]), int(rowd["HARe"]), har_muts
        )

        if not (did_single_pam or did_single_prot or did_double):
            work.at[idx, "warnings"] = _merge_warn(
                work.at[idx, "warnings"],
                [f"No synonymous mutation could be applied in target {target_arm} arm (kept original sequence)"],
            )

    out = df.copy()
    out.update(work)
    return out


# -----------------------------
# Primer helpers
# -----------------------------
def _count_trailing_gc(primer: str) -> int:
    """Count consecutive G/C bases from the 3' end of a primer."""
    c = 0
    for b in reversed(primer.upper()):
        if b in ("G", "C"):
            c += 1
        else:
            break
    return c


def _pick_single_primer(
    seq_template: str,
    which: str,
    rounds: list,
    max_end_gc: int = 3,
    show_reason: bool = False,
):
    """
    Try multiple parameter rounds to pick a single primer with primer3.

    which: 'left' (forward on template) or 'right' (reverse on template).

    Enforces:
      - Primer3 GC clamp (consecutive G/C at 3' end) via PRIMER_GC_CLAMP
      - Max consecutive G/C at 3' end via post-filter (<= max_end_gc)
    """
    which = which.lower()
    assert which in ("left", "right")

    last_res = None

    for ridx, rd in enumerate(rounds, start=1):
        p3_args = {
            "PRIMER_OPT_SIZE": rd["size_opt"],
            "PRIMER_MIN_SIZE": rd["size_min"],
            "PRIMER_MAX_SIZE": rd["size_max"],
            "PRIMER_OPT_TM": rd["tm_opt"],
            "PRIMER_MIN_TM": rd["tm_min"],
            "PRIMER_MAX_TM": rd["tm_max"],
            "PRIMER_MIN_GC": rd["gc_min"],
            "PRIMER_MAX_GC": rd["gc_max"],
            "PRIMER_MAX_HAIRPIN_TH": rd["hairpin_max"],
            "PRIMER_MAX_POLY_X": rd["poly_x_max"],
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_NUM_RETURN": 20,
            "PRIMER_EXPLAIN_FLAG": 1,
            "PRIMER_PICK_LEFT_PRIMER": 1 if which == "left" else 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1 if which == "right" else 0,
            "PRIMER_GC_CLAMP": rd["clamp_min"],
        }

        res = primer3.bindings.designPrimers(
            {"SEQUENCE_ID": "template", "SEQUENCE_TEMPLATE": seq_template},
            p3_args,
        )
        last_res = res

        cand = []
        key_prefix = "PRIMER_LEFT" if which == "left" else "PRIMER_RIGHT"
        nret = int(res.get(f"{key_prefix}_NUM_RETURNED", 0) or 0)

        for i in range(nret):
            seq = res[f"{key_prefix}_{i}_SEQUENCE"]
            tm = float(res.get(f"{key_prefix}_{i}_TM", float("nan")))
            gc = float(res.get(f"{key_prefix}_{i}_GC_PERCENT", float("nan")))
            ln = len(seq)

            run_gc = _count_trailing_gc(seq)
            if run_gc < rd["clamp_min"]:
                continue
            if run_gc > max_end_gc:
                continue

            cand.append((seq, tm, gc, ln))

        if cand:
            seq, tm, gc, ln = cand[0]  # primer3 sorts by objective
            note = f"round{ridx} ({which})"
            return seq, tm, gc, ln, note

    return None, None, None, None, ("no primer found" if not show_reason else str(last_res))


# -----------------------------
# Validation primers
# -----------------------------
def design_validation_primers(
    gRNA_df: pd.DataFrame,
    fasta_dict,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Design two gene-specific validation primers per row:
      PCR1: Forward primer 100–300 bp UPSTREAM of HAL (gene-oriented).
      PCR2: Reverse primer 100–300 bp DOWNSTREAM of HAR (gene-oriented).

    Requires:
      chromosome, feature, gene_strand, codon_start, codon_end, designable

    Adds:
      PCR1_forward_primer_seq / _Tm / _GC / _len
      PCR2_reverse_primer_seq / _Tm / _GC / _len
      warnings (appended)
    """
    rounds = [
        dict(tm_opt=60.0, tm_min=58.0, tm_max=62.0, gc_min=40.0, gc_max=60.0, size_opt=20, size_min=18, size_max=22, clamp_min=1, hairpin_max=24.0, poly_x_max=4),
        dict(tm_opt=60.0, tm_min=57.0, tm_max=63.0, gc_min=35.0, gc_max=65.0, size_opt=20, size_min=18, size_max=23, clamp_min=1, hairpin_max=48.0, poly_x_max=5),
        dict(tm_opt=60.0, tm_min=55.0, tm_max=65.0, gc_min=30.0, gc_max=70.0, size_opt=20, size_min=17, size_max=28, clamp_min=1, hairpin_max=48.0, poly_x_max=6),
    ]

    df = gRNA_df.copy()
    if "warnings" not in df.columns:
        df["warnings"] = "none"

    required = ["chromosome", "feature", "gene_strand", "codon_start", "codon_end", "designable"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"design_validation_primers missing required columns: {missing}")

    # Prepare output cols
    df["PCR1_forward_primer_seq"] = ""
    df["PCR1_forward_primer_Tm"] = np.nan
    df["PCR1_forward_primer_GC"] = np.nan
    df["PCR1_forward_primer_len"] = np.nan

    df["PCR2_reverse_primer_seq"] = ""
    df["PCR2_reverse_primer_Tm"] = np.nan
    df["PCR2_reverse_primer_GC"] = np.nan
    df["PCR2_reverse_primer_len"] = np.nan

    # Mark skipped rows
    mask_bad = ~df["designable"]
    df.loc[mask_bad, "warnings"] = df.loc[mask_bad, "warnings"].apply(lambda w: _merge_warn(w, ["skipped: prefilter_designable"]))

    iterator = df.iterrows()
    if show_progress:
        iterator = tqdm(iterator, total=len(df), desc="Designing validation primers", leave=False)

    for idx, row in iterator:
        if not bool(row["designable"]):
            continue

        chrom = row["chromosome"]
        feature = row["feature"]
        gene_strand = row["gene_strand"]
        cs = int(row["codon_start"])
        ce = int(row["codon_end"])
        clen = len(fasta_dict[chrom].seq)

        warn_msgs = [] if df.at[idx, "warnings"] in ("", "none") else [
            w.strip() for w in str(df.at[idx, "warnings"]).split(";") if w.strip()
        ]

        (HALs, HALe), (HARs, HARe) = _arm_ranges(feature, gene_strand, cs, ce)

        # Windows around arms (genomic coords)
        if gene_strand == "+":
            up_start, up_end = HALs - 300, HALs - 100
            dn_start, dn_end = HARe + 100, HARe + 300
        else:
            up_start, up_end = HALe + 100, HALe + 300
            dn_start, dn_end = HARs - 300, HARs - 100

        # prefilter_designable should guarantee these are in bounds
        if up_start < 1 or up_end > clen or dn_start < 1 or dn_end > clen or up_end < up_start or dn_end < dn_start:
            warn_msgs.append("Primer window out of contig despite designable=True (regression)")
            df.at[idx, "warnings"] = "; ".join(sorted(set(warn_msgs)))
            continue

        upstream_seq = get_sequence(fasta_dict, chrom, up_start, up_end, gene_strand)
        downstream_seq = get_sequence(fasta_dict, chrom, dn_start, dn_end, gene_strand)

        f_seq, f_tm, f_gc, f_len, _ = _pick_single_primer(upstream_seq, which="left", rounds=rounds, max_end_gc=3)
        if f_seq is None:
            warn_msgs.append("No PCR1 forward primer found after relaxation")
        else:
            df.at[idx, "PCR1_forward_primer_seq"] = f_seq
            df.at[idx, "PCR1_forward_primer_Tm"] = f_tm
            df.at[idx, "PCR1_forward_primer_GC"] = f_gc
            df.at[idx, "PCR1_forward_primer_len"] = f_len

        r_seq, r_tm, r_gc, r_len, _ = _pick_single_primer(downstream_seq, which="right", rounds=rounds, max_end_gc=3)
        if r_seq is None:
            warn_msgs.append("No PCR2 reverse primer found after relaxation")
        else:
            df.at[idx, "PCR2_reverse_primer_seq"] = r_seq
            df.at[idx, "PCR2_reverse_primer_Tm"] = r_tm
            df.at[idx, "PCR2_reverse_primer_GC"] = r_gc
            df.at[idx, "PCR2_reverse_primer_len"] = r_len

        df.at[idx, "warnings"] = "; ".join(sorted(set(warn_msgs))) if warn_msgs else "none"

    return df


def validation_primers(gRNA_df: pd.DataFrame, fasta_dict, show_progress: bool = True) -> pd.DataFrame:
    """Alias for pipeline compatibility."""
    return design_validation_primers(gRNA_df, fasta_dict, show_progress=show_progress)