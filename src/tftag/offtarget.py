"""
Cas-OFFinder I/O and parsing utilities.

Key invariant
-------------
Cas-OFFinder requires that every query sequence has exactly the same length
as the pattern sequence provided on line 2 of the input file.

For SpCas9 NGG searches, you typically run Cas-OFFinder on 23-mers:
  [20-nt protospacer][3-nt PAM]
with the pattern:
  NNNNNNNNNNNNNNNNNNNNNGG   (length 23)

This module therefore defaults to 23-mer queries and enforces length matching.
"""
from __future__ import annotations

import os
import re
import subprocess
from typing import Optional

import pandas as pd


# Matches either "... <strand> <mm>" or "... <strand><mm>" at end of line, e.g. "-0" or "+ 4"
_STRAND_MM_RE = re.compile(r"([+-])\s*([0-9]+)\s*$")


def write_cas_offinder_input(
    guides: pd.DataFrame,
    genome_fasta_path: str,
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
    outfile: str = "cas_offinder_input.txt",
    mismatches: int = 4,
    *,
    query_col: Optional[str] = None,
) -> str:
    """
    Write a Cas-OFFinder input file.

    Parameters
    ----------
    guides:
      Candidate guide DataFrame. Must contain either:
        - `grna_seq_23` (preferred), or
        - `spacer` + `pam_seq`
      unless `query_col` is explicitly set.

    pam_pattern:
      Pattern line written to Cas-OFFinder input. Length defines query length.

    query_col:
      If provided, use this column directly as the Cas-OFFinder query sequence column.

    Notes
    -----
    Cas-OFFinder requires: len(query_seq) == len(pam_pattern) for every query line.
    We enforce this and raise an informative error otherwise.
    """
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

    pat_len = len(pam_pattern)

    # Decide which column(s) provide the query sequences.
    if query_col is not None:
        if query_col not in guides.columns:
            raise KeyError(f"query_col='{query_col}' not found in guides columns.")
        seqs = guides[query_col].astype("string").str.upper()
    else:
        # Default behaviour: prefer explicit 23-mer if available, else build from spacer+pam_seq.
        if "grna_seq_23" in guides.columns:
            seqs = guides["grna_seq_23"].astype("string").str.upper()
        elif {"spacer", "pam_seq"}.issubset(guides.columns):
            seqs = (
                guides["spacer"].astype("string").str.upper()
                + guides["pam_seq"].astype("string").str.upper()
            )
        else:
            raise KeyError(
                "Need either 'grna_seq_23' or ('spacer' and 'pam_seq') to build Cas-OFFinder queries "
                "(or provide query_col=...)."
            )

    # Drop missing values *before* length checks to avoid 'nan' string artefacts.
    seqs = seqs.dropna()

    # Track observed lengths for better diagnostics if we fail length enforcement.
    lengths = seqs.str.len().value_counts(dropna=False).to_dict()

    # Enforce required length.
    seqs = seqs[seqs.str.len() == pat_len]
    if seqs.empty:
        raise ValueError(
            f"No valid query sequences of length {pat_len} were produced. "
            f"Observed sequence lengths: {lengths}. "
            f"Check pam_pattern and query sequence columns."
        )

    # Unique + sorted queries (Series.unique avoids extra pandas objects).
    queries = pd.unique(seqs).tolist()
    queries = [str(q) for q in queries if q]  # ensure pure strings, drop empty
    queries.sort()

    with open(outfile, "w") as f:
        f.write(genome_fasta_path + "\n")
        f.write(pam_pattern + "\n")
        mm = int(mismatches)
        for q in queries:
            f.write(f"{q} {mm}\n")

    return outfile


def run_cas_offinder(
    input_file: str,
    *,
    cas_offinder_bin: str = "cas-offinder",
    device_spec: str = "C",
    output_file: str = "cas_offinder_hits.txt",
) -> str:
    """
    Run Cas-OFFinder.

    Raises a RuntimeError with captured stdout/stderr on failure, which makes
    debugging installation/pattern/query issues much easier from CLI.
    """
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)

    cmd = [cas_offinder_bin, input_file, device_spec, output_file]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "Cas-OFFinder failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDERR: {e.stderr[-2000:] if e.stderr else '(none)'}\n"
            f"STDOUT: {e.stdout[-2000:] if e.stdout else '(none)'}"
        ) from e

    return output_file


def parse_cas_offinder_output(hit_file: str) -> pd.DataFrame:
    """
    Parse Cas-OFFinder output robustly across formats.

    Handles both:
      ... <strand> <mm>
      ... <strand><mm>   (e.g. -0, +4)

    Also handles chromosome descriptors that contain whitespace/metadata. We keep
    only the first token as the contig/chrom name (e.g. "2L", "3R", "4").

    Returns
    -------
    DataFrame with columns:
      spacer (20 nt),
      offtarget_target,
      offtarget_chrom,
      offtarget_pos,
      offtarget_strand,
      offtarget_mismatches
    """
    rows: list[dict] = []

    def _is_int(tok: str) -> bool:
        try:
            int(tok)
            return True
        except Exception:
            return False

    with open(hit_file, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            # 1) Extract strand/mismatch count from end of line.
            m = _STRAND_MM_RE.search(line)
            if not m:
                continue
            strand = m.group(1)
            mismatches = int(m.group(2))

            # 2) Strip trailing strand/mm and parse remaining tokens.
            core = line[: m.start()].strip()
            toks = core.split()
            if len(toks) < 3:
                continue

            query = toks[0].upper()
            target = toks[-1].upper()

            # 3) Genomic position is the last integer token before the target.
            pos_idx = None
            for i in range(len(toks) - 2, 0, -1):
                if _is_int(toks[i]):
                    pos_idx = i
                    break
            if pos_idx is None:
                continue
            pos = int(toks[pos_idx])

            # 4) Chrom descriptor spans toks[1:pos_idx]; contig is the first token.
            chrom_desc = " ".join(toks[1:pos_idx])
            chrom = chrom_desc.split()[0] if chrom_desc else ""

            rows.append(
                {
                    # Normalise key to 20-nt spacer for downstream merges.
                    "spacer": query[:20],
                    "offtarget_target": target,
                    "offtarget_chrom": chrom,
                    "offtarget_pos": pos,
                    "offtarget_strand": strand,
                    "offtarget_mismatches": mismatches,
                }
            )

    return pd.DataFrame(rows)


def summarize_specificity(offtargets: pd.DataFrame, *, max_mismatches: int = 4) -> pd.DataFrame:
    """
    Summarise off-target counts per spacer.

    Returns columns:
      spacer, n_hits, n_mm0..n_mm{max_mismatches}

    Notes
    -----
    - Cas-OFFinder reports mismatches over the full query length (e.g. 23-mer if you queried 23).
    - n_hits is the sum across 0..max_mismatches bins.
    """
    cols = ["spacer", "n_hits"] + [f"n_mm{k}" for k in range(max_mismatches + 1)]

    if offtargets.empty:
        return pd.DataFrame(columns=cols)

    agg = (
        offtargets.groupby(["spacer", "offtarget_mismatches"])
        .size()
        .unstack(fill_value=0)
    )

    for k in range(max_mismatches + 1):
        if k not in agg.columns:
            agg[k] = 0

    agg = agg[list(range(max_mismatches + 1))].rename(columns={k: f"n_mm{k}" for k in range(max_mismatches + 1)})
    agg["n_hits"] = agg.sum(axis=1)

    out = agg.reset_index()
    return out[cols]