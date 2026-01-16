"""
Off-target enumeration with Cas-OFFinder.
"""
from __future__ import annotations
import os, subprocess, pandas as pd
import re

_STRAND_MM_RE = re.compile(r"([+-])\s*([0-9]+)\s*$")

def write_cas_offinder_input(
    guides: pd.DataFrame,
    genome_fasta_path: str,
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNN",  # 20N (protospacer only)
    outfile: str = "cas_offinder_input.txt",
    mismatches: int = 4,
    seq_col: str = "spacer",
) -> str:
    """
    Write Cas-OFFinder input file.

    Notes:
      - Cas-OFFinder requires: len(pattern) == len(target_sequence) for every query line.
      - We use 20-nt protospacers by default (no PAM in the query).
    """
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)

    pattern = str(pam_pattern).strip().upper()
    patt_len = len(pattern)
    if patt_len == 0:
        raise ValueError("pam_pattern is empty.")

    if seq_col not in guides.columns:
        raise KeyError(f"Column '{seq_col}' not found in guides DataFrame.")

    seqs = (
        guides[seq_col]
        .dropna()
        .astype(str)
        .str.upper()
        .str.replace(r"\s+", "", regex=True)
        .unique()
        .tolist()
    )
    seqs = sorted(seqs)

    # Validate sequence lengths
    bad = [s for s in seqs if len(s) != patt_len]
    if bad:
        examples = ", ".join(bad[:5])
        raise ValueError(
            f"Cas-OFFinder input error: {len(bad)} sequences in '{seq_col}' do not match "
            f"pattern length {patt_len}. Examples: {examples}"
        )

    with open(outfile, "w") as f:
        f.write(genome_fasta_path + "\n")
        f.write(pattern + "\n")
        for seq in seqs:
            f.write(f"{seq} {int(mismatches)}\n")

    return outfile


def run_cas_offinder(input_file: str, cas_offinder_bin: str = "cas-offinder", device_spec: str = "C", output_file: str = "cas_offinder_hits.txt") -> str:
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    try:
        subprocess.run([cas_offinder_bin, input_file, device_spec, output_file], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "Cas-OFFinder failed.\n"
            f"Command: {cas_offinder_bin} {input_file} {device_spec} {output_file}\n"
            f"STDERR: {e.stderr[-2000:] if e.stderr else '(none)'}\n"
            f"STDOUT: {e.stdout[-2000:] if e.stdout else '(none)'}"
        ) from e
    return output_file


def parse_cas_offinder_output(hit_file: str) -> pd.DataFrame:
    """
    Parse Cas-OFFinder output robustly across formats, including lines ending with:
      ... <strand> <mm>
      ... <strand><mm>     (e.g. -0, +4)
    and where the chromosome field may contain whitespace/metadata.
    """
    rows = []

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

            # 1) Strand+mismatches are at the end, sometimes fused
            m = _STRAND_MM_RE.search(line)
            if not m:
                continue
            strand = m.group(1)
            mismatches = int(m.group(2))

            # 2) Remove trailing strand/mm and split remainder
            core = line[: m.start()].strip()
            toks = core.split()
            if len(toks) < 3:
                continue

            query = toks[0].upper()

            # target sequence is usually the last token in the "core"
            target = toks[-1].upper()

            # 3) Find the genomic position: last integer token before target
            pos_idx = None
            for i in range(len(toks) - 2, 0, -1):  # exclude query(0) and target(-1)
                if _is_int(toks[i]):
                    pos_idx = i
                    break
            if pos_idx is None:
                continue
            pos = int(toks[pos_idx])

            # 4) Chrom descriptor spans toks[1:pos_idx]; keep the first token as contig name
            chrom_desc = " ".join(toks[1:pos_idx])
            chrom = chrom_desc.split()[0] if chrom_desc else ""

            rows.append({
                "spacer": query[:20],
                "offtarget_target": target,
                "offtarget_chrom": chrom,
                "offtarget_pos": pos,
                "offtarget_strand": strand,
                "offtarget_mismatches": mismatches,
            })

    return pd.DataFrame(rows)


def summarize_specificity(offtargets: pd.DataFrame) -> pd.DataFrame:
    if offtargets.empty:
        return pd.DataFrame(columns=["spacer", "n_hits", "n_mm0", "n_mm1", "n_mm2", "n_mm3", "n_mm4"])
    agg = offtargets.groupby(["spacer", "offtarget_mismatches"]).size().unstack(fill_value=0)
    for k in range(5):
        if k not in agg.columns:
            agg[k] = 0
    agg = agg[[0, 1, 2, 3, 4]].rename(columns={i: f"n_mm{i}" for i in range(5)})
    agg["n_hits"] = agg.sum(axis=1)
    return agg.reset_index()[["spacer", "n_hits", "n_mm0", "n_mm1", "n_mm2", "n_mm3", "n_mm4"]]