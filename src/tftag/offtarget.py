"""
Off-target enumeration with Cas-OFFinder.
"""
from __future__ import annotations
import os, subprocess, pandas as pd


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
    rows = []

    def _is_int(tok: str) -> bool:
        try:
            int(tok)
            return True
        except Exception:
            return False

    with open(hit_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Split on any whitespace (handles tabs/spaces)
            toks = line.split()
            if len(toks) < 4:
                continue

            query = toks[0].upper()
            # Find first integer token after query => position
            pos_idx = None
            for i in range(1, len(toks)):
                if _is_int(toks[i]):
                    pos_idx = i
                    break
            if pos_idx is None:
                continue

            # Chromosome descriptor may contain many tokens; take first token for clean contig name
            chrom_field = " ".join(toks[1:pos_idx])
            chrom = chrom_field.split()[0] if chrom_field else ""

            pos = int(toks[pos_idx])

            # After pos: expect target, strand, mismatches (but strand+mismatch can be fused like '+4')
            if pos_idx + 1 >= len(toks):
                continue
            target = toks[pos_idx + 1].upper()

            strand = None
            mismatches = None

            # Case A: separate strand token
            if pos_idx + 2 < len(toks):
                s_tok = toks[pos_idx + 2]
                if s_tok in ("+", "-"):
                    strand = s_tok
                    # mismatches should be next token
                    if pos_idx + 3 < len(toks) and _is_int(toks[pos_idx + 3]):
                        mismatches = int(toks[pos_idx + 3])
                else:
                    # Case B: fused like '+4' or '-3'
                    if (s_tok.startswith("+") or s_tok.startswith("-")) and s_tok[1:].isdigit():
                        strand = s_tok[0]
                        mismatches = int(s_tok[1:])

            if strand is None or mismatches is None:
                # Could not parse; skip this line
                continue

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