"""
Off-target enumeration with Cas-OFFinder.
"""
from __future__ import annotations
import os, subprocess, pandas as pd


def write_cas_offinder_input(
    guides: pd.DataFrame,
    genome_fasta_path: str,
    pam_pattern: str = "NNNNNNNNNNNNNNNNNNNNNGG",
    outfile: str = "cas_offinder_input.txt",
    mismatches: int = 4,
) -> str:
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    with open(outfile, "w") as f:
        f.write(genome_fasta_path + "\n")
        f.write(pam_pattern + "\n")
        for seq in sorted(set(guides["spacer"].dropna().astype(str).str.upper())):
            f.write(f"{seq} {mismatches}\n")
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
    with open(hit_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                parts = line.strip().split()
            if len(parts) < 7:
                continue
            query_seq, target_seq, chrom, pos, strand, mismatches, patt = parts[:7]
            try:
                pos = int(pos)
                mismatches = int(mismatches)
            except Exception:
                continue
            rows.append(
                {
                    "spacer": query_seq[:20].upper(),
                    "offtarget_target": target_seq.upper(),
                    "offtarget_chrom": chrom,
                    "offtarget_pos": pos,
                    "offtarget_strand": strand,
                    "offtarget_mismatches": mismatches,
                }
            )
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