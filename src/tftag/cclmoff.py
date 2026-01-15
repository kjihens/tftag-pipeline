"""
CCLMoff adapter: build (on, off) pairs and aggregate predictions.
"""
from __future__ import annotations
import os, numpy as np, pandas as pd
from .coords import pam_coords, protospacer_coords
from .utils import get_sequence, revcomp


def _compute_onseq_23(row, fasta_dict):
    pam_s, pam_e = pam_coords(row)
    prot_s, prot_e = protospacer_coords(row)
    chrom = row["chromosome"]
    gstrand = row["grna_strand"]

    if pam_s is None or prot_s is None:
        return np.nan

    # In guide orientation, the contiguous 23-mer is protospacer + PAM
    start = prot_s
    end = pam_e
    return get_sequence(fasta_dict, chrom, start, end, gstrand)


def build_pairs_from_hits(guides_df: pd.DataFrame, offtargets_df: pd.DataFrame, fasta_dict, outfile: str = "out/cclmoff_pairs.tsv") -> pd.DataFrame:
    if guides_df.empty or offtargets_df.empty:
        raise ValueError("Need non-empty guides and offtargets to build pairs.")

    # pick one guide per spacer first (fast), then compute on_seq_23
    g = guides_df.sort_values(by=["cut_distance"]).drop_duplicates(subset=["spacer"], keep="first").copy()
    g["on_seq_23"] = [ _compute_onseq_23(r, fasta_dict) for r in g.to_dict(orient="records") ]
    g = g.dropna(subset=["on_seq_23"])

    hits = offtargets_df.copy()
    hits["offtarget_target"] = hits["offtarget_target"].astype(str).str.upper()

    merged = hits.merge(g[["spacer", "grna_strand", "on_seq_23"]], on="spacer", how="inner")

    # Orient off-target into guide orientation (if off-target strand differs from guide strand)
    same = merged["offtarget_strand"] == merged["grna_strand"]
    merged["off_seq_23"] = np.nan
    merged.loc[same, "off_seq_23"] = merged.loc[same, "offtarget_target"]

    opp = ~same
    merged.loc[opp, "off_seq_23"] = merged.loc[opp, "offtarget_target"].map(lambda s: revcomp(s) if isinstance(s, str) else np.nan)

    # Require exact 23-mer pairs
    merged = merged.dropna(subset=["off_seq_23", "on_seq_23"])
    merged = merged[merged["off_seq_23"].str.len() == 23]
    merged = merged[merged["on_seq_23"].str.len() == 23]

    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    merged[["on_seq_23", "off_seq_23"]].to_csv(outfile, sep="\t", header=False, index=False)

    return merged[["spacer", "on_seq_23", "off_seq_23"]].drop_duplicates()


def run_cclmoff_from_template(pairs_tsv: str, output_tsv: str, cmd_template: str) -> str:
    workdir = os.path.abspath(os.path.dirname(pairs_tsv) or ".")
    os.makedirs(workdir, exist_ok=True)

    # validate placeholders
    if "{pairs}" not in cmd_template or "{output}" not in cmd_template:
        raise ValueError("cclmoff_cmd template must contain {pairs} and {output} placeholders.")

    cmd = cmd_template.format(pairs=os.path.basename(pairs_tsv), output=os.path.basename(output_tsv), workdir=workdir)

    import subprocess
    try:
        subprocess.run(cmd, shell=True, check=True, cwd=workdir, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            "CCLMoff command failed.\n"
            f"CMD: {cmd}\n"
            f"STDERR: {e.stderr[-2000:] if e.stderr else '(none)'}\n"
            f"STDOUT: {e.stdout[-2000:] if e.stdout else '(none)'}"
        ) from e

    out_path = os.path.join(workdir, os.path.basename(output_tsv))
    if not os.path.exists(out_path):
        raise RuntimeError("CCLMoff did not produce the expected output file.")
    return out_path


def parse_cclmoff_output(preds_path: str) -> pd.DataFrame:
    def _try(path, sep, header):
        try:
            return pd.read_csv(path, sep=sep, header=header)
        except Exception:
            return None

    df = _try(preds_path, "\t", None)
    if df is None: df = _try(preds_path, ",", None)
    if df is None: df = _try(preds_path, "\t", "infer")
    if df is None: df = _try(preds_path, ",", "infer")

    if df is None or df.shape[1] < 3:
        raise ValueError("Unable to parse CCLMoff predictions (need 3+ columns).")

    df = df.iloc[:, :3].copy()
    df.columns = ["on_seq", "off_seq", "cclmoff_score"]
    df["on_seq"] = df["on_seq"].astype(str).str.upper()
    df["off_seq"] = df["off_seq"].astype(str).str.upper()
    df["cclmoff_score"] = pd.to_numeric(df["cclmoff_score"], errors="coerce")
    return df.dropna(subset=["cclmoff_score"])


def aggregate(pred_pairs: pd.DataFrame, preds: pd.DataFrame) -> pd.DataFrame:
    m = pred_pairs.merge(preds, left_on=["on_seq_23", "off_seq_23"], right_on=["on_seq", "off_seq"], how="inner")
    if m.empty:
        return pd.DataFrame(columns=["spacer", "cclmoff_max", "cclmoff_sum", "cclmoff_n"])
    grp = m.groupby("spacer")["cclmoff_score"]
    out = grp.agg(cclmoff_max="max", cclmoff_sum="sum", cclmoff_n="count").reset_index()
    return out