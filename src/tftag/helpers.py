from __future__ import annotations
import os
import pandas as pd


def parse_genes_arg(genes: str | None) -> list[str] | None:
    """
    Parse `genes` argument which may be:
      - None
      - path to a file (one gene per line)
      - comma-separated list
      - single gene id
    """
    if genes is None:
        return None
    genes = str(genes).strip()
    if not genes:
        return None
    if os.path.exists(genes):
        with open(genes) as fh:
            return [ln.strip() for ln in fh if ln.strip()]
    return [g.strip() for g in genes.split(",") if g.strip()]


