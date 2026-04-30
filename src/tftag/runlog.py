"""
Run logging utilities for TFTag.

Provides:
- human-readable text log
- machine-readable JSON log
- run metadata
- parameter capture
- summary statistics
- exception traceback logging
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any
import pandas as pd


def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _safe_git_commit() -> str | None:
    """Return current git commit hash if available."""
    try:
        res = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        return res.stdout.strip()
    except Exception:
        return None


def _json_safe(value: Any) -> Any:
    """Convert common non-JSON-safe values into JSON-safe forms."""
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(v) for v in value]
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    try:
        if not isinstance(value, (list, tuple, dict, set)) and pd.isna(value):
            return None
    except Exception:
        pass


@dataclass
class TFTagRunLogger:
    """
    Run logger for one TFTag pipeline execution.

    Usage
    -----
    logger = TFTagRunLogger(outdir="out", basename="tftag", run_id=run_id)
    logger.start(parameters=params)

    try:
        ...
        logger.add_summary("Input summary", {...})
        logger.success()
    except Exception:
        logger.failure()
        raise
    finally:
        logger.write()
    """

    outdir: str
    basename: str
    run_id: str
    command: list[str] | None = None

    start_time: float = field(default_factory=time.time)
    end_time: float | None = None
    status: str = "running"

    lines: list[str] = field(default_factory=list)
    data: dict[str, Any] = field(default_factory=dict)

    def start(self, parameters: dict[str, Any] | None = None) -> None:
        """Initialise the log."""
        os.makedirs(self.outdir, exist_ok=True)

        self.command = self.command or sys.argv[:]

        self.data = {
            "run_id": self.run_id,
            "status": self.status,
            "start_time": _now_iso(),
            "command": self.command,
            "working_directory": os.getcwd(),
            "python": sys.version.replace("\n", " "),
            "platform": platform.platform(),
            "git_commit": _safe_git_commit(),
            "parameters": parameters or {},
            "summaries": {},
            "errors": [],
        }

        self.log("=== TFTag run ===")
        self.log(f"Run ID: {self.run_id}")
        self.log(f"Start time: {self.data['start_time']}")
        self.log(f"Working directory: {self.data['working_directory']}")
        self.log(f"Command: {' '.join(self.command)}")

        if parameters:
            self.log("\n--- Parameters ---")
            for k, v in sorted(parameters.items()):
                self.log(f"{k}: {v}")

    def log(self, message: str = "") -> None:
        """Append a line to the text log and print it."""
        print(message)
        self.lines.append(str(message))

    def add_summary(self, title: str, summary: dict[str, Any]) -> None:
        """Add a named summary block to both text and JSON logs."""
        self.data["summaries"][title] = _json_safe(summary)

        self.log(f"\n--- {title} ---")
        for k, v in summary.items():
            self.log(f"{k}: {v}")

    def add_dataframe_summary(
        self,
        title: str,
        df: pd.DataFrame,
        *,
        score_col: str = "selection_score",
    ) -> None:
        """Add common useful statistics for a candidate/output dataframe."""
        summary: dict[str, Any] = {
            "rows": len(df),
        }

        if "gene_id" in df.columns:
            summary["genes"] = int(df["gene_id"].nunique())

        if {"gene_id", "tag"}.issubset(df.columns):
            summary["termini"] = int(df[["gene_id", "tag"]].drop_duplicates().shape[0])

        if "guide_found" in df.columns:
            guide_mask = df["guide_found"].fillna(True).astype(bool)
            summary["guide_rows"] = int(guide_mask.sum())
            summary["no_guide_rows"] = int((~guide_mask).sum())

        if "selection_tier" in df.columns:
            summary["selection_tier_counts"] = (
                df["selection_tier"]
                .value_counts(dropna=False)
                .sort_index()
                .to_dict()
            )

        if score_col in df.columns and df[score_col].notna().any():
            summary[f"{score_col}_min"] = float(df[score_col].min())
            summary[f"{score_col}_median"] = float(df[score_col].median())
            summary[f"{score_col}_max"] = float(df[score_col].max())

        for col in ["n_mm0", "n_mm1", "n_mm1_same_chr", "n_mm2_same_chr"]:
            if col in df.columns:
                summary[f"guides_with_{col}_gt0"] = int((pd.to_numeric(df[col], errors="coerce").fillna(0) > 0).sum())

        if "selection_warning" in df.columns:
            warnings = df["selection_warning"].fillna("none").astype(str)
            warnings = warnings[warnings != "none"]
            summary["rows_with_selection_warning"] = int(len(warnings))

        if "warnings" in df.columns:
            warnings = df["warnings"].fillna("none").astype(str)
            warnings = warnings[warnings != "none"]
            summary["rows_with_warning"] = int(len(warnings))

        self.add_summary(title, summary)

    def failure(self, exc: BaseException | None = None) -> None:
        """Record failure and traceback."""
        self.status = "failed"
        tb = traceback.format_exc()

        self.data["status"] = "failed"
        self.data["errors"].append(
            {
                "exception": repr(exc) if exc is not None else None,
                "traceback": tb,
            }
        )

        self.log("\n=== FAILED ===")
        if exc is not None:
            self.log(f"Exception: {repr(exc)}")
        self.log(tb)

    def success(self) -> None:
        """Record successful completion."""
        if self.status != "success":
            self.log("\n=== SUCCESS ===")
        self.status = "success"
        self.data["status"] = "success"

    def write(self) -> tuple[str, str]:
        """Write text and JSON logs to disk."""
        self.end_time = time.time()
        runtime_seconds = self.end_time - self.start_time

        self.data["end_time"] = _now_iso()
        self.data["runtime_seconds"] = runtime_seconds
        self.data["status"] = self.status

        self.log(f"End time: {self.data['end_time']}")
        self.log(f"Runtime seconds: {runtime_seconds:.2f}")
        self.log(f"Status: {self.status}")

        txt_path = os.path.join(self.outdir, f"{self.basename}_{self.run_id}.log")
        json_path = os.path.join(self.outdir, f"{self.basename}_{self.run_id}.json")

        with open(txt_path, "w", encoding="utf-8") as f:
            f.write("\n".join(self.lines) + "\n")

        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(_json_safe(self.data), f, indent=2)

        return txt_path, json_path
    
    def finish_early(self, reason: str) -> None:
        self.add_summary("Early stop", {
            "reason": reason,
            "early_termination": True
        })
        self.success()


