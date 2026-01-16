
from __future__ import annotations
import os, re, sqlite3
import pandas as pd


_VALID_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")



def _validate_ident(x: str, what: str) -> str:
    if not _VALID_IDENT.match(x):
        raise ValueError(f"Invalid {what}: {x!r}. Use letters/numbers/underscore; cannot start with number.")
    return x


def to_sqlite(
    df: pd.DataFrame,
    db_path: str,
    table: str,
    if_exists: str = "append",
    index: bool = False,
    create_indices: bool = False,
    chunksize: int = 5000,
) -> None:
    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)
    table = _validate_ident(table, "table name")

    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA journal_mode=WAL;")
        cur.execute("PRAGMA synchronous=NORMAL;")
        cur.execute("PRAGMA temp_store=MEMORY;")

        df.to_sql(table, conn, if_exists=if_exists, index=index, chunksize=chunksize, method="multi")

        if create_indices:
            # Create indices only for columns that exist.
            idx_specs = [
                ("gene_id", ["gene_id"]),
                ("gene_symbol", ["gene_symbol"]),
                ("tag", ["tag"]),
                ("spacer", ["spacer"]),
                ("coord", ["chromosome", "protospacer_start", "protospacer_end"]),
                ("rs3", ["rs3_score"]),
                ("cclmoff_max", ["cclmoff_max"]),
            ]
            for suffix, cols in idx_specs:
                if all(c in df.columns for c in cols):
                    idx_name = f"idx_{table}_{suffix}"
                    cols_sql = ", ".join([f'"{c}"' for c in cols])
                    cur.execute(f'CREATE INDEX IF NOT EXISTS "{idx_name}" ON "{table}"({cols_sql});')

        conn.commit()