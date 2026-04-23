
from __future__ import annotations
import os, re, sqlite3, math
import pandas as pd


_VALID_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")



def _validate_ident(x: str, what: str) -> str:
    if not _VALID_IDENT.match(x):
        raise ValueError(f"Invalid {what}: {x!r}. Use letters/numbers/underscore; cannot start with number.")
    return x


def to_sqlite(
    df,
    basename,
    outdir,
    table,
    if_exists="append",
    index=False,
    create_indices=True,
    chunksize=1000,
):
    os.makedirs(outdir, exist_ok=True)
    db_path = os.path.join(outdir, basename + ".sqlite")

    # SQLite has a limit on the number of variables in one statement.
    # With method="multi", total bound vars ~= n_columns * chunksize.
    ncols = len(df.columns) + (1 if index else 0)

    # Conservative default for SQLite builds that use 999 variables.
    max_sql_vars = 999

    # Leave a bit of margin.
    safe_chunksize = max(1, min(chunksize, max_sql_vars // max(1, ncols)))

    with sqlite3.connect(db_path) as conn:
        df.to_sql(
            table,
            conn,
            if_exists=if_exists,
            index=index,
            chunksize=safe_chunksize,
            method="multi",
        )

        if create_indices:
            cur = conn.cursor()

            # create indices only for columns that actually exist
            existing = set(df.columns)

            if "gene_symbol" in existing:
                cur.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_gene_symbol ON {table}(gene_symbol);")
            if "gene_id" in existing:
                cur.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_gene_id ON {table}(gene_id);")
            if "tag" in existing:
                cur.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_tag ON {table}(tag);")
            if {"chromosome", "grna_23_start", "grna_23_end"}.issubset(existing):
                cur.execute(
                    f"CREATE INDEX IF NOT EXISTS idx_{table}_coord "
                    f"ON {table}(chromosome, grna_23_start, grna_23_end);"
                )
            if "rs3_score" in existing:
                cur.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_rs3 ON {table}(rs3_score);")

            conn.commit()