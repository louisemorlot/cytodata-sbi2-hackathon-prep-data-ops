#!/usr/bin/env python
"""
Full Aconcagua aggregation to GENE and GUIDE level.

Run this on the hackathon VM (S3 egress is free there; locally it reads ~tens of GB over the network).
It aggregates the 31.2M-cell, 1,607-feature Aconcagua cell table (the Funk/Brieflow essential-gene OPS
screen in HeLa) by gene (perturbation_id) and by guide (barcode), using DuckDB out-of-core so it never
holds the 323 GB in memory. The mapping-QC filter keeps only confidently single-guide cells.

Outputs (Parquet, one row per gene / per guide, mean of every OrganelleProfiler feature + n_cells):
  aconcagua_gene.parquet    one row per gene  (non-targeting excluded)
  aconcagua_guide.parquet   one row per guide (keeps non-targeting guides as controls; `gene` column too)

Usage:
  pip install duckdb
  python aggregate_aconcagua.py            # both levels
  python aggregate_aconcagua.py gene       # one level
"""
import sys
import duckdb

URL = ("https://ops-explorer-public.s3.us-west-2.amazonaws.com/"
       "aconcagua/aconcagua/cell_data.parquet")
FEAT = "avg(COLUMNS('^(cell|cytoplasm|nucleus)_.*'))"   # all ~1,600 OrganelleProfiler features

LEVELS = {
    # name: (group key, extra select cols, extra WHERE)
    "gene":  ("perturbation_id", "",
              "AND perturbation_id NOT LIKE 'nontargeting%'"),
    "guide": ("barcode", ", any_value(perturbation_id) AS gene", ""),
}


def run(level: str) -> None:
    key, extra, where = LEVELS[level]
    out = f"aconcagua_{level}.parquet"
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs;")
    print(f"[{level}] aggregating by {key} -> {out} (this scans the full table once)", flush=True)
    con.sql(f"""
        COPY (
            SELECT {key} AS {level},
                   count(*) AS n_cells{extra},
                   {FEAT}
            FROM read_parquet('{URL}')
            WHERE mapped_single_gene {where}
            GROUP BY {key}
        ) TO '{out}' (FORMAT PARQUET)
    """)
    n = con.sql(f"SELECT count(*) FROM read_parquet('{out}')").fetchall()[0][0]
    print(f"[{level}] wrote {out}: {n} rows", flush=True)


if __name__ == "__main__":
    levels = sys.argv[1:] or ["gene", "guide"]
    for lv in levels:
        if lv not in LEVELS:
            sys.exit(f"unknown level {lv!r}; choose from {list(LEVELS)}")
        run(lv)
    print("DONE")
