#!/usr/bin/env python3
"""Build cohort-level summary tables from per-sample pipeline outputs.

Inputs are plain TSVs produced by earlier venus rules (TMB table, STAR QC
table, single-cell metrics). Rows are keyed by sample; the script outer-joins
them into one cohort TSV for the clinical report.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read_tsv(path: str, key: str) -> dict[str, dict[str, str]]:
    table: dict[str, dict[str, str]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            k = row.get(key, "")
            if k:
                table[k] = row
    return table


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tmb", help="tmb_from_maf output (key: sample)")
    ap.add_argument("--rna-qc", help="summarize_star_logs output (key: sample)")
    ap.add_argument("--scrna-metrics", help="scrna metrics TSV (key: sample)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sources = []
    if args.tmb:
        sources.append(("tmb_", read_tsv(args.tmb, "sample")))
    if args.rna_qc:
        sources.append(("rna_", read_tsv(args.rna_qc, "sample")))
    if args.scrna_metrics:
        sources.append(("sc_", read_tsv(args.scrna_metrics, "sample")))

    keys: list[str] = []
    seen: set[str] = set()
    for _, table in sources:
        for k in table:
            if k not in seen:
                seen.add(k)
                keys.append(k)

    # Build prefixed columns to avoid collisions between modules
    columns: list[str] = ["sample"]
    prefixed: dict[str, dict[str, str]] = {k: {} for k in keys}
    for prefix, table in sources:
        sample_cols = list(next(iter(table.values())).keys()) if table else []
        for col in sample_cols:
            if col == "sample":
                continue
            columns.append(prefix + col)
        for k in keys:
            for col in sample_cols:
                if col != "sample":
                    prefixed[k][prefix + col] = table.get(k, {}).get(col, "")

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for k in keys:
            w.writerow({"sample": k, **prefixed[k]})
    print(f"[cohort_tables] {len(keys)} samples, {len(columns)} columns -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
