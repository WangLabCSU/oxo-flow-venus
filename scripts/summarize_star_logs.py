#!/usr/bin/env python3
"""Summarize STAR *Log.final.out files into a cohort QC table.

STAR's Log.final.out is tab-separated key\tvalue pairs. We extract the
canonical QC metrics and emit one row per sample.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

METRICS = [
    ("Number of input reads", "total_reads"),
    ("Uniquely mapped reads %", "uniquely_mapped_pct"),
    ("Uniquely mapped reads number", "uniquely_mapped"),
    ("% of reads mapped to multiple loci", "multi_mapped_pct"),
    ("% of reads unmapped: too short", "unmapped_too_short_pct"),
    ("Average mapped length", "avg_mapped_length"),
    ("Number of splices: Total", "splices_total"),
    ("Mismatch rate per base, %", "mismatch_rate_pct"),
    ("Insertion rate per base", "insertion_rate"),
    ("Deletion rate per base", "deletion_rate"),
]


def parse_log(path: str) -> dict[str, str]:
    row: dict[str, str] = {"sample": Path(path).name.replace("Log.final.out", "").rstrip("._-")}
    text = Path(path).read_text()
    for key, column in METRICS:
        m = re.search(re.escape(key) + r"\s*\|\s*(.+)", text)
        if m:
            row[column] = m.group(1).strip().replace("%", "")
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--logs", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = [parse_log(p) for p in args.logs]
    columns = ["sample"] + [c for _, c in METRICS]
    with open(args.output, "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(row.get(c, "") for c in columns) + "\n")
    print(f"[summarize_star_logs] {len(rows)} samples -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
