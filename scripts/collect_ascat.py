#!/usr/bin/env python3
"""Concatenate per-pair ASCAT purity/ploidy tables into one cohort table.

Reads the expanded ascat_paired outputs (ascat/<pair>/purity_ploidy.tsv —
all share the same header) and writes a single header-once TSV so downstream
cohort_tables.py / clinical_report.py consume one file.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="per-pair purity_ploidy.tsv files")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    header: list[str] | None = None
    rows: list[str] = []
    for path in args.inputs:
        lines = Path(path).read_text().splitlines()
        if not lines:
            print(f"[collect_ascat] WARNING empty input {path}", file=sys.stderr)
            continue
        if header is None:
            header = lines[0]
        elif lines[0] != header:
            print(f"[collect_ascat] WARNING header mismatch in {path}", file=sys.stderr)
        rows.extend(lines[1:])

    if header is None:
        print("[collect_ascat] no input rows", file=sys.stderr)
        return 1
    Path(args.output).write_text("\n".join([header, *rows]) + "\n")
    print(f"[collect_ascat] {len(rows)} pairs -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
