#!/usr/bin/env python3
"""Extract per-sample single-cell metrics from dnbc4tools outputs.

dnbc4tools writes outs/metrics_summary.xls (tab-separated; key/value or a
headered table depending on version). We normalize whatever it emits into one
row per sample with the metrics that matter for cohort QC: estimated cells,
mean reads per cell, median UMI, sequencing saturation.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

KEY_ALIASES = {
    "sample": "sample",
    "cells": "estimated_cells",
    "estimated_cells": "estimated_cells",
    "mean_reads_per_cell": "mean_reads_per_cell",
    "median_umi_per_cell": "median_umi_per_cell",
    "median_genes_per_cell": "median_genes_per_cell",
    "sequencing_saturation": "sequencing_saturation",
    "saturation": "sequencing_saturation",
    "valid_barcodes": "valid_barcodes_pct",
    "valid_barcodes_pct": "valid_barcodes_pct",
    "total_reads": "total_reads",
}


def read_table(path: str) -> list[dict[str, str]]:
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        sample_lines = [ln for ln in fh if not ln.startswith("#")]
    if not sample_lines:
        return []
    try:
        reader = csv.DictReader(sample_lines, delimiter="\t")
        return [dict(r) for r in reader]
    except Exception:
        return []


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metrics", nargs="+", required=True,
                    help="one or more metrics_summary.xls paths")
    ap.add_argument("--samples", nargs="+", required=True,
                    help="sample names, same order as --metrics")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if len(args.metrics) != len(args.samples):
        print("[scrna_metrics] --metrics and --samples must have equal length", file=sys.stderr)
        return 1

    all_keys: list[str] = ["sample"]
    seen: set[str] = set()
    rows: list[dict[str, str]] = []
    for path, sample in zip(args.metrics, args.samples):
        recs = read_table(path)
        row: dict[str, str] = {"sample": sample}
        for rec in recs:
            for k, v in rec.items():
                norm = KEY_ALIASES.get((k or "").strip().lower().replace(" ", "_"))
                if not norm:
                    # keep unknown keys verbatim (snake_cased) for completeness
                    norm = "raw_" + (k or "col").strip().lower().replace(" ", "_")
                if norm == "sample":
                    continue
                row[norm] = (v or "").strip()
                if norm not in seen:
                    seen.add(norm)
                    all_keys.append(norm)
        rows.append(row)

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=all_keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"[scrna_metrics] {len(rows)} samples, {len(all_keys)} metrics -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
