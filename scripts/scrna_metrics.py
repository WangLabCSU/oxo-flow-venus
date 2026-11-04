#!/usr/bin/env python3
"""Extract per-sample single-cell metrics from dnbc4tools outputs.

dnbc4tools writes outs/metrics_summary.xls (tab-separated; key/value or a
headered table depending on version). We normalize whatever it emits into one
row per sample with the metrics that matter for cohort QC: estimated cells,
mean reads per cell, median UMI, sequencing saturation.

N13 (audit 20261104): dnbc4tools 3.0 has NO ambient-RNA or doublet metrics —
those columns are emitted as an explicit not_assessed marker (NOT blank, NOT
zero) so report consumers can tell "never measured" apart from "measured 0".
If a future dnbc4tools version starts reporting them, the KEY_ALIASES entries
promote the real values automatically.

Two input quirks handled here:
- oxo-flow config lists interpolate as one comma-joined string, so --samples
  entries are split on commas before pairing with --metrics paths.
- dnbc4tools 3.0 headers are human-readable ("Estimated number of cell") and
  values embed thousands separators ("7,390", "83.76%").
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

KEY_ALIASES = {
    "sample": "sample",
    "samplename": "sample",
    "cells": "estimated_cells",
    "estimated_cells": "estimated_cells",
    "estimated_number_of_cell": "estimated_cells",
    "mean_reads_per_cell": "mean_reads_per_cell",
    "mean_umi_counts_per_cell": "mean_umi_per_cell",
    "median_umi_per_cell": "median_umi_per_cell",
    "median_umi_counts_per_cell": "median_umi_per_cell",
    "median_genes_per_cell": "median_genes_per_cell",
    "mean_genes_per_cell": "mean_genes_per_cell",
    "total_genes_detected": "total_genes_detected",
    "sequencing_saturation": "sequencing_saturation",
    "saturation": "sequencing_saturation",
    "fraction_reads_in_cell": "fraction_reads_in_cell",
    "valid_barcodes": "valid_barcodes_pct",
    "valid_barcodes_pct": "valid_barcodes_pct",
    "cdna_valid_barcodes": "valid_barcodes_pct",
    "species": "species",
    "total_reads": "total_reads",
    # N13 (audit 20261104): ambient/doublet. dnbc4tools 3.0's
    # metrics_summary.xls carries NEITHER field — if a future version starts
    # emitting one, these aliases promote it to a real measurement instead of
    # letting it fall into the raw_* catch-all. Absent values stay an explicit
    # not_assessed marker (see NOT_ASSESSED_DEFAULTS) so downstream tables
    # never present a blank as if it were a measurement of zero.
    "doublet_pct": "doublet_pct",
    "doublet_rate": "doublet_pct",
    "doublet_fraction": "doublet_pct",
    "predicted_doublets": "doublet_pct",
    "doublet_score": "doublet_score",
    "ambient_pct": "ambient_pct",
    "ambient_rate": "ambient_pct",
    "ambient_fraction": "ambient_pct",
    "ambient_rna_fraction": "ambient_pct",
}

NOT_ASSESSED_DEFAULTS = {
    "doublet_pct": "not_assessed(tool_lacks)",
    "ambient_pct": "not_assessed(tool_lacks)",
}

# "7,390" / "1,234.5" / "83.76%" — commas used purely as thousands separators
NUMERIC_RE = re.compile(r"^-?\d{1,3}(,\d{3})+(\.\d+)?%?$")


def clean_value(v: str | None) -> str:
    v = (v or "").strip()
    if NUMERIC_RE.match(v):
        v = v.replace(",", "")
    return v


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
                    help="sample names, same order as --metrics; entries may be comma-joined")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    # config lists arrive as a single comma-joined string via oxo-flow's
    # {config.*} interpolation — split entries back into individual names.
    samples: list[str] = []
    for entry in args.samples:
        samples.extend(p.strip() for p in entry.split(",") if p.strip())

    if len(args.metrics) != len(samples):
        print(f"[scrna_metrics] --metrics ({len(args.metrics)}) and --samples ({len(samples)}) "
              "must have equal length", file=sys.stderr)
        return 1

    all_keys: list[str] = ["sample"]
    seen: set[str] = set()
    rows: list[dict[str, str]] = []
    for path, sample in zip(args.metrics, samples):
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
                row[norm] = clean_value(v)
                if norm not in seen:
                    seen.add(norm)
                    all_keys.append(norm)
        # N13: ambient/doublet are not assessed by dnbc4tools 3.0 — emit the
        # explicit marker unless the tool actually provided a value above.
        for col, marker in NOT_ASSESSED_DEFAULTS.items():
            if not row.get(col):
                row[col] = marker
                if col not in seen:
                    seen.add(col)
                    all_keys.append(col)
        rows.append(row)

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=all_keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"[scrna_metrics] {len(rows)} samples, {len(all_keys)} metrics -> {args.output}", file=sys.stderr)
    print("[scrna_metrics] ambient/doublet: dnbc4tools 3.0 lacks these fields; "
          "reported as not_assessed (P2: scDblFinder/DecontX on protected "
          "raw_matrix/)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
