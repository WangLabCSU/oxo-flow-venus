#!/usr/bin/env python3
"""Cohort-level single-cell integration summary.

Aggregates per-sample cluster tables (scrna_qc_cluster.py output) into a
cohort TSV: cells per sample, clusters detected, and shared marker signatures.
Pure-python (no scanpy/anndata) by design — heavy integration runs separately
when needed; this powers the report.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clusters", nargs="+", required=True,
                    help="scrna_qc_cluster outputs (one per sample)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    samples: dict[str, dict[str, str]] = {}
    marker_index: dict[str, set[str]] = {}
    for path in args.clusters:
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                s = row["sample"]
                rec = samples.setdefault(s, {"sample": s, "n_cells": "0",
                                             "n_clusters": "0", "clusters": "",
                                             "top_markers": ""})
                rec["n_cells"] = str(int(rec["n_cells"]) + int(row.get("n_cells", 0)))
                rec["clusters"] = (rec["clusters"] + "," + row["cluster"]).lstrip(",")
                for g in (row.get("top_markers") or "").split(","):
                    if g:
                        marker_index.setdefault(s, set()).add(g)
    for s, rec in samples.items():
        rec["n_clusters"] = str(len(rec["clusters"].split(",")))
        rec["top_markers"] = ",".join(sorted(marker_index.get(s, set()))[:15])

    cols = ["sample", "n_cells", "n_clusters", "clusters", "top_markers"]
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for s in sorted(samples):
            w.writerow(samples[s])
    print(f"[scrna_integrate] {len(samples)} samples -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
