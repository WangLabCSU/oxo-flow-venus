#!/usr/bin/env python3
"""Per-sample single-cell QC + clustering summary from dnbc4tools outputs.

Reads outs/analysis/QC_Cluster.h5ad (or falls back to cluster.csv) and emits
a compact TSV: cells per cluster, cluster fraction, and top marker genes per
cluster from marker.csv. Deliberately dependency-light (no scanpy import) —
the .h5ad path uses the cluster.csv that dnbc4tools already wrote.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cluster", required=True, help="outs/analysis/cluster.csv")
    ap.add_argument("--markers", required=True, help="outs/analysis/marker.csv")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    # cluster.csv: usually barcode,cluster (column names vary by version)
    with open(args.cluster) as fh:
        reader = csv.reader(fh)
        header = next(reader)
        col_c = next((i for i, h in enumerate(header)
                      if h.lower() in ("cluster", "louvain", "leiden", "seurat_clusters")),
                     len(header) - 1)
        counts: Counter[str] = Counter()
        for row in reader:
            if len(row) > col_c:
                counts[row[col_c]] += 1

    total = sum(counts.values()) or 1

    # markers: rows are (cluster, gene, score, ...) — take top 5 per cluster
    top_markers: dict[str, list[str]] = {}
    with open(args.markers) as fh:
        reader = csv.reader(fh)
        header = next(reader)
        cols = {h.lower(): i for i, h in enumerate(header)}
        c_cluster = cols.get("cluster", 0)
        c_gene = cols.get("gene", cols.get("gene_short_name", 1))
        c_score = cols.get("score", cols.get("avg_diff", 2))
        for row in reader:
            if len(row) <= max(c_cluster, c_gene, c_score):
                continue
            cl = row[c_cluster]
            try:
                score = float(row[c_score])
            except ValueError:
                score = 0.0
            top_markers.setdefault(cl, []).append((score, row[c_gene]))
    for cl in top_markers:
        top_markers[cl] = [g for _, g in sorted(top_markers[cl], reverse=True)[:5]]

    with open(args.output, "w") as out:
        out.write("sample\tcluster\tn_cells\tfraction\ttop_markers\n")
        for cl in sorted(counts, key=lambda c: -counts[c]):
            markers = ",".join(top_markers.get(cl, []))
            out.write(f"{args.sample}\t{cl}\t{counts[cl]}\t{counts[cl] / total:.3f}\t{markers}\n")
    print(f"[scrna_qc_cluster] {len(counts)} clusters for {args.sample} -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
