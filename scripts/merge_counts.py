#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs into a cohort gene-count matrix.

featureCounts writes a two-header TSV: line 1 is '# Geneid Chr Start End ...'
with sample path as the last column; the first data row repeats column names
(as of subread 2.x the second line is 'Geneid ...'). This script extracts
Geneid + the per-sample Counts column and outer-joins on gene id.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def sample_name(path: str) -> str:
    """Pt01-TR/rna/counts/Pt01-TR.featurecounts.txt -> Pt01-TR."""
    stem = Path(path).name
    for suffix in (".featurecounts.txt", ".txt"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem


def read_counts(path: str) -> tuple[dict[str, int], bool]:
    """Return (gene -> count, long-format-flag)."""
    counts: dict[str, int] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    start = 0
    if lines and lines[0].startswith("#"):
        start = 1
    header = lines[start].split("\t")
    # subread 2.x writes 'Geneid  Chr  Start ...' as the FIRST data row
    if header[0] == "Geneid" and start == 0:
        # no '#' prefix: first line IS the header
        pass
    col = header.index("Counts") if "Counts" in header else len(header) - 1
    gene_col = header.index("Geneid") if "Geneid" in header else 0
    for line in lines[start + 1:]:
        if not line:
            continue
        f = line.split("\t")
        if f[gene_col] == "Geneid":  # repeated header row emitted by 2.x
            continue
        try:
            counts[f[gene_col]] = int(f[col])
        except (ValueError, IndexError):
            continue
    return counts, True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    merged: dict[str, dict[str, int]] = {}
    order: list[str] = []
    for path in args.inputs:
        name = sample_name(path)
        order.append(name)
        counts, _ = read_counts(path)
        for gene, c in counts.items():
            merged.setdefault(gene, {})[name] = c
    samples = order
    with open(args.output, "w") as out:
        out.write("gene_id\t" + "\t".join(samples) + "\n")
        for gene in sorted(merged):
            out.write(gene + "\t" + "\t".join(str(merged[gene].get(s, 0)) for s in samples) + "\n")
    print(f"[merge_counts] {len(args.inputs)} samples, {len(merged)} genes -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
