#!/usr/bin/env python3
"""Merge per-sample featureCounts outputs into a cohort gene-count matrix.

featureCounts writes a two-header TSV: line 1 is '# Geneid Chr Start End ...'
with sample path as the last column; the first data row repeats column names
(as of subread 2.x the second line is 'Geneid ...'). This script extracts
Geneid + the per-sample Counts column and outer-joins on gene id.

A GTF (--gtf, optional but wired in the workflow) supplies the HGNC symbol
for each ENSEMBL gene id, emitted as the second column ``gene_name`` so the
matrix is directly readable without a separate id-mapping step.
"""
from __future__ import annotations

import argparse
import gzip
import re
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


def strip_version(gene_id: str) -> str:
    """ENSG00000223972.5 -> ENSG00000223972 (PAR_Y suffix preserved)."""
    return re.sub(r"\.\d+", "", gene_id)


def read_gene_names(gtf: str) -> dict[str, str]:
    """gene_id (versionless) -> gene_name, from GTF 'gene' feature rows."""
    names: dict[str, str] = {}
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = gn = None
            for attr in f[8].split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id "):
                    gid = attr[8:].strip().strip('"')
                elif attr.startswith("gene_name "):
                    gn = attr[10:].strip().strip('"')
                if gid and gn:
                    break
            if gid:
                names[strip_version(gid)] = gn or ""
    return names


def read_counts(path: str) -> dict[str, int]:
    """Return gene -> count."""
    counts: dict[str, int] = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    start = 0
    if lines and lines[0].startswith("#"):
        start = 1
    header = lines[start].split("\t")
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
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--gtf", default=None,
                    help="GENCODE GTF for gene_name mapping (optional)")
    args = ap.parse_args()

    gene_names = read_gene_names(args.gtf) if args.gtf else {}
    merged: dict[str, dict[str, int]] = {}
    order: list[str] = []
    for path in args.inputs:
        name = sample_name(path)
        order.append(name)
        for gene, c in read_counts(path).items():
            merged.setdefault(gene, {})[name] = c
    samples = order
    with open(args.output, "w") as out:
        out.write("gene_id\tgene_name\t" + "\t".join(samples) + "\n")
        for gene in sorted(merged):
            symbol = gene_names.get(strip_version(gene), "")
            out.write(gene + "\t" + symbol + "\t" +
                      "\t".join(str(merged[gene].get(s, 0)) for s in samples) + "\n")
    mapped = sum(1 for g in merged if gene_names.get(strip_version(g)))
    print(f"[merge_counts] {len(args.inputs)} samples, {len(merged)} genes "
          f"({mapped} with gene_name) -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
