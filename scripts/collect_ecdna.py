#!/usr/bin/env python3
"""Collect per-pair AmpliconSuite classification outputs into cohort tables.

Inputs arrive as ONE list mixing two file families (partitioned by name,
same convention as collect_msi.py):
  ecdna/<pair>/<pair>_ecDNA_counts.tsv  - one row per sample: #sample, ecDNA_count
  ecdna/<pair>/<pair>_result_table.tsv  - one row per amplicon feature

Outputs:
  --counts-output  cohort_ecdna.tsv: `sample<TAB>ecDNA_count` (the leading
                   `#` of the suite's `#sample` column is stripped so
                   cohort_tables.py can join on `sample` like every other
                   module table)
  --genes-output   ecdna_genes.tsv: the suite's gene_list flattened to one
                   row per (sample, amplicon, gene), with the amplicon
                   classification joined from the result_table via Feature ID
                   = sample_amplicon_feature (the suite's own join, cf.
                   make_results_table.py read_amplicon_gene_list). Kept
                   header-only when no pair has amplicons.

0-amplicon samples are safe: AmpliconClassifier still writes a counts TSV
(count 0) and an NA-classification result_table row, so the pair keeps its
row in the counts table.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

COUNTS_SUFFIX = "_ecDNA_counts.tsv"
RESULT_SUFFIX = "_result_table.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # The suite prefixes the counts header with '#' (`#sample<TAB>ecDNA_count`);
    # strip it so the column is named `sample` like every other module table.
    for row in rows:
        for key in list(row):
            if key.startswith("#"):
                row[key.lstrip("#")] = row.pop(key)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="per-pair <pair>_ecDNA_counts.tsv and "
                         "<pair>_result_table.tsv files (one list, any order)")
    ap.add_argument("--counts-output", required=True)
    ap.add_argument("--genes-output", required=True)
    args = ap.parse_args()

    counts_files: dict[str, Path] = {}
    result_files: dict[str, Path] = {}
    for p in args.inputs:
        path = Path(p)
        name = path.name
        if name.endswith(COUNTS_SUFFIX):
            counts_files[name.removesuffix(COUNTS_SUFFIX)] = path
        elif name.endswith(RESULT_SUFFIX):
            result_files[name.removesuffix(RESULT_SUFFIX)] = path
        else:
            print(f"[collect_ecdna] WARNING unrecognized input {p}",
                  file=sys.stderr)

    missing = sorted(set(counts_files) ^ set(result_files))
    for m in missing:
        print(f"[collect_ecdna] WARNING unmatched input for pair {m}",
              file=sys.stderr)

    counts_header = ["sample", "ecDNA_count"]
    count_rows: list[list[str]] = []
    genes_header = [
        "sample", "amplicon", "classification", "gene", "gene_cn",
        "truncated", "is_canonical_oncogene",
    ]
    gene_rows: list[list[str]] = []

    for pair in sorted(counts_files):
        for row in read_tsv(counts_files[pair]):
            count_rows.append([row.get("sample", ""), row.get("ecDNA_count", "")])

        # result_table: one row per amplicon feature; `Classification` carries
        # the ecDNA/BFB/linear verdict, keyed by Feature ID =
        # <sample>_<amplicon>_<feature> — exactly the suite's own join key
        # (make_results_table.py builds Feature ID from the first three
        # gene_list fields: sample_name, amplicon_number, feature). With 0
        # amplicons the result_table has only the `..._NA` row (skipped) and
        # the gene_list is header-only; the sample still keeps its counts row.
        classes: dict[str, str] = {}
        rt = result_files.get(pair)
        if rt:
            for row in read_tsv(rt):
                feat = row.get("Feature ID", "")
                if feat and feat != "NA":
                    classes[feat] = row.get("Classification", "")

        gene_path = counts_files[pair].with_name(f"{pair}_gene_list.tsv")
        if gene_path.exists():
            for row in read_tsv(gene_path):
                sample = row.get("sample_name", "")
                amp = row.get("amplicon_number", "")
                feat = row.get("feature", "")
                gene_rows.append([
                    sample, amp,
                    classes.get("_".join([sample, amp, feat]), ""),
                    row.get("gene", ""), row.get("gene_cn", ""),
                    row.get("truncated", ""),
                    row.get("is_canonical_oncogene", ""),
                ])

    with open(args.counts_output, "w") as out:
        out.write("\t".join(counts_header) + "\n")
        for row in count_rows:
            out.write("\t".join(row) + "\n")
    with open(args.genes_output, "w") as out:
        out.write("\t".join(genes_header) + "\n")
        for row in gene_rows:
            out.write("\t".join(row) + "\n")
    print(f"[collect_ecdna] {len(count_rows)} pairs, {len(gene_rows)} gene hits"
          f" -> {args.counts_output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
