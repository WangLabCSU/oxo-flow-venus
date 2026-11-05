#!/usr/bin/env python3
"""Collect per-sample neoantigen summary metrics into one cohort table.

Reads the long-format ``neoantigen/{pair}.summary.tsv`` files written by
neoantigen_report (venus#31) and pivots the metric rows into one wide row
per sample, so the clinical report can render a cohort-level neoantigen
burden table (2026-11-05: the final4 numbers previously lived only in the
per-sample files and never reached venus_clinical_report.md).

Missing metric rows become honest empty cells (never silently dropped);
the ``predictors`` / ``predictors_note`` rows carry the MHC-I/MHC-II
predictor provenance (MHC-II is license-gated and absent in P0).

Usage:
  collect_neoantigen.py --inputs neoantigen/GC01.summary.tsv ... \
      --output report/neoantigen_cohort.tsv
"""
import argparse
import os
import sys

COLUMNS = [
    "sample",
    "hla_alleles_typed",
    "hla_alleles_predictor_usable",
    "unique_candidates",
    "burden_strong_500nM",
    "burden_highaff_100nM",
    "expr_tpm_pass",
    "frameshift_neoantigens",
    "predictors",
    "predictors_note",
]


def read_tsv(path):
    import csv
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sample_name(path):
    return os.path.basename(path).split(".")[0]


def collect_row(path):
    row = {col: "" for col in COLUMNS}
    row["sample"] = sample_name(path)
    for line in read_tsv(path):
        metric = (line.get("metric") or "").strip()
        value = (line.get("value") or "").strip()
        if metric == "predictors":
            row["predictors"] = value
            row["predictors_note"] = (line.get("detail") or "").strip()
        elif metric in row:
            row[metric] = value
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="neoantigen/{pair}.summary.tsv files")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = [collect_row(p) for p in sorted(args.inputs)]
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[collect_neoantigen] wrote {len(rows)} sample(s) -> {args.output}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
