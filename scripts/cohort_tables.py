#!/usr/bin/env python3
"""Build the patient-level cohort summary from per-module venus tables.

Inputs arrive as ONE list (--inputs) mixing module tables; each file is
routed by its basename so argument order does not matter:
  rna_qc_summary.tsv  - STAR bulk-RNA QC
  scrna_metrics.tsv   - single-cell metrics
  cohort_ascat.tsv    - ASCAT purity/ploidy cohort table
  msi_cohort.tsv      - MSI cohort table (BH-corrected)
  anything else       - per-pair TMB tables (concatenated)

Rows are normalized to PATIENT level: assay suffixes (-TD tumor / -B blood
normal / -TR bulk RNA) are stripped so Pt01-TD, Pt01-TR and Pt01 all land on
one patient row. The original per-module sample id is preserved as
<module>sample_id for provenance. A module the patient lacks leaves its
columns empty.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


def read_tsv(path: str) -> list[dict[str, str]]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def route(path: str) -> str:
    name = Path(path).name
    if name == "rna_qc_summary.tsv":
        return "rna"
    if name == "scrna_metrics.tsv":
        return "sc"
    if name == "cohort_ascat.tsv":
        return "ascat"
    if name == "msi_cohort.tsv":
        return "msi"
    return "tmb"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="module tables: RNA QC, scRNA metrics, ASCAT/MSI cohort tables, "
        "per-pair TMB tables",
    )
    ap.add_argument(
        "--suffixes",
        default="TD,B,TR",
        help="comma-separated assay suffixes stripped from sample ids to get "
        "the patient key (default: TD,B,TR)",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    groups: dict[str, list[str]] = {"rna": [], "sc": [], "ascat": [], "msi": [], "tmb": []}
    for path in args.inputs:
        groups[route(path)].append(path)
    if not groups["rna"] or not groups["sc"]:
        print(
            "[cohort_tables] need at least RNA QC + scRNA metrics inputs",
            file=sys.stderr,
        )
        return 1

    suffix_re = re.compile(r"-(" + "|".join(args.suffixes.split(",")) + r")$")

    def patient_of(sample: str) -> str:
        return suffix_re.sub("", sample)

    # Module tables in a stable order; the TMB source is the concatenation
    # of every expanded per-pair table.
    sources = [
        ("rna_", read_tsv(groups["rna"][0])),
        ("sc_", read_tsv(groups["sc"][0])),
        ("tmb_", [row for path in groups["tmb"] for row in read_tsv(path)]),
    ]
    if groups["ascat"]:
        sources.append(("ascat_", read_tsv(groups["ascat"][0])))
    if groups["msi"]:
        sources.append(("msi_", read_tsv(groups["msi"][0])))

    # patient -> {column: value}; insertion order keeps cohort sample order.
    patients: dict[str, dict[str, str]] = {}
    columns: list[str] = ["sample"]
    seen_cols: set[str] = {"sample"}
    for prefix, rows in sources:
        for row in rows:
            sample = row.get("sample", "")
            if not sample:
                continue
            rec = patients.setdefault(patient_of(sample), {"sample": patient_of(sample)})
            for col, val in row.items():
                out_col = prefix + ("sample_id" if col == "sample" else col)
                if out_col not in seen_cols:
                    seen_cols.add(out_col)
                    columns.append(out_col)
                rec[out_col] = val

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for rec in patients.values():
            w.writerow(rec)
    print(
        f"[cohort_tables] {len(patients)} patients, {len(columns)} columns -> {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
