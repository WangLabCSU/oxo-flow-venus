#!/usr/bin/env python3
"""Build the patient-level cohort summary from per-module venus tables.

Inputs arrive as ONE ordered list (--inputs): [0] STAR bulk-RNA QC summary,
[1] single-cell metrics, [2:] per-pair TMB tables (the expanded tmb_per_pair
outputs — {input} is declared inputs followed by expanded ones).

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


def read_tsv(path: str) -> list[dict[str, str]]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="ordered module tables: [0] RNA QC, [1] scRNA metrics, [2:] TMB tables",
    )
    ap.add_argument(
        "--suffixes",
        default="TD,B,TR",
        help="comma-separated assay suffixes stripped from sample ids to get "
        "the patient key (default: TD,B,TR)",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if len(args.inputs) < 2:
        print(
            "[cohort_tables] need at least RNA QC + scRNA metrics inputs",
            file=sys.stderr,
        )
        return 1

    suffix_re = re.compile(r"-(" + "|".join(args.suffixes.split(",")) + r")$")

    def patient_of(sample: str) -> str:
        return suffix_re.sub("", sample)

    # Module tables in declared-input order; the TMB source is the
    # concatenation of every expanded per-pair table.
    sources = [
        ("rna_", read_tsv(args.inputs[0])),
        ("sc_", read_tsv(args.inputs[1])),
        ("tmb_", [row for path in args.inputs[2:] for row in read_tsv(path)]),
    ]

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
