#!/usr/bin/env python3
"""Parse SpecHLA typing results into the canonical venus HLA allele table.

Input : {outdir}/hla.result.txt   (SpecHLA wide table: Sample, HLA_A_1, ...)
        {outdir}/hla.result.g.group.txt (optional G-group resolution)
Output: long TSV — sample, gene, allele_1, allele_2 (2-field resolution,
        "HLA-A*02:01" style, pVACtools-compatible) + full-resolution strings
        kept in a separate column for provenance.

Tumor-only/organoid caveat is carried verbatim as a column so downstream
reports must render it (venus#31: organoid DNA typing lacks a germline
reference; LOH and purity distort tumor HLA — SpecHLA's built-in LOH output
should be cross-checked when available).
Stdlib only.
"""
import argparse
import csv
import re

GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
ALLELE_RE = re.compile(r"([A-Z0-9]+)\*(\d+):(\d+)")


def two_field(allele):
    """HLA-A*02:01:01:01 / A*02:01 / HLA_A*02:01 -> HLA-A*02:01 ('' if odd)."""
    s = str(allele or "").strip().replace("HLA_", "").replace("HLA-", "")
    m = ALLELE_RE.search(s)
    return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}" if m else ""


def parse_result_txt(path):
    """hla.result.txt -> {gene: (allele1_full, allele2_full)}."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    # line 0 is the SpecHLA version banner, line 1 the header
    if len(lines) < 2:
        raise SystemExit(f"[hla_parse] {path}: no result rows")
    header = lines[1].split("\t")
    data = lines[2].split("\t") if len(lines) > 2 else [""] * len(header)
    if len(data) < len(header):
        data += [""] * (len(header) - len(data))
    row = dict(zip(header, data))
    out = {}
    for g in GENES:
        out[g] = (row.get(f"HLA_{g}_1", "").strip(),
                  row.get(f"HLA_{g}_2", "").strip())
    return out


def parse_g_group(path):
    """hla.result.g.group.txt -> {gene: (g1, g2)}; tolerant to absence.
    The file is a wide table with the same HLA_<gene>_1/2 columns as
    hla.result.txt (audit M4)."""
    if not path:
        return {}
    try:
        with open(path) as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    except OSError:
        return {}
    if len(lines) < 2:
        return {}
    header = lines[0].split("\t")
    # skip a version/banner row if the first line is not the header
    if "HLA_A_1" not in header:
        if len(lines) < 3:
            return {}
        header = lines[1].split("\t")
        lines = lines[1:]
    data = lines[1].split("\t") if len(lines) > 1 else [""] * len(header)
    if len(data) < len(header):
        data += [""] * (len(header) - len(data))
    row = dict(zip(header, data))
    return {g: (row.get(f"HLA_{g}_1", "").strip(),
                row.get(f"HLA_{g}_2", "").strip()) for g in GENES}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--result", required=True)
    ap.add_argument("--g-group", default="")
    ap.add_argument("--spec-version", default="")
    ap.add_argument("--mode", default="tumor_only",
                    choices=["tumor_only", "paired"],
                    help="provenance caveat for the report layer")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    res = parse_result_txt(args.result)
    ggrp = parse_g_group(args.g_group)
    rows = []
    for g in GENES:
        a1, a2 = res[g]
        gg1, gg2 = ggrp.get(g, ("", ""))
        rows.append({
            "sample": args.sample, "gene": g,
            "allele_1": two_field(a1), "allele_2": two_field(a2),
            "allele_1_full": a1, "allele_2_full": a2,
            "allele_1_g_group": gg1, "allele_2_g_group": gg2,
            "typed": "yes" if (two_field(a1) or two_field(a2)) else "no",
            "spec_version": args.spec_version, "mode": args.mode,
            "caveat": ("organoid/tumor DNA typing — no matched germline; "
                       "HLA-LOH and purity may distort calls"
                       if args.mode == "tumor_only" else ""),
        })
    cols = ["sample", "gene", "allele_1", "allele_2", "allele_1_full",
            "allele_2_full", "allele_1_g_group", "allele_2_g_group",
            "typed", "spec_version", "mode", "caveat"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_typed = sum(1 for r in rows if r["typed"] == "yes")
    import sys
    sys.stderr.write(f"[hla_parse] {args.sample}: {n_typed}/8 genes typed "
                     f"-> {args.out}\n")


if __name__ == "__main__":
    main()
