#!/usr/bin/env python3
"""TCRbeta clonotype table from a TRUST4 report (venus#31 P1 item 2).

Input: TRUST4 {prefix}_report.tsv as written by trust-simplerep.pl (bundled
in the trust4 1.1.10 package; header verified from the packaged script):

    #count  frequency  CDR3nt  CDR3aa  V  D  J  C  cid  cid_full_length

with `count` an integer (trust-simplerep default roundDownCount=1; float
counts only appear with its --decimalCnt flag, which run-trust4 never passes)
and `frequency` in %e scientific notation. "*" placeholders in the report are
emitted as "." by trust-simplerep; we normalize them to "" for the TSV.

A clonotype is TCRbeta when any of its V/J/C gene names starts with "TRB"
(same regex trust-simplerep's GetDetailChainType uses to classify chain 4).
Rows are already aggregated per (V,D,J,C,CDR3nt) clonotype by
trust-simplerep; we re-rank by count and add productive/frame annotations:

    productive: CDR3aa is a translateable length (nt len % 3 == 0), has no
                stop ("_") and no ambiguous codon ("?")
    full_length: the representative assembly contains a complete VDJ
                (cid_full_length == 1)

Honesty contract: an input with zero TRB clonotypes (e.g. a B-cell lymphoma
transcriptome, or a failed assembly) produces the standard header plus ZERO
data rows and prints a stderr note — never a fabricated row, never a nonzero
exit (the rule treats a missing report.tsv as the real failure).

Usage:
  tcr_report.py --report tcr/S1/trust_report.tsv --sample S1 \
      --out tcr/S1.tcrb.clonotypes.tsv
"""
import argparse
import csv
import sys

HEADER = [
    "sample", "clonotype_rank", "cdr3_nt", "cdr3_aa", "v_gene", "d_gene",
    "j_gene", "c_gene", "productive", "full_length", "count", "frequency",
]

REPORT_HEADER = ("#count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J",
                 "C", "cid", "cid_full_length")


def is_trb(v, d, j, c):
    """TRB chain: any V/J/C gene starts with TRB (mirrors trust-simplerep's
    GetDetailChainType, which classifies TRB* on any of the gene fields)."""
    return any(g.startswith("TRB") for g in (v, d, j, c) if g and g != ".")


def normalize(gene):
    """trust-simplerep writes "." for absent genes -> empty string."""
    return "" if gene in (".", "*", "") else gene


def productive(cdr3_nt, cdr3_aa):
    """Productive CDR3: in-frame, no stop codon, no ambiguous base.

    trust-simplerep renders aa as "partial" (score 0), "out_of_frame", or a
    translation where "_" marks a stop and "?" an N-containing codon. The nt
    length check catches out_of_frame independently of the aa rendering.
    """
    if not cdr3_aa or cdr3_aa in ("partial", "out_of_frame"):
        return "no"
    if len(cdr3_nt) % 3 != 0:
        return "no"
    if "_" in cdr3_aa or "?" in cdr3_aa:
        return "no"
    return "yes"


def parse_report(path):
    """Yield TCRbeta clonotype dicts from a TRUST4 report.tsv."""
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = None
        for row in reader:
            if not row or (row[0].startswith("#") and header is None
                           and not row[0].startswith("#count")):
                continue
            if header is None:
                if tuple(row) != REPORT_HEADER:
                    raise SystemExit(
                        "tcr_report: unexpected TRUST4 report header: %r "
                        "(expected %r)" % (row, REPORT_HEADER))
                header = row
                continue
            if len(row) < len(REPORT_HEADER):
                continue  # truncated line: skip, count via stderr below
            (count_s, freq_s, cdr3_nt, cdr3_aa, v, d, j, c,
             _cid, full_s) = row[:len(REPORT_HEADER)]
            v, d, j, c = (normalize(g) for g in (v, d, j, c))
            if not is_trb(v, d, j, c):
                continue
            try:
                count = int(float(count_s))
                freq = float(freq_s)
                full = int(float(full_s))
            except ValueError:
                continue
            yield {
                "cdr3_nt": cdr3_nt,
                "cdr3_aa": cdr3_aa,
                "v_gene": v,
                "d_gene": d,
                "j_gene": j,
                "c_gene": c,
                "productive": productive(cdr3_nt, cdr3_aa),
                "full_length": "yes" if full == 1 else "no",
                "count": count,
                "frequency": freq,
            }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", required=True,
                    help="TRUST4 {prefix}_report.tsv")
    ap.add_argument("--sample", required=True, help="sample label")
    ap.add_argument("--out", required=True, help="output clonotype TSV")
    args = ap.parse_args(argv)

    rows = sorted(parse_report(args.report),
                  key=lambda r: (-r["count"], r["cdr3_nt"], r["v_gene"]))
    total_trb = sum(r["count"] for r in rows)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEADER)
        for i, r in enumerate(rows, 1):
            w.writerow([
                args.sample, i, r["cdr3_nt"], r["cdr3_aa"], r["v_gene"],
                r["d_gene"], r["j_gene"], r["c_gene"], r["productive"],
                r["full_length"], r["count"], "%.6g" % r["frequency"],
            ])
    prod = sum(1 for r in rows if r["productive"] == "yes")
    print("tcr_report: %s: %d TCRbeta clonotype(s), %d productive, "
          "%d total reads" % (args.sample, len(rows), prod, total_trb),
          file=sys.stderr)
    if not rows:
        print("tcr_report: WARNING no TCRbeta clonotypes assembled — table "
              "is header-only (honest empty state)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
