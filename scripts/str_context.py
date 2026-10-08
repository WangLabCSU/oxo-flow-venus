#!/usr/bin/env python3
"""Annotate MAF indels with reference STR (microsatellite) repeat-unit context.

For each indel row (Variant_Type DEL or INS) in a MAF, a reference window
[start-15, end+15] (1-based inclusive) is fetched and scanned for the
shortest tandem-repeat unit whose periodic run covers the locus:

  unit length 1, run >= 4 identical bases   -> mononucleotide (homopolymer)
  unit length 2, run >= 4 units (8 bases)   -> dinucleotide
  unit lengths 3-6, run >= 3 complete units -> longer STR
  none of the above                         -> non-STR (unit_len = 0)

The "periodic run" is the longest stretch of consecutive period-u repeats
found anywhere in the window (all offsets scanned). Shortest qualifying unit
wins, so homopolymers are never re-counted as longer repeats.

Reference access is pure Python via the FASTA .fai index (no samtools needed;
random seeks against the uncompressed FASTA are fast enough for MAF-scale
variant counts).

Standalone utility: not referenced by any rule in the workflow — used for
manual annotation of delivered MAFs (see notes/interpretation.md), so it is
not shipped into run workdirs.

Usage:
  python3 str_context.py --fasta hg38.primary.fa \
      --maf deliver/Pt09.maf [--maf deliver/Pt05.maf ...] \
      [--out-dir report/str]

Writes <out-dir>/<sample>.str.tsv (one row per indel) when --out-dir is set,
and always prints a per-MAF summary: total indels, annotated STR count,
mono-/di-/longer-unit counts and shares.
"""
import argparse
import csv
import os
from collections import Counter

WINDOWS = {1: 4, 2: 8}  # unit length -> minimum run in bases; u>=3 uses 3 units


def best_repeat(w):
    """Return (unit_len, unit_seq, run_bases) for the shortest qualifying unit."""
    n = len(w)
    for u in range(1, 7):
        need = WINDOWS.get(u, 3 * u)
        top, topu = 0, ""
        for i in range(0, n - 2 * u + 1):
            unit = w[i:i + u]
            j = i
            while j + u <= n and w[j:j + u] == unit:
                j += u
            if j - i > top:
                top, topu = j - i, unit
        if top >= need:
            return u, topu, top
    return 0, "", 0


class Fasta:
    """Random-access reader for an uncompressed FASTA via its .fai index."""

    def __init__(self, path):
        self.path = path
        self.idx = {}
        with open(path + ".fai") as fh:
            for l in fh:
                p = l.rstrip("\n").split("\t")
                self.idx[p[0]] = tuple(int(x) for x in p[1:5])
        self.fh = open(path, "rb")

    def fetch(self, chrom, start, end):
        """Return uppercase sequence for 1-based inclusive [start, end]."""
        if chrom not in self.idx:
            return ""
        length, offset, lb, lw = self.idx[chrom]
        start, end = max(1, start), min(end, length)
        if start > end:
            return ""
        n = end - start + 1
        self.fh.seek(offset + (start - 1) // lb * lw + (start - 1) % lb)
        raw = self.fh.read(n + (n // lb + 2) * (lw - lb) + 16)
        return raw.decode("ascii", "ignore").replace("\n", "")[:n].upper()


def maf_indels(path):
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        if r["Variant_Type"] not in ("DEL", "INS"):
            continue
        yield r


def annotate(path, fa):
    counts = Counter()
    out_rows = []
    for r in maf_indels(path):
        s, e = int(r["Start_Position"]), int(r["End_Position"])
        w = fa.fetch(r["Chromosome"], s - 15, e + 15)
        ul, us, rb = (0, "", 0) if len(w) < 8 else best_repeat(w)
        counts[ul] += 1
        out_rows.append([r["Chromosome"], s, e, r["Variant_Type"], ul, us, rb])
    return os.path.basename(path).replace(".maf", ""), counts, out_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="reference FASTA (with .fai)")
    ap.add_argument("--maf", action="append", required=True, help="MAF file (repeatable)")
    ap.add_argument("--out-dir", help="write per-sample TSVs here")
    args = ap.parse_args()

    fa = Fasta(args.fasta)
    for m in args.maf:
        sample, counts, rows = annotate(m, fa)
        total = sum(counts.values())
        mono, di, longer = counts[1], counts[2], sum(v for k, v in counts.items() if k >= 3)
        ann = mono + di + longer
        print(f"[{sample}] indels={total} STR_annotated={ann} "
              f"({100 * ann / total:.1f}%) mono={mono} di={di} longer={longer} "
              f"nonSTR={total - ann}")
        if ann:
            print(f"[{sample}] mono share of annotated: {100 * mono / ann:.1f}%")
        if args.out_dir:
            os.makedirs(args.out_dir, exist_ok=True)
            out = os.path.join(args.out_dir, f"{sample}.str.tsv")
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["Chromosome", "Start_Position", "End_Position",
                            "Variant_Type", "Unit_Length", "Unit_Seq", "Run_Bases"])
                w.writerows(rows)
            print(f"[{sample}] wrote {out}")


if __name__ == "__main__":
    main()
