#!/usr/bin/env python3
"""Fit COSMIC v3.1 SBS96 signatures to a venus MAF (SNPs only).

Steps:
  1. Count single-base substitutions (Variant_Type == SNP, single-base ACGT
     ref/alt) into the 96 SBS96 channels. Each variant gets the trinucleotide
     context centred on its reference base, read from the reference FASTA via
     .fai random access (same reader pattern as str_context.py), then
     pyrimidine-normalized: when the reference base is A or G the triplet and
     alleles are reverse-complemented so the altered pyrimidine (C or T) sits
     in the middle. Channels with non-ACGT context (reference N, IUPAC) are
     skipped and reported.
  2. NNLS-fit the count vector against the COSMIC v3.1 (GRCh38) signature
     matrix (96 rows x N signatures). Signatures are probability vectors
     summing to 1, so exposures are in units of attributed mutations.
  3. Write the 96-channel counts and the fitted exposures (fraction of the
     reconstruction per signature, sorted descending) plus the cosine
     similarity between counts and reconstruction.

The reconstruction may include negative-control-style flat signatures (SBS84
etc.) present in COSMIC v3.1; NNLS never assigns them negative weight, so a
spurious small exposure can appear on a weakly-identified signature — always
read exposures against the cosine similarity.

Usage:
  python3 sbs96_fit.py --maf deliver/Pt09.maf \
      --genome-fa /path/GRCh38.fa --signatures resources/cosmic_v3.1_sbs96_grch38.txt \
      --sample Pt09 --counts-out report/signatures/Pt09.sbs96.counts.tsv \
      --exposures-out report/signatures/Pt09.exposures.tsv
"""
import argparse
import csv
import sys
from collections import Counter

import numpy as np
from scipy.optimize import nnls

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
BASES = ("A", "C", "G", "T")
ALT_FROM_C = ("A", "G", "T")   # C>A, C>G, C>T
ALT_FROM_T = ("A", "C", "G")   # T>A, T>C, T>G

CHANNELS = [f"{x}[{c}>{a}]{z}"
            for x in BASES for c, alts in (("C", ALT_FROM_C), ("T", ALT_FROM_T))
            for a in alts for z in BASES]


class Fasta:
    """Random-access reader for an uncompressed FASTA via its .fai index."""

    def __init__(self, path):
        self.idx = {}
        with open(path + ".fai") as fh:
            for l in fh:
                p = l.rstrip("\n").split("\t")
                self.idx[p[0]] = tuple(int(x) for x in p[1:5])
        self.fh = open(path, "rb")

    def _key(self, chrom):
        """Resolve MAF chromosome naming against .fai keys (chr-prefix drift)."""
        if chrom in self.idx:
            return chrom
        if ("chr" + chrom) in self.idx:
            return "chr" + chrom
        if chrom.startswith("chr") and chrom[3:] in self.idx:
            return chrom[3:]
        return None

    def fetch(self, chrom, start, end):
        """Return uppercase sequence for 1-based inclusive [start, end]."""
        key = self._key(chrom)
        if key is None:
            return ""
        length, offset, lb, lw = self.idx[key]
        start, end = max(1, start), min(end, length)
        if start > end:
            return ""
        n = end - start + 1
        self.fh.seek(offset + (start - 1) // lb * lw + (start - 1) % lb)
        raw = self.fh.read(n + (n // lb + 2) * (lw - lb) + 16)
        return raw.decode("ascii", "ignore").replace("\n", "")[:n].upper()


def load_signatures(path):
    """Return (channel-ordered matrix [96 x n], signature names).

    Rows are matched to CHANNELS by their Type label; the matrix keeps the
    canonical CHANNELS order so count vectors index into it directly.
    """
    table = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        names = header[1:]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                table[parts[0]] = [float(x) for x in parts[1:]]
    missing = [c for c in CHANNELS if c not in table]
    if missing:
        sys.exit(f"[sbs96_fit] signature matrix missing channels: {missing[:3]}")
    mat = np.array([table[c] for c in CHANNELS], dtype=float)
    return mat, names


def sbs96_counts(maf_path, fa):
    """Count SBS96 channels from a MAF; return (counts, skipped, badref)."""
    counts = Counter()
    skipped = badref = 0
    with open(maf_path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    for r in csv.DictReader(lines, delimiter="\t"):
        if r.get("Variant_Type") != "SNP":
            continue
        ref = (r.get("Reference_Allele") or "").upper()
        alt = (r.get("Tumor_Seq_Allele2") or "").upper()
        if len(ref) != 1 or len(alt) != 1 or ref not in COMPLEMENT or alt not in COMPLEMENT:
            continue
        pos = int(r["Start_Position"])
        ctx = fa.fetch(r["Chromosome"], pos - 1, pos + 1)
        if len(ctx) != 3 or any(b not in COMPLEMENT for b in ctx):
            skipped += 1
            continue
        if ctx[1] != ref:
            badref += 1
            continue
        if ref in ("A", "G"):  # pyrimidine normalization: reverse-complement
            ctx = "".join(COMPLEMENT[b] for b in reversed(ctx))
            ref, alt = COMPLEMENT[ref], COMPLEMENT[alt]
        counts[f"{ctx[0]}[{ref}>{alt}]{ctx[2]}"] += 1
    return counts, skipped, badref


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maf", required=True)
    ap.add_argument("--genome-fa", required=True, help="reference FASTA (with .fai)")
    ap.add_argument("--signatures", required=True, help="COSMIC SBS96 matrix TSV")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--counts-out", required=True)
    ap.add_argument("--exposures-out", required=True)
    args = ap.parse_args()

    fa = Fasta(args.genome_fa)
    counts, skipped, badref = sbs96_counts(args.maf, fa)
    b = np.array([counts[c] for c in CHANNELS], dtype=float)
    total = int(b.sum())
    print(f"[{args.sample}] SNPs counted={total} skipped_nonACGT_context={skipped} "
          f"ref_mismatch={badref}")

    with open(args.counts_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Channel", "Count"])
        w.writerows(zip(CHANNELS, b.astype(int)))

    mat, names = load_signatures(args.signatures)
    if total == 0:
        sys.exit(f"[{args.sample}] no countable SNPs — exposures not fitted")
    x, _ = nnls(mat, b)
    recon = mat @ x
    cos = float(np.dot(b, recon) / (np.linalg.norm(b) * np.linalg.norm(recon)))

    order = np.argsort(-x)
    with open(args.exposures_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Signature", "Attributed_Mutations", "Exposure_Fraction"])
        for i in order:
            if x[i] <= 0:
                continue
            w.writerow([names[i], f"{x[i]:.1f}", f"{x[i] / recon.sum():.4f}"])
    print(f"[{args.sample}] signatures>0={int((x > 0).sum())} "
          f"cosine={cos:.4f} -> {args.exposures_out}")


if __name__ == "__main__":
    main()
