#!/usr/bin/env python3
"""Offline tests for sbs96_fit (issue #33 N9).

Covers the .fai random-access FASTA reader (incl. chr-prefix drift and
contig-edge clamping), pyrimidine normalization (A/G ref reverse-complement),
ref-mismatch / non-ACGT-context skip counters, signature-matrix loading
(real COSMIC v3.1 file), and a full CLI run where NNLS must recover two
orthogonal synthetic signatures exactly.

Run: python3 tests/test_sbs96_fit.py   (needs numpy+scipy, no network)
"""
import csv
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "scripts"))

import sbs96_fit as sf  # noqa: E402

# 200 bases, period 10: position p (1-based) has SEQ[(p-1) % 10]
# CTTCAGGATC -> 1=C 2=T 3=T 4=C 5=A 6=G 7=G 8=A 9=T 10=C (then repeats)
SEQ = "CTTCAGGATC" * 20
# contig carrying non-ACGT bases (tests the skipped counter)
SEQ_N = "ACGTACGTNN"


def write_fa(dirpath, seqs):
    """Write a FASTA (60 bases/line) + matching .fai; return the fa path."""
    fa = os.path.join(dirpath, "ref.fa")
    entries = []
    with open(fa, "w") as fh:
        for chrom, seq in seqs.items():
            fh.write(f">{chrom}\n")
            offset = fh.tell()  # .fai offset = first base AFTER the header
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
            entries.append((chrom, len(seq), offset))
    with open(fa + ".fai", "w") as fh:
        for chrom, length, offset in entries:
            fh.write(f"{chrom}\t{length}\t{offset}\t60\t61\n")
    return fa


MAF_COLS = ("Hugo_Symbol\tChromosome\tStart_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\tVariant_Type\n")


def maf_row(gene, chrom, pos, ref, alt, vtype="SNP"):
    return f"{gene}\t{chrom}\t{pos}\t{ref}\t{alt}\t{vtype}\n"


def write_maf(path, rows):
    with open(path, "w") as fh:
        fh.write("## venus MAF\n")
        fh.write(MAF_COLS)
        fh.writelines(rows)


class FastaFetch(unittest.TestCase):
    def test_fetch_line_crossing_and_clamps(self):
        with tempfile.TemporaryDirectory() as td:
            fa = write_fa(td, {"chr1": SEQ})
            fao = sf.Fasta(fa)
            self.assertEqual(fao.fetch("chr1", 1, 3), "CTT")
            # spans the 60-base line break: pos 59..62 = T C C T
            self.assertEqual(fao.fetch("chr1", 59, 62), "TCCT")
            # 1-based inclusive single base
            self.assertEqual(fao.fetch("chr1", 5, 5), "A")
            # end clamped to contig length; start beyond -> empty
            self.assertEqual(fao.fetch("chr1", len(SEQ), len(SEQ) + 50), "C")
            self.assertEqual(fao.fetch("chr1", len(SEQ) + 10, len(SEQ) + 60),
                             "")
            # unknown contig -> empty string, no exception
            self.assertEqual(fao.fetch("chrUn", 1, 3), "")

    def test_chr_prefix_drift_both_directions(self):
        with tempfile.TemporaryDirectory() as td:
            # index has "1", query uses "chr1" -> resolved
            fa_plain = write_fa(td, {"1": SEQ})
            self.assertEqual(sf.Fasta(fa_plain).fetch("chr1", 1, 3), "CTT")
        with tempfile.TemporaryDirectory() as td:
            # index has "chr2", query uses "2" -> resolved
            fa_chr = write_fa(td, {"chr2": SEQ})
            self.assertEqual(sf.Fasta(fa_chr).fetch("2", 1, 3), "CTT")


class Signatures(unittest.TestCase):
    COSMIC = os.path.join(REPO, "resources", "cosmic_v3.1_sbs96_grch38.txt")

    def test_load_real_cosmic_matrix(self):
        mat, names = sf.load_signatures(self.COSMIC)
        self.assertEqual(mat.shape, (96, 72))
        self.assertEqual(names[0], "SBS1")
        # signature columns are probability vectors summing to 1
        self.assertTrue(np.allclose(mat.sum(axis=0), 1.0, atol=1e-4))

    def test_missing_channel_exits(self):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, "sig.tsv")
            with open(path, "w") as fh:
                fh.write("Type\tSigX\n")
                for c in sf.CHANNELS[:-1]:  # drop the last channel
                    fh.write(f"{c}\t0.5\n")
            with self.assertRaises(SystemExit):
                sf.load_signatures(path)


class Counts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.fa = sf.Fasta(write_fa(self.tmp.name, {"chr1": SEQ,
                                                    "chr2": SEQ_N}))

    def tearDown(self):
        self.tmp.cleanup()

    def test_channels_normalization_and_counters(self):
        rows = [
            maf_row("EGFR", "1", 2, "T", "A"),    # ctx CTT  -> C[T>A]T
            maf_row("KRAS", "1", 3, "T", "C"),    # ctx TTC  -> T[T>C]C
            maf_row("BRAF", "1", 5, "A", "G"),    # ctx CAG  -> rc CTG -> C[T>C]G
            maf_row("MET", "1", 8, "A", "T"),     # ctx GAT  -> rc ATC -> A[T>A]C
            maf_row("TP53", "1", 10, "C", "T"),   # ctx TCC  -> T[C>T]C
            maf_row("ERBB2", "1", 11, "C", "A"),  # ctx CCT  -> C[C>A]T
            maf_row("BAD", "1", 4, "A", "G"),     # MAF ref A != seq C -> badref
            maf_row("NSD1", "2", 8, "T", "A"),    # ctx GTN -> non-ACGT skip
            maf_row("EDGE", "1", 1, "C", "T"),    # contig edge -> 2-base ctx skip
            maf_row("INDEL", "1", 4, "CT", "-", "DEL"),  # not a SNP
        ]
        maf = os.path.join(self.tmp.name, "x.maf")
        write_maf(maf, rows)
        counts, skipped, badref = sf.sbs96_counts(maf, self.fa)
        self.assertEqual(badref, 1)
        self.assertEqual(skipped, 2)
        self.assertEqual(counts["C[T>A]T"], 1)
        self.assertEqual(counts["T[T>C]C"], 1)
        self.assertEqual(counts["C[T>C]G"], 1)
        self.assertEqual(counts["A[T>A]C"], 1)
        self.assertEqual(counts["T[C>T]C"], 1)
        self.assertEqual(counts["C[C>A]T"], 1)
        self.assertEqual(sum(counts.values()), 6)

    def test_non_snp_and_multibase_ignored(self):
        maf = os.path.join(self.tmp.name, "x.maf")
        write_maf(maf, [
            maf_row("G1", "1", 4, "CT", "-", "DEL"),
            maf_row("G2", "1", 4, "C", "G", "INS"),
        ])
        counts, skipped, badref = sf.sbs96_counts(maf, self.fa)
        self.assertEqual((sum(counts.values()), skipped, badref), (0, 0, 0))


def write_synth_matrix(path):
    """Two orthogonal unit-vector signatures over the full 96 channels."""
    with open(path, "w") as fh:
        fh.write("Type\tSigA\tSigB\n")
        for c in sf.CHANNELS:
            a = 1.0 if c == "T[C>T]A" else 0.0
            b = 1.0 if c == "T[T>C]C" else 0.0
            fh.write(f"{c}\t{a}\t{b}\n")


class Cli(unittest.TestCase):
    def test_exact_recovery_of_two_signatures(self):
        with tempfile.TemporaryDirectory() as td:
            fa = write_fa(td, {"chr1": SEQ})
            maf = os.path.join(td, "x.maf")
            # pos 4 = C (ctx TCA -> T[C>T]A), pos 3 = T (ctx TTC -> T[T>C]C)
            write_maf(maf, [maf_row("G", "1", 4, "C", "T")] * 5
                      + [maf_row("H", "1", 3, "T", "C")] * 3)
            sig = os.path.join(td, "sig.tsv")
            write_synth_matrix(sig)
            counts_out = os.path.join(td, "counts.tsv")
            expo_out = os.path.join(td, "expo.tsv")
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts", "sbs96_fit.py"),
                 "--maf", maf, "--genome-fa", fa, "--signatures", sig,
                 "--sample", "S1", "--counts-out", counts_out,
                 "--exposures-out", expo_out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            self.assertIn("SNPs counted=8", proc.stdout)
            self.assertIn("cosine=1.0000", proc.stdout)
            with open(counts_out) as fh:
                crows = {r["Channel"]: r["Count"] for r in
                         csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(crows["T[C>T]A"], "5")
            self.assertEqual(crows["T[T>C]C"], "3")
            self.assertEqual(len(crows), 96)
            with open(expo_out) as fh:
                erows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual([r["Signature"] for r in erows],
                             ["SigA", "SigB"])  # sorted desc by exposure
            self.assertEqual(erows[0]["Attributed_Mutations"], "5.0")
            self.assertEqual(erows[0]["Exposure_Fraction"], "0.6250")
            self.assertEqual(erows[1]["Attributed_Mutations"], "3.0")
            self.assertEqual(erows[1]["Exposure_Fraction"], "0.3750")

    def test_zero_snp_maf_fails_loud(self):
        with tempfile.TemporaryDirectory() as td:
            fa = write_fa(td, {"chr1": SEQ})
            maf = os.path.join(td, "x.maf")
            write_maf(maf, [])  # header only
            sig = os.path.join(td, "sig.tsv")
            write_synth_matrix(sig)
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts", "sbs96_fit.py"),
                 "--maf", maf, "--genome-fa", fa, "--signatures", sig,
                 "--sample", "S1", "--counts-out", os.path.join(td, "c.tsv"),
                 "--exposures-out", os.path.join(td, "e.tsv")],
                capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("no countable SNPs", proc.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
