#!/usr/bin/env python3
"""Tests for the WGS alignment-QC collector (N10/N11, audit 20261104).

Covers: flagstat parsing (modern + legacy pct), samtools stats SN/GC parsing,
mosdepth summary/thresholds parsing (plain + gzip), autosome-normalised sex
calls with disclosed thresholds, suffix routing (order-independent), empty
input -> header-only table.

Run: python3 tests/test_alignqc.py   (stdlib only)
"""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "scripts"))

import collect_alignqc as aq  # noqa: E402

FLAGSTAT = (
    "20000 + 0 in total (QC-passed reads + QC-failed reads)\n"
    "0 + 0 secondary\n"
    "0 + 0 supplementary\n"
    "0 + 0 duplicates\n"
    "19800 + 0 mapped (99.00% : N/A)\n"
    "20000 + 0 primary\n")

STATS = (
    "This file was produced by samtools stats (1.21+htslib-1.21)\n"
    "SN\traw total sequences:\t20000\n"
    "SN\tinsert size: average\t345.67\n"
    "SN\tinsert size: median\t330.0\n"
    "SN\tinsert size: standard deviation\t17.89\n"
    "GC\t10.0\t100\n"
    "GC\t20.0\t900\n")

SUMMARY = (
    "chrom\tlength\tbases\tmean\tmin\tmax\n"
    "chr1\t100000\t150000\t1.50\t0\t30\n"
    "chrX\t10000\t8000\t0.40\t0\t12\n"
    "chrY\t5000\t300\t0.03\t0\t9\n"
    "total\t115000\t158300\t1.38\t0\t30\n")

THRESH = ("# chrom\tstart\tend\tregion\t10X\t20X\n"
          "chr1\t0\t100000\tunknown\t10000\t2000\n"
          "chrX\t0\t10000\tunknown\t100\t10\n")


def write_set(d, sample="Pt01-TD", gz_thresholds=False):
    """Write the four QC fixtures into dir d; return the four paths."""
    def p(ext):
        return os.path.join(d, sample + ext)
    for ext, payload in ((".flagstat", FLAGSTAT), (".stats", STATS),
                         (".mosdepth.summary.txt", SUMMARY)):
        with open(p(ext), "w") as fh:
            fh.write(payload)
    thresh_path = p(".mosdepth.thresholds.bed.gz" if gz_thresholds
                    else ".mosdepth.thresholds.bed")
    if gz_thresholds:
        with gzip.open(thresh_path, "wt") as fh:
            fh.write(THRESH)
    else:
        with open(thresh_path, "w") as fh:
            fh.write(THRESH)
    return [p(".flagstat"), p(".stats"), p(".mosdepth.summary.txt"),
            thresh_path]


def run_cli(paths, out):
    proc = subprocess.run(
        [sys.executable, os.path.join(REPO, "scripts", "collect_alignqc.py"),
         "--inputs", *paths, "--output", out],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return proc.stderr


class Flagstat(unittest.TestCase):
    def test_modern_pct(self):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "S1.flagstat")
            with open(f, "w") as fh:
                fh.write(FLAGSTAT)
            row = aq.parse_flagstat(f)
        self.assertEqual(row["total_reads"], "20000")
        self.assertEqual(row["mapped_pct"], "99.00")
        self.assertEqual(row["duplicate_pct"], "0.0000")

    def test_legacy_pct_computed(self):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "S1.flagstat")
            with open(f, "w") as fh:
                fh.write("10 + 0 in total ...\n5 + 0 mapped\n")
            row = aq.parse_flagstat(f)
        self.assertEqual(row["mapped_pct"], "50.0000")


class Stats(unittest.TestCase):
    def test_sn_and_gc(self):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "S1.stats")
            with open(f, "w") as fh:
                fh.write(STATS)
            row = aq.parse_stats(f)
        self.assertEqual(row["insert_size_median"], "330.0")
        self.assertEqual(row["insert_size_average"], "345.67")
        self.assertEqual(row["insert_size_stdev"], "17.89")
        # weighted: (10*100 + 20*900) / 1000
        self.assertEqual(row["gc_pct"], "19.00")


class SexCall(unittest.TestCase):
    def test_male(self):
        self.assertEqual(aq.call_sex(0.5, 0.3)[0], "M")

    def test_female(self):
        self.assertEqual(aq.call_sex(1.0, 0.001)[0], "F")

    def test_ambiguous_bands(self):
        # X between the M/F bands -> ambiguous regardless of Y
        self.assertEqual(aq.call_sex(0.8, 0.001)[0], "ambiguous")
        # X low but Y below the male floor -> ambiguous (contamination class)
        self.assertEqual(aq.call_sex(0.4, 0.02)[0], "ambiguous")

    def test_missing_contigs(self):
        self.assertEqual(aq.call_sex(0.9, None)[0], "F(no_chrY_contig)")
        self.assertEqual(aq.call_sex(0.5, None)[0], "ambiguous(no_chrY)")
        self.assertEqual(aq.call_sex(None, None)[0], "ambiguous(no_chrX)")


class Mosdepth(unittest.TestCase):
    def test_summary_autosome_mean_excludes_xy(self):
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "S1.mosdepth.summary.txt")
            with open(f, "w") as fh:
                fh.write(SUMMARY)
            rows = aq.parse_mosdepth_summary(f)
        am = aq.autosome_mean(rows)
        self.assertAlmostEqual(am, 1.5)  # only chr1 counts

    def test_thresholds_plain_and_gzip(self):
        for gz in (False, True):
            with tempfile.TemporaryDirectory() as d:
                paths = write_set(d, gz_thresholds=gz)
                thr = [p for p in paths if "thresholds" in p][0]
                row = aq.parse_thresholds(thr)
                # (10000+100)/110000 and (2000+10)/110000
                self.assertAlmostEqual(float(row["coverage_10x_pct"]),
                                       9.1818, places=3)
                self.assertAlmostEqual(float(row["coverage_20x_pct"]),
                                       1.8273, places=3)

    def test_route_and_sample_name(self):
        self.assertEqual(aq.route("a/S1.mosdepth.thresholds.bed.gz"),
                         ".mosdepth.thresholds.bed.gz")
        self.assertEqual(aq.sample_of("a/Pt01-B.mosdepth.summary.txt"),
                         "Pt01-B")
        self.assertEqual(aq.sample_of("a/Pt01-TD.flagstat"), "Pt01-TD")


class Cli(unittest.TestCase):
    def test_full_row_and_order_independence(self):
        with tempfile.TemporaryDirectory() as d:
            paths = write_set(d)
            out = os.path.join(d, "summary.tsv")
            run_cli(list(reversed(paths)), out)  # shuffled order
            import csv
            rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["sample"], "Pt01-TD")
        self.assertEqual(r["mapped_pct"], "99.00")
        self.assertEqual(r["insert_size_median"], "330.0")
        self.assertEqual(r["gc_pct"], "19.00")
        self.assertEqual(r["mean_depth"], "1.38")
        self.assertEqual(r["predicted_sex"], "ambiguous")  # 0.267/0.02 fixture
        self.assertEqual(r["contamination_freemix"],
                         "not_assessed(no_genotype_panel_p2)")

    def test_empty_inputs_header_only(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "summary.tsv")
            run_cli([], out)
            with open(out) as fh:
                header = fh.readline().strip().split("\t")
            self.assertIn("predicted_sex", header)
            self.assertIn("contamination_freemix", header)

    def test_unknown_suffix_fails_loud(self):
        with tempfile.TemporaryDirectory() as d:
            junk = os.path.join(d, "S1.unknown")
            with open(junk, "w") as fh:
                fh.write("x")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "collect_alignqc.py"),
                 "--inputs", junk, "--output", os.path.join(d, "o.tsv")],
                capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("unrecognised", proc.stderr + proc.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
