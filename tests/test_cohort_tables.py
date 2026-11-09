#!/usr/bin/env python3
"""Offline tests for cohort_tables (issue #33 N9).

Covers basename routing (order-independent --inputs), patient-level
normalization by assay-suffix stripping (-TD/-B/-TR), per-module column
prefixing with sample_id provenance, optional module sources
(wgs/ascat/ecdna/msi), TMB concatenation over multiple per-pair tables, and
the missing rna/sc exit. Sample-id shapes mirror the real pipeline: pair ids
are bare (Pt01), DNA/ecDNA rows carry the tumor sample id (Pt01-TD), bulk
RNA rows carry Pt01-TR, scRNA rows carry the bare patient id.

Run: python3 tests/test_cohort_tables.py   (stdlib only)
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCRIPT = os.path.join(REPO, "scripts", "cohort_tables.py")


def write_tsv(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return path


def run_cohort(tmp, inputs, extra=()):
    out = os.path.join(tmp, "cohort.tsv")
    proc = subprocess.run(
        [sys.executable, SCRIPT, "--inputs", *inputs, "--output", out,
         *extra],
        capture_output=True, text=True)
    return proc, out


class CohortTables(unittest.TestCase):
    def test_full_merge_suffix_stripping_and_prefixing(self):
        with tempfile.TemporaryDirectory() as td:
            rna = write_tsv(os.path.join(td, "rna_qc_summary.tsv"),
                            "sample\tpct_dup\tmapped_pct\n"
                            "Pt01-TR\t55.1\t99.2\n"
                            "Pt02-TR\t40.0\t98.0\n")
            sc = write_tsv(os.path.join(td, "scrna_metrics.tsv"),
                           "sample\tcells\tdoublet_pct\n"
                           "Pt01\t8100\t3.1\n")
            tmb1 = write_tsv(os.path.join(td, "Pt01.tmb.tsv"),
                             "sample\ttmb_mut_per_mb\tcoding_non_silent\n"
                             "Pt01\t6.10\t88\n")
            tmb2 = write_tsv(os.path.join(td, "Pt02.tmb.tsv"),
                             "sample\ttmb_mut_per_mb\tcoding_non_silent\n"
                             "Pt02\t2.50\t31\n")
            ascat = write_tsv(os.path.join(td, "cohort_ascat.tsv"),
                              "pair\tsample\tpurity\tploidy\n"
                              "Pt01\tPt01-TD\t0.85\t2.9\n")
            msi = write_tsv(os.path.join(td, "msi_cohort.tsv"),
                            "sample\tbh_pct\tmsi_status_wgs\n"
                            "Pt01\t1.00\tMSS\n")
            ecdna = write_tsv(os.path.join(td, "cohort_ecdna.tsv"),
                              "sample\tecDNA_count\nPt01-TD\t3\n")
            wgs = write_tsv(os.path.join(td, "wgs_alignqc_summary.tsv"),
                            "sample\tmean_depth\tpredicted_sex\n"
                            "Pt01-TD\t38.2\tM\n")
            # deliberately shuffled argument order: routing is by basename
            proc, out = run_cohort(td, [msi, tmb2, ecdna, wgs, sc, ascat,
                                        tmb1, rna])
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                rows = {ln.split("\t")[0]: dict(zip(header, ln.rstrip("\n").split("\t")))
                        for ln in fh}
            self.assertEqual(sorted(rows), ["Pt01", "Pt02"])
            # provenance: original per-module ids preserved as <mod>sample_id
            p1 = rows["Pt01"]
            self.assertEqual(p1["rna_sample_id"], "Pt01-TR")
            self.assertEqual(p1["sc_sample_id"], "Pt01")
            self.assertEqual(p1["tmb_sample_id"], "Pt01")
            self.assertEqual(p1["ascat_sample_id"], "Pt01-TD")
            self.assertEqual(p1["msi_sample_id"], "Pt01")
            self.assertEqual(p1["ecdna_sample_id"], "Pt01-TD")
            self.assertEqual(p1["wgs_sample_id"], "Pt01-TD")
            # suffixed rows (-TD/-TR) all land on the same patient row
            self.assertEqual(p1["sample"], "Pt01")
            self.assertEqual(p1["rna_pct_dup"], "55.1")
            self.assertEqual(p1["sc_cells"], "8100")
            self.assertEqual(p1["tmb_tmb_mut_per_mb"], "6.10")
            self.assertEqual(p1["ascat_purity"], "0.85")
            self.assertEqual(p1["msi_bh_pct"], "1.00")
            self.assertEqual(p1["ecdna_ecDNA_count"], "3")
            self.assertEqual(p1["wgs_mean_depth"], "38.2")
            # Pt02 lacks sc/ascat/ecdna/wgs rows -> those columns empty
            p2 = rows["Pt02"]
            self.assertEqual(p2["tmb_tmb_mut_per_mb"], "2.50")
            self.assertEqual(p2["ascat_purity"], "")
            self.assertEqual(p2["sc_cells"], "")
            self.assertEqual(p2["ecdna_ecDNA_count"], "")
            # Pt02-TR RNA row also folded into Pt02
            self.assertEqual(p2["rna_pct_dup"], "40.0")

    def test_custom_suffixes_flag(self):
        with tempfile.TemporaryDirectory() as td:
            rna = write_tsv(os.path.join(td, "rna_qc_summary.tsv"),
                            "sample\tmapped_pct\nPt01-TR\t99.0\n")
            sc = write_tsv(os.path.join(td, "scrna_metrics.tsv"),
                           "sample\tcells\nPt01\t7000\n")
            # default suffixes strip -TR already: one patient row
            proc, out = run_cohort(td, [rna, sc])
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out) as fh:
                body = fh.read().splitlines()
            self.assertEqual(len(body), 2)  # header + one Pt01 row
            # --suffixes replaces the default list: only -X stripped now, so
            # the -TR row becomes its own patient
            td2 = tempfile.mkdtemp(dir=td)
            rna2 = write_tsv(os.path.join(td2, "rna_qc_summary.tsv"),
                             "sample\tmapped_pct\nPt01-TR\t99.0\n")
            sc2 = write_tsv(os.path.join(td2, "scrna_metrics.tsv"),
                            "sample\tcells\nPt01\t7000\n")
            proc2, out2 = run_cohort(td, [rna2, sc2],
                                     extra=["--suffixes", "X"])
            self.assertEqual(proc2.returncode, 0, proc2.stderr)
            with open(out2) as fh:
                body2 = fh.read().splitlines()
            self.assertEqual(len(body2), 3)  # header + Pt01-TR + Pt01 rows

    def test_missing_rna_or_sc_exits_1(self):
        with tempfile.TemporaryDirectory() as td:
            sc = write_tsv(os.path.join(td, "scrna_metrics.tsv"),
                           "sample\tcells\nPt01\t7000\n")
            tmb = write_tsv(os.path.join(td, "x.tmb.tsv"),
                            "sample\ttmb\nPt01\t5\n")
            proc, _ = run_cohort(td, [sc, tmb])  # no rna
            self.assertEqual(proc.returncode, 1)
            self.assertIn("need at least RNA QC + scRNA metrics", proc.stderr)
            rna = write_tsv(os.path.join(td, "rna_qc_summary.tsv"),
                            "sample\tmapped_pct\nPt01\t99\n")
            proc2, _ = run_cohort(td, [rna, tmb])  # no sc
            self.assertEqual(proc2.returncode, 1)
            self.assertIn("need at least RNA QC + scRNA metrics", proc2.stderr)

    def test_unknown_basename_routes_to_tmb(self):
        with tempfile.TemporaryDirectory() as td:
            rna = write_tsv(os.path.join(td, "rna_qc_summary.tsv"),
                            "sample\tmapped_pct\nPt01\t99\n")
            sc = write_tsv(os.path.join(td, "scrna_metrics.tsv"),
                           "sample\tcells\nPt01\t7000\n")
            odd = write_tsv(os.path.join(td, "Pt01.some.unknown.table.tsv"),
                            "sample\ttotal_pass_variants\nPt01\t120\n")
            proc, out = run_cohort(td, [rna, sc, odd])
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                row = dict(zip(header, fh.readline().rstrip("\n").split("\t")))
            self.assertIn("tmb_total_pass_variants", header)
            self.assertEqual(row["tmb_total_pass_variants"], "120")
            self.assertEqual(row["tmb_sample_id"], "Pt01")


if __name__ == "__main__":
    unittest.main(verbosity=2)
