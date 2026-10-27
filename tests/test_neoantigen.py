#!/usr/bin/env python3
"""Tests for the neoantigen module scripts (venus#31 P0).

Covers: hla_parse (SpecHLA wide tables incl. G-group, audit M4), pvac_alleles
(supported-set intersection + min-allele QC), neoantigen_report (pVACtools
7.x REAL column names, burden math, expression filter, NA robustness,
fixed-schema honesty — audit B6/m3/m10).

Run: python3 tests/test_neoantigen.py   (stdlib only)
"""
import csv
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "scripts"))

import hla_parse  # noqa: E402
import neoantigen_report as nr  # noqa: E402
import pvac_alleles  # noqa: E402

SPEC_RESULT = (
    "SpecHLA v1.0.12\n"
    "Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\t"
    "HLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\t"
    "HLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2\n"
    "T1\tA*02:01:01\tA*11:01:01\tB*40:01:02\tB*58:01:01\tC*07:01:01\t"
    "C*03:04:01\tDPA1*01:03:01\tDPA1*02:01:02\tDPB1*04:01:01\tDPB1*14:01:01\t"
    "DQA1*01:02:01\tDQA1*05:05:01\tDQB1*06:01:01\tDQB1*03:01:01\t"
    "DRB1*13:02:01\tDRB1*04:03:01\n")

# audit M4: the G-group file is a WIDE table with the same column layout
SPEC_GGROUP = (
    "Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\t"
    "HLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\t"
    "HLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2\n"
    "T1\tA*02:01:01G\tA*11:01:01G\tB*40:01:01G\tB*58:01:01G\tC*07:01:01G\t"
    "C*03:04:01G\tDPA1*01:03:01G\tDPA1*02:01:02G\tDPB1*04:01:01G\t"
    "DPB1*14:01:01G\tDQA1*01:02:01G\tDQA1*05:05:01G\tDQB1*06:01:01G\t"
    "DQB1*03:01:01G\tDRB1*13:02:01G\tDRB1*04:03:01G\n")

# pVACtools 7.x filtered.tsv REAL column names (audit B6)
PVAC_7X_COLS = ("Gene Name\tVariant Type\tHLA Allele\tMT Epitope Seq\t"
                "Best MT IC50 Score\tBest MT Percentile\tTumor DNA VAF\t"
                "Mutation\tProtein Position\tTranscript\n")


def write_tsv(path, text):
    with open(path, "w") as fh:
        fh.write(text)


class HlaParse(unittest.TestCase):
    def test_wide_table_to_long(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt",
                                         delete=False) as fh:
            fh.write(SPEC_RESULT)
            path = fh.name
        res = hla_parse.parse_result_txt(path)
        self.assertEqual(res["A"][0], "A*02:01:01")
        self.assertEqual(hla_parse.two_field("A*02:01:01"), "HLA-A*02:01")
        self.assertEqual(hla_parse.two_field("HLA-A*02:01"), "HLA-A*02:01")
        self.assertEqual(hla_parse.two_field(""), "")
        os.unlink(path)

    def test_g_group_wide_table(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt",
                                         delete=False) as fh:
            fh.write(SPEC_GGROUP)
            path = fh.name
        ggrp = hla_parse.parse_g_group(path)
        self.assertEqual(ggrp["A"], ("A*02:01:01G", "A*11:01:01G"))
        self.assertEqual(ggrp["DRB1"][0], "DRB1*13:02:01G")
        os.unlink(path)

    def test_cli_writes_12_col_schema(self):
        with tempfile.TemporaryDirectory() as d:
            res = os.path.join(d, "hla.result.txt")
            gg = os.path.join(d, "hla.result.g.group.txt")
            write_tsv(res, SPEC_RESULT)
            write_tsv(gg, SPEC_GGROUP)
            out = os.path.join(d, "alleles.tsv")
            subprocess.run([sys.executable,
                            os.path.join(REPO, "scripts", "hla_parse.py"),
                            "--sample", "T1", "--result", res,
                            "--g-group", gg,
                            "--mode", "tumor_only", "--out", out], check=True)
            rows = list(csv.DictReader(open(out), delimiter="\t"))
            self.assertEqual(len(rows), 8)
            self.assertTrue(all(c in rows[0] for c in
                                ("sample", "gene", "allele_1", "allele_2",
                                 "allele_1_g_group", "typed", "mode",
                                 "caveat")))
            a = next(r for r in rows if r["gene"] == "A")
            self.assertEqual(a["allele_1"], "HLA-A*02:01")
            self.assertEqual(a["allele_1_g_group"], "A*02:01:01G")
            self.assertIn("LOH", a["caveat"])


class PvacAlleles(unittest.TestCase):
    def test_intersection_and_min_qc(self):
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, "a.tsv")
            write_tsv(tsv,
                      "sample\tgene\tallele_1\tallele_2\n"
                      "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n"
                      "T1\tB\tHLA-B*07:02\t\n")
            out = os.path.join(d, "alleles.txt")
            orig = pvac_alleles.mhcflurry_supported
            pvac_alleles.mhcflurry_supported = (
                lambda: {"HLA-A*02:01", "HLA-B*07:02"})
            try:
                sys.argv = ["pvac_alleles.py", "--alleles-tsv", tsv,
                            "--out", out]
                pvac_alleles.main()
            finally:
                pvac_alleles.mhcflurry_supported = orig
            with open(out) as fh:
                self.assertEqual(fh.read().strip(),
                                 "HLA-A*02:01,HLA-B*07:02")

            # collapsed typing must fail loud (audit m9 escape hatch documented)
            tsv2 = os.path.join(d, "b.tsv")
            write_tsv(tsv2,
                      "sample\tgene\tallele_1\tallele_2\n"
                      "T1\tA\tHLA-A*02:01\t\n")
            out2 = os.path.join(d, "alleles2.txt")
            sys.argv = ["pvac_alleles.py", "--alleles-tsv", tsv2,
                        "--min-alleles", "2", "--out", out2]
            pvac_alleles.mhcflurry_supported = (
                lambda: {"HLA-A*02:01", "HLA-B*07:02"})
            with self.assertRaises(SystemExit):
                pvac_alleles.main()


def run_report(tmp, pvac_rows, hla_rows, counts_lines=None, extra=None,
               usable=None):
    pvac = os.path.join(tmp, "filtered.tsv")
    write_tsv(pvac, PVAC_7X_COLS + pvac_rows)
    hla = os.path.join(tmp, "alleles.tsv")
    write_tsv(hla, "sample\tgene\tallele_1\tallele_2\n" + hla_rows)
    args = ["--pvac-tsv", pvac, "--hla-tsv", hla, "--patient", "T1",
            "--predictors", "MHCflurry",
            "--out-candidates", os.path.join(tmp, "c.tsv"),
            "--out-summary", os.path.join(tmp, "s.tsv"),
            "--out-md", os.path.join(tmp, "r.md")]
    if usable is not None:
        ua = os.path.join(tmp, "alleles.txt")
        write_tsv(ua, usable)
        args += ["--usable-alleles", ua]
    if counts_lines is not None:
        counts = os.path.join(tmp, "counts.txt")
        write_tsv(counts,
                  "# Geneid\tChr\tStart\tEnd\tStrand\tLength\tT1-TR\n"
                  + counts_lines)
        gmap = os.path.join(tmp, "gmap.tsv")
        write_tsv(gmap, "gene_id\tsymbol\nENSG1\tEGFR\nENSG2\tKRAS\n")
        args += ["--counts", counts, "--gene-map", gmap]
    if extra:
        args += extra
    proc = subprocess.run(
        [sys.executable,
         os.path.join(REPO, "scripts", "neoantigen_report.py"), *args],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    with open(os.path.join(tmp, "s.tsv")) as fh:
        summary = {(r["metric"], r["value"]) for r in
                   csv.DictReader(fh, delimiter="\t")}
    with open(os.path.join(tmp, "c.tsv")) as fh:
        cands = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(tmp, "r.md")) as fh:
        md = fh.read()
    return summary, cands, md


class NeoantigenReport(unittest.TestCase):
    def test_burden_and_expr_filter_with_7x_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = ("EGFR\tmissense_variant\tHLA-A*02:01\tGVLVGVAPE\t2.5\t0.1\t"
                    "0.4\tp.E12D\t12\tENST1\n"
                    "KRAS\tmissense_variant\tHLA-A*02:01\tVLVVGAVGV\t80\t1.2\t"
                    "0.3\tp.G12C\t12\tENST2\n"
                    "BRCA1\tframeshift_variant\tHLA-B*07:02\tYTYQTHTIL\t450\t"
                    "2\t0.2\tp.R1fs\t1\tENST3\n"
                    "ERBB2\tinframe_insertion\tHLA-B*07:02\tNLVLVGVA\t900\t5\t"
                    "0.1\tp.V2insI\t2\tENST4\n")
            counts = ("ENSG1\tchr7\t1\t100\t+\t2000\t5000\n"
                      "ENSG2\tchr12\t1\t100\t+\t2000\t0\n")
            summary, cands, md = run_report(
                tmp, rows, "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n",
                counts_lines=counts, usable="HLA-A*02:01,HLA-A*24:02")
            # strong (<=500): 2.5, 80, 450 -> 3; high-aff (<=100): 2.5, 80 -> 2
            self.assertIn(("burden_strong_500nM", "3"), summary)
            self.assertIn(("burden_highaff_100nM", "2"), summary)
            self.assertIn(("frameshift_neoantigens", "1"), summary)
            self.assertIn(("hla_alleles_typed", "2"), summary)
            self.assertIn(("hla_alleles_predictor_usable", "2"), summary)
            egfr = next(c for c in cands if c["gene"] == "EGFR")
            self.assertEqual(egfr["expr_pass"], "yes")
            kras = next(c for c in cands if c["gene"] == "KRAS")
            self.assertEqual(kras["expr_pass"], "no")

    def test_na_score_dedupe_robust(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = ("EGFR\tmissense_variant\tHLA-A*02:01\tGVLVGVAPE\tNA\tNA\t"
                    "0.4\tp.E12D\t12\tENST1\n"
                    "EGFR\tmissense_variant\tHLA-A*02:01\tGVLVGVAPE\t2\t0.1\t"
                    "0.4\tp.E12D\t12\tENST1\n")
            summary, cands, _ = run_report(
                tmp, rows, "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
            self.assertEqual(len(cands), 1)
            self.assertEqual(cands[0]["mt_score_nm"], "2")

    def test_empty_candidates_fixed_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary, cands, md = run_report(
                tmp, "", "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
            self.assertIn(("burden_strong_500nM", "0"), summary)
            self.assertEqual(len(cands), 0)

    def test_tumor_only_caveat_rendered(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, md = run_report(
                tmp, "", "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
            self.assertIn("germline", md)
            self.assertIn("not run in P0", md)
            self.assertNotIn("5-14%", md)  # audit M6: unsourced number gone


if __name__ == "__main__":
    unittest.main(verbosity=2)
