#!/usr/bin/env python3
"""Tests for the neoantigen module scripts (venus#31 P0).

Covers: hla_parse (SpecHLA wide table -> canonical alleles), pvac_alleles
(supported-set intersection + min-allele QC), neoantigen_report (pVACseq
parsing, burden math, TESLA TPM filter, fixed-schema honesty).

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

SPEC_RESULT = (
    "SpecHLA v1.0.12\n"
    "Sample\tHLA_A_1\tHLA_A_2\tHLA_B_1\tHLA_B_2\tHLA_C_1\tHLA_C_2\t"
    "HLA_DPA1_1\tHLA_DPA1_2\tHLA_DPB1_1\tHLA_DPB1_2\tHLA_DQA1_1\tHLA_DQA1_2\t"
    "HLA_DQB1_1\tHLA_DQB1_2\tHLA_DRB1_1\tHLA_DRB1_2\n"
    "T1\tA*02:01:01\tA*11:01:01\tB*40:01:02\tB*58:01:01\tC*07:01:01\t"
    "C*03:04:01\tDPA1*01:03:01\tDPA1*02:01:02\tDPB1*04:01:01\tDPB1*14:01:01\t"
    "DQA1*01:02:01\tDQA1*05:05:01\tDQB1*06:01:01\tDQB1*03:01:01\t"
    "DRB1*13:02:01\tDRB1*04:03:01\n")


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

    def test_cli_writes_12_col_schema(self):
        with tempfile.NamedTemporaryFile("w", suffix=".txt",
                                         delete=False) as fh:
            fh.write(SPEC_RESULT)
            path = fh.name
        out = path + ".alleles.tsv"
        subprocess.run([sys.executable, os.path.join(REPO, "scripts",
                                                     "hla_parse.py"),
                        "--sample", "T1", "--result", path,
                        "--mode", "tumor_only", "--out", out], check=True)
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(c in rows[0] for c in
                            ("sample", "gene", "allele_1", "allele_2",
                             "typed", "mode", "caveat")))
        a = next(r for r in rows if r["gene"] == "A")
        self.assertEqual(a["allele_1"], "HLA-A*02:01")
        self.assertIn("LOH", a["caveat"])
        os.unlink(path)
        os.unlink(out)


class PvacAlleles(unittest.TestCase):
    def test_intersection_and_min_qc(self):
        script = os.path.join(REPO, "scripts", "pvac_alleles.py")
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, "a.tsv")
            with open(tsv, "w") as fh:
                fh.write("sample\tgene\tallele_1\tallele_2\n")
                fh.write("T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
                fh.write("T1\tB\tHLA-B*07:02\t\n")
            # fake supported list: only A*02:01 and B*07:02 are supported
            out = os.path.join(d, "alleles.txt")
            # patch the supported() helper instead of the CLI
            import pvac_alleles
            orig = pvac_alleles.mhcflurry_supported
            pvac_alleles.mhcflurry_supported = (
                lambda: {"HLA-A*02:01", "HLA-B*07:02"})
            try:
                sys.argv = ["pvac_alleles.py", "--alleles-tsv", tsv,
                            "--out", out]
                pvac_alleles.main()
            finally:
                pvac_alleles.mhcflurry_supported = orig
            self.assertEqual(open(out).read().strip(),
                             "HLA-A*02:01,HLA-B*07:02")

            # an HLA typing that collapses below min-alleles must fail loud
            tsv2 = os.path.join(d, "b.tsv")
            with open(tsv2, "w") as fh:
                fh.write("sample\tgene\tallele_1\tallele_2\n")
                fh.write("T1\tA\tHLA-A*02:01\t\n")
            out2 = os.path.join(d, "alleles2.txt")
            sys.argv = ["pvac_alleles.py", "--alleles-tsv", tsv2,
                        "--min-alleles", "2", "--out", out2]
            with self.assertRaises(SystemExit):
                pvac_alleles.mhcflurry_supported = (
                    lambda: {"HLA-A*02:01", "HLA-B*07:02"})
                pvac_alleles.main()


PVAC_FIXTURE_COLS = ("Gene\tVariant Type\tHLA Allele\tMT Epitope Seq\t"
                     "MT Score\tMT Percentile\tPresentation Score\t"
                     "Tumor DNA VAF\tTranscript\n")


def run_report(tmp, pvac_rows, hla_rows, counts_lines=None, extra=None):
    pvac = os.path.join(tmp, "filtered.tsv")
    with open(pvac, "w") as fh:
        fh.write(PVAC_FIXTURE_COLS)
        fh.write(pvac_rows)
    hla = os.path.join(tmp, "alleles.tsv")
    with open(hla, "w") as fh:
        fh.write("sample\tgene\tallele_1\tallele_2\n")
        fh.write(hla_rows)
    args = ["--pvac-tsv", pvac, "--hla-tsv", hla, "--patient", "T1",
            "--predictors", "MHCflurry",
            "--out-candidates", os.path.join(tmp, "c.tsv"),
            "--out-summary", os.path.join(tmp, "s.tsv"),
            "--out-md", os.path.join(tmp, "r.md")]
    if counts_lines is not None:
        counts = os.path.join(tmp, "counts.txt")
        with open(counts, "w") as fh:
            fh.write("# Geneid\tChr\tStart\tEnd\tStrand\tLength\tT1-TR\n")
            fh.write(counts_lines)
        gmap = os.path.join(tmp, "gmap.tsv")
        with open(gmap, "w") as fh:
            fh.write("gene_id\tsymbol\n")
            fh.write("ENSG1\tEGFR\n")
            fh.write("ENSG2\tKRAS\n")
        args += ["--counts", counts, "--gene-map", gmap]
    if extra:
        args += extra
    subprocess.run([sys.executable,
                    os.path.join(REPO, "scripts", "neoantigen_report.py"),
                    *args], check=True, capture_output=True)
    with open(os.path.join(tmp, "s.tsv")) as fh:
        summary = {(r["metric"], r["value"]) for r in
                   csv.DictReader(fh, delimiter="\t")}
    with open(os.path.join(tmp, "c.tsv")) as fh:
        cands = list(csv.DictReader(fh, delimiter="\t"))
    md = open(os.path.join(tmp, "r.md")).read()
    return summary, cands, md


class NeoantigenReport(unittest.TestCase):
    def test_burden_and_tesla_filter(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = ("EGFR\tmissense\tHLA-A*02:01\tGVLVGVAPE\t2.5\t0.1\t0.9\t"
                    "0.4\tENST1\n"
                    "KRAS\tmissense\tHLA-A*02:01\tVLVVGAVGV\t80\t1.2\t0.5\t"
                    "0.3\tENST2\n"
                    "BRCA1\tframeshift\tHLA-B*07:02\tYTYQTHTIL\t450\t2\t0.3\t"
                    "0.2\tENST3\n"
                    "ERBB2\tinframe\tHLA-B*07:02\tNLVLVGVA\t900\t5\t0.1\t"
                    "0.1\tENST4\n")
            counts = ("ENSG1\tchr7\t1\t100\t+\t2000\t5000\n"
                      "ENSG2\tchr12\t1\t100\t+\t2000\t0\n")
            summary, cands, md = run_report(
                tmp, rows, "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n",
                counts_lines=counts)
            # strong (<=500): 2.5, 80, 450 -> 3; tesla (<=100): 2.5, 80 -> 2
            self.assertIn(("burden_strong_500nM", "3"), summary)
            self.assertIn(("burden_tesla_100nM", "2"), summary)
            self.assertIn(("frameshift_neoantigens", "1"), summary)
            # TPM: EGFR high (pass), KRAS low (fail), BRCA1 absent
            egfr = next(c for c in cands if c["gene"] == "EGFR")
            self.assertEqual(egfr["tesla_tpm_pass"], "yes")
            kras = next(c for c in cands if c["gene"] == "KRAS")
            self.assertEqual(kras["tesla_tpm_pass"], "no")
            self.assertIn(("hla_alleles_typed", "2"), summary)

    def test_dedup_best_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = ("EGFR\tmissense\tHLA-A*02:01\tGVLVGVAPE\t9\t1\t0.9\t"
                    "0.4\tENST1\n"
                    "EGFR\tmissense\tHLA-A*02:01\tGVLVGVAPE\t2\t1\t0.9\t"
                    "0.4\tENST1\n")
            summary, cands, _ = run_report(
                tmp, rows, "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
            self.assertEqual(len(cands), 1)
            self.assertEqual(cands[0]["mt_score_nm"], "2")

    def test_empty_candidates_fixed_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary, cands, md = run_report(
                tmp, "", "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
            metrics = dict(summary)
            self.assertIn(("burden_strong_500nM", "0"), summary)
            self.assertEqual(len(cands), 0)

    def test_tumor_only_caveat_rendered(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, md = run_report(
                tmp, "", "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n")
            self.assertIn("5-14%", md)
            self.assertIn("not run in P0", md)  # class-II honesty note


if __name__ == "__main__":
    unittest.main(verbosity=2)
