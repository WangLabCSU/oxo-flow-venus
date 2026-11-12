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

import collect_neoantigen  # noqa: E402
import hla_parse  # noqa: E402
import icb_linkage as il  # noqa: E402
import manta_to_arriba as mta  # noqa: E402
import msi_frameshift_join as mfj  # noqa: E402
import neoantigen_report as nr  # noqa: E402
import pvac_alleles  # noqa: E402
import tcr_report  # noqa: E402

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
        # real featureCounts layout: '# Program:...' line THEN the Geneid header
        write_tsv(counts,
                  "# Program:featureCounts v2.1.1; Command:\"featureCounts\"\n"
                  "Geneid\tChr\tStart\tEnd\tStrand\tLength\tT1-TR\n"
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


class CollectNeoantigenTest(unittest.TestCase):
    """venus#34: per-sample neoantigen summaries -> one wide cohort table."""

    SUMMARY_A = (
        "class\tmetric\tvalue\tdetail\n"
        "neoantigen\thla_alleles_typed\t13\tHLA-A*02:07;HLA-A*11:01\n"
        "neoantigen\thla_alleles_predictor_usable\t6\t\n"
        "neoantigen\tpredictors\tMHCflurry\tMHC-II predictions not run in P0\n"
        "neoantigen\tcandidates_reported\t263\t\n"
        "neoantigen\tunique_candidates\t261\t\n"
        "neoantigen\tburden_strong_500nM\t231\t\n"
        "neoantigen\tburden_highaff_100nM\t139\t\n"
        "neoantigen\tframeshift_neoantigens\t0\tframeshift-derived among strong binders\n"
        "neoantigen\texpr_tpm_pass\t138\t\n")

    SUMMARY_B = (
        "class\tmetric\tvalue\tdetail\n"
        "neoantigen\thla_alleles_typed\t14\t\n"
        "neoantigen\thla_alleles_predictor_usable\t6\t\n"
        "neoantigen\tpredictors\tMHCflurry\tMHC-II predictions not run in P0\n"
        "neoantigen\tunique_candidates\t262\t\n"
        "neoantigen\tburden_strong_500nM\t232\t\n"
        "neoantigen\tburden_highaff_100nM\t161\t\n"
        "neoantigen\tframeshift_neoantigens\t0\t\n"
        "neoantigen\texpr_tpm_pass\t125\t\n")

    def test_wide_pivot_and_honest_missing(self):
        with tempfile.TemporaryDirectory() as td:
            a = os.path.join(td, "GC01.summary.tsv")
            b = os.path.join(td, "GC02.summary.tsv")
            with open(a, "w") as fh:
                fh.write(self.SUMMARY_A)
            with open(b, "w") as fh:
                fh.write(self.SUMMARY_B)
            out = os.path.join(td, "cohort.tsv")
            sys.argv = ["collect_neoantigen.py", "--inputs", a, b,
                        "--output", out]
            self.assertEqual(collect_neoantigen.main(), 0)
            with open(out) as fh:
                lines = fh.read().splitlines()
            header = lines[0].split("\t")
            self.assertEqual(header[0], "sample")
            self.assertIn("predictors_note", header)
            by_a = dict(zip(header, lines[1].split("\t")))
            by_b = dict(zip(header, lines[2].split("\t")))
            self.assertEqual(by_a["sample"], "GC01")
            self.assertEqual(by_a["burden_strong_500nM"], "231")
            self.assertEqual(by_a["expr_tpm_pass"], "138")
            self.assertEqual(by_a["predictors"], "MHCflurry")
            self.assertIn("MHC-II", by_a["predictors_note"])
            # honest missing: candidates_reported absent in B -> empty cell
            self.assertEqual(by_b["candidates_reported"], "")
            self.assertEqual(by_b["burden_highaff_100nM"], "161")
            # deterministic sample order (sorted by input path)
            self.assertEqual([by_a["sample"], by_b["sample"]], ["GC01", "GC02"])


class TcrReportTest(unittest.TestCase):
    """venus#31 P1 item 2: TRUST4 report.tsv -> TCRbeta clonotype TSV.

    Fixture mirrors trust-simplerep.pl output (packaged trust4 1.1.10):
    header `#count frequency CDR3nt CDR3aa V D J C cid cid_full_length`,
    counts as integers, frequency in %e scientific notation, absent genes
    rendered as ".".
    """

    REPORT = (
        "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\t"
        "cid_full_length\n"
        # productive TRB clonotype, full-length assembly (24 nt = 8 aa)
        "42\t1.000000e+00\tTGTGCCAGCAGTCTGGGAGAATGG\tCASSLGEW\t"
        "TRBV5-1\t.\tTRBJ2-1\tTRBC1\tasm1\t1\n"
        # non-productive: out_of_frame rendering
        "12\t2.857143e-01\tTGTGCCAGCAGTCTG\tout_of_frame\tTRBV12-3\t.\t"
        "TRBJ1-1\tTRBC2\tasm2\t1\n"
        # TRA clonotype (alpha, not beta) -> filtered out
        "100\t6.800000e-01\tGCCTGAGGAC\tCAASSTGNQFYF\tTRAJ45\t.\t"
        "TRAC\t.\tasm3\t1\n"
        # stop codon in aa
        "7\t1.666667e-01\tTGTGGCCCTTTG\tC*G_L\tTRBV20-1\t.\tTRBJ2-7\t"
        "TRBC2\tasm4\t0\n"
    )

    def run_cli(self, tmp, report_text):
        rep = os.path.join(tmp, "trust_report.tsv")
        write_tsv(rep, report_text)
        out = os.path.join(tmp, "clonotypes.tsv")
        proc = subprocess.run(
            [sys.executable, os.path.join(REPO, "scripts", "tcr_report.py"),
             "--report", rep, "--sample", "S1", "--out", out],
            capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        with open(out) as fh:
            return list(csv.DictReader(fh, delimiter="\t")), proc.stderr

    def test_trb_filter_ranking_and_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, err = self.run_cli(tmp, self.REPORT)
            self.assertEqual(len(rows), 3)
            # ranked by descending count
            self.assertEqual([r["clonotype_rank"] for r in rows],
                             ["1", "2", "3"])
            self.assertEqual([r["count"] for r in rows], ["42", "12", "7"])
            top = rows[0]
            self.assertEqual(top["sample"], "S1")
            self.assertEqual(top["v_gene"], "TRBV5-1")
            self.assertEqual(top["c_gene"], "TRBC1")
            self.assertEqual(top["productive"], "yes")
            self.assertEqual(top["full_length"], "yes")
            of = rows[1]
            self.assertEqual(of["productive"], "no")
            self.assertEqual(of["cdr3_aa"], "out_of_frame")
            stop = rows[2]
            self.assertEqual(stop["productive"], "no")
            # honest stderr note about filtering is allowed but not the
            # no-clonotype warning
            self.assertNotIn("WARNING no TCRbeta", err)
            # no TRA gene leaked in
            self.assertFalse(any("TRA" in r["v_gene"] for r in rows))

    def test_empty_trb_writes_header_only(self):
        bcell = (
            "#count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\tC\tcid\t"
            "cid_full_length\n"
            "99\t1.000000e+00\tATGGAGCTG\tMIYSGI\tIGHV1-2\t.\tIGHJ4\t"
            "IGHM\tasm9\t1\n")
        with tempfile.TemporaryDirectory() as tmp:
            rows, err = self.run_cli(tmp, bcell)
            self.assertEqual(rows, [])
            self.assertIn("WARNING no TCRbeta clonotypes", err)
            out = os.path.join(tmp, "clonotypes.tsv")
            with open(out) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(lines[0].split("\t"), tcr_report.HEADER)
            self.assertEqual(len(lines), 1)

    def test_bad_header_fails_loud(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(AssertionError):
                self.run_cli(tmp, "count\tfreq\n1\t1e0\n")


class MantaToArribaTest(unittest.TestCase):
    """venus#31 P1 item 1: Manta BND -> Arriba-format TSV for pVACfuse.

    Tiny synthetic genome + GTF + VCF; fetch is injected so the conversion
    core stays stdlib-only. Donor: GENEA (+ strand, chr1, CDS 51-150, ATG at
    51..53, AAA codons, TAA stop at 141..143); acceptor: GENEB (- strand,
    chr2, CDS 51-150, ATG at 51..53, CCC codons, TTT stop at 141..143).
    """

    GTF = (
        "chr1\tg\tCDS\t51\t150\t.\t+\t.\t"
        "gene_id \"T1\"; transcript_id \"T1\"; gene_name \"GENEA\";\n"
        "chr2\tg\tCDS\t51\t150\t.\t-\t.\t"
        "gene_id \"T2\"; transcript_id \"T2\"; gene_name \"GENEB\";\n")

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        tmp = self._tmp.name
        write_tsv(os.path.join(tmp, "g.gtf"), self.GTF)
        self.gtf = os.path.join(tmp, "g.gtf")

    def tearDown(self):
        self._tmp.cleanup()

    @staticmethod
    def fetch(chrom, s, e):
        genome = {
            "chr1": ("N" * 50 + "ATG" + "AAA" * 29 + "TAA" + "GGG" * 2),
            "chr2": ("N" * 50 + "ATG" + "CCC" * 29 + "TTT" + "GGG" * 2),
        }
        genome = {c: s2.ljust(200, "N") for c, s2 in genome.items()}
        return genome[chrom][s:e]

    def make_vcf(self, alt_a, alt_b, samples=("TUMOR",), pr="1,12",
                 sr="0,7"):
        lines = [
            "##fileformat=VCFv4.2",
            "##INFO=<ID=MATEID,Number=.,Type=String>",
            "##INFO=<ID=SVTYPE,Number=1,Type=String>",
            "##FORMAT=<ID=PR,Number=1,Type=Integer>",
            "##FORMAT=<ID=SR,Number=1,Type=Integer>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(samples),
            "chr1\t100\tBND_A\tN\t%s\t.\tPASS\tSVTYPE=BND;MATEID=BND_B\t"
            "PR:SR\t%s" % (alt_a, ":".join([pr, sr] * len(samples))),
            "chr2\t120\tBND_B\tN\t%s\t.\tPASS\tSVTYPE=BND;MATEID=BND_A\t"
            "PR:SR\t%s" % (alt_b, ":".join([pr, sr] * len(samples))),
        ]
        path = os.path.join(self._tmp.name, "bnd.vcf")
        write_tsv(path, "\n".join(lines) + "\n")
        return path

    def run_convert(self, alt_a, alt_b, samples=("TUMOR",)):
        rows, skipped = mta.convert(self.make_vcf(alt_a, alt_b, samples),
                                    self.gtf, self.fetch)
        return rows, skipped

    def test_valid_pair_yields_row_with_junction_and_support(self):
        # BND_A keeps LEFT of chr1:100 (five' for the + gene), BND_B keeps
        # RIGHT of chr2:120 (three' for the - gene) -> one 5'/3' fusion.
        rows, skipped = self.run_convert("N]chr2:120]", "]chr1:100]N")
        self.assertEqual(skipped, {})
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["#gene1"], "GENEA")
        self.assertEqual(r["gene2"], "GENEB")
        self.assertEqual(r["breakpoint1"], "chr1:100")
        self.assertEqual(r["breakpoint2"], "chr2:120")
        self.assertIn("|", r["peptide_sequence"])
        # donor side = MKKKKK...; acceptor side = revcomp of chr2 piece
        self.assertTrue(r["peptide_sequence"].startswith("MK"))
        self.assertEqual(r["reading_frame"], "out-of-frame")
        self.assertEqual(r["split_reads1"], 7)
        self.assertEqual(r["split_reads2"], 7)
        self.assertEqual(r["discordant_mates"], 12)
        self.assertEqual(r["FUSION_ID"], "BND_A|BND_B")

    def test_head_to_head_has_no_53_model_and_is_skipped(self):
        # Both sides keep a transcript-UPSTREAM piece -> no 5'/3' model.
        rows, skipped = self.run_convert("N[chr2:120[", "]chr1:100]N")
        self.assertEqual(rows, [])
        self.assertIn("no 5'/3' transcript model (head/tail-to-tail or "
                      "strand clash)", skipped)

    def test_two_sample_vcf_picks_highest_support_sample(self):
        # NORMAL carries more PR+SR alt support than TUMOR -> picked, with a
        # stderr warning (auto-detection is honest, never silent).
        rows, _ = self.run_convert("N]chr2:120]", "]chr1:100]N",
                                   samples=("NORMAL", "TUMOR"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["split_reads1"], 7)

    def test_zero_fusions_writes_header_only_tsv(self):
        out = os.path.join(self._tmp.name, "out.tsv")
        vcf = self.make_vcf("N[chr2:120[", "]chr1:100]N")
        rc = mta.MantaToArriba(self.gtf, self.fetch).run(vcf, out)
        self.assertEqual(rc, 0)
        with open(out) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[0].split("\t"), mta.ARRIBA_COLUMNS)
        self.assertEqual(len(lines), 1)  # header only


class HlaConsistencyTest(unittest.TestCase):
    """SpecHLA (DNA) vs arcasHLA (RNA) per-gene QC (venus#31 P1 item 3)."""

    DNA_TSV = (
        "sample\tgene\tallele_1\tallele_2\tallele_1_full\tallele_2_full\t"
        "allele_1_g_group\tallele_2_g_group\ttyped\tspec_version\tmode\t"
        "caveat\n"
        "Pt01-TD\tA\tHLA-A*02:01\tHLA-A*29:02\t\t\t\t\tTRUE\t1.0.12\t"
        "tumor_only\tcaveat\n"
        "Pt01-TD\tB\tHLA-B*07:02\tHLA-B*44:03\t\t\t\t\tTRUE\t1.0.12\t"
        "tumor_only\tcaveat\n"
        "Pt01-TD\tC\tHLA-C*07:01\tHLA-C*01:02\t\t\t\t\tTRUE\t1.0.12\t"
        "tumor_only\tcaveat\n")

    def run_cli(self, rna_json_text, dna_text=None, dna_sample="Pt01-TD",
                rna_sample="Pt01-TR", patient="Pt01"):
        with tempfile.TemporaryDirectory() as tmp:
            dna = os.path.join(tmp, "dna.tsv")
            rna = os.path.join(tmp, "rna.json")
            out = os.path.join(tmp, "out.tsv")
            write_tsv(dna, dna_text or self.DNA_TSV)
            write_tsv(rna, rna_json_text)
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts",
                                              "hla_consistency.py"),
                 "--dna-tsv", dna, "--dna-sample", dna_sample,
                 "--rna-json", rna, "--rna-sample", rna_sample,
                 "--patient", patient, "--out", out],
                capture_output=True, text=True)
            if proc.returncode != 0:
                return proc.returncode, [], proc.stderr
            with open(out) as fh:
                return 0, list(csv.DictReader(fh, delimiter="\t")), proc.stderr

    def test_match_mismatch_absent_not_run(self):
        # A: match (order-insensitive, 3-field RNA normalizes to 2-field);
        # B: arcasHLA omitted the locus (below min_count, absent from the
        # JSON) -> absent_rna;
        # C: arcasHLA called only one allele (half-call) -> mismatch;
        # D-class genes: not typed anywhere -> not_run.
        rna = ('{"A": ["HLA-A*29:02:01", "HLA-A*02:01:01"], '
               '"C": ["HLA-C*07:01"]}')
        rc, rows, err = self.run_cli(rna)
        self.assertEqual(rc, 0)
        self.assertEqual(len(rows), 8)
        by_gene = {r["gene"]: r for r in rows}
        self.assertEqual(by_gene["A"]["concordance"], "match")
        self.assertEqual(by_gene["A"]["rna_alleles"],
                         "HLA-A*02:01;HLA-A*29:02")
        self.assertEqual(by_gene["B"]["concordance"], "absent_rna")
        self.assertIn("min_count", by_gene["B"]["detail"])
        self.assertEqual(by_gene["C"]["concordance"], "mismatch")
        self.assertEqual(by_gene["DRB1"]["concordance"], "not_run")
        self.assertIn("1/2", err)

    def test_absent_dna_when_spechla_skips_locus(self):
        # SpecHLA leaves a locus untyped (tumor DNA gap / LOH) while arcas
        # typed it -> absent_dna, with the RNA call carried in the row.
        dna = (
            "sample\tgene\tallele_1\tallele_2\tallele_1_full\tallele_2_full\t"
            "allele_1_g_group\tallele_2_g_group\ttyped\tspec_version\tmode\t"
            "caveat\n"
            "Pt01-TD\tA\tHLA-A*02:01\tHLA-A*29:02\t\t\t\t\tTRUE\t1.0.12\t"
            "tumor_only\tcaveat\n")
        rna = ('{"A": ["HLA-A*02:01:01", "HLA-A*29:02:01"], '
               '"B": ["HLA-B*07:02:01", "HLA-B*44:03:01"]}')
        rc, rows, _ = self.run_cli(rna, dna_text=dna)
        self.assertEqual(rc, 0)
        by_gene = {r["gene"]: r for r in rows}
        self.assertEqual(by_gene["A"]["concordance"], "match")
        self.assertEqual(by_gene["B"]["concordance"], "absent_dna")
        self.assertEqual(by_gene["B"]["rna_alleles"],
                         "HLA-B*07:02;HLA-B*44:03")
        self.assertEqual(by_gene["C"]["concordance"], "not_run")

    def test_heterozygote_orientation_does_not_fake_mismatch(self):
        # allele1/allele2 swap between layers must not be a mismatch.
        rna = ('{"B": ["HLA-B*44:03:01:02", "HLA-B*07:02:01:01"], '
               '"C": ["HLA-C*01:02:01", "HLA-C*07:01:01:01"]}')
        rc, rows, _ = self.run_cli(rna)
        self.assertEqual(rc, 0)
        by_gene = {r["gene"]: r for r in rows}
        self.assertEqual(by_gene["B"]["concordance"], "match")
        self.assertEqual(by_gene["C"]["concordance"], "match")

    def test_unreadable_rna_json_fails_loud(self):
        rc, rows, err = self.run_cli('["not", "an", "object"]')
        self.assertNotEqual(rc, 0)
        self.assertIn("unreadable", err)


class HlaLohTest(unittest.TestCase):
    """spechla-loh driver (venus#31 P1 item 4) — fake-spechla-loh harness."""

    PURITY_TSV = (
        "pair\tsample\tpurity\tploidy\tpsi\taberrant_cell_fraction\t"
        "goodness_of_fit\tnon_aberrant\tfailed\tn_segments\n"
        "Pt01\tPt01-TD\t0.72\t2.9\t2\t0.71\t0.91\tFALSE\tFALSE\t41\n")

    # merge.hla.copy.txt as cal.hla.copy.pl writes it (verified header;
    # the tool's tag==0 homogeneous row puts a real allele in KeptHLA and
    # the literal "homogeneous" in LossHLA)
    MERGE = (
        "Sample\tHLA\tAllele1\tAllele2\tcopyratio\tKeptHLA\tLossHLA\tFreq1\t"
        "Freq2\tPurity\tHet_num\tLOH\n"
        "Pt01-TD\tA\tHLA-A*02:01\tHLA-A*29:02\t2:1\tHLA-A*29:02\t"
        "HLA-A*02:01\t0.62\t0.38\t0.72\t17\tY\n"
        "Pt01-TD\tB\tHLA-B*07:02\tHLA-B*44:03\t1:1\tHLA-B*44:03\t"
        "HLA-B*07:02\t0.51\t0.49\t0.72\t23\tN\n"
        "Pt01-TD\tC\tHLA-C*07:01\tHLA-C*01:02\t3:3\tHLA-C*07:01\t"
        "homogeneous\t1.0\t0.0\t0.72\t9\tN\n")

    FREQ_SAMPLE = (
        "# Haplotype\tFrequency\n"
        "hla.allele.1.HLA_A.fasta\t0.62\n"
        "hla.allele.2.HLA_A.fasta\t0.38\n"
        "# The number of heterozygous variant is 17\n")

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        tmp = self._tmp.name
        # fake spechla-loh: echoes args to a log, writes merge.hla.copy.txt
        # into the -O dir, exits 0 — mimics cal.hla.copy.pl's contract.
        self.bin_dir = os.path.join(tmp, "bin")
        os.mkdir(self.bin_dir)
        fake = os.path.join(self.bin_dir, "spechla-loh")
        write_tsv(fake, "#!/bin/bash\n")
        with open(fake, "a") as fh:
            fh.write(
                'args=("$@"); outdir=""; while [ $# -gt 0 ]; do case "$1" in '
                '-O) outdir="$2";; esac; shift; done\n'
                'printf \'%s\\n\' "${args[@]}" > "$outdir/spechla-loh.args"\n'
                'cat > "$outdir/merge.hla.copy.txt" <<\'MERGE_EOF\'\n'
                + self.MERGE + "MERGE_EOF\n")
        os.chmod(fake, 0o755)
        # inputs
        self.purity = os.path.join(tmp, "purity.tsv")
        write_tsv(self.purity, self.PURITY_TSV)
        self.result = os.path.join(tmp, "hla.result.txt")
        write_tsv(self.result, "SpecHLA v1.0.12\n")
        self.freqs = []
        for gene in ("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"):
            f = os.path.join(tmp, f"HLA_{gene}_freq.txt")
            write_tsv(f, self.FREQ_SAMPLE)
            self.freqs.append(f)

    def tearDown(self):
        self._tmp.cleanup()

    def run_driver(self, outdir, out, extra_env=None):
        env = dict(os.environ)
        env["PATH"] = self.bin_dir + os.pathsep + env["PATH"]
        if extra_env:
            env.update(extra_env)
        return subprocess.run(
            [sys.executable, os.path.join(REPO, "scripts", "hla_loh.py"),
             "--pair", "Pt01", "--sample", "Pt01-TD",
             "--result", self.result,
             "--freq", *self.freqs,
             "--purity-tsv", self.purity,
             "--het-cutoff", "5",
             "--outdir", outdir, "--out", out],
            capture_output=True, text=True, env=env)

    def test_loh_table_and_homogeneous_note(self):
        outdir = os.path.join(self._tmp.name, "loh_work")
        out = os.path.join(self._tmp.name, "loh.tsv")
        proc = self.run_driver(outdir, out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["pair"], "Pt01")
        self.assertEqual(rows[0]["LOH"], "Y")
        self.assertEqual(rows[1]["LOH"], "N")
        # tool's own homogeneous tag -> note, never re-interpreted as LOH
        self.assertEqual(rows[2]["LOH"], "N")
        self.assertIn("homogeneous", rows[2]["note"])
        self.assertIn("1/3", proc.stderr)
        # purity/ploidy actually forwarded to spechla-loh
        with open(os.path.join(outdir, "spechla-loh.args")) as fh:
            args = fh.read().splitlines()
        self.assertIn("0.720000", args)
        self.assertIn("2.900000", args)

    def test_missing_ascat_fit_fails_loud(self):
        empty = os.path.join(self._tmp.name, "empty.tsv")
        write_tsv(empty, "pair\tsample\tpurity\tploidy\n")
        out = os.path.join(self._tmp.name, "loh.tsv")
        proc = subprocess.run(
            [sys.executable, os.path.join(REPO, "scripts", "hla_loh.py"),
             "--pair", "Pt01", "--sample", "Pt01-TD",
             "--result", self.result, "--freq", *self.freqs,
             "--purity-tsv", empty, "--outdir", os.path.join(
                 self._tmp.name, "w2"), "--out", out],
            capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("REQUIRES purity and ploidy", proc.stderr)

    def test_missing_freq_file_fails_loud(self):
        os.remove(self.freqs[3])  # DPA1
        out = os.path.join(self._tmp.name, "loh.tsv")
        proc = self.run_driver(os.path.join(self._tmp.name, "w3"), out)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("missing freq file", proc.stderr)


class ScrnaCorroborationTest(unittest.TestCase):
    """venus#31 P1 item 5: scRNA %positive corroboration per candidate gene.

    h5ad fixtures are built IN-TEST with h5py (importable where these tests
    run; envs/scrna is the only venus env carrying it) and mirror the anndata
    0.1.0 on-disk schema the h5py-only reader supports: var/obs _index
    (nullable-string Group or flat Dataset) and dense X / csr / csc_matrix
    Groups. Covers: dense vs CSR vs CSC equivalence, explicit-zero handling,
    gene-absent note (never a silent zero), the '-1' barcode-suffix fallback,
    and the corroboration merge in neoantigen_report (--scrna-corr optional
    input; base candidate schema stays stable without it).
    """

    GENES = ["EGFR", "KRAS", "TP53", "ABCF1"]
    # 6 cells x 4 genes; EGFR positive in cells 0,2,4 -> 3/6, mean (5+3+1)/3
    M = [[5.0, 0.0, 1.0, 0.0],
         [0.0, 2.0, 0.0, 0.0],
         [3.0, 0.0, 0.0, 4.0],
         [0.0, 0.0, 2.0, 0.0],
         [1.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0]]
    BARCODES = ["AAAC-1", "AAAD-1", "AAAE-1", "AAAF-1", "AAAG-1", "AAAH-1"]
    CLUSTER_CSV = ("barcode,cluster\nAAAC-1,c0\nAAAD-1,c1\nAAAE-1,c0\n"
                   "AAAF-1,c1\nAAAG-1,c0\nAAAH-1,c0\n")
    PVAC = ("Gene Name\tVariant Type\tHLA Allele\tBest MT IC50 Score\n"
            "EGFR\tsnv\tHLA-A*02:01\t12.5\n"
            "KRAS\tsnv\tHLA-A*02:01\t88.0\n"
            "NOTINMATRIX\tframeshift\tHLA-B*07:02\t50.0\n"
            "TP53\tsnv\tHLA-C*01:02\t6.2\n")

    def setUp(self):
        try:
            import h5py  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("h5py unavailable on this runner")
        self._tmp = tempfile.TemporaryDirectory()
        tmp = self._tmp.name
        self.pvac = os.path.join(tmp, "filtered.tsv")
        write_tsv(self.pvac, self.PVAC)
        self.cluster_csv = os.path.join(tmp, "cluster.csv")
        write_tsv(self.cluster_csv, self.CLUSTER_CSV)

    def tearDown(self):
        self._tmp.cleanup()

    def _write_h5ad(self, kind, barcodes=None, name=None):
        """Build a minimal anndata-schema h5ad with X in the given encoding."""
        import h5py
        import numpy as np
        barcodes = barcodes or self.BARCODES
        hp = os.path.join(self._tmp.name,
                          name or f"{kind}_{len(barcodes)}.h5ad")
        m = np.array(self.M, dtype=np.float64)
        with h5py.File(hp, "w") as f:
            # anndata 0.1.0 nullable-string encoding: <frame>/_index/Group{values}
            for frame, names in (("var", self.GENES), ("obs", barcodes)):
                fr = f.create_group(frame)
                ix = fr.create_group("_index")
                ix.create_dataset("values",
                                  data=[n.encode() for n in names])
            r, c = np.nonzero(m)
            if kind == "dense":
                x = f.create_dataset("X", data=m)
                x.attrs["encoding-type"] = "array"
            else:
                g = f.create_group("X")
                g.attrs["encoding-type"] = f"{kind}_matrix"
                g.attrs["shape"] = m.shape
                if kind == "csr":
                    order = np.argsort(r, kind="stable")
                    vals, idx = m[r[order], c[order]], c[order]
                    counts = np.bincount(r, minlength=m.shape[0])
                else:
                    order = np.argsort(c, kind="stable")
                    vals, idx = m[r[order], c[order]], r[order]
                    counts = np.bincount(c, minlength=m.shape[1])
                g.create_dataset("data", data=vals)
                g.create_dataset("indices", data=idx)
                g.create_dataset("indptr", data=np.concatenate(
                    ([0], np.cumsum(counts))).astype(np.int64))
        return hp

    def run_corr(self, h5ad_path, with_cluster=False, out=None):
        out = out or os.path.join(self._tmp.name, "corr.tsv")
        cmd = [sys.executable,
               os.path.join(REPO, "scripts", "scrna_corroboration.py"),
               "--h5ad", h5ad_path, "--pvac-tsv", self.pvac, "--out", out]
        if with_cluster:
            cmd += ["--cluster", self.cluster_csv]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        with open(out) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows, proc.stderr

    def test_dense_csr_csc_agree_and_missing_gene_note(self):
        for kind in ("dense", "csr", "csc"):
            rows, err = self.run_corr(self._write_h5ad(kind))
            by_gene = {r["gene"]: r for r in rows}
            self.assertEqual(by_gene["EGFR"]["n_cells_total"], "6")
            self.assertEqual(by_gene["EGFR"]["n_expressing"], "3")
            self.assertEqual(by_gene["EGFR"]["pct_expressing"], "50.0")
            self.assertEqual(by_gene["EGFR"]["mean_expr_positive"], "3.00")
            self.assertEqual(by_gene["KRAS"]["n_expressing"], "1")
            self.assertEqual(by_gene["KRAS"]["pct_expressing"], "16.7")
            # TP53 column [1,0,0,2,0,0]: 2/6 positive, mean (1+2)/2
            self.assertEqual(by_gene["TP53"]["n_expressing"], "2")
            self.assertEqual(by_gene["TP53"]["mean_expr_positive"], "1.50")
            # absent gene: explicit note, never a silent zero row
            self.assertEqual(by_gene["NOTINMATRIX"]["note"],
                             "gene_not_in_scRNA_matrix")
            self.assertEqual(by_gene["NOTINMATRIX"]["pct_expressing"], "")
            self.assertIn("no malignant-cell annotation", err)

    def test_per_cluster_breakdown(self):
        rows, _ = self.run_corr(self._write_h5ad("dense"), with_cluster=True)
        by_gene = {r["gene"]: r for r in rows}
        # c0 = cells 0,2,4,5; c1 = cells 1,3
        self.assertEqual(by_gene["EGFR"]["per_cluster_pct"],
                         "c0:3/4; c1:0/2")
        self.assertEqual(by_gene["KRAS"]["per_cluster_pct"], "c0:0/4; c1:1/2")
        self.assertEqual(by_gene["TP53"]["per_cluster_pct"],
                         "c0:1/4; c1:1/2")

    def test_barcode_suffix_fallback(self):
        # h5ad barcodes bare (dnbc4tools style), cluster.csv carries '-1'
        # (10x style): exact match matches nothing -> suffix-strip retry.
        rows, err = self.run_corr(
            self._write_h5ad("dense", barcodes=[b.split("-")[0]
                                                for b in self.BARCODES]),
            with_cluster=True)
        by_gene = {r["gene"]: r for r in rows}
        self.assertEqual(by_gene["EGFR"]["per_cluster_pct"], "c0:3/4; c1:0/2")

    def test_bad_header_fails_loud(self):
        bad = os.path.join(self._tmp.name, "bad.tsv")
        write_tsv(bad, "gene\tpct\nEGFR\t50\n")
        with self.assertRaises(SystemExit):
            nr.read_scrna_corr(bad)

    def test_report_merges_corroboration_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            corr = os.path.join(tmp, "corr.tsv")
            self.run_corr(self._write_h5ad("dense"), with_cluster=True,
                          out=corr)
            rows_text = ("EGFR\tmissense_variant\tHLA-A*02:01\tGVLVGVAPE\t"
                         "2.5\t0.1\t0.4\tp.E12D\t12\tENST1\n"
                         "ZZZ9\tframeshift_variant\tHLA-B*07:02\tYYYYYYYYL\t"
                         "100\t2\t0.2\tp.R1fs\t1\tENST3\n")
            hla = "T1\tA\tHLA-A*02:01\tHLA-A*24:02\n"
            summary, cands, md = run_report(
                tmp, rows_text, hla, extra=["--scrna-corr", corr])
            # one candidate gene (EGFR) carries a %positive value
            self.assertIn(("scrna_pct_expressing", "1"), summary)
            egfr = next(c for c in cands if c["gene"] == "EGFR")
            self.assertEqual(egfr["scrna_pct_expressing"], "50.0")
            self.assertEqual(egfr["scrna_per_cluster_pct"], "c0:3/4; c1:0/2")
            # gene in pVAC output but not in the corroboration TSV at all
            zz = next(c for c in cands if c["gene"] == "ZZZ9")
            self.assertEqual(zz["scrna_note"], "no_scrna_corroboration_run")
            # absent-from-matrix gene keeps the explicit note through merge
            # (NOTINMATRIX is not among pVAC rows here; use a corr row check)
            # schema stability: a run WITHOUT --scrna-corr keeps base columns
            _, cands2, _ = run_report(tmp, rows_text, hla)
            self.assertNotIn("scrna_pct_expressing", cands2[0])
            self.assertNotIn("scrna_note", cands2[0])
            self.assertEqual(sorted(cands2[0].keys()),
                             sorted(nr.CAND_COLS))
            # md: caveat bullet + extended top-N table
            self.assertIn("over ALL filtered cells", md)
            self.assertIn("no malignant-cell annotation", md)


class MsiFrameshiftJoinTest(unittest.TestCase):
    """venus#31 P1 item 7: frameshift_neoantigens x MSI status cohort page.

    Covers the three honesty-critical shapes of the join:
      both sides present (MSI-H first, frameshift share computed),
      tumor-only MSI (header-only msi_cohort -> every MSI cell
      "not_available", never a silent blank), and paired-only neoantigen
      (header-only neoantigen_cohort -> frameshift cells not_available).
    """

    NEO = ("sample\thla_alleles_typed\thla_alleles_predictor_usable\t"
           "candidates_reported\tunique_candidates\tburden_strong_500nM\t"
           "burden_highaff_100nM\texpr_tpm_pass\tframeshift_neoantigens\t"
           "predictors\tpredictors_note\n"
           "GC01\t13\t6\t263\t261\t20\t139\t138\t2\tMHCflurry\tnote\n"
           "GC02\t14\t6\t200\t200\t10\t100\t125\t0\tMHCflurry\tnote\n")
    MSI = ("sample\traw_total_sites\tsummary_total_sites\treported_unstable\t"
           "reported_pct\tbh_unstable\tbh_pct\tunstable_delta\t"
           "uint16_saturated\tmsi_status_wgs\tmsi_status_wes\n"
           "GC01\t5000\t5000\t3\t0.06\t4\t0.08\t1\tfalse\tMSI-H\tMSS\n")
    MSI_HEADER_ONLY = (
        "sample\traw_total_sites\tsummary_total_sites\treported_unstable\t"
        "reported_pct\tbh_unstable\tbh_pct\tunstable_delta\t"
        "uint16_saturated\tmsi_status_wgs\tmsi_status_wes\n")

    def run_join(self, tmp, neo_text, msi_text, extra=()):
        neo = os.path.join(tmp, "neo.tsv")
        msi = os.path.join(tmp, "msi.tsv")
        tsv = os.path.join(tmp, "joined.tsv")
        md = os.path.join(tmp, "joined.md")
        write_tsv(neo, neo_text)
        write_tsv(msi, msi_text)
        rc = mfj.main_with_args(
            ["--neoantigen-cohort", neo, "--msi-cohort", msi,
             "--strong-nm", "500", "--out-tsv", tsv, "--out-md", md,
             *extra])
        assert rc == 0
        with open(tsv) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        with open(md) as fh:
            md_text = fh.read()
        return rows, md_text

    def test_join_msih_first_and_frameshift_share(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(tmp, self.NEO, self.MSI)
            self.assertEqual([r["sample"] for r in rows], ["GC01", "GC02"])
            gc01 = rows[0]
            self.assertEqual(gc01["frameshift_neoantigens"], "2")
            self.assertEqual(gc01["msi_status_wgs"], "MSI-H")
            self.assertEqual(gc01["msi_status_wes"], "MSS")
            self.assertEqual(gc01["msi_bh_pct"], "0.08")
            # frameshift share of strong binders: 2/20 = 10.0%
            self.assertEqual(gc01["frameshift_pct_of_strong"], "10.0")
            gc02 = rows[1]
            self.assertEqual(gc02["frameshift_neoantigens"], "0")
            self.assertEqual(gc02["msi_status_wgs"], "not_available")
            self.assertEqual(gc02["frameshift_pct_of_strong"], "0.0")
            self.assertIn("msi_not_available", gc02["note"])
            # md states cutoffs and the sentinel
            self.assertIn("3.5%", md)
            self.assertIn("15%", md)
            self.assertIn("not_available", md)
            self.assertIn("1 sample(s) have both sides joined", md)

    def test_tumor_only_msi_header_only_honest_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(tmp, self.NEO, self.MSI_HEADER_ONLY)
            self.assertEqual(len(rows), 2)
            for r in rows:
                self.assertEqual(r["msi_status_wgs"], "not_available")
                self.assertEqual(r["msi_status_wes"], "not_available")
                self.assertEqual(r["msi_bh_pct"], "not_available")
                self.assertIn("msi_not_available", r["note"])
            self.assertIn("header-only", md)
            self.assertIn("matched normal", md)
            self.assertIn("frameshift_neoantigens", md)

    def test_paired_neo_header_only_honest_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp,
                "sample\thla_alleles_typed\tframeshift_neoantigens\n",
                self.MSI)
            self.assertEqual(rows[0]["sample"], "GC01")
            self.assertEqual(rows[0]["frameshift_neoantigens"],
                             "not_available")
            self.assertEqual(rows[0]["msi_status_wgs"], "MSI-H")
            self.assertEqual(rows[0]["frameshift_pct_of_strong"],
                             "not_available")
            self.assertIn("frameshift_not_available", rows[0]["note"])

    def test_sample_only_on_one_side(self):
        msi = self.MSI.replace("GC01", "GC09")
        with tempfile.TemporaryDirectory() as tmp:
            rows, _ = self.run_join(tmp, self.NEO, msi)
            # the unmatched MSI sample keeps its row (honest presence);
            # sort order: MSI-H rows first, then by sample id
            self.assertEqual([r["sample"] for r in rows],
                             ["GC09", "GC01", "GC02"])
            gc09 = rows[0]
            self.assertEqual(gc09["frameshift_neoantigens"],
                             "not_available")
            self.assertEqual(gc09["frameshift_pct_of_strong"],
                             "not_available")
            self.assertIn("frameshift_not_available", gc09["note"])
            gc01 = rows[1]
            self.assertEqual(gc01["msi_status_wgs"], "not_available")
            self.assertIn("msi_not_available", gc01["note"])

    def test_empty_inputs_render_honest_page(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp,
                "sample\tframeshift_neoantigens\n",
                "sample\tmsi_status_wgs\tmsi_status_wes\tbh_pct\n")
            self.assertEqual(rows, [])
            self.assertIn("no rows", md)


class IcbLinkageTest(unittest.TestCase):
    """venus#31: TMB x MSI x neoantigen burden cross-evidence page.

    Covers the honesty-critical shapes: full join with both ICB flag
    conventions (TMB-H FDA >=13 kept separate from 10-<13 legacy; MSI
    wgs/wes side by side), patient-key normalization (-TD stripped),
    header-only MSI (tumor-only) and header-only neoantigen (paired),
    a cohort summary without the TMB column, and fully empty inputs.
    """

    COHORT = (
        "sample\ttmb_sample_id\ttmb_tmb_mut_per_mb\ttmb_all_somatic_mut_per_mb\n"
        "Pt01\tPt01-TD\t13.5\t20.1\n"
        "Pt02\tPt02-TD\t11.2\t15.0\n"
        "Pt03\tPt03-TD\t3.0\t4.5\n")
    MSI = (
        "sample\traw_total_sites\tsummary_total_sites\treported_unstable\t"
        "reported_pct\tbh_unstable\tbh_pct\tunstable_delta\t"
        "uint16_saturated\tmsi_status_wgs\tmsi_status_wes\n"
        "Pt01-TD\t5000\t5000\t3\t0.06\t4\t0.08\t1\tfalse\tMSI-H\tMSS\n"
        "Pt02-TD\t5000\t5000\t1\t0.02\t2\t0.04\t1\tfalse\tMSS\tMSS\n")
    NEO = (
        "sample\thla_alleles_typed\thla_alleles_predictor_usable\t"
        "candidates_reported\tunique_candidates\tburden_strong_500nM\t"
        "burden_highaff_100nM\texpr_tpm_pass\tframeshift_neoantigens\t"
        "predictors\tpredictors_note\n"
        "Pt01\t13\t6\t263\t261\t20\t139\t138\t2\tMHCflurry\tnote\n"
        "Pt02\t14\t6\t200\t200\t10\t100\t125\t0\tMHCflurry\tnote\n")
    MSI_HEADER_ONLY = (
        "sample\traw_total_sites\tsummary_total_sites\treported_unstable\t"
        "reported_pct\tbh_unstable\tbh_pct\tunstable_delta\t"
        "uint16_saturated\tmsi_status_wgs\tmsi_status_wes\n")

    def run_join(self, tmp, cohort_text, msi_text, neo_text, extra=()):
        cohort = os.path.join(tmp, "cohort.tsv")
        msi = os.path.join(tmp, "msi.tsv")
        neo = os.path.join(tmp, "neo.tsv")
        tsv = os.path.join(tmp, "icb.tsv")
        md = os.path.join(tmp, "icb.md")
        write_tsv(cohort, cohort_text)
        write_tsv(msi, msi_text)
        write_tsv(neo, neo_text)
        rc = il.main_with_args(
            ["--cohort-summary", cohort, "--msi-cohort", msi,
             "--neoantigen-cohort", neo, "--out-tsv", tsv, "--out-md", md,
             *extra])
        assert rc == 0
        with open(tsv) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        with open(md) as fh:
            md_text = fh.read()
        return rows, md_text

    def test_full_join_flags_and_patient_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp, self.COHORT, self.MSI, self.NEO)
            self.assertEqual(
                [r["sample"] for r in rows], ["Pt01", "Pt02", "Pt03"])
            pt01 = rows[0]
            self.assertEqual(pt01["tmb_coding_mut_per_mb"], "13.5")
            self.assertEqual(pt01["tmb_status"],
                             "TMB-H (FDA >=13 mut/Mb coding)")
            self.assertEqual(pt01["msi_status_wgs"], "MSI-H")
            self.assertEqual(pt01["msi_status_wes"], "MSS")
            self.assertEqual(pt01["msi_bh_pct"], "0.08")
            self.assertEqual(pt01["neo_burden_strong_500nM"], "20")
            self.assertEqual(pt01["neo_frameshift_neoantigens"], "2")
            self.assertEqual(pt01["note"], "")
            pt02 = rows[1]
            # 11.2 mut/Mb coding: legacy context wording, never a label
            self.assertEqual(pt02["tmb_status"],
                             "TMB 10-<13 (legacy high-TMB context, not FDA label)")
            self.assertEqual(pt02["msi_status_wgs"], "MSS")
            pt03 = rows[2]
            self.assertEqual(pt03["tmb_status"], "not_high")
            for f in ("msi_status_wgs", "neo_burden_strong_500nM"):
                self.assertEqual(pt03[f], "not_available")
            self.assertIn("msi_not_available", pt03["note"])
            self.assertIn("neoantigen_not_available", pt03["note"])
            # md: conventions + counts
            self.assertIn(">= 13 mut/Mb coding", md)
            self.assertIn("3.5%", md)
            self.assertIn("15%", md)
            self.assertIn("never collapsed", md)
            self.assertIn("1 flagged TMB-H", md)
            self.assertIn("1 flagged MSI-H", md)

    def test_tumor_only_msi_header_only_honest_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp, self.COHORT, self.MSI_HEADER_ONLY, self.NEO)
            for r in rows:
                self.assertEqual(r["msi_status_wgs"], "not_available")
                self.assertEqual(r["msi_status_wes"], "not_available")
                self.assertIn("msi_not_available", r["note"])
            self.assertIn("header-only", md)
            self.assertIn("matched normal", md)

    def test_paired_neo_header_only_honest_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp, self.COHORT, self.MSI,
                "sample\thla_alleles_typed\tframeshift_neoantigens\n")
            for r in rows:
                self.assertEqual(r["neo_burden_strong_500nM"],
                                 "not_available")
                self.assertIn("neoantigen_not_available", r["note"])
            self.assertIn("header-only", md)

    def test_cohort_without_tmb_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp,
                "sample\trna_sample_id\trna_total_reads\n"
                "Pt01\tPt01-TD\t1000000\n",
                self.MSI, self.NEO)
            for r in rows:
                self.assertEqual(r["tmb_coding_mut_per_mb"],
                                 "not_available")
                self.assertEqual(r["tmb_status"], "not_available")
                self.assertIn("tmb_not_available", r["note"])
            self.assertIn("carries no", md)

    def test_empty_inputs_render_honest_page(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, md = self.run_join(
                tmp,
                "sample\ttmb_sample_id\ttmb_tmb_mut_per_mb"
                "\ttmb_all_somatic_mut_per_mb\n",
                "sample\tmsi_status_wgs\tmsi_status_wes\tbh_pct\n",
                "sample\tburden_strong_500nM\n")
            self.assertEqual(rows, [])
            self.assertIn("no cohort-summary row matched", md)
            self.assertIn("header-only", md)


if __name__ == "__main__":
    unittest.main(verbosity=2)
