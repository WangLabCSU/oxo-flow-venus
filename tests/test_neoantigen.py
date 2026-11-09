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
import manta_to_arriba as mta  # noqa: E402
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
