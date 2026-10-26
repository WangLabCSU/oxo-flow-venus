#!/usr/bin/env python3
"""Tests for clinical_interpret.py — audit M7 (therapy-matching layer).

Covers: NMPA reachability contract (every curated row must be reachable from
at least one code key-construction path), ERBB2 amp wiring (M2), NTRK folding
(M3), BRCA1/2 somatic-loss context row (M4), per-gene CNV tiers (M5), and
cohort-mode caveats (m1).

Run: python3 tests/test_clinical_interpret.py   (needs pandas)
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

import clinical_interpret as ci  # noqa: E402

NMPA = os.path.join(REPO, "resources", "clinical", "nmpa_oncology_approvals.tsv")
FUSIONS = os.path.join(REPO, "resources", "clinical", "actionable_fusion_pairs.tsv")
CNVGENES = os.path.join(REPO, "resources", "clinical", "actionable_cnv_genes.tsv")

SNV_COLS = ["chrom", "pos", "ref", "alt", "gene", "hgvsp", "consequence",
            "tumor_vaf", "tumor_dp", "gnomad_af", "clinvar_id", "clinvar_sig",
            "kb_tag", "civic_variant_id", "civic_variant",
            "civic_evidence_levels", "civic_evidence_types", "civic_therapies",
            "civic_diseases", "amp_tier", "amp_level", "evidence_source",
            "tier_basis"]


def write_tsv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def run_interpret(tmp, snv_rows, cnv_rows, sv_rows, tmb_rows=None,
                  cohort_mode="tumor_only"):
    d = tmp
    write_tsv(os.path.join(d, "snv.tsv"), SNV_COLS, snv_rows)
    write_tsv(os.path.join(d, "cnv.tsv"),
              ["gene", "call", "segment_log2", "segment_cn", "segment",
               "ecDNA_evidence", "disease_context", "therapy_context", "notes",
               "amp_tier"], cnv_rows)
    write_tsv(os.path.join(d, "sv.tsv"),
              ["fusion_candidate", "gene_5_call", "gene_3_call", "breakpoint1",
               "breakpoint2", "scoring"], sv_rows)
    if tmb_rows is None:
        tmb_path = ""
    else:
        write_tsv(os.path.join(d, "tmb.tsv"),
                  ["pair_id", "coding_mb", "tmb_mut_per_mb",
                   "all_somatic_mut_per_mb"], tmb_rows)
        tmb_path = os.path.join(d, "tmb.tsv")
    proc = subprocess.run(
        [sys.executable, os.path.join(REPO, "scripts", "clinical_interpret.py"),
         "--patient", "T1", "--snv-evidence", os.path.join(d, "snv.tsv"),
         "--cnv-evidence", os.path.join(d, "cnv.tsv"),
         "--sv-evidence", os.path.join(d, "sv.tsv"),
         "--tmb", tmb_path, "--nmpa", NMPA,
         "--cohort-mode", cohort_mode,
         "--out-biomarkers", os.path.join(d, "bm.tsv"),
         "--out-therapies", os.path.join(d, "th.tsv")],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    with open(os.path.join(d, "bm.tsv")) as fh:
        bm = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(d, "th.tsv")) as fh:
        th = list(csv.DictReader(fh, delimiter="\t"))
    return bm, th


def empty_snv():
    return ["7", "1", "A", "T", "TEST1", "A1T", "missense_variant",
            "0.5", "100", "0.0001", "", "", "", "", "", "", "", "", "", "", "",
            "", "no_kb_match"]


class Reachability(unittest.TestCase):
    def test_every_nmpa_row_is_reachable(self):
        """Audit M7: static contract — each DNA/status-keyed curated row must
        be constructible by at least one code path. Protein-keyed rows
        (IHC-defined: HER2-low RC48, CLDN18.2 zolbetuximab) are exempt by
        design: WGS cannot map them automatically; they belong to the
        orthogonal-testing recommendation."""
        table = list(csv.DictReader(open(NMPA), delimiter="\t"))
        table_keys = {r["biomarker_key"] for r in table
                      if r["keying_modality"] != "protein"}
        keys = set(ci.ALLELE_NMPA_KEYS.values())
        keys |= {"MSI-H positive", "TMB-H (FDA >=13 mut/Mb)",
                 "TMB borderline (10-<13 mut/Mb)",
                 "BRCA1/2 loss (somatic, context)"}
        for row in csv.DictReader(open(CNVGENES), delimiter="\t"):
            if row["call_type"] == "amp":
                keys.add(f"{row['gene']} amp")
        for row in csv.DictReader(open(FUSIONS), delimiter="\t"):
            keys.add(ci.fusion_nmpa_key(row["gene_5"]))
        self.assertEqual(set(), table_keys - keys,
                         "unreachable NMPA rows (dead keys)")
        protein = [r for r in table if r["keying_modality"] == "protein"]
        for r in protein:
            self.assertEqual("false", r["context_only"],
                             "protein rows are IHC label realities; they must "
                             "stay visible in curation docs")

    def test_ntrk_folding(self):
        self.assertEqual(ci.fusion_nmpa_key("NTRK1"), "NTRK fusion resolved")
        self.assertEqual(ci.fusion_nmpa_key("NTRK3"), "NTRK fusion resolved")
        self.assertEqual(ci.fusion_nmpa_key("ALK"), "ALK fusion resolved")


class Behaviors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_erbb2_amp_wires_her2_rows(self):
        bm, th = run_interpret(self.tmp.name, [empty_snv()], [
            ["ERBB2", "focal_amplification", "3.1", "12", "17:35-38(0.5Mb)",
             "yes", "gastric/GEJ", "anti-HER2", "", "II"]], [])
        bio = [b for b in bm if b["biomarker"].startswith("ERBB2")]
        self.assertTrue(bio and bio[0]["amp_tier"] == "II")
        ths = {t["therapy"] for t in th if t["biomarker"].startswith("ERBB2")}
        self.assertIn("trastuzumab + chemo (FL/XELOX)", ths)
        self.assertIn("trastuzumab deruxtecan (T-DXd)", ths)
        req = next(t["requires"] for t in th
                   if t["therapy"].startswith("trastuzumab + chemo"))
        self.assertIn("IHC3+", req)

    def test_ntrk1_fusion_maps_curated_row(self):
        bm, th = run_interpret(self.tmp.name, [empty_snv()], [], [
            ["TPM3--NTRK1", "TPM3", "NTRK1", "chr1:1", "chr1:2", "paired-both"]])
        bio = [b for b in bm if "fusion" in b["biomarker"] and "screening" not in b["class"]]
        self.assertTrue(bio and bio[0]["amp_tier"] == "II")
        ths = {t["therapy"] for t in th}
        self.assertTrue(any("larotrectinib" in t for t in ths))

    def test_brca1_del_is_context_only(self):
        bm, th = run_interpret(self.tmp.name, [empty_snv()], [
            ["BRCA1", "deletion", "-1.5", "0", "17:1-2(1.0Mb)",
             "", "gastric", "PARP context", "", "III"]], [])
        bio = [b for b in bm if b["biomarker"].startswith("BRCA1")][0]
        self.assertEqual(bio["nmpa_key"], "BRCA1/2 loss (somatic, context)")
        row = [t for t in th if "olaparib" in t["therapy"]]
        self.assertTrue(row and "context only" in row[0]["region_or_approval"])

    def test_cnv_tier_from_table_not_blanket_ii(self):
        bm, th = run_interpret(self.tmp.name, [empty_snv()], [
            ["CCNE1", "focal_amplification", "2.0", "8", "19:1-2(0.3Mb)",
             "yes", "gastric/ovarian", "trial", "", "III"]], [])
        bio = [b for b in bm if b["biomarker"].startswith("CCNE1")][0]
        self.assertEqual(bio["amp_tier"], "III")

    def test_cohort_mode_paired_caveat(self):
        row = list(empty_snv())
        row[13] = "v9"  # civic match -> variant is retained by `keep`
        bm, _ = run_interpret(self.tmp.name, [row], [], [],
                              cohort_mode="paired")
        snv = [b for b in bm if b["class"] == "SNV/indel"
               and b["biomarker"].startswith("TEST1")]
        self.assertTrue(snv, bm)
        self.assertIn("matched-normal filtered", snv[0]["caveat"])

    def test_germline_screen_counts_clinvar(self):
        row = list(empty_snv())
        row[11] = "Pathogenic"  # clinvar_sig
        bm, _ = run_interpret(self.tmp.name, [row], [], [],
                              tmb_rows=[["T1", "30", "15.0", "20.0"]])
        tmb = [b for b in bm if b["class"] == "TMB"][0]
        self.assertIn("germline screen: 1 ClinVar-flagged", tmb["detail"])
        self.assertIn("TMB-H (FDA >=13", tmb["biomarker"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
