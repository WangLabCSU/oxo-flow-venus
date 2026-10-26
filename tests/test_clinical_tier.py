#!/usr/bin/env python3
"""Unit + integration tests for clinical_tier.py (issue #30 acceptance).

Run: python3 tests/test_clinical_tier.py   (needs pandas + pysam)
Integration fixtures mirror the CIViC monthly-dump TSV column contract used
by load_civic() and a minimal VEP-style CSQ VCF.
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

import clinical_tier as ct  # noqa: E402

VARIANT_COLS = ["variant_id", "feature_name", "variant", "hgvs_descriptions",
                "vicc_compliant_name", "single_variant_molecular_profile_id",
                "clinvar_ids", "reference_build", "chromosome", "start",
                "stop", "reference_bases", "variant_bases"]
EVIDENCE_COLS = ["molecular_profile_id", "therapies", "disease",
                 "evidence_level", "evidence_type"]
ASSERTION_COLS = ["molecular_profile_id", "therapies", "disease",
                  "amp_category", "assertion_type", "regulatory_approval",
                  "nccn_guideline", "significance"]


def tsv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


# variant_id, gene, name, hgvs, vicc, mp
VARIANTS = [
    ("v1", "EGFR", "E746_A750del", "", "", "mp1"),
    ("v2", "EGFR", "L858R", "", "", "mp2"),
    ("v3", "EGFR", "L747_P753delinsS", "", "", "mp3"),
    ("v4", "EGFR", "V769_D770insASV", "NP_005219.2:p.Val769_Asp770insAlaSerVal", "", "mp4"),
    ("v5", "EGFR", "H773_V774insH", "", "", "mp5"),
    ("v6", "MET", "EX14 skipping", "", "", "mp6"),
    ("v7", "EGFR", "E746_T751del", "", "", "mp7"),
]
EVIDENCE = [
    ("mp1", "erlotinib", "lung non-small cell carcinoma", "A", "Predictive"),
    ("mp2", "osimertinib", "lung non-small cell carcinoma", "A", "Predictive"),
    ("mp3", "afatinib", "lung non-small cell carcinoma", "B", "Predictive"),
    ("mp4", "alectinib", "lung non-small cell carcinoma", "A", "Predictive"),
    ("mp5", "osimertinib", "lung non-small cell carcinoma", "B", "Predictive"),
    ("mp6", "tepotinib", "lung non-small cell carcinoma", "B", "Predictive"),
    ("mp7", "", "lung non-small cell carcinoma", "C", "Predictive"),
]
ASSERTIONS = [("mp2", "", "", "Tier I - Level A", "Predictive", "", "", "")]


class IntervalParsing(unittest.TestCase):
    def test_del_range(self):
        self.assertEqual(ct.parse_p_interval("ENSP:x:p.Glu746_Ala750del"),
                         ("del", 746, 750))

    def test_del_1letter(self):
        self.assertEqual(ct.parse_p_interval("E746_A750del"), ("del", 746, 750))

    def test_delins(self):
        self.assertEqual(ct.parse_p_interval("p.Val600_Lys601delinsGluAla"),
                         ("delins", 600, 601))

    def test_ins(self):
        self.assertEqual(ct.parse_p_interval("p.Met751_Ile752insArg"),
                         ("ins", 751, 752))

    def test_dup_normalized_to_ins(self):
        self.assertEqual(ct.parse_p_interval("p.His773dup"), ("ins", 773, 773))

    def test_single_residue_del(self):
        self.assertEqual(ct.parse_p_interval("p.Lys7del"), ("del", 7, 7))

    def test_substitution_rejected(self):
        self.assertIsNone(ct.parse_p_interval("p.Thr315Ile"))

    def test_frameshift_rejected(self):
        self.assertIsNone(ct.parse_p_interval("p.Gly746ValfsTer12"))

    def test_nonsense_rejected(self):
        self.assertIsNone(ct.parse_p_interval("p.Glu384Ter"))

    def test_nucleotide_form_rejected(self):
        self.assertIsNone(ct.parse_p_interval("c.2235_2236delGG"))
        self.assertIsNone(ct.parse_p_interval("g.55181319_55181320insGGT"))

    def test_prose_rejected(self):
        self.assertIsNone(ct.parse_p_interval("EX14 skipping"))
        self.assertIsNone(ct.parse_p_interval(""))


class CivicIndex(unittest.TestCase):
    def test_name_interval_key(self):
        self.assertEqual(ct.civic_interval_keys("E746_A750del", "", ""),
                         {("del", 746, 750)})

    def test_p_hgvs_key(self):
        self.assertEqual(ct.civic_interval_keys(
            "V769_D770insASV", "NP_005219.2:p.Val769_Asp770insAlaSerVal", ""),
            {("ins", 769, 770)})

    def test_named_event_yields_nothing(self):
        self.assertEqual(ct.civic_interval_keys("EX14 skipping", "", ""), set())

    def test_vicc_source(self):
        self.assertEqual(ct.civic_interval_keys("", "", "p.Glu746_Ala750del"),
                         {("del", 746, 750)})


class Integration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.tmp.name
        tsv(os.path.join(d, "var.tsv"), VARIANT_COLS,
            [(*v[:5], *v[5:], "", "", "", "", "", "", "") for v in VARIANTS])
        tsv(os.path.join(d, "ev.tsv"), EVIDENCE_COLS, EVIDENCE)
        tsv(os.path.join(d, "as.tsv"), ASSERTION_COLS, ASSERTIONS)

        csq_def = "SYMBOL|Consequence|HGVSp"
        recs = [
            # chrom pos id ref alt -> CSQ, extra INFO ; expected (civic, tier, basis)
            ("chr7", 55174721, ".", "A", "T", "EGFR|inframe_deletion|ENSP00000275493.2:p.Glu746_Ala750del", "", ("E746_A750del", "I", "interval:del:746-750")),
            ("chr7", 55191822, ".", "C", "T", "EGFR|missense_variant|ENSP00000275493.2:p.Leu858Arg", ";CLNSIG=Pathogenic", ("L858R", "I", "protein:EGFR:L858R")),
            ("chr7", 55174715, ".", "AC", "T", "EGFR|inframe_deletion|ENSP00000275493.2:p.Leu747_Pro753delinsSer", "", ("L747_P753delinsS", "I", "interval:delins:747-753")),
            ("chr7", 55181296, ".", "A", "AGGT", "EGFR|inframe_insertion|ENSP00000275493.2:p.Val769_Asp770insAlaSerVal", "", ("V769_D770insASV", "I", "interval:ins:769-770")),
            ("chr7", 55181290, ".", "A", "AGGT", "EGFR|inframe_insertion|ENSP00000275493.2:p.His773dup", "", ("H773_V774insH", "I", "interval:ins:773-774")),
            ("chr7", 55174718, ".", "AGT", "A", "EGFR|inframe_deletion|ENSP00000275493.2:p.Glu749_Ala750del", "", ("", "III", "no_kb_match")),
            ("chr7", 55181378, ".", "A", "ACGG", "MET|inframe_insertion|ENSP00000286626.4:p.Met751_Ile752insArg", "", ("", "III", "no_kb_match")),
            ("chr7", 55249071, ".", "G", "A", "EGFR|missense_variant|ENSP00000275493.2:p.Thr790Met", "", ("", "III", "no_kb_match")),
        ]
        vcf = ["##fileformat=VCFv4.2", "##contig=<ID=chr7>",
               "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Format: " + csq_def + "\">",
               "##INFO=<ID=AF,Number=A,Type=Float,Description=\"gnomAD AF\">",
               "##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"ClinVar significance\">",
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"G\">",
               "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">",
               "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"]
        for c, p, i, r, a, csq, extra, _ in recs:
            vcf.append(f"{c}\t{p}\t{i}\t{r}\t{a}\t100\tPASS\tCSQ={csq}{extra}\tGT:AD:DP\t0/1:30,70:100")
        cls.vcf_path = os.path.join(d, "in.vcf")
        with open(cls.vcf_path, "w") as fh:
            fh.write("\n".join(vcf) + "\n")

        cls.out_path = os.path.join(d, "out.tsv")
        import subprocess
        script = os.path.join(os.path.dirname(__file__), "..", "scripts",
                              "clinical_tier.py")
        proc = subprocess.run([sys.executable, script, "--vcf", cls.vcf_path,
                               "--civic-variant-tsv", os.path.join(d, "var.tsv"),
                               "--civic-evidence-tsv", os.path.join(d, "ev.tsv"),
                               "--civic-assertion-tsv", os.path.join(d, "as.tsv"),
                               "--out", cls.out_path],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise AssertionError(
                f"clinical_tier.py failed rc={proc.returncode}\n"
                f"stderr:\n{proc.stderr}")
        import csv
        with open(cls.out_path) as fh:
            cls.rows = list(csv.DictReader(fh, delimiter="\t"))
        cls.by_pos = {r["pos"]: r for r in cls.rows}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _check(self, pos, civic, tier, basis):
        r = self.by_pos[pos]
        self.assertEqual(r["civic_variant"], civic, f"pos {pos}")
        self.assertEqual(r["amp_tier"], tier, f"pos {pos}")
        self.assertEqual(r["tier_basis"], basis, f"pos {pos}")

    def test_egfr_exon19_del_matches(self):
        self._check("55174721", "E746_A750del", "I", "interval:del:746-750")

    def test_substitution_exact_and_assertion_amp(self):
        self._check("55191822", "L858R", "I", "protein:EGFR:L858R")

    def test_delins_matches(self):
        self._check("55174715", "L747_P753delinsS", "I", "interval:delins:747-753")

    def test_insertion_matches(self):
        self._check("55181296", "V769_D770insASV", "I", "interval:ins:769-770")

    def test_dup_query_unifies_with_ins_variant(self):
        self._check("55181290", "H773_V774insH", "I", "interval:ins:773-774")

    def test_unindexed_named_event_no_match(self):
        self._check("55181378", "", "III", "no_kb_match")
        self.assertEqual(self.by_pos["55181378"]["hgvsp"], "Met751_Ile752insArg")

    def test_no_match_still_tier3(self):
        self._check("55249071", "", "III", "no_kb_match")

    def test_edge_grazing_interval_rejected(self):
        # audit M6: query del 749-750 overlaps CIViC 746-750 by only 2/5
        # residues — must NOT inherit the hotspot's tier/therapy mapping
        self._check("55174718", "", "III", "no_kb_match")

    def test_clnsig_scalar_not_charjoined(self):
        # audit M1: pysam returns a scalar str for Number=. String INFO;
        # char-joining it used to corrupt the value and zero the germline screen
        self.assertEqual(self.by_pos["55191822"]["clinvar_sig"], "Pathogenic")


if __name__ == "__main__":
    unittest.main(verbosity=2)
