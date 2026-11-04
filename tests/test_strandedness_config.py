#!/usr/bin/env python3
"""Tests for N14 (featureCounts strandedness config) + workflow TOML hygiene.

Covers: every .oxoflow fragment and envs/*/pixi.toml parses with tomllib; the
`featurecounts_strandedness` config key exists in venus.oxoflow with the
unstranded default "0"; the featurecounts rule actually references the key;
star_align carries the chimera flags required by rna_fusion_screen (N12);
the qc.oxoflow wgs_* rules and the multiqc anchors/scan roots for qc/alignqc
are wired (N10).

Run: python3 tests/test_strandedness_config.py   (stdlib only)
"""
import glob
import os
import tomllib
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)


def load(path):
    with open(os.path.join(REPO, path), "rb") as fh:
        return tomllib.load(fh)


def rules_of(fragment):
    data = load(os.path.join("rules", fragment))
    return {r["name"]: r for r in data.get("rules", [])}


class TomlValid(unittest.TestCase):
    def test_all_oxoflow_parse(self):
        for path in glob.glob(os.path.join(REPO, "rules", "*.oxoflow")) + \
                [os.path.join(REPO, "venus.oxoflow")]:
            with self.subTest(path=os.path.basename(path)):
                tomllib.load(open(path, "rb"))

    def test_all_pixi_toml_parse(self):
        for path in glob.glob(os.path.join(REPO, "envs", "*", "pixi.toml")):
            with self.subTest(path=os.path.basename(os.path.dirname(path))):
                tomllib.load(open(path, "rb"))


class N14Strandedness(unittest.TestCase):
    def test_config_key_exists_with_unstranded_default(self):
        cfg = load("venus.oxoflow")["config"]
        self.assertIn("featurecounts_strandedness", cfg)
        self.assertEqual(str(cfg["featurecounts_strandedness"]), "0")

    def test_rule_references_config_key(self):
        rule = rules_of("rna.oxoflow")["featurecounts"]
        shell = rule["shell"]
        self.assertIn("{config.featurecounts_strandedness}", shell)
        self.assertIn("-s {config.featurecounts_strandedness}", shell)

    def test_documented_in_readme(self):
        readme = open(os.path.join(REPO, "README.md")).read()
        self.assertIn("featurecounts_strandedness", readme)


class N12ChimeraWiring(unittest.TestCase):
    def test_star_align_emits_chimeric_junction(self):
        rule = rules_of("rna.oxoflow")["star_align"]
        self.assertIn("Chimeric.out.junction", " ".join(rule["output"]))
        self.assertIn("Chimeric.out.junction", " ".join(rule["protected_output"]))

    def test_star_align_chimera_flags(self):
        shell = rules_of("rna.oxoflow")["star_align"]["shell"]
        self.assertIn("--chimSegmentMin 12", shell)
        # STAR writes Chimeric.out.junction ONLY when "Junctions" is in
        # --chimOutType (default WithinBAM HardClip writes NONE)
        self.assertIn("--chimOutType WithinBAM SoftClip Junctions", shell)
        self.assertIn("--chimJunctionOverhangMin 12", shell)

    def test_fusion_rules_exist(self):
        rna = rules_of("rna.oxoflow")
        self.assertIn("rna_fusion_screen", rna)
        self.assertIn("rna_fusion_cohort", rna)


class N10N11Wiring(unittest.TestCase):
    def test_wgs_qc_rules_exist_with_gates(self):
        qc = rules_of("qc.oxoflow")
        for name in ("wgs_flagstat", "wgs_bamstats", "wgs_mosdepth"):
            self.assertIn(name, qc, name)
            self.assertEqual(qc[name].get("when"),
                             "wildcard.assay == 'wgs'", name)
        # the collector is cohort-level: no per-sample gate, but it must fan
        # out over the DNA group and produce a header-only table when empty
        self.assertIn("wgs_alignqc_summary", qc)
        self.assertIsNone(qc["wgs_alignqc_summary"].get("when"))
        expand = {e["pattern"] for e in qc["wgs_alignqc_summary"]["expand_inputs"]}
        self.assertIn("qc/alignqc/{sample}.mosdepth.summary.txt", expand)

    def test_mosdepth_no_whole_genome_perbase(self):
        # windowed --by keeps the ~100GB-class per-base output off disk
        shell = rules_of("qc.oxoflow")["wgs_mosdepth"]["shell"]
        self.assertIn("--by 1000000", shell)
        self.assertIn("--thresholds 10,20", shell)
        self.assertNotIn("per-base", shell)

    def test_multiqc_scans_alignqc_with_anchors(self):
        qc = rules_of("qc.oxoflow")
        expand = {e["pattern"] for e in qc["multiqc"]["expand_inputs"]}
        for pat in ("qc/alignqc/{sample}.flagstat",
                    "qc/alignqc/{sample}.stats",
                    "qc/alignqc/{sample}.mosdepth.summary.txt"):
            self.assertIn(pat, expand, pat)
        self.assertIn("qc/alignqc", qc["multiqc"]["shell"])

    def test_report_chain_consumes_new_tables(self):
        report = rules_of("report.oxoflow")
        self.assertIn("qc/alignqc/wgs_alignqc_summary.tsv",
                      report["cohort_tables"]["input"])
        self.assertIn("report/fusion_screen_cohort.tsv",
                      report["clinical_report"]["input"])
        self.assertIn("--alignqc {input[3]}", report["clinical_report"]["shell"])
        self.assertIn("--fusion-screen {input[4]}",
                      report["clinical_report"]["shell"])

    def test_n11_env_has_mosdepth(self):
        deps = load("envs/qc/pixi.toml")["dependencies"]
        self.assertIn("mosdepth", deps)
        self.assertIn("samtools", deps)


if __name__ == "__main__":
    unittest.main(verbosity=2)
