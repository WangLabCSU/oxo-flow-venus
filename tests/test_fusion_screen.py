#!/usr/bin/env python3
"""Tests for the RNA fusion screen (N12, audit 20261104).

Covers: Chimeric.out.junction per-read parsing (14/15-col), mirror-junction
collapsing, unique-read aggregation, anchor math from p-extended CIGARs,
GENCODE gene annotation (+10kb promoter window), panel matching (pair +
wildcard single-gene), anchor/split-read floors, empty junction -> header-only,
cohort concat ordering, and the screening disclosure column.

Run: python3 tests/test_fusion_screen.py   (stdlib only)
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

import collect_fusion  # noqa: E402
import fusion_screen as fs  # noqa: E402

# 15-column rows (with RG) — one row per split READ.
# cols: d_chr d_junc d_strand a_chr a_junc a_strand motif hom_l hom_r read
#       d_start d_cigar a_start a_cigar RG
JUNCTIONS = (
    # EGFR--DGKB-like intra-chromosomal: 3 unique reads (one mirrored,
    # one duplicate read with a better CIGAR)
    "chr7\t55086714\t+\tchr7\t55141787\t-\t1\t0\t0\tr1\t55086702\t"
    "50M251p20M\t55141775\t12M100p30M\tRGX\n"
    "chr7\t55141787\t-\tchr7\t55086714\t+\t1\t0\t0\tr2\t55141775\t"
    "15M30p40M\t55086702\t25M50p50M\tRGX\n"
    "chr7\t55086714\t+\tchr7\t55141787\t-\t1\t0\t0\tr3\t55086702\t"
    "12M99p12M\t55141775\t20M10p20M\tRGX\n"
    # TET1--MSI2 inter-chromosomal panel pair, 1 read (kept: panel hit)
    "chr14\t90000000\t+\tchr17\t40000000\t+\t0\t0\t0\tx1\t89999990\t"
    "30M20p20M\t39999990\t20M10p20M\tRGX\n"
    # anchor below floor (5/5) — must be dropped
    "chr1\t1000000\t+\tchr2\t2000000\t+\t1\t0\t0\ty1\t999990\t"
    "5M10p5M\t1999990\t5M5p5M\tRGX\n"
    # no gene within 10kb on either side, 2 reads, no panel hit -> kept but
    # intergenic-labelled? no: no gene pair and <3 reads -> dropped; make it
    # 3 reads so it survives as intergenic
    "chr3\t1000000\t+\tchr4\t2000000\t+\t1\t0\t0\tz1\t999990\t"
    "30M10p30M\t1999990\t30M10p30M\tRGX\n"
    "chr3\t1000000\t+\tchr4\t2000000\t+\t1\t0\t0\tz2\t999990\t"
    "30M10p30M\t1999990\t30M10p30M\tRGX\n"
    "chr3\t1000000\t+\tchr4\t2000000\t+\t1\t0\t0\tz3\t999990\t"
    "30M10p30M\t1999990\t30M10p30M\tRGX\n")

GTF = (
    "#!genome-build GRCh38\n"
    "chr7\tensembl\tgene\t55086000\t55087000\t.\t+\t.\t"
    'gene_id "ENSGA"; gene_name "EGFR";\n'
    "chr7\tensembl\tgene\t55141000\t55142000\t.\t-\t.\t"
    'gene_id "ENSGD"; gene_name "DGKB";\n'
    "chr14\tensembl\tgene\t89999000\t90010000\t.\t+\t.\t"
    'gene_id "ENSGT"; gene_name "TET1";\n'
    "chr17\tensembl\tgene\t39990000\t40010000\t.\t+\t.\t"
    'gene_id "ENSGM"; gene_name "MSI2";\n')

PANEL = ("gene_5\tgene_3\tdisease_context\tevidence_context\n"
         "EGFR\t*\tgastric\tNCCN\n"
         "TET1\tMSI2\tglioma\tcase_report\n")

HEADERONLY_HEADER = fs.COLUMNS


def write_fixtures(d):
    j = os.path.join(d, "S.junction")
    g = os.path.join(d, "mini.gtf")
    p = os.path.join(d, "panel.tsv")
    for path, payload in ((j, JUNCTIONS), (g, GTF), (p, PANEL)):
        with open(path, "w") as fh:
            fh.write(payload)
    return j, g, p


def run_screen(tmp, junction="S.junction", min_anchor=None,
               min_split=None):
    j, g, p = write_fixtures(tmp)
    out = os.path.join(tmp, "cand.tsv")
    args = [sys.executable, os.path.join(REPO, "scripts", "fusion_screen.py"),
            "--junction", os.path.join(tmp, junction), "--gtf", g,
            "--panel", p, "--sample", "Pt01-TR", "--output", out]
    if min_anchor is not None:
        args += ["--min-anchor", str(min_anchor)]
    if min_split is not None:
        args += ["--min-split-reads", str(min_split)]
    proc = subprocess.run(args, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Anchors(unittest.TestCase):
    def test_anchor_len_from_cigarp(self):
        self.assertEqual(fs.anchor_len("50M251p20M"), 70)
        self.assertEqual(fs.anchor_len("30M10p30M"), 60)
        self.assertEqual(fs.anchor_len("20M"), 20)
        self.assertEqual(fs.anchor_len("10S30M200N20M5S"), 50)  # S/N excluded
        self.assertEqual(fs.anchor_len(""), 0)
        self.assertEqual(fs.anchor_len("*"), 0)


class PanelMatch(unittest.TestCase):
    PANEL = [{"gene_5": "EGFR", "gene_3": "*"},
             {"gene_5": "TET1", "gene_3": "MSI2"}]

    def test_exact_pair(self):
        hit, label = fs.panel_match(self.PANEL, "TET1", "MSI2")
        self.assertTrue(hit)
        self.assertEqual(label, "TET1--MSI2")

    def test_wildcard_with_partner(self):
        hit, label = fs.panel_match(self.PANEL, "EGFR", "DGKB")
        self.assertTrue(hit)
        self.assertIn("wildcard", label)

    def test_wildcard_unresolved_partner(self):
        hit, label = fs.panel_match(self.PANEL, "EGFR", "")
        self.assertTrue(hit)
        self.assertIn("partner unresolved", label)

    def test_no_match(self):
        hit, _ = fs.panel_match(self.PANEL, "KRAS", "EGFR2")
        self.assertFalse(hit)

    def test_wildcard_never_matches_both_wildcards(self):
        # both sides '*' rows cannot produce a gene-specific call
        hit, _ = fs.panel_match([{"gene_5": "*", "gene_3": "*"}], "A", "B")
        self.assertFalse(hit)


class Screen(unittest.TestCase):
    def test_aggregation_annotation_filters(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = run_screen(tmp)
        self.assertEqual(len(rows), 3)
        # mirror junction collapsed: 3 unique reads (r1/r2/r3)
        egfr = next(r for r in rows if r["gene_a"] == "EGFR")
        self.assertEqual(int(egfr["split_reads"]), 3)
        self.assertEqual(egfr["motif"], "GT/AG")
        self.assertEqual(int(egfr["anchor_a"]), 75)   # r2's better CIGAR
        self.assertEqual(int(egfr["anchor_b"]), 55)
        self.assertEqual(egfr["panel_hit"], "yes")
        self.assertEqual(egfr["evidence_class"], "panel_hit")
        # 1-read inter-chromosomal panel pair kept as single-read panel hit
        tet1 = next(r for r in rows if r["gene_a"] == "TET1")
        self.assertEqual(tet1["evidence_class"], "panel_hit_single_read")
        self.assertEqual(tet1["panel_label"], "TET1--MSI2")
        # intergenic 3-read junction kept with intergenic label
        ig = [r for r in rows if r["gene_a"] == "" and r["gene_b"] == ""]
        self.assertEqual(len(ig), 1)
        self.assertTrue(ig[0]["fusion_pair"].startswith("intergenic_"))
        self.assertIn("screening_level", ig[0]["disclosure"])
        # every row carries the screening disclosure
        self.assertTrue(all("screening_level" in r["disclosure"] for r in rows))
        # panel hits sort first
        self.assertEqual(rows[0]["panel_hit"], "yes")

    def test_stricter_floor_applies_to_generic_candidates(self):
        # min-split-reads governs NON-panel candidates: at floor 5 the 3-read
        # generic junction is DROPPED, while actionable single-read panel hits
        # survive (kept + disclosed via evidence_class, never silently gone)
        with tempfile.TemporaryDirectory() as tmp:
            strict = run_screen(tmp, min_split=5)
            default = run_screen(tmp)
        self.assertFalse([r for r in strict if r["panel_hit"] != "yes"])
        self.assertTrue(any(r["evidence_class"] == "panel_hit_single_read"
                            for r in strict))
        # default floor keeps the 3-read generic junction as supported
        self.assertTrue(any(r["evidence_class"] == "supported_3reads"
                            for r in default))

    def test_starfusion_format_header_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            j, g, p = write_fixtures(tmp)
            with open(os.path.join(tmp, "hdr.junction"), "w") as fh:
                fh.write("chr_donorA\tbrkpt_donorA\tstrand_donorA\tmore\n")
                fh.write(JUNCTIONS)
            out = os.path.join(tmp, "o.tsv")
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts",
                                              "fusion_screen.py"),
                 "--junction", os.path.join(tmp, "hdr.junction"),
                 "--gtf", g, "--panel", p, "--sample", "S",
                 "--output", out], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(len(rows), 3)

    def test_empty_junction_header_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            j, g, p = write_fixtures(tmp)
            with open(os.path.join(tmp, "empty.junction"), "w") as fh:
                fh.write("")
            out = os.path.join(tmp, "o.tsv")
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts",
                                              "fusion_screen.py"),
                 "--junction", os.path.join(tmp, "empty.junction"),
                 "--gtf", g, "--panel", p, "--sample", "S",
                 "--output", out], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out) as fh:
                lines = fh.read().strip().split("\n")
        self.assertEqual(lines[0].split("\t"), HEADERONLY_HEADER)
        self.assertEqual(len(lines), 1)


class Collect(unittest.TestCase):
    def test_concat_order_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = os.path.join(tmp, "A.fusion_candidates.tsv")
            b = os.path.join(tmp, "B.fusion_candidates.tsv")
            with open(a, "w") as fh:
                fh.write("sample\tfusion_pair\tsplit_reads\tpanel_hit\t"
                         "evidence_class\n"
                         "A\tX--Y\t3\tno\tsupported_3reads\n"
                         "A\tP1--P2\t1\tyes\tpanel_hit_single_read\n")
            with open(b, "w") as fh:
                fh.write("sample\tfusion_pair\tsplit_reads\tpanel_hit\t"
                         "evidence_class\n"
                         "B\tQ\t5\tyes\tpanel_hit\n")
            out = os.path.join(tmp, "cohort.tsv")
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts",
                                              "collect_fusion.py"),
                 "--inputs", a, b, "--output", out],
                capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(len(rows), 3)
        # panel hits first, then by descending split_reads
        self.assertEqual([r["sample"] for r in rows], ["B", "A", "A"])
        self.assertEqual(int(rows[0]["split_reads"]), 5)

    def test_empty_input_list_header_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "cohort.tsv")
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, "scripts",
                                              "collect_fusion.py"),
                 "--inputs", "--output", out],
                capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out) as fh:
                self.assertTrue(fh.readline().startswith("sample"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
