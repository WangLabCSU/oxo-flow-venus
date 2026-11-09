#!/usr/bin/env python3
"""Offline tests for the MAF-domain tools (issue #33 N9).

Covers vcf2maf_light (CSQ header parsing, indel anchor trimming and MAF
coordinates, classification incl. frameshift by type, sample-column
resolution, ref==alt skip) and tmb_from_maf (dedupe key, coding vs total
denominators, BED merging, gzip input, malformed-row skip).

Run: python3 tests/test_maf_tools.py   (stdlib only)
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

import tmb_from_maf  # noqa: E402
import vcf2maf_light as v2m  # noqa: E402

# vep-rs style CSQ header (Format: inside Description)
CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations'
    ' from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature|HGVSc|'
    'HGVSp">\n')
CHROM_LINE = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n")


def write_vcf(path, body, gz=False):
    text = CSQ_HEADER + CHROM_LINE + body
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)


def run_vcf2maf(tmp, body, gz=False, tumor="TUMOR", normal="NORMAL"):
    vcf = os.path.join(tmp, "in.vcf.gz" if gz else "in.vcf")
    write_vcf(vcf, body, gz=gz)
    out = os.path.join(tmp, "out.maf")
    proc = subprocess.run(
        [sys.executable, os.path.join(REPO, "scripts", "vcf2maf_light.py"),
         "--vcf", vcf, "--tumor-id", tumor, "--normal-id", normal,
         "--output", out],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    with open(out) as fh:
        lines = fh.read().splitlines()
    header = lines[2].split("\t")
    rows = [dict(zip(header, ln.split("\t"))) for ln in lines[3:]]
    return rows, proc


class Vcf2mafHelpers(unittest.TestCase):
    def test_csq_format_header_variants(self):
        f1 = v2m.parse_csq_format(
            '##INFO=<ID=CSQ,Type=String,Description="x. Format: A|B|C">')
        self.assertEqual(f1, ["A", "B", "C"])
        f2 = v2m.parse_csq_format(
            '##INFO=<ID=CSQ,Number=.,Type=String,Format=A|B|C">')
        self.assertEqual(f2, ["A", "B", "C"])
        self.assertEqual(v2m.parse_csq_format("##INFO=<ID=DP>"), [])

    def test_trim_alleles_and_types(self):
        # 100 AGT->A (deletion): trim common suffix G? no — suffix 'T' vs 'A'
        # differ, so prefix trim: A|AGT -> GT vs '' with prefix 1
        ref, alt, prefix = v2m.trim_alleles("A", "AGT")
        self.assertEqual((ref, alt, prefix), ("", "GT", 1))
        # 100 AGT->A the other order
        ref, alt, prefix = v2m.trim_alleles("AGT", "A")
        self.assertEqual((ref, alt, prefix), ("GT", "", 1))
        self.assertEqual(v2m.variant_type("GT", ""), "DEL")
        self.assertEqual(v2m.variant_type("", "GT"), "INS")
        self.assertEqual(v2m.variant_type("A", "T"), "SNP")
        self.assertEqual(v2m.variant_type("AT", "GC"), "DNP")
        self.assertEqual(v2m.variant_type("ATG", "CCA"), "TNP")
        self.assertEqual(v2m.variant_type("ATGC", "GGGG"), "ONP")

    def test_classify_frameshift_by_type(self):
        self.assertEqual(v2m.classify("frameshift_variant", "DEL", "GT", ""),
                         "Frame_Shift_Del")
        self.assertEqual(v2m.classify("frameshift_variant", "INS", "", "GT"),
                         "Frame_Shift_Ins")
        # compound terms resolve to the first classed term
        self.assertEqual(
            v2m.classify("frameshift_variant&splice_region_variant", "DEL",
                         "GT", ""),
            "Frame_Shift_Del")
        self.assertEqual(v2m.classify("missense_variant", "SNP", "A", "T"),
                         "Missense_Mutation")
        self.assertEqual(v2m.classify("intergenic_variant", "SNP", "A", "T"),
                         "IGR")
        self.assertEqual(v2m.classify("totally_unknown_term", "SNP", "A", "T"),
                         "Targeted_Region")
        self.assertEqual(v2m.classify("", "SNP", "A", "T"), "IGR")

    def test_resolve_sample_columns(self):
        # tumor/normal resolved BY NAME, not position
        self.assertEqual(
            v2m.resolve_sample_columns(["NORMAL", "TUMOR"], "TUMOR", "NORMAL"),
            (1, 0))
        self.assertEqual(
            v2m.resolve_sample_columns(["TUMOR", "NORMAL"], "TUMOR", "NORMAL"),
            (0, 1))
        # tumor-only: one column, normal -1
        self.assertEqual(
            v2m.resolve_sample_columns(["TUMOR"], "TUMOR", ""), (0, -1))


class Vcf2mafEndToEnd(unittest.TestCase):
    def test_snp_row_schema_and_vaf(self):
        with tempfile.TemporaryDirectory() as tmp:
            body = (
                "7\t55191822\trs123\tG\tA\t.\tPASS\t"
                "CSQ=A|missense_variant|EGFR|ENST1|c.12A>C|p.T4A;"
                "DP=100;AD=70,30\tGT:DP:AD\t"
                "0/1:100:70,30\t0/1:100:70,30\n")
            rows, proc = run_vcf2maf(tmp, body)
            self.assertEqual(len(rows), 1)
            r = rows[0]
            self.assertEqual(r["Hugo_Symbol"], "EGFR")
            self.assertEqual(r["Chromosome"], "7")
            self.assertEqual(r["Start_Position"], "55191822")
            self.assertEqual(r["End_Position"], "55191822")
            self.assertEqual(r["Reference_Allele"], "G")
            self.assertEqual(r["Tumor_Seq_Allele1"], "G")
            self.assertEqual(r["Tumor_Seq_Allele2"], "A")
            self.assertEqual(r["Variant_Classification"], "Missense_Mutation")
            self.assertEqual(r["Variant_Type"], "SNP")
            self.assertEqual(r["Tumor_Sample_Barcode"], "TUMOR")
            self.assertEqual(r["Matched_Norm_Sample_Barcode"], "NORMAL")
            self.assertEqual(r["dbSNP_RS"], "rs123")
            self.assertEqual(r["t_depth"], "100")
            self.assertEqual(r["t_ref_count"], "70")
            self.assertEqual(r["t_alt_count"], "30")
            self.assertEqual(r["t_vaf"], "0.3000")
            self.assertEqual(r["HGVSc"], "c.12A>C")
            self.assertEqual(r["HGVSp"], "p.T4A")
            self.assertEqual(r["Transcript_ID"], "ENST1")

    def test_deletion_coordinates_and_indel_padding(self):
        with tempfile.TemporaryDirectory() as tmp:
            # VCF del at 100: REF=AGT ALT=A -> MAF Start=101 End=102, GT vs -
            body = (
                "1\t100\t.\tAGT\tA\t.\tPASS\tCSQ=A|frameshift_variant|KRAS|"
                "ENST2|c.1del|p.R2fs;DP=50;AD=30,20\tGT:DP:AD\t"
                "0/1:50:30,20\t0/1:50:30,20\n")
            rows, _ = run_vcf2maf(tmp, body)
            r = rows[0]
            self.assertEqual(r["Start_Position"], "101")
            self.assertEqual(r["End_Position"], "102")
            self.assertEqual(r["Reference_Allele"], "GT")
            self.assertEqual(r["Tumor_Seq_Allele2"], "-")
            self.assertEqual(r["Variant_Type"], "DEL")
            self.assertEqual(r["Variant_Classification"], "Frame_Shift_Del")

    def test_insertion_and_gz_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            # VCF ins at 100: A -> AGT -> MAF Start=End=101, - vs GT
            body = (
                "2\t100\t.\tA\tAGT\t.\tPASS\tCSQ=G|frameshift_variant|MET|"
                "ENST3|c.1ins|p.V2insI;DP=40;AD=20,20\tGT:DP:AD\t"
                "0/1:40:20,20\t0/1:40:20,20\n")
            rows, proc = run_vcf2maf(tmp, body, gz=True)
            self.assertIn("wrote 1 variants", proc.stderr)
            r = rows[0]
            self.assertEqual(r["Start_Position"], "101")
            self.assertEqual(r["End_Position"], "101")
            self.assertEqual(r["Reference_Allele"], "-")
            self.assertEqual(r["Tumor_Seq_Allele2"], "GT")
            self.assertEqual(r["Variant_Classification"], "Frame_Shift_Ins")

    def test_ref_equals_alt_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            body = ("1\t100\t.\tA\tA\t.\tPASS\t;\tGT:DP:AD\t0/1:10:5,5\n")
            rows, proc = run_vcf2maf(tmp, body)
            self.assertEqual(rows, [])
            self.assertIn("wrote 0 variants", proc.stderr)


class TmbHelpers(unittest.TestCase):
    def test_merge_bed_overlaps_and_bookends(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = os.path.join(tmp, "t.bed")
            with open(bed, "w") as fh:
                # chr1: [100,200) + [150,250) + bookend [250,300) -> 200
                # chr2: [0,100); the degenerate [500,500) is excluded from both
                fh.write("chr1\t100\t200\n"
                         "chr1\t150\t250\n"
                         "chr1\t250\t300\n"
                         "chr2\t0\t100\n"
                         "chr2\t500\t500\n")
            merged, raw = tmb_from_maf.merge_bed(bed)
            self.assertEqual(raw, 350)  # 100+100+50+100 (degenerate dropped)
            self.assertEqual(merged, 300)

    def test_read_maf_skips_comments_and_malformed(self):
        with tempfile.TemporaryDirectory() as tmp:
            maf = os.path.join(tmp, "x.maf")
            with open(maf, "w") as fh:
                fh.write("## banner\n#version 2.4.1\n")
                fh.write("Chromosome\tStart_Position\tEnd_Position\t"
                         "Reference_Allele\tTumor_Seq_Allele2\t"
                         "Tumor_Sample_Barcode\tVariant_Type\t"
                         "Variant_Classification\n")
                fh.write("1\t10\t10\tA\tG\tT1\tSNP\tMissense_Mutation\n")
                fh.write("short-but-malformed\n")
            rows = tmb_from_maf.read_maf(maf)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["Tumor_Sample_Barcode"], "T1")


def maf_text(rows):
    cols = ("Chromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
            "Tumor_Seq_Allele2\tTumor_Sample_Barcode\tVariant_Type\t"
            "Variant_Classification\n")
    return cols + "".join("\t".join(r) + "\n" for r in rows)


class TmbEndToEnd(unittest.TestCase):
    def test_counts_denominators_and_dedupe(self):
        with tempfile.TemporaryDirectory() as tmp:
            m1 = os.path.join(tmp, "a.maf")
            with open(m1, "w") as fh:
                fh.write(maf_text([
                    # duplicate key must be counted once
                    ("1", "100", "100", "A", "G", "Pt01-TD", "SNP",
                     "Missense_Mutation"),
                    ("1", "100", "100", "A", "G", "Pt01-TD", "SNP",
                     "Missense_Mutation"),
                    # silent -> not coding, not SNP-typed? SNP yes
                    ("2", "200", "200", "C", "T", "Pt01-TD", "SNP", "Silent"),
                    ("3", "300", "301", "AT", "-", "Pt01-TD", "DEL",
                     "Frame_Shift_Del"),
                ]))
            m2 = os.path.join(tmp, "b.maf")
            with open(m2, "w") as fh:
                fh.write(maf_text([
                    ("7", "55191822", "55191822", "G", "A", "Pt02-TD", "SNP",
                     "Nonsense_Mutation"),
                ]))
            out = os.path.join(tmp, "tmb.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "tmb_from_maf.py"),
                 "--mafs", m1, m2,
                 "--coding-mb", "34", "--target-mb", "3000",
                 "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            import csv
            with open(out) as fh:
                rows = {r["sample"]: r for r in csv.DictReader(fh,
                                                               delimiter="\t")}
            p1 = rows["Pt01-TD"]
            self.assertEqual(p1["total_pass_variants"], "3")  # dedupe dropped 1
            self.assertEqual(p1["coding_non_silent"], "2")  # missense + fs
            self.assertEqual(p1["snps"], "2")
            self.assertEqual(p1["indels"], "1")
            self.assertEqual(p1["tmb_mut_per_mb"], f"{2/34:.3f}")
            self.assertEqual(p1["all_somatic_mut_per_mb"], f"{3/3000:.3f}")
            self.assertIn("Pt02-TD", rows)

    def test_target_bed_shrinks_coding_mb(self):
        with tempfile.TemporaryDirectory() as tmp:
            maf = os.path.join(tmp, "a.maf")
            with open(maf, "w") as fh:
                fh.write(maf_text([
                    ("1", "100", "100", "A", "G", "T1", "SNP",
                     "Missense_Mutation"),
                ]))
            bed = os.path.join(tmp, "targets.bed")
            with open(bed, "w") as fh:
                fh.write("1\t0\t1000000\n")  # 1 Mb < 34 Mb coding default
            out = os.path.join(tmp, "tmb.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "tmb_from_maf.py"),
                 "--mafs", maf, "--target-bed", bed, "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            with open(out) as fh:
                body = [ln.split("\t") for ln in fh.read().splitlines()[1:]]
            # coding territory clamped to the captured 1.0 Mb; target also 1.0
            self.assertEqual(body[0][5], "1.00")   # coding_mb
            self.assertEqual(body[0][6], "1.00")   # target_mb
            self.assertEqual(body[0][7], "1.000")  # 1 mut / 1 Mb

    def test_empty_maf_zero_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            maf = os.path.join(tmp, "a.maf")
            with open(maf, "w") as fh:
                fh.write(maf_text([]))
            out = os.path.join(tmp, "tmb.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "tmb_from_maf.py"),
                 "--mafs", maf, "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            import csv
            with open(out) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            # no samples -> header-only table (no fabricated zero rows)
            self.assertEqual(rows, [])
            with open(out) as fh:
                self.assertTrue(fh.readline().startswith("sample\t"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
