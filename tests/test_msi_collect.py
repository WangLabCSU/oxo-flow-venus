#!/usr/bin/env python3
"""Offline tests for collect_msi (issue #33 N9).

Covers the BH recompute from per-site p-values, the uint16 saturation flag
(MSIsensor-msi v1.3.0 rank wrap at 65535), pair matching by file name,
unmatched-input warnings, cutoff classification (WGS 3.5% / WES 15%), and
malformed per-site row handling.

Run: python3 tests/test_msi_collect.py   (stdlib only)
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "scripts"))

import collect_msi  # noqa: E402


class BhHelpers(unittest.TestCase):
    def test_bh_count_sorted_input(self):
        bh = collect_msi.bh_unstable_count
        # 5 tests: q = [0.05, 0.05, 0.05, 0.125, 0.20] -> first 3 significant
        self.assertEqual(bh([0.01, 0.02, 0.03, 0.10, 0.20]), 3)
        # single tiny p is significant; single large p is not
        self.assertEqual(bh([0.001]), 1)
        self.assertEqual(bh([0.5]), 0)
        self.assertEqual(bh([]), 0)

    def test_bh_ties_can_wipe_significance(self):
        bh = collect_msi.bh_unstable_count
        # 10 tied p=0.001 among n=1000: q_1 = 0.001*1000/1 = 1.0 -> none pass
        self.assertEqual(bh([0.001] * 10 + [0.9] * 990), 0)
        # 10x rarer events do pass: q_1 = 0.00001*1000 = 0.01
        self.assertEqual(bh([0.00001] * 10 + [0.9] * 990), 10)

    def test_bh_monotone_boundary(self):
        bh = collect_msi.bh_unstable_count
        # 20 sites, exactly 1 at p=0.002: q = 0.002*20/1 = 0.04 <= 0.05
        self.assertEqual(bh([0.002] + [0.9] * 19), 1)
        # borderline: p=0.0026 -> q = 0.052 > 0.05 -> nothing significant
        self.assertEqual(bh([0.0026] + [0.9] * 19), 0)


class EndToEnd(unittest.TestCase):
    SUMMARY_SAT = "100000\t65535\t65.53\n"      # uint16-saturated report
    SUMMARY_OK = "1000\t10\t1.0\n"              # honest small report

    RAW_HEAD = "chr\tpos\tscore\tread_depth_tumor\tread_depth_normal\t"
    RAW_2 = "repeat\tP_value\n"

    @staticmethod
    def raw_rows(n, p):
        body = "".join(
            f"chr1\t{i}\t0.9\t40\t40\t-\tA.15A.16\t{p}\n" for i in range(n))
        return "chr\tstart\tend\tscore\tdepth_t\tdepth_n\trepeat\tp\n" + body

    def test_pairing_saturation_and_status(self):
        with tempfile.TemporaryDirectory() as td:
            # saturated pair: reported 65535, BH finds 70000 unstable
            summ_sat = os.path.join(td, "GC01-TD_vs_GC01-B.msi.tsv")
            raw_sat = os.path.join(td, "GC01-TD_vs_GC01-B_all")
            with open(summ_sat, "w") as fh:
                fh.write(self.SUMMARY_SAT)
            with open(raw_sat, "w") as fh:
                # 70000 sites at p=0 -> bh=70000 > 65535 and reported==65535
                fh.write("header\tline\n")
                fh.write("chr1\t1\t0\t40\t40\t-\tA.15A.16\t0\n" * 70000)
            # clean pair: 1000 sites, 10 unstable by BH (p=1e-5 ties),
            # summary reports 10 as well
            summ_ok = os.path.join(td, "GC02-TD_vs_GC02-B.msi.tsv")
            raw_ok = os.path.join(td, "GC02-TD_vs_GC02-B_all")
            with open(summ_ok, "w") as fh:
                fh.write(self.SUMMARY_OK)
            with open(raw_ok, "w") as fh:
                ps = [0.00001] * 10 + [0.9] * 990
                fh.write("chr\tstart\tscore\tdepth_t\tdepth_n\t"
                         "repeat2\trepeat\tp\n")
                for i, p in enumerate(ps):
                    fh.write(f"chr1\t{i}\t0\t40\t40\t-\tA.15A.16\t{p}\n")
            out = os.path.join(td, "msi_cohort.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "collect_msi.py"),
                 "--inputs", summ_sat, raw_sat, summ_ok, raw_ok,
                 "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            import csv
            with open(out) as fh:
                rows = {r["sample"]: r for r in
                        csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(sorted(rows), ["GC01-TD_vs_GC01-B",
                                            "GC02-TD_vs_GC02-B"])
            sat = rows["GC01-TD_vs_GC01-B"]
            self.assertEqual(sat["raw_total_sites"], "70000")
            self.assertEqual(sat["reported_unstable"], "65535")
            self.assertEqual(sat["bh_unstable"], "70000")
            self.assertEqual(sat["uint16_saturated"], "true")
            self.assertEqual(sat["unstable_delta"], "4465")
            self.assertEqual(sat["msi_status_wgs"], "MSI-H")
            self.assertEqual(sat["msi_status_wes"], "MSI-H")
            ok = rows["GC02-TD_vs_GC02-B"]
            # 10/1000 = 1.0% -> MSS at both cutoffs
            self.assertEqual(ok["bh_unstable"], "10")
            self.assertEqual(ok["bh_pct"], "1.00")
            self.assertEqual(ok["uint16_saturated"], "false")
            self.assertEqual(ok["msi_status_wgs"], "MSS")
            self.assertEqual(ok["msi_status_wes"], "MSS")
            self.assertEqual(ok["summary_total_sites"], "1000")
            self.assertEqual(ok["reported_pct"], "1.0")

    def test_wes_cutoff_msi_h(self):
        with tempfile.TemporaryDirectory() as td:
            summ = os.path.join(td, "P1.msi.tsv")
            raw = os.path.join(td, "P1_all")
            with open(summ, "w") as fh:
                fh.write("200\t40\t20.0\n")
            with open(raw, "w") as fh:
                fh.write("chr\tstart\tscore\tdepth_t\tdepth_n\t"
                         "repeat2\trepeat\tp\n")
                for i in range(40):
                    fh.write(f"chr1\t{i}\t0\t40\t40\t-\tA.15A.16\t0.0001\n")
                for i in range(160):
                    fh.write(f"chr1\t{1000+i}\t0\t40\t40\t-\tA.15A.16\t0.9\n")
            out = os.path.join(td, "cohort.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "collect_msi.py"),
                 "--inputs", summ, raw, "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            import csv
            with open(out) as fh:
                r = list(csv.DictReader(fh, delimiter="\t"))[0]
            self.assertEqual(r["bh_unstable"], "40")
            self.assertEqual(r["bh_pct"], "20.00")
            # 20% >= 15% (WES) and >= 3.5% (WGS)
            self.assertEqual(r["msi_status_wgs"], "MSI-H")
            self.assertEqual(r["msi_status_wes"], "MSI-H")

    def test_unmatched_and_unrecognized_inputs(self):
        with tempfile.TemporaryDirectory() as td:
            summ = os.path.join(td, "P1.msi.tsv")
            raw = os.path.join(td, "P2_all")
            orphan = os.path.join(td, "stray.txt")
            for p, text in ((summ, "10\t1\t10.0\n"), (raw, "h\tp\n"),
                            (orphan, "x")):
                with open(p, "w") as fh:
                    fh.write(text)
            out = os.path.join(td, "cohort.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "collect_msi.py"),
                 "--inputs", summ, raw, orphan, "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            # P1 has no _all, P2 has no summary, stray.txt unrecognized
            self.assertIn("P1", proc.stderr)
            self.assertIn("P2", proc.stderr)
            self.assertIn("stray.txt", proc.stderr)
            with open(out) as fh:
                lines = fh.read().splitlines()
            self.assertEqual(len(lines), 1)  # header only, no matched pairs

    def test_malformed_raw_rows_skipped(self):
        with tempfile.TemporaryDirectory() as td:
            summ = os.path.join(td, "P1.msi.tsv")
            raw = os.path.join(td, "P1_all")
            with open(summ, "w") as fh:
                fh.write("2\t1\t50.0\n")
            with open(raw, "w") as fh:
                # 8-column rows where column 8 (index 7) is the p-value;
                # short rows / empty p / non-float p must be skipped
                fh.write("chr\tstart\tscore\tdepth_t\tdepth_n\t"
                         "repeat2\trepeat\tP_value\n")
                fh.write("chr1\t1\t0\t40\t40\t-\tA.15A.16\n")        # 7 cols
                fh.write("chr1\t2\t0\t40\t40\t-\tA.15A.16\t\n")      # empty p
                fh.write("chr1\t3\t0\t40\t40\t-\tA.15A.16\tnope\n")  # non-float
                fh.write("chr1\t4\t0\t40\t40\t-\tA.15A.16\t0.01\n")
                fh.write("chr1\t5\t0\t40\t40\t-\tA.15A.16\t0.00001\n")
            out = os.path.join(td, "cohort.tsv")
            proc = subprocess.run(
                [sys.executable,
                 os.path.join(REPO, "scripts", "collect_msi.py"),
                 "--inputs", summ, raw, "--output", out],
                capture_output=True, text=True)
            assert proc.returncode == 0, proc.stderr
            import csv
            with open(out) as fh:
                r = list(csv.DictReader(fh, delimiter="\t"))[0]
            # only 2 rows parsed (the float ones with >=8 fields)
            self.assertEqual(r["raw_total_sites"], "2")
            self.assertEqual(r["bh_unstable"], "2")


if __name__ == "__main__":
    unittest.main(verbosity=2)
