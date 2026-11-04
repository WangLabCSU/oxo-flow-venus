#!/usr/bin/env python3
"""Tests for scrna_metrics N13 disclosure (audit 20261104).

dnbc4tools 3.0 metrics_summary.xls carries no ambient-RNA / doublet fields:
the collector must emit explicit not_assessed markers (never blank, never
zero), and must promote REAL values if a future tool version starts reporting
them (KEY_ALIASES entries).

Run: python3 tests/test_scrna_metrics.py   (stdlib only)
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

import scrna_metrics as sm  # noqa: E402

# dnbc4tools 3.0-style headered metrics (no doublet/ambient anywhere)
DNBC3_XLS = (
    "Sample\tEstimated number of cell\tMean reads per cell\t"
    "Median UMI counts per cell\tSequencing saturation\tFraction reads in cell\n"
    "Pt01\t7,390\t45,000\t28,000\t83.76%\t97.25%\n")

# hypothetical future version that reports doublets/ambient
DNBC4_XLS = (
    "Sample\tEstimated number of cell\tDoublet rate\tAmbient rate\n"
    "Pt01\t7390\t1.2%\t0.8%\n")


def run_cli(metrics_payload, out):
    with tempfile.TemporaryDirectory() as d:
        m = os.path.join(d, "metrics_summary.xls")
        with open(m, "w") as fh:
            fh.write(metrics_payload)
        proc = subprocess.run(
            [sys.executable, os.path.join(REPO, "scripts", "scrna_metrics.py"),
             "--metrics", m, "--samples", "Pt01", "--output", out],
            capture_output=True, text=True)
        assert proc.returncode == 0, proc.stderr
        with open(out) as fh:
            return list(csv.DictReader(fh, delimiter="\t")), proc.stderr


class Disclosure(unittest.TestCase):
    def test_dnbc3_emits_not_assessed_markers(self):
        with tempfile.TemporaryDirectory() as d:
            rows, stderr = run_cli(
                DNBC3_XLS, os.path.join(d, "out.tsv"))
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual(r["doublet_pct"], "not_assessed(tool_lacks)")
        self.assertEqual(r["ambient_pct"], "not_assessed(tool_lacks)")
        # real fields unaffected (thousands separators cleaned)
        self.assertEqual(r["estimated_cells"], "7390")
        self.assertEqual(r["sequencing_saturation"], "83.76%")
        self.assertIn("not_assessed", stderr)

    def test_future_tool_values_promoted(self):
        with tempfile.TemporaryDirectory() as d:
            rows, _ = run_cli(DNBC4_XLS, os.path.join(d, "out.tsv"))
        r = rows[0]
        self.assertEqual(r["doublet_pct"], "1.2%")
        self.assertEqual(r["ambient_pct"], "0.8%")
        self.assertNotIn("not_assessed", r.values())

    def test_marker_module_constant(self):
        self.assertEqual(sm.NOT_ASSESSED_DEFAULTS["doublet_pct"],
                         "not_assessed(tool_lacks)")
        self.assertEqual(sm.NOT_ASSESSED_DEFAULTS["ambient_pct"],
                         "not_assessed(tool_lacks)")
        # alias wiring: a future doublet-rate key lands on the promoted column
        self.assertEqual(sm.KEY_ALIASES["doublet_rate"], "doublet_pct")
        self.assertEqual(sm.KEY_ALIASES["ambient_fraction"], "ambient_pct")


if __name__ == "__main__":
    unittest.main(verbosity=2)
