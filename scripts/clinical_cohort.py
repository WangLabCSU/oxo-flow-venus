#!/usr/bin/env python3
"""Cohort collector for venus clinical interpretation.

Gathers per-pair biomarkers.tsv files into report/clinical_cohort.tsv with a
patient column, and writes a Tier I/II highlight table used by
scripts/clinical_report.py (--clinical). Missing pair files are skipped with a
stderr note (honest degradation, same read_tsv philosophy).
"""
import argparse
import os
import sys

import pandas as pd

TIER_ORDER = {"I": 0, "II": 1, "III": 2, "IV": 3, "": 4}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True, help="config/pairs.tsv")
    ap.add_argument("--dir", required=True, help="report/clinical directory")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t", dtype=str).fillna("")
    frames, skipped = [], []
    for pid in pairs["pair_id"]:
        path = os.path.join(args.dir, f"{pid}.biomarkers.tsv")
        if not os.path.exists(path):
            skipped.append(pid)
            continue
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        df.insert(0, "pair_id", pid)
        frames.append(df)
    if skipped:
        sys.stderr.write(f"[clinical_cohort] no biomarkers for: {', '.join(skipped)}\n")
    cols = ["pair_id", "class", "biomarker", "amp_tier", "detail",
            "evidence", "therapies_civic", "clinvar_sig", "caveat"]
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    if not out.empty:
        out["_o"] = out["amp_tier"].map(TIER_ORDER)
        out = out.sort_values(["_o", "pair_id", "class"]).drop(columns="_o")
    out.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[clinical_cohort] wrote {len(out)} rows over "
                     f"{len(frames)}/{len(pairs)} pairs -> {args.out}\n")


if __name__ == "__main__":
    main()
