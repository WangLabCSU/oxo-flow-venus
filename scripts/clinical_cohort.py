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
    ap.add_argument("--out", required=True, help="cohort biomarkers TSV")
    ap.add_argument("--therapies-out", default="",
                    help="optional cohort therapies TSV")
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t", dtype=str).fillna("")
    frames, th_frames, skipped = [], [], []
    for pid in pairs["pair_id"]:
        path = os.path.join(args.dir, f"{pid}.biomarkers.tsv")
        if not os.path.exists(path):
            skipped.append(pid)
            continue
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        df.insert(0, "pair_id", pid)
        frames.append(df)
        th_path = os.path.join(args.dir, f"{pid}.therapies.tsv")
        if args.therapies_out and os.path.exists(th_path):
            th = pd.read_csv(th_path, sep="\t", dtype=str).fillna("")
            th.insert(0, "pair_id", pid)
            th_frames.append(th)
    if skipped:
        sys.stderr.write(f"[clinical_cohort] no biomarkers for: {', '.join(skipped)}\n")
    cols = ["pair_id", "class", "biomarker", "amp_tier", "detail",
            "evidence", "therapies_civic", "clinvar_sig", "caveat",
            "civic_diseases", "evidence_types"]
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    if not out.empty:
        out["_o"] = out["amp_tier"].map(TIER_ORDER)
        out = out.sort_values(["_o", "pair_id", "class"]).drop(columns="_o")
    out.to_csv(args.out, sep="\t", index=False)
    sys.stderr.write(f"[clinical_cohort] wrote {len(out)} rows over "
                     f"{len(frames)}/{len(pairs)} pairs -> {args.out}\n")
    if args.therapies_out:
        thcols = ["pair_id", "biomarker", "therapy", "target_or_context",
                  "region_or_approval", "line_or_level", "requires", "source",
                  "kb_tag"]
        th = pd.concat(th_frames, ignore_index=True) if th_frames \
            else pd.DataFrame(columns=thcols)
        th.to_csv(args.therapies_out, sep="\t", index=False)
        sys.stderr.write(f"[clinical_cohort] wrote {len(th)} therapy rows "
                         f"-> {args.therapies_out}\n")


if __name__ == "__main__":
    main()
