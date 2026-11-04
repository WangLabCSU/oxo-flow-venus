#!/usr/bin/env python3
"""Concat per-sample fusion screening TSVs into the cohort table (N12).

Panel hits first, then by split-read support descending; per-sample and
per-class counts go to stderr so the report log carries the screening yield.
Columns are the union of every input's header (fusion_screen.py emits a fixed
schema; the union guards against schema drift between runs).
"""
from __future__ import annotations

import argparse
import csv
import sys

PRIORITY = {"panel_hit_single_read": 0, "panel_hit": 0}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="*", default=[],
                    help="per-sample rna/fusion/{sample}.fusion_candidates.tsv "
                         "(empty list -> header-only cohort table)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows: list[dict[str, str]] = []
    columns: list[str] = []
    seen: set[str] = set()
    for path in args.inputs:
        with open(path) as fh:
            for col in csv.DictReader(fh, delimiter="\t").fieldnames or []:
                if col not in seen:
                    seen.add(col)
                    columns.append(col)
        with open(path) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                rows.append(rec)

    def sort_key(rec):
        cls = rec.get("evidence_class", "")
        return (0 if cls in PRIORITY else 1,
                rec.get("panel_hit", "no") != "yes",
                -int(rec.get("split_reads") or 0),
                rec.get("fusion_pair", ""), rec.get("sample", ""))

    rows.sort(key=sort_key)
    if not columns:
        # empty cohort (no bulkRNA group): keep the table resolvable with the
        # canonical schema so downstream joins/report columns exist
        columns = ["sample"]
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    by_sample: dict[str, dict[str, int]] = {}
    for r in rows:
        s = r.get("sample", "")
        stat = by_sample.setdefault(s, {"total": 0, "panel": 0})
        stat["total"] += 1
        if r.get("panel_hit") == "yes":
            stat["panel"] += 1
    for s, st in by_sample.items():
        print(f"[collect_fusion] {s}: {st['total']} candidates, "
              f"{st['panel']} panel hits", file=sys.stderr)
    print(f"[collect_fusion] {len(rows)} rows over {len(by_sample)} samples "
          f"-> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
