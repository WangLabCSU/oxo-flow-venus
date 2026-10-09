#!/usr/bin/env python3
"""Collect MSI calls into one cohort table, correcting the MSIsensor uint16 bug.

MSIsensor-msi (v1.3.0) stores each site's rank in a uint16 inside
calculateFDR; once a pair has >65535 ranked sites the rank wraps and the
reported "Number_of_Unstable_Sites" is a truncated artifact (e.g. exactly
65,535 = 2^16-1). The fix is to recompute a proper Benjamini-Hochberg FDR
from the raw per-site output (<pair>_all: P_value is column 8), which is
p-sorted, so unstable_bh = #{ q_i <= 0.05 } with
q_i = min_{j >= i} min(1, p_j * N / j).

Inputs arrive as ONE list mixing two file families (partitioned by name):
  deliver/<pair>.msi.tsv  - the msisensor summary as delivered
  msi/<pair>_all          - raw per-site table (header + p-sorted rows)

Output columns: sample, raw_total_sites, summary_total_sites,
reported_unstable, reported_pct, bh_unstable, bh_pct, unstable_delta,
uint16_saturated, msi_status_wgs (3.5% cutoff, Niu et al. 2014),
msi_status_wes (15% cutoff, WES/panel convention).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

WGS_CUTOFF_PCT = 3.5   # Niu et al. 2014, Sci Rep (WGS msisensor legacy)
WES_CUTOFF_PCT = 15.0  # WES/panel convention
SATURATION = 65535     # 2^16 - 1


def bh_unstable_count(p_values: list[float]) -> int:
    """Count sites with BH-adjusted p <= 0.05 (input assumed sorted asc)."""
    n = len(p_values)
    q = [0.0] * n
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, p_values[i] * n / (i + 1))
        q[i] = min(1.0, running)
    return sum(1 for v in q if v <= 0.05)


def read_summary(path: Path) -> tuple[str, str, str]:
    """Return (total, unstable, pct) from the 3-column msisensor summary."""
    lines = path.read_text().splitlines()
    if not lines:
        return "", "", ""
    rows = lines if not lines[0].startswith("Total_Number_of_Sites") else lines[1:]
    if not rows:
        return "", "", ""
    f = rows[0].split("\t")
    return f[0], f[1] if len(f) > 1 else "", f[2] if len(f) > 2 else ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="deliver/<pair>.msi.tsv and msi/<pair>_all files")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    summaries: dict[str, Path] = {}
    raws: dict[str, Path] = {}
    for p in args.inputs:
        path = Path(p)
        name = path.name
        if name.endswith(".msi.tsv"):
            summaries[name.removesuffix(".msi.tsv")] = path
        elif name.endswith("_all"):
            raws[name.removesuffix("_all")] = path
        else:
            print(f"[collect_msi] WARNING unrecognized input {p}", file=sys.stderr)

    pairs = sorted(set(summaries) & set(raws))
    missing = (set(summaries) ^ set(raws)) - set(pairs)
    for m in sorted(missing):
        print(f"[collect_msi] WARNING unmatched input for pair {m}", file=sys.stderr)

    fields = ["sample", "raw_total_sites", "summary_total_sites",
              "reported_unstable", "reported_pct", "bh_unstable", "bh_pct",
              "unstable_delta", "uint16_saturated", "msi_status_wgs",
              "msi_status_wes"]
    n_rows = 0
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for pair in pairs:
            st, su, sp = read_summary(summaries[pair])
            ps: list[float] = []
            with open(raws[pair]) as fh:
                header = fh.readline()  # skip header
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 8 or not f[7]:
                        continue
                    try:
                        ps.append(float(f[7]))
                    except ValueError:
                        continue
            ps.sort()
            n = len(ps)
            bh = bh_unstable_count(ps) if n else 0
            pct = 100.0 * bh / n if n else 0.0
            reported = int(su) if su.isdigit() else 0
            sat = reported == SATURATION and bh > SATURATION
            w.writerow({
                "sample": pair,
                "raw_total_sites": n,
                "summary_total_sites": st,
                "reported_unstable": su,
                "reported_pct": sp,
                "bh_unstable": bh,
                "bh_pct": f"{pct:.2f}",
                "unstable_delta": bh - reported,
                "uint16_saturated": str(sat).lower(),
                "msi_status_wgs": "MSI-H" if pct >= WGS_CUTOFF_PCT else "MSS",
                "msi_status_wes": "MSI-H" if pct >= WES_CUTOFF_PCT else "MSS",
            })
            n_rows += 1
    print(f"[collect_msi] {n_rows} pairs -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
