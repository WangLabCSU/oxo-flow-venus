#!/usr/bin/env python3
"""Join cohort frameshift-neoantigen counts with MSI status (venus#31 P1 item 7).

Builds the MSI-H x frameshift cohort page from two report-module tables:

  report/neoantigen_cohort.tsv - per-sample neoantigen burden pivot
                                 (collect_neoantigen.py; the
                                 frameshift_neoantigens column counts
                                 frameshift-derived candidates among
                                 strong binders, computed per pair by
                                 neoantigen_report.py)
  report/msi_cohort.tsv        - per-pair MSI calls (collect_msi.py;
                                 BH-corrected unstable-site fraction with
                                 msi_status_wgs at the 3.5% Niu et al. 2014
                                 WGS cutoff and msi_status_wes at the 15%
                                 WES/panel convention)

Every joined cell is either a real value or the explicit sentinel
"not_available" -- never a silent blank. The two sides rarely co-occur
today: MSIsensor-pro needs a matched normal (paired cohorts only) while
the per-pair neoantigen report is gated to tumor-only cohorts
(rules/neoantigen.oxoflow), so the page states which side carries data
instead of implying a join that never happened. A cohort that ran the
neoantigen module AND has matched normals populates both sides.

Usage:
  msi_frameshift_join.py \
      --neoantigen-cohort report/neoantigen_cohort.tsv \
      --msi-cohort report/msi_cohort.tsv \
      --strong-nm 500 \
      --out-tsv report/msi_frameshift_cohort.tsv \
      --out-md report/msi_frameshift_cohort.md
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

NOT_AVAILABLE = "not_available"

TSV_FIELDS = [
    "sample",
    "frameshift_neoantigens",
    "burden_strong_500nM",
    "frameshift_pct_of_strong",
    "msi_status_wgs",
    "msi_status_wes",
    "msi_bh_pct",
    "note",
]


def read_rows(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cell(row: dict | None, key: str) -> str | None:
    """Whitespace-stripped cell value; None when absent or empty."""
    if not row:
        return None
    value = (row.get(key) or "").strip()
    return value or None


def pct_of_strong(fs: str | None, burden: str | None) -> str | None:
    """frameshift share of strong binders in percent, 1 decimal."""
    try:
        f = float(fs)   # type: ignore[arg-type]
        s = float(burden)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if s <= 0:
        return None
    return f"{100.0 * f / s:.1f}"


def join(neo_rows: list[dict], msi_rows: list[dict]) -> list[dict]:
    neo = {(r.get("sample") or "").strip(): r for r in neo_rows}
    msi = {(r.get("sample") or "").strip(): r for r in msi_rows}
    neo.pop("", None)
    msi.pop("", None)

    out = []
    for sample in sorted(set(neo) | set(msi)):
        n, m = neo.get(sample), msi.get(sample)
        fs = cell(n, "frameshift_neoantigens")
        burden = cell(n, "burden_strong_500nM")
        fpct = pct_of_strong(fs, burden)
        swgs = cell(m, "msi_status_wgs")
        swes = cell(m, "msi_status_wes")
        bpct = cell(m, "bh_pct")
        notes = []
        if fs is None:
            notes.append("frameshift_not_available")
        if swgs is None:
            notes.append("msi_not_available")
        out.append({
            "sample": sample,
            "frameshift_neoantigens": fs or NOT_AVAILABLE,
            "burden_strong_500nM": burden or NOT_AVAILABLE,
            "frameshift_pct_of_strong": fpct or NOT_AVAILABLE,
            "msi_status_wgs": swgs or NOT_AVAILABLE,
            "msi_status_wes": swes or NOT_AVAILABLE,
            "msi_bh_pct": bpct or NOT_AVAILABLE,
            "note": "; ".join(notes),
        })
    # Complete rows first, MSI-H on top (the page exists for them), then
    # sample order.
    out.sort(key=lambda r: (
        0 if (r["frameshift_neoantigens"] != NOT_AVAILABLE
              and r["msi_status_wgs"] != NOT_AVAILABLE) else 1,
        0 if "MSI-H" in (r["msi_status_wgs"], r["msi_status_wes"]) else 1,
        r["sample"],
    ))
    return out


def write_tsv(path: str, rows: list[dict]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_FIELDS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_md(path: str, rows: list[dict], strong_nm: str,
             neo_header_only: bool, msi_header_only: bool) -> None:
    n_msi = sum(1 for r in rows if r["msi_status_wgs"] != NOT_AVAILABLE)
    n_fs = sum(1 for r in rows if r["frameshift_neoantigens"] != NOT_AVAILABLE)
    n_both = sum(1 for r in rows
                 if r["msi_status_wgs"] != NOT_AVAILABLE
                 and r["frameshift_neoantigens"] != NOT_AVAILABLE)
    lines = [
        "# MSI-H x frameshift neoantigen cohort page",
        "",
        "Generated by `scripts/msi_frameshift_join.py` (venus#31 P1 item 7)",
        "from `report/neoantigen_cohort.tsv` and `report/msi_cohort.tsv`.",
        "",
        "## Method",
        "",
        f"- `frameshift_neoantigens`: frameshift-derived candidates among"
        f" strong binders (<= {strong_nm} nM), computed per pair by"
        " `neoantigen_report.py` and pivoted cohort-wide by"
        " `collect_neoantigen.py`.",
        "- `msi_status_wgs` applies the 3.5% BH-corrected unstable-site"
        " cutoff (Niu et al. 2014, WGS); `msi_status_wes` the 15% WES/panel"
        " convention. Both are recomputed from raw per-site p-values by"
        " `collect_msi.py`, which also corrects the MSIsensor uint16 rank"
        " saturation (>65,535 sites).",
        f"Missing values are the explicit sentinel `{NOT_AVAILABLE}` --"
        " never a silent blank.",
        "",
        "## Data availability",
        "",
    ]
    if msi_header_only:
        lines += [
            "- MSI: no rows. MSIsensor-pro requires a matched normal for"
            " both calling and the baseline (a panel of normals cannot be"
            " built from a single tumor-only cohort), so a tumor-only"
            " cohort carries the header-only `report/msi_cohort.tsv`"
            " companion table and every MSI cell below is honestly"
            f" `{NOT_AVAILABLE}`.",
        ]
    elif n_msi == 0:
        lines += [
            "- MSI: `report/msi_cohort.tsv` has rows but none matched a"
            " sample in the neoantigen cohort table.",
        ]
    else:
        lines += [
            f"- MSI: {n_msi} sample(s) with MSI status",
            "  (BH-corrected, paired-normal MSIsensor-pro).",
        ]
    if neo_header_only:
        lines += [
            "- Frameshift neoantigens: the per-pair neoantigen report is"
            " gated to tumor-only cohorts (rules/neoantigen.oxoflow), so a"
            " paired cohort carries the header-only"
            " `report/neoantigen_cohort.tsv` companion and every"
            " frameshift cell below is honestly"
            f" `{NOT_AVAILABLE}`.",
        ]
    elif n_fs == 0:
        lines += [
            "- Frameshift neoantigens: no sample in"
            " `report/neoantigen_cohort.tsv` carries the metric.",
        ]
    else:
        lines += [
            f"- Frameshift neoantigens: {n_fs} sample(s) with the metric;",
            f" {n_both} sample(s) have both sides joined.",
        ]
    lines += [
        "",
        "## Cohort table",
        "",
        "| " + " | ".join(TSV_FIELDS) + " |",
        "|" + "|".join(["---"] * len(TSV_FIELDS)) + "|",
    ]
    for r in rows:
        lines.append("| " + " | ".join(r[f] for f in TSV_FIELDS) + " |")
    lines += ["", "## Interpretation caveats", "",
              "- A frameshift count is a lower bound: it covers only"
              " candidates that passed pVACseq binding prediction, so"
              " MSI-H samples with 0 frameshift neoantigens may simply"
              " lack expressed frameshift indels in predicted genes.",
              "- `msi_status_wgs` (3.5%) and `msi_status_wes` (15%) can"
              " disagree for cohorts whose BH-corrected unstable-site"
              " fraction falls between the two cutoffs; both are shown,"
              " never collapsed.",
              ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def main_with_args(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--neoantigen-cohort", required=True,
                    help="report/neoantigen_cohort.tsv (collect_neoantigen)")
    ap.add_argument("--msi-cohort", required=True,
                    help="report/msi_cohort.tsv (collect_msi)")
    ap.add_argument("--strong-nm", default="500",
                    help="strong-binder threshold used upstream (label only)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args(argv)

    neo_rows = read_rows(args.neoantigen_cohort)
    msi_rows = read_rows(args.msi_cohort)
    joined = join(neo_rows, msi_rows)
    write_tsv(args.out_tsv, joined)
    write_md(args.out_md, joined, str(args.strong_nm),
             neo_header_only=not neo_rows, msi_header_only=not msi_rows)
    print(f"[msi_frameshift_join] {len(joined)} sample(s) -> {args.out_tsv}",
          file=sys.stderr)
    return 0


def main() -> int:
    return main_with_args(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main())
