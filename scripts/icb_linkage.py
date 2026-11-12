#!/usr/bin/env python3
"""ICB-relevant linkage page: TMB x MSI x neoantigen burden (venus#31 P0 gap).

Joins the three independently-computed immunotherapy biomarker tables into
one patient-level cross-evidence page:

  report/cohort_summary.tsv   - patient-level summary (cohort_tables.py);
                                the TMB side reads the tmb_tmb_mut_per_mb
                                (coding, mut/Mb) and tmb_all_somatic_mut_per_mb
                                columns produced from per-pair tmb_from_maf
                                tables
  report/msi_cohort.tsv       - per-pair MSI calls (collect_msi.py;
                                BH-corrected unstable-site fraction with
                                msi_status_wgs at the 3.5% Niu et al. 2014
                                WGS cutoff and msi_status_wes at the 15%
                                WES/panel convention; header-only for
                                tumor-only cohorts)
  report/neoantigen_cohort.tsv - per-sample neoantigen burden pivot
                                (collect_neoantigen.py; header-only for
                                paired cohorts, see rules/report.oxoflow)

Status flags reuse the exact conventions already rendered by
clinical_report.py -- TMB-H at the FDA label threshold (>= 13 mut/Mb
coding) with 10-<13 kept as legacy high-TMB context, never upgraded to a
label -- and both MSI conventions side by side, never collapsed.

Sample keys are normalized to a patient id (trailing -TD/-B/-TR by
default, --suffixes) so the patient-level TMB table joins the pair-level
MSI / neoantigen tables. Every joined cell is either a real value or the
explicit sentinel "not_available" -- never a silent blank. The page is
descriptive cross-evidence for review; it is not an ICB response
prediction.

Usage:
  icb_linkage.py \
      --cohort-summary report/cohort_summary.tsv \
      --msi-cohort report/msi_cohort.tsv \
      --neoantigen-cohort report/neoantigen_cohort.tsv \
      --out-tsv report/icb_linkage_cohort.tsv \
      --out-md report/icb_linkage_cohort.md
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

NOT_AVAILABLE = "not_available"

TSV_FIELDS = [
    "sample",
    "tmb_coding_mut_per_mb",
    "tmb_all_somatic_mut_per_mb",
    "tmb_status",
    "msi_status_wgs",
    "msi_status_wes",
    "msi_bh_pct",
    "neo_burden_strong_500nM",
    "neo_burden_highaff_100nM",
    "neo_frameshift_neoantigens",
    "note",
]

TMB_COL_CODING = "tmb_tmb_mut_per_mb"
TMB_COL_ALL = "tmb_all_somatic_mut_per_mb"
# Same wording split as clinical_report.py interpretation_markdown: the FDA
# label threshold vs the legacy 10-<13 context that must not read as a label.
TMB_STATUS_H = "TMB-H (FDA >=13 mut/Mb coding)"
TMB_STATUS_LEGACY = "TMB 10-<13 (legacy high-TMB context, not FDA label)"
TMB_STATUS_NOT_HIGH = "not_high"


def read_rows(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cell(row: dict | None, key: str) -> str | None:
    """Whitespace-stripped cell value; None when absent or empty."""
    if not row:
        return None
    value = (row.get(key) or "").strip()
    return value or None


def tmb_status_of(tmb: str | None) -> str | None:
    """FDA-label vs legacy wording, identical to clinical_report.py."""
    if tmb is None:
        return None
    try:
        value = float(tmb)
    except ValueError:
        return None
    if value >= 13:
        return TMB_STATUS_H
    if value >= 10:
        return TMB_STATUS_LEGACY
    return TMB_STATUS_NOT_HIGH


def msi_flagged(status: str | None) -> bool:
    return status is not None and "MSI-H" in status


def join(cohort_rows: list[dict], msi_rows: list[dict],
         neo_rows: list[dict], suffixes: str) -> list[dict]:
    suffix_re = re.compile(r"-(" + "|".join(suffixes.split(",")) + r")$")

    def patient_of(sample: str) -> str:
        return suffix_re.sub("", sample)

    # First-wins on duplicate patients everywhere, matching
    # cohort_tables.py patients.setdefault: a repeated id is a data bug,
    # and keeping the first row is deterministic, unlike last-wins.
    cohort: dict[str, dict] = {}
    for row in cohort_rows:
        sample = (row.get("tmb_sample_id") or "").strip()
        if sample:
            cohort.setdefault(patient_of(sample), row)
    msi: dict[str, dict] = {}
    for row in msi_rows:
        sample = (row.get("sample") or "").strip()
        if sample:
            msi.setdefault(patient_of(sample), row)
    neo: dict[str, dict] = {}
    for row in neo_rows:
        sample = (row.get("sample") or "").strip()
        if sample:
            neo.setdefault(patient_of(sample), row)

    out = []
    for patient in sorted(set(cohort) | set(msi) | set(neo)):
        c, m, n = cohort.get(patient), msi.get(patient), neo.get(patient)
        tmb = cell(c, TMB_COL_CODING)
        tmb_all = cell(c, TMB_COL_ALL)
        status = tmb_status_of(tmb)
        swgs = cell(m, "msi_status_wgs")
        swes = cell(m, "msi_status_wes")
        bpct = cell(m, "bh_pct")
        burden = cell(n, "burden_strong_500nM")
        burden_ha = cell(n, "burden_highaff_100nM")
        fs = cell(n, "frameshift_neoantigens")
        notes = []
        if c is None:
            notes.append("tmb_not_available")
        elif tmb is None:
            notes.append("tmb_column_missing")
        if m is None:
            notes.append("msi_not_available")
        if n is None:
            notes.append("neoantigen_not_available")
        out.append({
            "sample": patient,
            "tmb_coding_mut_per_mb": tmb or NOT_AVAILABLE,
            "tmb_all_somatic_mut_per_mb": tmb_all or NOT_AVAILABLE,
            "tmb_status": status or NOT_AVAILABLE,
            "msi_status_wgs": swgs or NOT_AVAILABLE,
            "msi_status_wes": swes or NOT_AVAILABLE,
            "msi_bh_pct": bpct or NOT_AVAILABLE,
            "neo_burden_strong_500nM": burden or NOT_AVAILABLE,
            "neo_burden_highaff_100nM": burden_ha or NOT_AVAILABLE,
            "neo_frameshift_neoantigens": fs or NOT_AVAILABLE,
            "note": "; ".join(notes),
        })
    # ICB-flagged rows first (TMB-H or MSI-H under either convention),
    # then rows carrying more available sides, then sample order.
    def sides(r: dict) -> int:
        return sum(1 for f in ("tmb_coding_mut_per_mb", "msi_status_wgs",
                               "neo_burden_strong_500nM")
                   if r[f] != NOT_AVAILABLE)

    out.sort(key=lambda r: (
        0 if (r["tmb_status"] == TMB_STATUS_H
              or msi_flagged(r["msi_status_wgs"])
              or msi_flagged(r["msi_status_wes"])) else 1,
        -sides(r),
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


def write_md(path: str, rows: list[dict], suffixes: str,
             cohort_has_tmb_col: bool, msi_header_only: bool,
             neo_header_only: bool) -> None:
    n_tmb = sum(1 for r in rows
                if r["tmb_coding_mut_per_mb"] != NOT_AVAILABLE)
    n_msi = sum(1 for r in rows if r["msi_status_wgs"] != NOT_AVAILABLE)
    n_neo = sum(1 for r in rows
                if r["neo_burden_strong_500nM"] != NOT_AVAILABLE)
    n_tmbh = sum(1 for r in rows if r["tmb_status"] == TMB_STATUS_H)
    n_msih = sum(1 for r in rows
                 if msi_flagged(r["msi_status_wgs"])
                 or msi_flagged(r["msi_status_wes"]))
    lines = [
        "# ICB-relevant linkage: TMB x MSI x neoantigen burden",
        "",
        "Generated by `scripts/icb_linkage.py` (venus#31) from"
        " `report/cohort_summary.tsv`, `report/msi_cohort.tsv` and"
        " `report/neoantigen_cohort.tsv`.",
        "",
        "## Method",
        "",
        f"- TMB side: `tmb_tmb_mut_per_mb` (coding mut/Mb) and"
        f" `tmb_all_somatic_mut_per_mb` from the patient-level cohort"
        f" summary; TMB-H flag applies the FDA label threshold"
        f" (>= 13 mut/Mb coding), 10-<13 is legacy high-TMB context and"
        f" never shown as a label (same wording as the clinical report).",
        "- MSI side: `msi_status_wgs` applies the 3.5% BH-corrected"
        " unstable-site cutoff (Niu et al. 2014, WGS); `msi_status_wes`"
        " the 15% WES/panel convention. Both are shown, never collapsed.",
        "- Neoantigen side: pVACseq candidates at the strong (500 nM)"
        " and high-affinity (100 nM) binding thresholds, pivoted"
        " cohort-wide by `collect_neoantigen.py`.",
        f"- Sample keys are normalized to patient level (trailing"
        f" -{suffixes.replace(',', '/-')} stripped) so the patient-level"
        " TMB table joins the pair-level MSI / neoantigen tables.",
        f"Missing values are the explicit sentinel `{NOT_AVAILABLE}` --"
        " never a silent blank.",
        "",
        "## Data availability",
        "",
    ]
    if not cohort_has_tmb_col:
        lines += [
            "- TMB: `report/cohort_summary.tsv` carries no"
            f" `{TMB_COL_CODING}` column; every TMB cell below is honestly"
            f" `{NOT_AVAILABLE}`.",
        ]
    elif n_tmb == 0:
        lines += [
            "- TMB: no cohort-summary row matched a sample in the MSI or"
            " neoantigen tables.",
        ]
    else:
        lines += [f"- TMB: {n_tmb} patient(s) with coding TMB;"
                  f" {n_tmbh} flagged TMB-H (FDA >=13)."]
    if msi_header_only:
        lines += [
            "- MSI: header-only `report/msi_cohort.tsv` (MSIsensor-pro"
            " needs a matched normal; tumor-only cohort), so every MSI"
            f" cell below is honestly `{NOT_AVAILABLE}`.",
        ]
    elif n_msi == 0:
        lines += [
            "- MSI: `report/msi_cohort.tsv` has rows but none matched a"
            " patient in the other tables.",
        ]
    else:
        lines += [
            f"- MSI: {n_msi} sample(s) with status"
            "  (BH-corrected, paired-normal MSIsensor-pro);"
            f" {n_msih} flagged MSI-H under either convention.",
        ]
    if neo_header_only:
        lines += [
            "- Neoantigen burden: header-only"
            " `report/neoantigen_cohort.tsv` (the neoantigen module is"
            " gated to tumor-only cohorts, venus#31), so every neoantigen"
            f" cell below is honestly `{NOT_AVAILABLE}`.",
        ]
    elif n_neo == 0:
        lines += [
            "- Neoantigen burden: no sample in"
            " `report/neoantigen_cohort.tsv` matched a patient here.",
        ]
    else:
        lines += [f"- Neoantigen burden: {n_neo} patient(s) with the"
                  " strong-binder burden metric."]
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
              "- This page is descriptive cross-evidence for manual"
              " review, not an ICB response prediction: TMB, MSI and"
              " neoantigen burden are correlated but partially independent"
              " axes, and response depends on additional factors (HLA"
              " presentation breadth, tumor microenvironment, prior"
              " therapy).",
              "- The TMB-H threshold applies to coding mut/Mb over the"
              " cohort's callable territory (panel or WGS); values are"
              " not directly comparable across capture designs.",
              "- `msi_status_wgs` (3.5%) and `msi_status_wes` (15%) can"
              " disagree for cohorts whose BH-corrected unstable-site"
              " fraction falls between the two cutoffs; both are shown,"
              " never collapsed.",
              "- A high neoantigen burden with no ICB flag is expected in"
              " MSS/low-TMB samples; burden alone does not imply"
              " immunogenicity (binding prediction only, expression not"
              " required for the burden columns).",
              ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def main_with_args(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-summary", required=True,
                    help="report/cohort_summary.tsv (cohort_tables)")
    ap.add_argument("--msi-cohort", required=True,
                    help="report/msi_cohort.tsv (collect_msi or the"
                         " header-only tumor-only companion)")
    ap.add_argument("--neoantigen-cohort", required=True,
                    help="report/neoantigen_cohort.tsv (collect_neoantigen"
                         " or the header-only paired companion)")
    ap.add_argument("--suffixes", default="TD,B,TR",
                    help="comma-separated assay suffixes stripped from"
                         " sample ids to get the patient key"
                         " (default: TD,B,TR; mirrors cohort_tables)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args(argv)

    cohort_rows = read_rows(args.cohort_summary)
    msi_rows = read_rows(args.msi_cohort)
    neo_rows = read_rows(args.neoantigen_cohort)
    # Column presence must come from the header line, not the rows: a
    # header-only table yields no DictReader rows at all.
    with open(args.cohort_summary, newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"), [])
    cohort_has_tmb_col = TMB_COL_CODING in header
    joined = join(cohort_rows, msi_rows, neo_rows, args.suffixes)
    write_tsv(args.out_tsv, joined)
    write_md(args.out_md, joined, str(args.suffixes),
             cohort_has_tmb_col=cohort_has_tmb_col,
             msi_header_only=not msi_rows,
             neo_header_only=not neo_rows)
    print(f"[icb_linkage] {len(joined)} patient(s) -> {args.out_tsv}",
          file=sys.stderr)
    return 0


def main() -> int:
    return main_with_args(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main())
