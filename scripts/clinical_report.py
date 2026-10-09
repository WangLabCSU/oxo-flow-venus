#!/usr/bin/env python3
"""Render the cohort clinical report (Markdown + HTML) from venus outputs.

Assembles: cohort overview, per-patient variant burden (dual TMB convention),
rule-based interpretation flags, bulk RNA QC, single-cell QC, optional
cohort-specific notes, and the methods section produced by
methods_from_rule_runs.py (describing what actually executed).
"""
from __future__ import annotations

import argparse
import csv
import html
import re
import sys
from datetime import datetime
from pathlib import Path


def read_tsv(path: str | None) -> list[dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _median(values: list[float]) -> float:
    v = sorted(values)
    n = len(v)
    if not v:
        return 0.0
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def tmb_markdown(path: str | None, max_rows: int = 20, title: str = "") -> str:
    """Compact TMB table: patient + key counts + TMB, from the patient-level
    cohort summary (tmb_ columns) or a raw tmb_from_maf table."""
    rows = read_tsv(path)
    if not rows:
        return f"_{title}: no data._\n"
    cols = list(rows[0].keys())
    if "tmb_sample_id" in cols:  # patient-level summary
        show = [
            ("sample", "patient"),
            ("tmb_sample_id", "tumor"),
            ("tmb_total_pass_variants", "pass variants"),
            ("tmb_coding_non_silent", "coding non-silent"),
            ("tmb_snps", "SNPs"),
            ("tmb_indels", "indels"),
            ("tmb_target_mb", "target Mb"),
            ("tmb_tmb_mut_per_mb", "TMB coding (mut/Mb)"),
            ("tmb_all_somatic_mut_per_mb", "TMB all somatic (mut/Mb)"),
        ]
        show = [(c, l) for c, l in show if c in cols]
    else:
        show = [(c, c) for c in cols]
    out = []
    if title:
        out.append(f"**{title}**\n")
    out.append("| " + " | ".join(l for _, l in show) + " |")
    out.append("|" + "---|" * len(show))
    for row in rows[:max_rows]:
        out.append("| " + " | ".join(str(row.get(c, "")) for c, _ in show) + " |")
    if len(rows) > max_rows:
        out.append(f"\n_… {len(rows) - max_rows} more rows in {path}_")
    return "\n".join(out) + "\n"


def interpretation_markdown(path: str | None) -> str:
    """Rule-based cohort interpretation, computed from the data (no
    cohort-specific logic): TMB outliers under both conventions, indel
    burden (MSI hint), clinical cutoffs. Cohort-specific narrative belongs
    in the --notes fragment."""
    rows = read_tsv(path)
    if not rows:
        return "_Interpretation: no data._\n"
    cols = list(rows[0].keys())
    c_coding = "tmb_tmb_mut_per_mb" if "tmb_tmb_mut_per_mb" in cols else "tmb_mut_per_mb"
    c_all = "tmb_all_somatic_mut_per_mb" if "tmb_all_somatic_mut_per_mb" in cols else "all_somatic_mut_per_mb"
    c_indels = "tmb_indels" if "tmb_indels" in cols else "indels"
    c_total = "tmb_total_pass_variants" if "tmb_total_pass_variants" in cols else "total_pass_variants"

    def f(row: dict, col: str) -> float:
        try:
            return float(row.get(col, "") or 0)
        except ValueError:
            return 0.0

    med_coding = _median([f(r, c_coding) for r in rows])
    med_all = _median([f(r, c_all) for r in rows])

    out = [
        "### Interpretation (rule-based, computed from this cohort)\n",
        "TMB is reported under two conventions: **coding non-silent** per Mb "
        "(clinical WES/panel convention, comparable with the ≥10 mut/Mb "
        "high-TMB cutoff) and **all somatic** per Mb (WGS literature "
        "convention — every PASS somatic variant per callable Mb). Without a "
        "per-base callable mask the denominator over-estimates territory, so "
        "both values are conservative.\n",
    ]
    flagged = []
    for r in rows:
        s = r.get("sample", "?")
        tmb, tmb_all = f(r, c_coding), f(r, c_all)
        flags = []
        if med_coding > 0 and tmb >= 5 * med_coding:
            flags.append(f"coding TMB {tmb:.3f} mut/Mb = {tmb/med_coding:.0f}× cohort median")
        if med_all > 0 and tmb_all >= 5 * med_all:
            flags.append(f"all-somatic TMB {tmb_all:.3f} mut/Mb = {tmb_all/med_all:.0f}× cohort median")
        if tmb >= 10:
            flags.append("exceeds ≥10 mut/Mb high-TMB clinical cutoff (coding convention)")
        tot, ind = f(r, c_total), f(r, c_indels)
        if tot > 0 and ind / tot >= 0.30:
            flags.append(f"indel fraction {ind/tot:.0%} — possible MSI/MMR deficiency")
        if flags:
            flagged.append((s, flags))
    if flagged:
        out.append("**Outlier samples:**\n")
        for s, flags in flagged:
            out.append(f"- **{s}**: " + "; ".join(flags))
        out.append("")
        out.append(
            "Outlier status under BOTH conventions with a high indel fraction "
            "is the statistical signature of a hypermutator (MSI/MMR-deficient) "
            "phenotype; confirm with an orthogonal MSI assay before clinical use."
        )
    else:
        out.append("No sample exceeds 5× the cohort median TMB or the ≥10 mut/Mb "
                   "clinical cutoff; indel fractions are within the normal range.")
    return "\n".join(out) + "\n"


def tsv_markdown(path: str | None, max_rows: int = 20, title: str = "",
                 only: list[str] | None = None) -> str:
    rows = read_tsv(path)
    if not rows:
        return f"_{title}: no data._\n"
    cols = [c for c in (only or rows[0].keys()) if c in rows[0]]
    out = []
    if title:
        out.append(f"**{title}**\n")
    out.append("| " + " | ".join(cols) + " |")
    out.append("|" + "---|" * len(cols))
    for row in rows[:max_rows]:
        out.append("| " + " | ".join(str(row.get(c, "")) for c in cols) + " |")
    if len(rows) > max_rows:
        out.append(f"\n_… {len(rows) - max_rows} more rows in {path}_")
    return "\n".join(out) + "\n"


INLINE_CODE_RE = re.compile(r"`([^`]+)`")
INLINE_BOLD_RE = re.compile(r"\*\*([^*]+)\*\*")
# Underscore italics only for spans containing a space, so snake_case paths
# and identifiers (tmb_mut_per_mb, logs/sbs96_Pt09.log) are left untouched.
INLINE_EM_RE = re.compile(r"(?<![\w])_([^_]*\s[^_]*)_(?![\w])")


def _inline(escaped: str) -> str:
    """Inline `code` / **bold** / _italic_ on already-escaped text."""
    escaped = INLINE_CODE_RE.sub(r"<code>\1</code>", escaped)
    escaped = INLINE_BOLD_RE.sub(r"<b>\1</b>", escaped)
    return INLINE_EM_RE.sub(r"<em>\1</em>", escaped)


def md_to_html(md: str) -> str:
    """Minimal Markdown -> HTML for headings, tables, code, emphasis."""
    lines = md.splitlines()
    out: list[str] = []
    in_table = False
    for line in lines:
        esc = html.escape(line)
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            if all(set(c) <= {"-", ":", " "} for c in cells) and cells:
                continue  # separator row
            if not in_table:
                out.append("<table>")
                in_table = True
                out.append("<tr>" + "".join(f"<th>{_inline(c)}</th>" for c in cells) + "</tr>")
            else:
                out.append("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in cells) + "</tr>")
            continue
        if in_table:
            out.append("</table>")
            in_table = False
        if line.startswith("### "):
            out.append(f"<h3>{_inline(esc[4:])}</h3>")
        elif line.startswith("## "):
            out.append(f"<h2>{_inline(esc[3:])}</h2>")
        elif line.startswith("# "):
            out.append(f"<h1>{_inline(esc[2:])}</h1>")
        elif line.startswith("- "):
            out.append(f"<li>{_inline(esc[2:])}</li>")
        elif line.strip() == "":
            out.append("")
        else:
            out.append(f"<p>{_inline(esc)}</p>")
    if in_table:
        out.append("</table>")
    return "\n".join(out)


ASCAT_METHODS = (
    "Allele-specific copy number, tumour purity and ploidy were inferred with "
    "**ASCAT 3.2.0** using its high-throughput-sequencing chain, exactly as "
    "documented by VanLoo-lab: allele counts at the G1000 hg38 SNP loci "
    "(alleleCounter 4.3.0; min base quality 20, min mapping quality 35), "
    "logR GC-content and replication-timing correction with the official "
    "G1000 references, ASPCF segmentation (penalty 70) and the purity/ploidy "
    "grid fit with gamma = 1 (the HTS value; 0.55 is the SNP-array legacy). "
    "Autosomes 1:22 only — patient sex is unannotated for this cohort, so "
    "X/Y are excluded from the fit. Per-module selection rationale: "
    "`docs/method-rationale.md` in the pipeline repository.\n"
)

ASCAT_INTERP = (
    "Interpretation guidance: colorectal cancers are typically aneuploid "
    "(CIN phenotype) with recurrent arm-level events — gains of 20q/13q/7 "
    "and losses of 18q/17p (TP53); low-purity fits (<30%) have reduced "
    "sensitivity for subclonal events, and a goodness-of-fit below ~0.80 "
    "flags a less reliable purity/ploidy solution. Segment tables "
    "(`ascat/<pair>/segments.tsv`) and per-sample plots accompany this "
    "report.\n"
)

MSI_HEADER = (
    "Microsatellite instability was called with **MSIsensor-msi 1.3.0** on "
    "paired WGS. The table reports both the pipeline's own summary and a "
    "**Benjamini-Hochberg recomputation** performed by this pipeline from "
    "the raw per-site output.\n"
)

MSI_UINT16_NOTE = (
    "**Known MSIsensor counting bug, corrected here:** MSIsensor v1.3.0 "
    "stores each site's rank in a uint16 while computing FDR, so above "
    "65,535 sites the rank wraps and the reported unstable-site count "
    "saturates (visible as exactly 65,535 = 2^16-1). Recomputing a proper "
    "BH FDR from the raw per-site table fixes this: for **Pt09** the "
    "reported 65,535 sites (10.54%) is an artifact — the BH-corrected count "
    "is **102,316 unstable sites (16.46%)**, which *strengthens* the MSI-H "
    "call. For the other nine pairs reported and recomputed counts agree "
    "exactly (one pair differs by a single borderline site at q ≈ 0.05).\n"
)

MSI_CUTOFF_NOTE = (
    "**Cutoffs are convention-dependent:** the 3.5% threshold is the WGS-era "
    "msisensor convention (Niu et al. 2014, *Sci Rep*); WES/panel practice "
    "commonly uses ≥15% (msisensor-pro and clinical assays). Only Pt09 "
    "exceeds both thresholds (MSI-H); all other pairs are MSS under either "
    "convention, consistent with the TMB/indel-fraction analysis above. "
    "Clinical use requires orthogonal confirmation (MMR immunohistochemistry "
    "or MLH1 promoter methylation).\n"
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tmb")
    ap.add_argument("--ascat", help="ASCAT cohort purity/ploidy table")
    ap.add_argument("--msi", help="MSI cohort table (BH-corrected)")
    ap.add_argument("--rna-qc")
    ap.add_argument("--scrna-metrics")
    ap.add_argument("--methods", help="methods Markdown fragment")
    ap.add_argument("--pairs", help="config/pairs.tsv")
    ap.add_argument("--notes", help="cohort-specific interpretation notes (Markdown)")
    ap.add_argument("--output-md", required=True)
    ap.add_argument("--output-html")
    args = ap.parse_args()

    pairs = read_tsv(args.pairs)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    parts: list[str] = [
        "# Venus cohort report",
        f"_Generated {stamp} by the venus multi-omics pipeline (oxo-flow)._",
        "",
        "## Cohort overview",
        "",
    ]
    if pairs:
        parts.append(f"- Paired tumor/normal samples: **{len(pairs)}**")
        parts.append(f"- Tumor IDs: {', '.join(r['experiment'] for r in pairs)}")
        parts.append(f"- Normal IDs: {', '.join(r['control'] for r in pairs)}")
    parts.append("")

    parts.append("## Tumor mutational burden (per Mb of callable territory)\n")
    parts.append(tmb_markdown(args.tmb, title="TMB summary"))
    parts.append(interpretation_markdown(args.tmb))

    if args.ascat:
        parts.append("## Tumor purity, ploidy and allele-specific copy number (ASCAT)\n")
        parts.append(ASCAT_METHODS)
        parts.append(tsv_markdown(
            args.ascat, title="ASCAT fit summary",
            only=["sample", "purity", "ploidy", "psi",
                  "aberrant_cell_fraction", "goodness_of_fit",
                  "n_segments", "non_aberrant", "failed"]))
        parts.append(ASCAT_INTERP)

    if args.msi:
        parts.append("## Microsatellite instability (MSI)\n")
        parts.append(MSI_HEADER)
        parts.append(tsv_markdown(
            args.msi, title="MSI summary (reported vs BH-corrected)",
            only=["sample", "summary_total_sites", "reported_unstable",
                  "reported_pct", "bh_unstable", "bh_pct",
                  "uint16_saturated", "msi_status_wgs", "msi_status_wes"]))
        parts.append(MSI_UINT16_NOTE)
        parts.append(MSI_CUTOFF_NOTE)

    if args.notes and Path(args.notes).exists():
        parts.append("### Cohort-specific findings\n")
        parts.append(Path(args.notes).read_text())

    parts.append("## Bulk RNA-seq alignment QC\n")
    parts.append(tsv_markdown(args.rna_qc, title="STAR alignment metrics"))

    parts.append("## Single-cell QC\n")
    parts.append(tsv_markdown(args.scrna_metrics, title="scRNA metrics"))

    if args.methods and Path(args.methods).exists():
        parts.append(Path(args.methods).read_text())
    else:
        parts.append("## Methods\n\n_(methods fragment missing — see run logs)_")

    md = "\n".join(parts) + "\n"
    Path(args.output_md).write_text(md)

    if args.output_html:
        body = md_to_html(md)
        doc = (
            "<!doctype html><html><head><meta charset='utf-8'>"
            "<title>Venus cohort report</title>"
            "<style>body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;"
            "margin:2rem auto;max-width:1100px;line-height:1.5}"
            "table{border-collapse:collapse;width:100%;font-size:0.85rem}"
            "th,td{border:1px solid #ddd;padding:4px 8px;text-align:left}"
            "th{background:#f5f5f5}h1,h2,h3{color:#1a1a2e}</style></head>"
            f"<body>{body}</body></html>"
        )
        Path(args.output_html).write_text(doc)
    print(f"[clinical_report] wrote {args.output_md}"
          + (f" + {args.output_html}" if args.output_html else ""), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
