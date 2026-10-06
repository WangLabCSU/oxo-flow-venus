#!/usr/bin/env python3
"""Render the cohort clinical report (Markdown + HTML) from venus outputs.

Assembles: cohort overview, per-patient variant burden, CNV summaries, bulk
RNA QC, single-cell QC, and the methods section produced by
methods_from_rule_runs.py (describing what actually executed).
"""
from __future__ import annotations

import argparse
import csv
import html
import sys
from datetime import datetime
from pathlib import Path


def read_tsv(path: str | None) -> list[dict[str, str]]:
    if not path or not Path(path).exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))




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
            ("tmb_tmb_mut_per_mb", "TMB (mut/Mb)"),
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


def tsv_markdown(path: str | None, max_rows: int = 20, title: str = "") -> str:
    rows = read_tsv(path)
    if not rows:
        return f"_{title}: no data._\n"
    cols = list(rows[0].keys())
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
                out.append("<tr>" + "".join(f"<th>{c}</th>" for c in cells) + "</tr>")
            else:
                out.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
            continue
        if in_table:
            out.append("</table>")
            in_table = False
        if line.startswith("### "):
            out.append(f"<h3>{esc[4:]}</h3>")
        elif line.startswith("## "):
            out.append(f"<h2>{esc[3:]}</h2>")
        elif line.startswith("# "):
            out.append(f"<h1>{esc[2:]}</h1>")
        elif line.startswith("- "):
            out.append(f"<li>{esc[2:]}</li>")
        elif line.strip() == "":
            out.append("")
        else:
            out.append(f"<p>{esc}</p>")
    if in_table:
        out.append("</table>")
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tmb")
    ap.add_argument("--rna-qc")
    ap.add_argument("--scrna-metrics")
    ap.add_argument("--methods", help="methods Markdown fragment")
    ap.add_argument("--pairs", help="config/pairs.tsv")
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
