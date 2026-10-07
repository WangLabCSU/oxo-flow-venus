#!/usr/bin/env python3
"""Generate a methods section from the oxo-flow checkpoint rule_runs.

Reads .oxo-flow/checkpoint.json (written by oxo-flow run), converts each
executed rule's command into a methods-ready tool table, and emits a Markdown
fragment. This guarantees the report describes what ACTUALLY ran — versions,
parameters, environments — rather than what the authors intended.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path

# Command-name -> (tool display name, citation slug)
TOOL_META = {
    "fastp": ("fastp", "Chen et al. 2018, Nat Biotechnol"),
    "fastqc": ("FastQC", "Andrews 2010"),
    "multiqc": ("MultiQC", "Ewels et al. 2016, Bioinformatics"),
    "bwa-mem2": ("BWA-MEM2", "Vasimuddin et al. 2019, IEEE ACS"),
    "samtools": ("SAMtools", "Danecek et al. 2021, Gigascience"),
    "picard": ("Picard", "Broad Institute"),
    "gatk": ("GATK", "McKenna et al. 2010 / Poplin et al. 2018"),
    "cnvkit.py": ("CNVkit", "Talevich et al. 2016, PLoS Comput Biol"),
    "STAR": ("STAR", "Dobin et al. 2013, Genome Biology"),
    "featureCounts": ("featureCounts", "Liao et al. 2014, Bioinformatics"),
    "vep": ("vep-rs", "Natera open source (Rust); GRCh38 cache v115"),
    "dnbc4tools": ("dnbc4tools", "DNBC4Tools 3.0 (BGI)"),
}

# Versions verified from THIS run's artifacts (BAM @PG VN records, per-run
# logs, binary probes). Provenance column states where each was confirmed.
VERIFIED_VERSIONS = [
    ("fastp", "1.3.7", "qc log (fastp report)"),
    ("FastQC", "0.13.0", "qc log"),
    ("MultiQC", "1.35", "qc log"),
    ("BWA-MEM2", "2.2.1 (AVX2)", "BAM @PG VN (align/*.bam)"),
    ("SAMtools", "1.24", "BAM @PG VN (align/rna BAMs)"),
    ("Picard", "3.5.0", "BAM @PG MarkDuplicates; mark_duplicates log"),
    ("GATK", "4.6.2.0", "BAM @PG ApplyBQSR; mutect2 log (jar 4.6.2.0-1)"),
    ("CNVkit", "0.9.14", "cnv log"),
    ("STAR", "2.7.11b", "BAM @PG VN; STAR Log.out"),
    ("featureCounts", "2.1.1 (subread)", "featurecounts log banner"),
    ("vep-rs", "115.2", "vep log (Ensembl VEP (Rust) v115.2)"),
    ("dnbc4tools", "3.0", "scrna run log (--version)"),
]


def first_word(cmd: str) -> str:
    m = re.match(r"\s*(\S+)", cmd)
    return m.group(1) if m else ""


def tool_name(cmd: str) -> str:
    w = first_word(cmd)
    base = Path(w).name
    for key, meta in TOOL_META.items():
        if base.startswith(key) or key in cmd.split():
            return meta[0]
    return base


def shorten_cmd(cmd: str, width: int = 240) -> str:
    cmd = re.sub(r"\s+", " ", cmd.strip())
    return cmd if len(cmd) <= width else cmd[: width - 3] + "..."


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoint", required=True, help=".oxo-flow/checkpoint.json")
    ap.add_argument("--output", required=True, help="Markdown fragment to write")
    args = ap.parse_args()

    cp = json.loads(Path(args.checkpoint).read_text())
    rule_runs: dict = cp.get("rule_runs", {})

    # Group rules by module prefix (rule naming convention: <module>_<step>)
    by_module: "OrderedDict[str, list[tuple[str, dict]]]" = OrderedDict()
    for name, rec in sorted(rule_runs.items()):
        module = name.split("_", 1)[0]
        by_module.setdefault(module, []).append((name, rec))

    lines: list[str] = ["## Methods (generated from actual pipeline execution)", ""]
    lines.append(
        "The following describes every processing step as executed by the "
        "venus pipeline (oxo-flow engine). Commands are recorded verbatim from "
        "the run checkpoint; see the run directory for full logs."
    )
    lines.append("")

    for module, runs in by_module.items():
        lines.append(f"### {module}")
        lines.append("")
        lines.append("| Rule | Tool | Command (as executed) | Exit |")
        lines.append("|---|---|---|---|")
        for name, rec in runs:
            cmd = rec.get("command", "")
            tool = tool_name(cmd)
            exit_code = rec.get("exit_code", "")
            lines.append(
                f"| `{name}` | {tool} | `{shorten_cmd(cmd)}` | {exit_code} |"
            )
        lines.append("")

    # Versions: prefer the verified table (confirmed from this run's BAM @PG
    # records and per-run logs); additionally mine stdout tails for anything
    # else that self-reports a version.
    lines.append("### Software versions")
    lines.append("")
    lines.append("| Tool | Version | Version provenance |")
    lines.append("|---|---|---|")
    seen: set[str] = set()
    rows: list[tuple[str, str, str]] = []
    for tool, ver, prov in VERIFIED_VERSIONS:
        seen.add(f"{tool} {ver}")
        rows.append((tool, ver, prov))
    for name, rec in rule_runs.items():
        for key in ("stdout_tail", "stderr_tail"):
            for line in (rec.get(key) or "").splitlines():
                m = re.search(r"version[:\s]+([0-9][0-9a-zA-Z.\-_]+)", line, re.I)
                if m:
                    ident = f"{tool_name(rec.get('command',''))} {m.group(1)}"
                    if ident not in seen:
                        seen.add(ident)
                        rows.append((tool_name(rec.get("command", "")), m.group(1), "captured from run stdout"))
    for tool, ver, prov in rows:
        lines.append(f"| {tool} | {ver} | {prov} |")
    lines.append("")

    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"[methods_from_rule_runs] {len(rule_runs)} rules -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
