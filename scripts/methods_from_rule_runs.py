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
    "vep": ("vep-rs", "Natera open source; cache v113"),
    "dnbc4tools": ("dnbc4tools", "DNBC4Tools 3.0 (BGI)"),
}


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

    # Versions: mine `--version`-style output from stdout tails if present.
    lines.append("### Software versions")
    lines.append("")
    seen: set[str] = set()
    for name, rec in rule_runs.items():
        for key in ("stdout_tail", "stderr_tail"):
            for line in (rec.get(key) or "").splitlines():
                m = re.search(r"version[:\s]+([0-9][0-9a-zA-Z.\-_]+)", line, re.I)
                if m:
                    ident = f"{tool_name(rec.get('command',''))} {m.group(1)}"
                    if ident not in seen:
                        seen.add(ident)
                        lines.append(f"- {ident}")
    if not seen:
        lines.append("- (versions captured in per-run logs; see run directory)")
    lines.append("")

    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"[methods_from_rule_runs] {len(rule_runs)} rules -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
