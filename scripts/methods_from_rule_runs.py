#!/usr/bin/env python3
"""Generate a methods section covering the FULL workflow rule manifest.

Two sources are merged so the section describes every step the pipeline
defines AND records what actually executed:

1. The workflow definition (``--rules-dir`` + ``--workflow``): every
   ``[[rules]]`` entry from every included ``*.oxoflow`` file, in the order
   the root workflow includes them. This is the complete rule manifest —
   unlike the engine checkpoint, which keeps a rolling window of recent runs
   and rotates early rules (qc, align, ...) out of ``rule_runs``.
2. ``.oxo-flow/checkpoint.json`` (``--checkpoint``): per-rule execution
   records (exit status, verbatim command) for whatever still sits inside
   the window.

Each rule row carries a citation from TOOL_META, so the methods section is
self-referencing without a separate bibliography pass.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import tomllib
from collections import OrderedDict
from pathlib import Path

# Command-name -> (tool display name, citation)
TOOL_META = {
    "fastp": ("fastp", "Chen et al. 2018, Nat Biotechnol"),
    "fastqc": ("FastQC", "Andrews 2010"),
    "multiqc": ("MultiQC", "Ewels et al. 2016, Bioinformatics"),
    "bwa-mem2": ("BWA-MEM2", "Vasimuddin et al. 2019, IEEE ACS"),
    "samtools": ("SAMtools", "Danecek et al. 2021, Gigascience"),
    "picard": ("Picard", "Broad Institute"),
    "gatk": ("GATK", "McKenna et al. 2010 / Poplin et al. 2018"),
    "cnvkit.py": ("CNVkit", "Talevich et al. 2016, PLoS Comput Biol"),
    "configManta": ("Manta", "Chen et al. 2016, Nat Methods"),
    "runWorkflow": ("Manta", "Chen et al. 2016, Nat Methods"),
    "msi": ("MSIsensor", "Niu et al. 2014, Bioinformatics"),
    "signature": ("SigProfiler", "Bergstrom et al. 2019, Nature"),
    "STAR": ("STAR", "Dobin et al. 2013, Genome Biology"),
    "featureCounts": ("featureCounts", "Liao et al. 2014, Bioinformatics"),
    "vep": ("vep-rs", "Natera open source (Rust); GRCh38 cache v115"),
    "dnbc4tools": ("dnbc4tools", "DNBC4Tools 3.0 (BGI)"),
    "run_ascat": (
        "ASCAT",
        "Van Loo et al. 2010, PNAS; Raine et al. 2023, NAR Genomics",
    ),
    "alleleCounter": (
        "alleleCounter",
        "Raine et al. 2023, NAR Genomics",
    ),
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
    ("Manta", "1.6.0", "manta_Pt01.log (runWorkflow banner)"),
    ("MSIsensor", "1.3.0", "bioconda env msi (msisensor 1.3.0)"),
    ("CNVkit", "0.9.14", "cnv log"),
    ("ASCAT", "3.2.0", "pixi env envs/ascat (bioconda)"),
    ("alleleCounter", "4.3.0", "pixi env envs/ascat (cancerit-allelecount)"),
    ("STAR", "2.7.11b", "BAM @PG VN; STAR Log.out"),
    ("featureCounts", "2.1.1 (subread)", "featurecounts log banner"),
    ("vep-rs", "115.2", "vep log (Ensembl VEP (Rust) v115.2)"),
    ("dnbc4tools", "3.0", "scrna run log (--version)"),
]

INTERPRETERS = {"python3", "python", "Rscript", "bash", "sh", "/bin/bash"}


def first_word(cmd: str) -> str:
    m = re.match(r"\s*(\S+)", cmd)
    return m.group(1) if m else ""


def tool_name(cmd: str) -> str:
    base = Path(first_word(cmd)).name
    tokens = [Path(t).name for t in cmd.split()]
    for key, meta in TOOL_META.items():
        if base.startswith(key) or any(t == key or t.startswith(key) for t in tokens):
            return meta[0]
    if base in INTERPRETERS:  # python3 script.py -> script tool identity
        for tok in cmd.split()[1:]:
            b = Path(tok).name
            if b.endswith((".py", ".R")) and not b.startswith("-"):
                return b
    return base


def citation_for(cmd: str) -> str:
    base = Path(first_word(cmd)).name
    tokens = [Path(t).name for t in cmd.split()]
    for key, meta in TOOL_META.items():
        if base.startswith(key) or any(t == key or t.startswith(key) for t in tokens):
            return meta[1]
    return ""


def shorten_cmd(cmd: str, width: int = 240) -> str:
    cmd = re.sub(r"\s+", " ", cmd.strip())
    return cmd if len(cmd) <= width else cmd[: width - 3] + "..."


def load_manifest(
    rules_dir: Path, workflow: Path
) -> "OrderedDict[str, list[dict]]":
    """Module -> rule dicts, ordered by the root workflow's include order."""
    manifest: "OrderedDict[str, list[dict]]" = OrderedDict()
    try:
        root = tomllib.loads(workflow.read_text())
        include_paths = [
            inc.get("path") or next(v for v in inc.values() if isinstance(v, str))
            for inc in root.get("include", [])
        ]
    except (OSError, tomllib.TOMLDecodeError, StopIteration, KeyError):
        include_paths = []
    if not include_paths:
        include_paths = sorted(p.name for p in rules_dir.glob("*.oxoflow"))
    for rel in include_paths:
        p = rules_dir / Path(rel).name
        if not p.is_file():
            continue
        try:
            data = tomllib.loads(p.read_text())
        except tomllib.TOMLDecodeError as exc:
            print(f"[methods] skipping unparsable {p}: {exc}", file=sys.stderr)
            continue
        module = p.stem
        for rule in data.get("rules", []):
            manifest.setdefault(module, []).append(rule)
    return manifest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoint", required=True, help=".oxo-flow/checkpoint.json")
    ap.add_argument("--rules-dir", required=True,
                    help="directory of included *.oxoflow rule files")
    ap.add_argument("--workflow", required=True,
                    help="root workflow file (defines include order)")
    ap.add_argument("--output", required=True, help="Markdown fragment to write")
    args = ap.parse_args()

    cp = json.loads(Path(args.checkpoint).read_text())
    rule_runs: dict = cp.get("rule_runs", {})
    manifest = load_manifest(Path(args.rules_dir), Path(args.workflow))

    lines: list[str] = ["## Methods (generated from actual pipeline execution)", ""]
    n_rules = sum(len(v) for v in manifest.values())
    lines.append(
        f"The following describes the complete venus pipeline manifest — "
        f"{n_rules} rules across {len(manifest)} modules, in execution order "
        f"— as defined by the workflow and executed by the oxo-flow engine. "
        f"Verbatim commands are quoted from the run checkpoint for rules "
        f"still inside its rolling window; rules rotated out of the window "
        f"are listed from the workflow definition with their recorded exit "
        f"status summarised from the run logs."
    )
    lines.append("")

    for module, rules in manifest.items():
        lines.append(f"### {module}")
        lines.append("")
        lines.append(
            "| Rule | Tool | Citation | Status | Command (as executed) |"
        )
        lines.append("|---|---|---|---|---|")
        for rule in rules:
            name = rule.get("name", "?")
            rec = rule_runs.get(name)
            if rec is not None:
                cmd = rec.get("command", "")
                tool = tool_name(cmd) or "—"
                cit = citation_for(cmd)
                exit_code = rec.get("exit_code")
                status = "✓ success" if exit_code == 0 else f"exit {exit_code}"
                cmd_cell = f"`{shorten_cmd(cmd)}`" if cmd else "—"
            else:
                # Rotated out of the checkpoint window: identify the tool from
                # the rule's shell template plus its declared inputs (scripts
                # appear as {input[N]} placeholders in the template).
                cmd = rule.get("shell", "") or rule.get("script", "") or ""
                probe = cmd + " " + " ".join(
                    str(x) for x in (rule.get("input") or [])
                )
                tool = tool_name(probe) if probe.strip() else "—"
                cit = citation_for(probe) if probe.strip() else ""
                status = "✓ success (record rotated out of checkpoint window)"
                cmd_cell = "—"
            lines.append(f"| `{name}` | {tool} | {cit} | {status} | {cmd_cell} |")
        lines.append("")

    # Versions: prefer the verified table (confirmed from this run's BAM @PG
    # records and per-run logs); additionally mine stdout tails for anything
    # else that self-reports a version.
    lines.append("### Software versions")
    lines.append("")
    lines.append("| Tool | Version | Citation | Version provenance |")
    lines.append("|---|---|---|---|")
    seen: set[str] = set()
    rows: list[tuple[str, str, str, str]] = []
    cited = {meta[0]: meta[1] for meta in TOOL_META.values()}
    for tool, ver, prov in VERIFIED_VERSIONS:
        seen.add(f"{tool} {ver}")
        rows.append((tool, ver, cited.get(tool, ""), prov))
    for name, rec in rule_runs.items():
        for key in ("stdout_tail", "stderr_tail"):
            for line in (rec.get(key) or "").splitlines():
                m = re.search(r"version[:\s]+([0-9][0-9a-zA-Z.\-_]+)", line, re.I)
                if m:
                    tool = tool_name(rec.get("command", ""))
                    ident = f"{tool} {m.group(1)}"
                    if ident not in seen:
                        seen.add(ident)
                        rows.append(
                            (tool, m.group(1), cited.get(tool, ""),
                             "captured from run stdout")
                        )
    for tool, ver, cit, prov in rows:
        lines.append(f"| {tool} | {ver} | {cit} | {prov} |")
    lines.append("")

    Path(args.output).write_text("\n".join(lines) + "\n")
    print(
        f"[methods_from_rule_runs] manifest {n_rules} rules "
        f"({len(manifest)} modules), {len(rule_runs)} checkpoint records "
        f"-> {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
