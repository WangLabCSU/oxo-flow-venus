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
    "AmpliconSuite-pipeline.py": (
        "AmpliconSuite",
        "Deshpande et al. 2022, Nat Biotechnol",
    ),
}

# Versions verified from pipeline-run artifacts (BAM @PG VN records, per-run
# logs, binary probes). Provenance column states the KIND of artifact each was
# confirmed from — kept cohort-agnostic so the table stays truthful in any
# cohort's report.
VERIFIED_VERSIONS = [
    ("fastp", "1.3.7", "fastp report (qc log)"),
    ("FastQC", "0.13.0", "qc log"),
    ("MultiQC", "1.35", "qc log"),
    ("BWA-MEM2", "2.2.1 (AVX2)", "BAM @PG VN records"),
    ("SAMtools", "1.24", "BAM @PG VN records"),
    ("Picard", "3.5.0", "BAM @PG record (MarkDuplicates); run log"),
    ("GATK", "4.6.2.0", "BAM @PG record (ApplyBQSR); run log jar banner"),
    ("Manta", "1.6.0", "run log (runWorkflow banner)"),
    ("MSIsensor", "1.3.0", "bioconda env probe"),
    ("CNVkit", "0.9.14", "cnv run log"),
    ("ASCAT", "3.2.0", "pixi env (bioconda)"),
    ("alleleCounter", "4.3.0", "pixi env (cancerit-allelecount)"),
    ("STAR", "2.7.11b", "BAM @PG VN record; STAR Log.out"),
    ("featureCounts", "2.1.1 (subread)", "featurecounts log banner"),
    ("vep-rs", "115.2", "vep run log (Ensembl VEP (Rust) banner)"),
    ("dnbc4tools", "3.0", "scrna run log (--version)"),
]

INTERPRETERS = {"python3", "python", "Rscript", "bash", "sh", "/bin/bash"}


def first_word(cmd: str) -> str:
    m = re.match(r"\s*(\S+)", cmd)
    return m.group(1) if m else ""


def unwrap_bash_c(cmd: str) -> str | None:
    """Inner script of a `pixi run --manifest-path X bash -c '<script>'` wrapper.

    Engine-recorded commands are pixi-wrapped; the analytical tool identity
    lives inside the -c payload, not in the wrapper tokens. The shell prelude
    (`set -euo pipefail`) is stripped so the first token is a real command.
    """
    m = re.search(r"\bbash\s+-c\s+'(.*)'", cmd, re.S)
    if not m:
        return None
    return re.sub(r"^\s*set\s+-[euo\w]+\s+pipefail\s*", "", m.group(1))


def _meta_scan(cmd: str) -> tuple[str, str] | None:
    base = Path(first_word(cmd)).name
    tokens = [Path(t).name for t in cmd.split()]
    for key, meta in TOOL_META.items():
        if base.startswith(key) or any(t == key or t.startswith(key) for t in tokens):
            return meta
    return None


def _scripts_tool(cmd: str) -> str:
    base = Path(first_word(cmd)).name
    scan = [t for t in cmd.split()[1:] if not t.startswith("-")]
    if base in INTERPRETERS:  # python3 script.py -> script tool identity
        for tok in scan:
            b = Path(tok).name
            if b.endswith((".py", ".R")):
                return b
    # Multi-statement inner scripts (post bash -c unwrap): the first command
    # may be mkdir/export housekeeping while the analytical identity is a
    # later script token — scan for it before falling back to the base.
    for tok in scan:
        b = Path(tok).name
        if b.endswith((".py", ".R")):
            return b
    return base


def tool_name(cmd: str) -> str:
    if first_word(cmd) == "pixi":
        inner = unwrap_bash_c(cmd)
        if inner:
            return tool_name(inner)
    meta = _meta_scan(cmd)
    if meta:
        return meta[0]
    return _scripts_tool(cmd)


def citation_for(cmd: str) -> str:
    if first_word(cmd) == "pixi":
        inner = unwrap_bash_c(cmd)
        if inner:
            return citation_for(inner)
    meta = _meta_scan(cmd)
    return meta[1] if meta else ""


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
        # No usable [[include]] list: fall back to scanning the rules dir
        # (legacy workdir layout where rule files are copied under rules/).
        candidates = sorted(rules_dir.glob("*.oxoflow"))
    else:
        # Mirror the engine's [[include]] semantics: absolute paths pass
        # through verbatim, relative paths resolve against the workflow
        # file's OWN directory. Rewriting every include to
        # `rules_dir/<basename>` silently skipped all includes in
        # include-based workdirs (rule files live in the reference repo,
        # not under rules/), yielding an empty manifest and an empty
        # methods section (#5).
        wf_dir = workflow.parent
        candidates = []
        for rel in include_paths:
            p = Path(rel)
            candidates.append(p if p.is_absolute() else wf_dir / p)
    for p in candidates:
        if not p.is_file():
            print(f"[methods] WARNING: include target not found: {p}",
                  file=sys.stderr)
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


def attribute_rule_runs(
    rule_runs: dict, manifest: "OrderedDict[str, list[dict]]"
) -> dict[str, dict[str, dict]]:
    """Map manifest rule name -> {checkpoint key: record}.

    Engine checkpoint keys are `<rule>_<instance>` for wildcard-instantiated
    rules (e.g. ``ampliconsuite_GC01``) and the bare rule name for
    singleton rules. Exact-name records were previously found, but every
    per-instance record was missed (``rule_runs.get(name)`` on key
    ``ampliconsuite_GC01``), so fresh per-sample runs were misreported as
    "rotated out of the checkpoint window". Attribution here walks manifest
    rule names LONGEST-FIRST so ``ascat_cohort_tumoronly`` is never
    mis-attributed to its prefix-sibling ``ascat_cohort``.
    """
    names = sorted(
        (r.get("name", "?") for rules in manifest.values() for r in rules),
        key=len, reverse=True,
    )
    attributed: dict[str, dict[str, dict]] = {}
    for key, rec in rule_runs.items():
        if not isinstance(rec, dict):
            continue
        for name in names:
            if key == name or key.startswith(name + "_"):
                attributed.setdefault(name, {})[key] = rec
                break
    return attributed


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
    if not manifest:
        # An empty manifest means the workflow definition was not resolved at
        # all — the section would be silently wrong (a methods section that
        # describes nothing, in a clinical deliverable). Fail loudly instead.
        print(
            "[methods_from_rule_runs] FATAL: resolved 0 rules from the "
            f"workflow manifest (workflow={args.workflow}, "
            f"rules-dir={args.rules_dir}); refusing to write an empty "
            "methods section. Check [[include]] paths and file existence.",
            file=sys.stderr,
        )
        return 2
    attributed = attribute_rule_runs(rule_runs, manifest)

    lines: list[str] = ["## Methods (generated from actual pipeline execution)", ""]
    n_rules = sum(len(v) for v in manifest.values())
    lines.append(
        f"The following describes the complete venus pipeline manifest — "
        f"{n_rules} rules across {len(manifest)} modules, in execution order "
        f"— as defined by the workflow and executed by the oxo-flow engine. "
        f"Verbatim commands are quoted from the run checkpoint for rules "
        f"with checkpoint records (per-instance rules show the instance "
        f"count and one representative command); rules rotated out of the "
        f"window are listed from the workflow definition only."
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
            recs = attributed.get(name) or {}
            if recs:
                keys = sorted(recs)
                rep = recs[keys[0]]
                cmd = rep.get("command", "")
                tool = tool_name(cmd) or "—"
                cit = citation_for(cmd)
                failed = [k for k in keys if recs[k].get("exit_code") != 0]
                if failed:
                    code = recs[failed[0]].get("exit_code")
                    status = (f"exit {code} "
                              f"({len(failed)}/{len(keys)} instance(s) failed)")
                else:
                    status = f"✓ success ({len(keys)} instance(s))"
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
                # No record: rotated out of the engine's rolling window, or
                # never instantiated in this cohort (paired-only rules in a
                # tumor-only run) — do NOT claim success for either case.
                status = ("no checkpoint record (rotated out of the window, "
                          "or not instantiated in this cohort)")
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
        if not isinstance(rec, dict):
            continue
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
