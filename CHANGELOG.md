# Changelog

All notable changes to Venus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Intermediate-file cleanup support: pure-scratch rules (fastp, bwa_mem2,
  mutect2 scatter/merge/filter chain) carry `temporary = true`; delivery,
  report, ASCAT, CNVkit, RNA-count, scRNA-analysis and QC-report rules declare
  `protected_output` so `oxo-flow clean` preserves them.
- `scripts/clean_intermediates.sh` — one-click cleanup wrapper: dry-run
  preview by default, `--mode all --apply` backs up `deliver/ report/ ascat/
  msi/` to a timestamped tarball then runs `oxo-flow clean --force`.
- `docs/cleanup.md` — what each mode keeps/deletes, the recompute-on-demand
  trade-off, the MSI caveat, and undeclared leftovers (`markdup.bam.bai`,
  manta workspaces, CNVkit helper `.cnn` files).
- Known engine gap documented: tombstone cleanup for wildcard-rule outputs is
  not yet supported, so `temporary` cleanup is dormant until the engine fix.

## [0.3.0] - 2026-10-05

### Changed
- Complete pipeline rebuild on oxo-flow ≥ 0.23: tumor multi-omics (WGS/WES somatic
  SNV/indel + CNV + bulk RNA + single-cell RNA + cohort report).
- Somatic chain: Mutect2 per-chromosome scatter → GatherVcfs → MergeMutectStats →
  LearnReadOrientationModel → FilterMutectCalls → PASS selection; both paired and
  tumor-only modes gated per pair by `wildcard.control`.
- Annotation: vep-rs (JSON cache) → annotated VCF + MAF (zero-Perl conversion).
- Report: clinical report whose methods section is generated from the commands
  that actually ran (rule_runs in `.oxo-flow/checkpoint.json`).
- Environments: pixi tomls under `envs/`, one per module.
- Removed the historical Rust binary + crates.io release workflow; CI now
  validates workflow TOML, helper scripts, and config tables.
