# Changelog

All notable changes to Venus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
