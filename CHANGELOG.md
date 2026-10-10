# Changelog

All notable changes to Venus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Intermediate-file cleanup support: pure-scratch rules (fastp, bwa_mem2,
  mutect2 scatter/merge chain) carry `temporary = true`; delivery, report,
  ASCAT, CNVkit, RNA-count, scRNA-analysis and QC-report rules declare
  `protected_output` so `oxo-flow clean` preserves them.
- Full filtered somatic call sets preserved: `filter_mutect` outputs
  (`vcf.filtered/{pair_id}.vcf.gz` + `.tbi`) and Manta scored calls
  (`manta/{pair_id}/results/variants/somaticSV.vcf.gz`/`candidateSV.vcf.gz`
  + `.tbi`) are protected — the PASS-only delivery VCFs discard the
  filtered-out records, which would otherwise be unrecoverable without full
  realignment.
- MultiQC aggregation inputs protected (`fastp.json`/`.fastp.html`,
  fastqc HTML reports) so the protected cohort QC report stays re-renderable
  after a clean.
- MSI products protected: `deliver/{pair_id}.msi.tsv` (per-pair clinical
  score) and `msi/reference.list` (cohort-shared msisensor scan, 1–3 h) —
  detection results survive every clean.
- Second protection wave (independent bioinformatics audit): Picard
  `dup_metrics.txt`, Mutect2 `.filteringStats.tsv`, and scRNA FastQC HTMLs are
  protected — all are measurement results and/or ingestion inputs of the
  protected MultiQC report (an unprotected input would silently vanish from a
  re-rendered report).
- `docs/artifacts.md` — complete output inventory: every declared output
  (1,553 files on the reference cohort) classified A/B/C by cleanup handling,
  per-module tables, the never-tracked (D) side-product list, and what each
  cleanup strategy (S1–S5) does to each class.
- `clean_intermediates.sh --apply` tarball extended to capture the
  evidence-grade undeclared side products before deletion: `msi/{pair_id}_all`
  (backs the MSI-cohort BH correction), Manta germline/candidate VCFs, scRNA
  per-cell results (`QC_Cluster.h5ad`, `singlecell.csv`, report HTMLs), and
  STAR `SJ.out.tab` junction tables.
- `scripts/clean_intermediates.sh` — one-click cleanup wrapper: dry-run
  preview by default, `--mode all --apply` backs up `deliver/ report/ ascat/
  msi/` to a timestamped tarball then runs `oxo-flow clean --force`.
- `docs/cleanup.md` — what each mode keeps/deletes, the recompute-on-demand
  trade-off, the MSI caveat, and undeclared leftovers (`markdup.bai`, manta
  workspaces, CNVkit BEDs, dnbc4tools side products).
- Known engine gap documented: tombstone cleanup for wildcard-rule outputs is
  not yet supported ([Traitome/oxo-flow#844](https://github.com/Traitome/oxo-flow/issues/844)),
  so `temporary` cleanup is dormant until the engine fix.
- Third protection wave (archive-grade intermediates): BQSR analysis-ready
  BAMs + indexes + recal tables (`bqsr/{sample}/{sample}.bqsr.bam`/`.bai`/
  `.recal.table`, ~1 TB cohort-wide), STAR aligned BAMs (~31 GB), scRNA
  matrix trios (`raw_matrix/` for SoupX/DecontX/scDblFinder and
  `filter_matrix/` — the 10x MEX Seurat ingests) + per-sample report HTMLs,
  and MSI per-pair locus detail (`msi/{pair_id}_all`, backs the
  `report/msi_cohort.tsv` BH correction) + `_unstable` summaries are
  protected. Any re-analysis after a clean now needs **no** realignment and
  no recount; the undeclared-path protections are pattern-based, so no
  completed instance goes stale. Census: 693 protected / 840 temporary /
  20 declared-deletable (only `markdup.bam`) of 1,553; a full mode-B clean
  reclaims ≈1.7 TB.

### Fixed
- `star_align` no longer hardcodes the read-group platform: the SAM `PL` tag
  now comes from `config.read_group_platform` (the same key the DNA-side
  `bwa_mem2` params use, default `"BGI"`), so DNA and RNA BAM headers always
  agree. Rendered commands are unchanged on the current config, but the
  engine fingerprints the shell *template*, so merging this stales the 10
  completed `star_align` instances — merge and sync the workdir only in a
  planned re-run window
  ([#2](https://github.com/WangLabCSU/oxo-flow-venus/issues/2)).
- Output census arithmetic: the declared-output total is **1,553** (protected
  **593**), matching the per-module tables; the headline previously
  undercounted by 10.
- S4 backup tarball now captures STAR junctions: the glob matched the wrong
  path (`*_STAR/SJ.out.tab`); STAR writes
  `rna/star/{sample}/{sample}.SJ.out.tab` because `--outFileNamePrefix`
  embeds the sample name. Docs corrected to the real path.
- `rules/msi.oxoflow` comment no longer suggests passing a BED to msisensor-pro:
  its `-b` flag sets the thread count (the rule binds it to `{threads}`), so
  `config.msi_targets_args` must stay empty; capture-restricted MSI is
  unsupported and the only WES/panel lever is `config.msi_coverage`.
- Inventory omissions filled: undeclared `vcf.raw/{pair_id}/{chr}.vcf.gz.tbi`
  indexes (×240, survive `clean`), ASCAT `{experiment}_normalBAF_rawBAF.txt`
  (30 germline-evidence files, not 20), the full msisensor trio
  (`{pair_id}_all` / `_dis` / `_unstable`), and STAR run logs + first-pass
  genome dirs.
- `clean_intermediates.sh --help` no longer prints the `set -euo pipefail`
  line (sed range trimmed).

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
