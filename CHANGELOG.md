# Changelog

All notable changes to Venus will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- ecDNA / focal-amplification module (issue #19) wrapping AmpliconSuite-
  pipeline 1.6 (bioconda `ampliconsuite`; CNVkit seed intervals ->
  AmpliconArchitect -> AmpliconClassifier). New rule `ampliconsuite`
  (`rules/ecdna.oxoflow`, env `envs/amplicon/pixi.toml`) takes the BQSR tumor
  BAM only — no matched-normal dependency, instantiates in BOTH paired and
  tumor-only cohorts. Outputs the five classification TSVs lifted to
  `ecdna/<pair>/` plus a `.suite_finished` sentinel (the suite's own finish
  flag, renamed — a partially-failed run never satisfies the outputs). New
  `[config]` keys: `aa_data_repo` (suite AA_DATA_REPO dir), `aa_solver`
  (default `clarabel` — no license needed; the rule pip-installs the mosek
  wheel only so the module-level import succeeds), `ecdna_args` (verbatim
  passthrough, e.g. `--normal_bam/--sv_vcf/--cnv_bed`). `scripts/
  verify_refdata.sh` checks the AA data repo (GRCh38/GRCh38_viral/mm10 +
  CNVkit ref + seed blacklist).
- Cohort ecDNA reporting: `ecdna_cohort` rule -> `report/cohort_ecdna.tsv` +
  `report/ecdna_genes.tsv` via `scripts/collect_ecdna.py` (strips the
  suite's `#sample` header prefix; joins amplicon classification from the
  result_table; 0-amplicon samples keep their counts row). `cohort_ecdna.tsv`
  joins into the patient-level summary (`ecdna_count` column) and
  `clinical_report.py` gains an `--ecdna` section (counts table + methods
  note).
- CNVkit gene annotation support, new `[config]` key `cnv_annotate_args`
  (default `""`): names bins with a gene model, e.g.
  `--annotate /path/gene.bed`, so the built reference `.cnn` and downstream
  `.cnn/.cnr/.cns` carry gene names instead of `-`. Valid in **paired mode
  only** (`cnvkit_paired`, where the reference is built): cnvkit ≥0.9 hard-
  rejects `--annotate` (and `-n/-t/-a`) alongside `--reference`, so
  `cnvkit_tumoronly` deliberately does not pass it — tumor-only gene naming
  comes from the reference `.cnn`'s own gene column.
- `scripts/annotate_cnv_reference.py`: re-annotate the gene column of an
  existing reference `.cnn` from a BED4 gene model without rebuilding the
  reference (numeric columns untouched; output verified bin-for-bin against
  `cnvkit.py target --annotate`, ERBB2 bins named on the chr17 test slice,
  end-to-end `batch --reference` run confirmed named `.cnr`/`.cns`).
- `docs/refdata.md` §7: two-level fix description (paired flag vs tumor-only
  reference re-annotation) with the cnvkit 0.9.14 flag-conflict evidence.

### Fixed
- Issue #17 (empty gene column in CNVkit outputs): root cause is the
  reference `.cnn` being built without `--annotate`; the gene column flows
  reference → target bed → coverage → fix → `.cnr/.cns` unchanged, so
  annotating the reference is both necessary and sufficient for tumor-only
  naming.

## [0.4.0] - 2026-10-21

### Added
- Tumor-only cohort gating for the report chain: new `[config]` flag
  `tumor_only` (string, default `"false"`; set with `--arg tumor_only=true`).
  `ascat_cohort` / `msi_cohort` carry `when = "config.tumor_only != 'true'"`,
  and companion rules `ascat_cohort_tumoronly` / `msi_cohort_tumoronly` share
  their output paths (`ascat/cohort_ascat.tsv`, `report/msi_cohort.tsv`) with
  the complementary gate, writing header-only tables (column names mirror
  `run_ascat.R` / `collect_msi.py` exactly). Without this, a tumor-only
  cohort could not finish: `ascat_paired` / `msi_paired` are never
  instantiated (per-pair `wildcard.control` gating), so the cohort rules'
  expanded inputs were unproducible and `cohort_tables` /
  `clinical_report` / `methods_from_rule_runs` all failed. Mixed cohorts
  (paired + tumor-only pairs in one run) remain unsupported for the ASCAT/MSI
  aggregations — documented in README.
- `docs/method-rationale.md`: why ASCAT/MSI are paired-only, and the
  header-only-companion gating pattern.
### Changed
- `mutect2_tumoronly` now consumes `config.pon` (`--panel-of-normals`), same
  as `mutect2_paired` — tumor-only calls previously skipped PoN filtering
  entirely. Default stays the bundled 1000G panel; override with a cohort
  PoN via `--arg pon=...`.
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
- Disk-pressure guard armed: `[engine]` section in `venus.oxoflow`
  (`min_free_disk = "2T"`, `reclaim_free_disk = "5T"`) — under disk pressure
  the engine tombstones completed temporary outputs mid-run, holds new
  dispatches below the floor, and aborts with a checkpoint as a last resort
  (oxo-flow PR [#845](https://github.com/Traitome/oxo-flow/pull/845); requires
  engine commit `a456c796` or later, not yet in a release tag).
- `temporary = true` tombstone cleanup is live: the engine gap that silently
  skipped wildcard-rule outputs ([Traitome/oxo-flow#844](https://github.com/Traitome/oxo-flow/issues/844))
  is fixed in the same engine PR, so the seven scratch rules (fastp, bwa_mem2,
  mutect2 scatter/merge chain) are cleaned automatically at run completion.
  `docs/cleanup.md` mode A status updated accordingly.
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
- `scripts/verify_refdata.sh` — host-side one-shot verification of every
  reference path declared in `venus.oxoflow` [config]: genome FASTA + BWA-MEM2
  index, STAR index + `genomeParameters.txt`, GENCODE GTF, GATK resources
  (+`.tbi`), vep-rs binary + cache `info.json`, dnbc4tools reference/binary,
  CNVkit reference (WARN-only, refdata.md §7), and the full ASCAT table
  contract (22 plain alleles + 22 loci, 18-field GC/RT rows, RT chr22
  sentinel) re-validated without re-downloading. Resolves and reports compat
  symlinks, so a broken link into the canonical root fails the check before a
  run starts.

### Changed
- Reference layout unified on the host (2026-10-11): all assets physically
  live under the canonical root `/home/data/reference/` — the dnbc4tools
  bundle (35 GB) moved from `/home/data/data-resource/` and the ASCAT set
  (1.8 GB) from `/home/data/wsx/hnzl_data/refdata/` via same-filesystem
  moves, with compat symlinks left at the old locations. Config strings are
  deliberately unchanged (rule commands embed them; changing them would
  stale every fingerprint touching those references), so
  `scrna_reference`/`scrna_bin`/`ascat_refdir` now resolve through symlinks —
  documented in `docs/refdata.md` and checked by `scripts/verify_refdata.sh`.

### Fixed
- MultiQC aggregation no longer silently drops QC sources (regression of the
  issue-#9 scan whitelist): the shell now also scans `align/` (Picard
  duplicate metrics, 20 WGS samples) and `scrna/fastqc/` (raw scRNA FastQC,
  10 samples), and `expand_inputs` mirrors every scanned root — Picard
  dup metrics over DNA, scRNA FastQC over scRNA, STAR `Log.final.out` +
  `ReadsPerGene.out.tab` and featureCounts `.summary` over bulkRNA — so the
  engine schedules multiqc only after each source exists instead of letting
  the protected report freeze without those sections
  ([#9](https://github.com/WangLabCSU/oxo-flow-venus/issues/9)).
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
