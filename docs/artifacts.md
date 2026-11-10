# Output artifact inventory and cleanup classification

Every file the venus pipeline produces in a run workdir, classified by how the
cleanup tooling treats it. The inventory is generated from `rules/*.oxoflow`
and verified against the live cohort workdir (`processed-wsx`, 10 CRC pairs,
2026-10 ground truth). See [cleanup.md](cleanup.md) for the strategy rationale
and [../scripts/clean_intermediates.sh](../scripts/clean_intermediates.sh) for
the operator wrapper.

## Reference cohort

| Dimension | Value |
|---|---|
| DNA (WGS tumour/normal) | 10 pairs = 20 samples |
| Bulk RNA | 10 samples |
| Single-cell RNA | 10 samples |
| Chromosome scatter | 24 units per pair (Mutect2) |
| Completed rule instances | 561 (`completed_rules` in `.oxo-flow/checkpoint.json`) |
| Declared output files | 1,553 (sum of output patterns × instances) |
| Paired vs tumor-only | all 10 pairs matched → every tumor-only rule = 0 instances |

Reproduce the counts:

```bash
python3 - <<'EOF'
import tomllib, glob
for f in sorted(glob.glob('rules/*.oxoflow')):
    d = tomllib.load(open(f, 'rb'))
    for r in d['rules']:
        po, tmp = set(r.get('protected_output', [])), 'temporary' in r
        for o in r.get('output', []):
            cls = 'A' if o in po else ('B' if tmp else 'C')
            print(f"{f} {r['name']:22} [{cls}] {o}")
EOF
```

## Classification scheme

| Class | Meaning | Engine mechanism | Files (this cohort) |
|---|---|---|---|
| **A — protected** | Detection/measurement results, delivery artifacts, expensive report-facing products, archive-grade BAMs/matrices. Survive every clean. | `protected_output = [...]` | **693** (45%) |
| **B — temporary scratch** | Fully regenerable intermediates, only feed downstream rules. Slated for automatic post-success deletion. | `temporary = true` | **840** (54%) |
| **C — declared, deletable** | Intermediates without either marker. Deleted by `clean`, survive normal operation. | declared `output` only | **20** (1%) — only `align/{sample}/{sample}.markdup.bam` |
| **D — undeclared leftovers** | Side products the engine never tracks. `clean` cannot see them. | not declared | logs + per-tool extras |

## What each strategy does to each class

| Strategy | A protected | B temporary | C deletable | D undeclared |
|---|---|---|---|---|
| **S1** Auto-tombstone after successful run (Mode A) | keep | delete *(live on engine ≥ `a456c796`, PR [#845](https://github.com/Traitome/oxo-flow/pull/845))* | keep | keep |
| **S2** `oxo-flow clean` (dry-run) | nothing deleted, preview only | | | |
| **S3** `oxo-flow clean --force` (Mode B) | keep | delete | delete | keep |
| **S4** `scripts/clean_intermediates.sh --mode all --apply` | keep (tar backup of `deliver/ report/ ascat/ msi/` first) | delete | delete | keep (tar additionally captures the remaining D-class evidence: Manta variant VCFs, scRNA per-cell results, STAR junctions, MSI `_dis` distributions) |
| **S5** manual `rm` of a declared file | destroys the result or marks its producer **stale** → next `run` recomputes it and every ancestor | | | no bookkeeping impact, but no safety net |

Notes:

- **Mode A requires the cleanup-fixed engine**: oxo-flow commit `a456c796`
  (PR [#845](https://github.com/Traitome/oxo-flow/pull/845)) fixed the tombstone
  loop's wildcard skip ([#844](https://github.com/Traitome/oxo-flow/issues/844));
  on an older engine the seven B-class venus scratch rules are silently kept,
  and only S3/S4 reclaim scratch. Never work around this with manual `rm` —
  deleting declared outputs without checkpoint bookkeeping makes producers
  stale and would trigger recomputation on the next `oxo-flow run`.
- **S5 on an A-class file is a real loss**: regenerating it re-runs the rule
  *and all ancestors* (e.g. any `deliver/` VCF → full BWA-MEM2 realignment,
  days of compute). Treat the A column as load-bearing.
- After any clean, `oxo-flow plan` reports cleaned rules as stale — expected.
  The pre-clean invariant is `would run: 0 | skip: 561 | completed: 561`.
- `logs/**`, `.oxo-flow/` (checkpoint + env manifests) and `backups/` are
  never touched by any engine strategy.

## Inventory by module

Classes: **A** protected · **B** temporary · **C** deletable · **D** undeclared.
"×n" = instances on the reference cohort; a 24-chr pair produces 24 files per
chr-scattered pattern.

### QC (`rules/qc.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `fastp` ×30 (20 DNA + 10 bulkRNA) | `qc/trim/{sample}_R{1,2}.trim.fastq.gz` | B | ~1 TB cohort-wide; regenerates from `raw/` |
| `fastp` | `qc/trim/{sample}.fastp.json` / `.html` | A | inputs to protected MultiQC |
| `fastqc` ×30 | `qc/fastqc/{sample}_R{1,2}.trim_fastqc.html` | A | inputs to protected MultiQC |
| `multiqc` ×1 | `qc/multiqc/multiqc_report.html`, `.../multiqc_general_stats.txt` | A | cohort QC report |
| `wgs_flagstat` ×20 | `qc/alignqc/{sample}.flagstat` | **A** | mapping rate / duplicates (20261104 N10); MultiQC input |
| `wgs_bamstats` ×20 | `qc/alignqc/{sample}.stats` | **A** | insert-size distribution / GC (N10); MultiQC input |
| `wgs_mosdepth` ×20 | `qc/alignqc/{sample}.mosdepth.summary.txt` + `.mosdepth.thresholds.bed.gz` | **A** | depth + ≥10x/20x fractions + sex-check input (N10/N11); 1Mb windows — no whole-genome per-base file |
| `wgs_mosdepth` | `qc/alignqc/{sample}.mosdepth.regions.bed.gz` (+ side products `*.mosdepth.global.dist.txt`, `*.mosdepth.region.dist.txt`) | C | per-window means + depth histograms; regenerable in minutes |
| `wgs_alignqc_summary` ×1 | `qc/alignqc/wgs_alignqc_summary.tsv` | **A** | one row per DNA sample; feeds MultiQC anchors, cohort_summary, clinical report |

### Alignment (`rules/align.oxoflow`) — DNA

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `bwa_mem2` ×20 | `align/{sample}/{sample}.sorted.bam` | B | hundreds of GB; the biggest single scratch class |
| `mark_duplicates` ×20 | `align/{sample}/{sample}.markdup.bam` | C | regenerable from sorted.bam; **the only declared deletable class left** — its information is fully subsumed by the bqsr BAM |
| `mark_duplicates` | `align/{sample}/{sample}.dup_metrics.txt` | **A** | library QC measurement; feeds the protected MultiQC report |
| `bqsr` ×20 | `bqsr/{sample}/{sample}.recal.table` / `.bqsr.bam` / `.bqsr.bai` | **A** | analysis-ready BAMs every caller consumes (Mutect2/Manta/ASCAT/CNVkit/MSIsensor) — TCGA-style archive layer; ~1 TB cohort-wide, protected so re-analysis never means realignment |

### Somatic SNV/indel (`rules/varcall.oxoflow`, `rules/varcall_merge.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `mutect2_paired` ×240 (10 pairs × 24 chr) | `vcf.raw/{pair_id}/{chr}.vcf.gz` + `.vcf.gz.stats` + `{chr}.tar.gz` | B | 720 files; tens of GB |
| `mutect2_tumoronly` | (same patterns) | — | defined, 0 instances (all pairs matched) |
| `gather_chr_vcfs` ×10 | `vcf.raw/{pair_id}.vcf.gz` + `.tbi` | B | |
| `merge_mutect_stats` ×10 | `vcf.raw/{pair_id}.stats` | B | |
| `learn_orientation` ×10 | `vcf.filtered/{pair_id}.orientation-bias.tar.gz` | B | |
| `filter_mutect` ×10 | `vcf.filtered/{pair_id}.vcf.gz` + `.tbi` | **A** | only copy of the filtered-out records (germline/strand/PoN); PASS delivery discards them |
| `filter_mutect` | `vcf.filtered/{pair_id}.vcf.gz.filteringStats.tsv` | **A** | audit record of the filtering step (per-filter rejection counts) |
| `select_pass` ×10 | `deliver/{pair_id}.somatic.pass.vcf.gz` + `.tbi` | **A** | primary somatic delivery |

### Annotation (`rules/annotation.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `vep_annotate` ×10 | `deliver/{pair_id}.somatic.pass.vep.vcf.gz` + `.tbi` | **A** | annotated clinical call set |
| `vcf2maf` ×10 | `deliver/{pair_id}.maf` | **A** | MAF for downstream tools |

### CNV — ASCAT (`rules/ascat.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `ascat_paired` ×10 | `ascat/{pair_id}/{pair_id}.ASCAT.rds`, `purity_ploidy.tsv`, `segments.tsv`, `metrics.tsv`, `{experiment}_tumour{LogR,BAF}.txt`, `plots/{pair_id}.{experiment}.{tumour,germline,ASPCF}.png`, `.ASCATprofile.pdf` (10 files/pair) | **A** | most expensive CNV step; report-facing |
| `ascat_cohort` ×1 | `ascat/cohort_ascat.tsv` | **A** | |

### CNV — CNVkit (`rules/cnv.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `cnvkit_paired` ×10 | `cnv/{pair_id}/reference.cnn`, `{experiment}.bqsr.cnr` / `.call.cns` / `-scatter.png`, `{pair_id}.cnr` / `.call.cns` / `.cnv.png` | **A** | `{pair_id}` files symlink the `{experiment}` files; both sides protected |
| `cnvkit_tumoronly` | (same patterns) | — | defined, 0 instances |
| — | `cnv/{pair_id}/GRCh38*.bed` | D | target/antitarget BEDs copied in as inputs |

### SV — Manta (`rules/sv.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `manta_paired` ×10 | `manta/{pair_id}/results/variants/somaticSV.vcf.gz` + `.tbi` | **A** | full scored somatic set behind the PASS-only delivery |
| `manta_tumoronly` | `.../candidateSV.vcf.gz` + `.tbi` | — (would be A) | defined, 0 instances |
| `sv_deliver_paired` ×10 | `deliver/{pair_id}.sv.vcf.gz` + `.tbi` | **A** | |
| `sv_deliver_tumoronly` | `deliver/{pair_id}.sv.candidate.vcf.gz` + `.tbi` | — (would be A) | defined, 0 instances |
| — | `manta/{pair_id}/workspace/`, `results/{evidence,stats}/`, run metadata, undeclared `results/variants/{candidateSV,candidateSmallIndels,diploidSV}.vcf.gz(.tbi)` | D | `diploidSV` = germline calls, worth archiving before manual rm |

### MSI (`rules/msi.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `msisensor_scan` ×1 | `msi/reference.list` | **A** | cohort-shared microsatellite scan, 1–3 h on hg38 |
| `msi_paired` ×10 | `deliver/{pair_id}.msi.tsv` | **A** | per-pair clinical MSI score (detection result) |
| `msi_paired` | `msi/{pair_id}_all` + `msi/{pair_id}_unstable` (undeclared but protected) | **A** | `_all` is evidence-grade per-locus detail — the protected `report/msi_cohort.tsv` BH-corrects over its counts; `_unstable` (114 B) is the instability summary. Protection is pattern-based, valid without an output-list change |
| — | `msi/{pair_id}_dis` | D | per-locus distributions, ~1.1 GB/pair (~11 GB cohort-wide) — the only MSI leftover worth reclaiming |

### Signatures (`rules/signature.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `sbs96_fit` ×10 | `report/signatures/{pair_id}.sbs96.counts.tsv` + `.exposures.tsv` | **A** | |

### Bulk RNA (`rules/rna.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `star_align` ×10 | `rna/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam` | **A** | archive-grade RNA BAM (~31 GB cohort-wide): junction discovery / IGV / re-counting without re-alignment |
| `star_align` | `rna/star/{sample}/{sample}.Log.final.out` + `.ReadsPerGene.out.tab` + `.Chimeric.out.junction` | **A** | QC products; the junction file feeds `rna_fusion_screen` (20261104 N12) |
| `featurecounts` ×10 | `rna/counts/{sample}.featurecounts.txt` + `.summary` | **A** | strand mode via `config.featurecounts_strandedness` (N14, default 0 = unstranded) |
| `merge_counts` ×1 | `rna/counts/cohort.gene_counts.txt` | **A** | cohort matrix |
| `rna_qc_summary` ×1 | `rna/qc/rna_qc_summary.tsv` | **A** | |
| `rna_fusion_screen` ×10 | `rna/fusion/{sample}.fusion_candidates.tsv` | **A** | screening-level split-read candidates (N12); every row carries the orthogonal-confirmation disclosure |
| `rna_fusion_cohort` ×1 | `report/fusion_screen_cohort.tsv` | **A** | panel hits first; clinical report section |

### Single-cell RNA (`rules/scrna.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `scrna_fastqc` ×10 | `scrna/fastqc/{sample}_SC_R{1,2}_fastqc.html` | **A** | QC measurement; feeds the protected MultiQC report |
| `scrna_count` ×10 | `scrna/count/{sample}/outs/raw_matrix/{matrix.mtx,features.tsv,barcodes.tsv}.gz` | **A** | unfiltered (all-barcodes) count matrix — the layer ambient-RNA/doublet tools (SoupX, DecontX, scDblFinder) ingest; **re-running scrna_count is the most expensive recompute in the pipeline**, so the matrix trio is protected |
| `scrna_count` | `outs/filter_matrix/{matrix.mtx,features.tsv,barcodes.tsv}.gz` (undeclared but protected) + `outs/filter_feature.h5ad`, `outs/analysis/{cluster,marker}.csv`, `outs/metrics_summary.xls`, `outs/{sample}_scRNA_report.html` (undeclared but protected) | **A** | filter_matrix/ is the 10x MEX trio Seurat's `Read10X` ingests directly; h5ad/CSVs/xls are detection results; the 2.8 MB report HTML is the dnbc4tools QC deliverable. The protection is pattern-based, so these survive every clean even though `clean` itself only deletes declared outputs |
| `scrna_qc_cluster` ×10 | `scrna/qc/{sample}.clusters.tsv` | **A** | |
| `scrna_integrate` ×1 | `scrna/cohort_clusters.tsv` | **A** | |
| — | `outs/anno_decon_sorted.bam` + `.bai`, `outs/analysis/QC_Cluster.h5ad` (~1 GB/sample), `outs/singlecell.csv` (~98 MB/sample) | D | dnbc4tools side products: the BAM is scRNA's weakest asset (downstream reads matrices, not BAMs); the two per-cell results are terminal measurement outputs — S4's tarball archives them |

### ecDNA (`rules/ecdna.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `ampliconsuite` ×10 | `ecdna/{pair_id}/{pair_id}_{ecDNA_counts,gene_list,feature_basic_properties,amplicon_classification_profiles,result_table}.tsv` + `.suite_finished` | **A** | AmpliconClassifier products lifted from `<pair>_classification/`; the suite's internal CNVkit/AA intermediates stay D-class under their own `<pair>_*` dirs |
| `ecdna_cohort` ×1 | `report/cohort_ecdna.tsv` + `report/ecdna_genes.tsv` | **A** | collected by `collect_ecdna.py` |

### Neoantigen (`rules/neoantigen.oxoflow`, venus#31 — tumor_only == 'true' only)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `vep_reannotate_plugins` ×pairs | `neoantigen/{pair_id}/vep_plugins.vcf.gz` + `.tbi` | (undeclared-class: no marker) | plugin CSQ fields pVACtools requires; regenerable from the PASS VCF in minutes |
| `hla_type_dna` ×DNA-samples | `hla/{sample}/hla.result.txt`, `hla/{sample}.alleles.tsv` | (undeclared-class) | SpecHLA typing; extract FASTQs under `hla/{sample}/` are D-class |
| `gene_symbol_map` ×1 | `ref/gene_symbols.tsv` | (undeclared-class) | GTF-derived gene_id→symbol cache |
| `pvacseq_run` ×pairs | `neoantigen/{pair_id}/filtered.tsv`, `alleles.txt` | (undeclared-class) | `run/` dir (MHC_Class_I predictions) is D-class scratch |
| `neoantigen_report` ×pairs | `neoantigen/{pair_id}.{candidates.tsv,summary.tsv,report.md}` | (undeclared-class) | per-patient burden deliverables |
| `pvacfuse_run_paired` ×pairs | `neoantigen/{pair_id}/{fuse.filtered.tsv,fuse.alleles.txt}` | (undeclared-class) | pVACfuse on Manta PASS SV via `manta_to_arriba.py`; `fuse_run/` scratch is D-class |
| `pvacfuse_run_tumoronly` ×pairs | `neoantigen/{pair_id}/{fuse.candidate.filtered.tsv,fuse.alleles.txt}` | (undeclared-class) | complementary when-gate twin on candidateSV — distinct filename so unscored calls are never confused with PASS-scored ones |
| `scrna_corroboration` ×pairs | `scrna_corr/{pair_id}.tsv` (protected) | (undeclared-class) | opt-in `run_scrna_corr`; per-gene %positive from the pair's filter_feature.h5ad (no malignant annotation — caveat in report text) |
| `neoantigen_report_scrna` ×pairs | (same 3 paths as `neoantigen_report`) | (undeclared-class) | complementary when-gate twin adding `--scrna-corr`; exactly one of the two runs per pair |
| `msi_frameshift_page` ×1 | `report/msi_frameshift_cohort.{tsv,md}` | (undeclared-class) | MSI-H × frameshift join page; every missing cell is the explicit sentinel `not_available` (see `scripts/msi_frameshift_join.py`) |

### HLA QC (`rules/hla_qc.oxoflow`, venus#31 P1 items 3+4)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `hla_type_rna` ×RNA-samples | `hla_rna/{sample}/{sample}.genotype.json` + `logs/{extract,genotype}.log` | (undeclared-class) | opt-in `run_hla_rna`; arcasHLA extract+genotype from the STAR BAM |
| `hla_dna_rna_qc` ×pairs | `hla_qc/{pair_id}.consistency.tsv` (protected) | (undeclared-class) | SpecHLA-vs-arcasHLA 2-field consistency; fires only for pairs with both layers (expand_inputs self-gating) |
| `hla_loh` ×pairs | `hla_loh/{pair_id}.loh.tsv` (protected) | (undeclared-class) | opt-in `run_hla_loh` + paired-only; spechla-loh with ASCAT purity/ploidy (DASH has no bioconda/PyPI package — see docs/method-rationale.md) |

### TCR (`rules/tcr.oxoflow`, venus#31 P1 item 2)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `tcr_from_rna` ×RNA-samples | `tcr/{sample}.tcrb.clonotypes.tsv` + `tcr/{sample}/trust_report.tsv` (both protected) | (undeclared-class) | opt-in `run_tcr`; TRUST4 on the STAR BAM, TCRbeta-clonotype filtered by `tcr_report.py` |

### Reports & cohort tables (`rules/report.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `methods_from_rule_runs` ×1 | `report/methods_section.md` | **A** | generated from the commands that actually ran |
| `tmb_per_pair` ×10 | `report/tmb/{pair_id}.tmb.tsv` | **A** | detection result |
| `scrna_metrics_cohort` ×1 | `report/scrna_metrics.tsv` | **A** | |
| `msi_cohort` ×1 | `report/msi_cohort.tsv` | **A** | detection result |
| `msi_cohort_tumoronly` | (same path) | — | defined, 0 instances (paired cohort) |
| `ecdna_cohort` ×1 | `report/cohort_ecdna.tsv` + `report/ecdna_genes.tsv` | **A** | see ecDNA module table |
| `cohort_tables` ×1 | `report/cohort_summary.tsv` | **A** | |
| `clinical_report` ×1 | `report/venus_clinical_report.md` + `.html` | **A** | final deliverable |

### Clinical interpretation (`rules/clinical.oxoflow`, added post-0.5.0)

Every output here is cheap to regenerate from protected upstreams (VEP VCF,
CNVkit cns, Manta SV, TMB/MSI/signatures) plus the host KB dumps — class **B**
by regeneration cost, **A** by informativeness. Per the module-protection
convention the evidence/biomarker/therapy tables are declared
`protected_output`, so they survive all clean strategies regardless.

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `clinical_kb_prep` ×1 | `report/clinical/kb_versions.tsv` | **A** | KB provenance consumed by the report |
| `clinical_snv_annotate` ×pairs | `report/clinical/{pair_id}.snv_evidence.tsv` (protected); `{pair_id}.annotated.vcf.gz`(+`.tbi`), `{pair_id}.gnomad.vcf.gz` | B | gnomAD/ClinVar joins + CIViC tiering; VCFs regenerable in minutes |
| `clinical_oncokb` / `_skip` ×pairs | `report/clinical/{pair_id}.snv_evidence.onco.tsv` | B | token-gated enrichment or byte-copy passthrough |
| `clinical_cnv_genes` ×pairs | `report/clinical/{pair_id}.cnv_evidence.tsv` | B | panel calls from call.cns + ecDNA cross-flag |
| `clinical_sv_fusions_{paired,tumoronly}` ×pairs | `report/clinical/{pair_id}.sv_evidence.tsv` | B | complementary when-gates share the path |
| `clinical_biomarkers_{paired,tumoronly}` ×pairs | `report/clinical/{pair_id}.{biomarkers,therapies}.tsv` | B/A | per-patient deliverables |
| `clinical_cohort` ×1 | `report/clinical_cohort.tsv` | **A** | feeds `clinical_report --clinical` |

### Never-tracked files (all D)

| Path | Producer | Count (this cohort) |
|---|---|---|
| `align/{sample}/{sample}.markdup.bai` | Picard `--CREATE_INDEX` | 20 |
| `manta/{pair_id}/workspace/`, `results/{evidence,stats}/`, `runWorkflow.py`, `*.config.pickle`, `workflow.*.log.txt` | Manta | 10 dirs |
| `manta/{pair_id}/results/variants/{candidateSV,candidateSmallIndels,diploidSV}.vcf.gz(.tbi)` (paired) | Manta | ~5 files/pair |
| `vcf.raw/{pair_id}/{chr}.vcf.gz.tbi` | auto-created scatter-VCF index (undeclared) — survives `clean` | 240 |
| `cnv/{pair_id}/GRCh38*.bed` | CNVkit inputs | ~2/pair |
| `ascat/{pair_id}/{experiment}_normal{LogR,BAF}.txt`, `{experiment}_normalBAF_rawBAF.txt` | ASCAT germline allele evidence (written by `run_ascat.R`) | 30 |
| `rna/star/{sample}/{sample}.SJ.out.tab` | STAR splice-junction inventory (`--outFileNamePrefix rna/star/{sample}/{sample}.`) | 10 |
| `rna/star/{sample}/{sample}.Log.out` / `.Log.progress.out`, `{sample}._STARgenome/`, `{sample}._STARpass1/`, one interrupted-run `*.bam.oxo-failed` | STAR run logs + 1st-pass genome dirs | ~20+ |
| `scrna/count/{sample}/outs/anno_decon_sorted.bam` + `.bai` | dnbc4tools (scRNA's weakest asset — downstream reads matrices, not BAMs) | 20 |
| `scrna/count/{sample}/outs/filter_matrix/` + `{sample}_scRNA_report.html` | dnbc4tools — never-tracked but **protected** (see scrna module table), so they survive every clean anyway | 30 |
| `scrna/count/{sample}/outs/analysis/QC_Cluster.h5ad`, `outs/singlecell.csv` | dnbc4tools per-cell results (terminal, not consumed by `scrna_qc_cluster`) | 20 |
| `msi/{pair_id}_dis` (~1.1 GB each ≈ 11 GB cohort-wide) | msisensor-pro per-pair distributions (`-o msi/{pair_id}` prefix) — the only deletable MSI leftover | 10 |
| `msi/{pair_id}_all` (~20 MB), `{pair_id}_unstable` (114 B) | msisensor-pro side products — never-tracked but **protected** (see MSI module table) | 10 each |
| `logs/**` | every rule | ~560 |
| `backups/cleanup-*.tar.gz` | `clean_intermediates.sh --apply` | per invocation |
| `.oxo-flow/` | engine state (checkpoint, env manifests) | — engine-internal, never delete |

## Deliberate losses (sanity-checked)

| Artifact | Why losing it is acceptable | Cost to regret it |
|---|---|---|
| trim FASTQs (B) | regenerate from `raw/` by re-running fastp | ~1 h/sample |
| sorted + markdup BAMs (B/C) | realign from trim FASTQs | **days** (full BWA-MEM2 over 20 WGS) |
| `vcf.raw` scatter chain (B) | re-run Mutect2 | hours–~day: the bqsr BAMs survive (A), so **no** realignment is needed |
| manta `workspace/`, `evidence/`, cnvkit BEDs (D) | tool-internal scratch | none for results |
| manta `workspace/`, `evidence/`, cnvkit BEDs (D) | tool-internal scratch | none for results (`{pair_id}_all` MSI locus detail is NOT on this list — it backs the protected `report/msi_cohort.tsv` BH correction and is tarred by S4) |

## Coverage statement

The census tables above were generated when the repo had 13 rule files
(40 rules); the manifest now holds **18 `rules/*.oxoflow` files (76 rules)**
after the ecDNA/clinical/neoantigen/hla_qc/tcr additions — those module tables
are included above, but the 561-instance / 1,553-file cohort census predates
them and counts only the rules that existed in the reference `processed-wsx` run
(`2026-10` ground truth). The per-instance counts for the newer modules depend
on cohort mode (neoantigen rules instantiate only with `tumor_only=true`;
tcr/hla_qc rules are opt-in via `run_tcr`/`run_hla_rna`/`run_hla_loh`).
If a rule is added, regenerate the census with the snippet at the top
and extend the matching module table.
