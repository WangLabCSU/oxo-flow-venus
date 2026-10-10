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
| Declared output files | 1,543 (sum of output patterns × instances) |
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
| **A — protected** | Detection/measurement results, delivery artifacts, expensive report-facing products. Survive every clean. | `protected_output = [...]` | **583** (38%) |
| **B — temporary scratch** | Fully regenerable intermediates, only feed downstream rules. Slated for automatic post-success deletion. | `temporary = true` | **840** (54%) |
| **C — declared, deletable** | Intermediates without either marker. Deleted by `clean`, survive normal operation. | declared `output` only | **120** (8%) |
| **D — undeclared leftovers** | Side products the engine never tracks. `clean` cannot see them. | not declared | logs + per-tool extras |

## What each strategy does to each class

| Strategy | A protected | B temporary | C deletable | D undeclared |
|---|---|---|---|---|
| **S1** Auto-tombstone after successful run (Mode A) | keep | delete *(dormant — pending engine [#844](https://github.com/Traitome/oxo-flow/issues/844))* | keep | keep |
| **S2** `oxo-flow clean` (dry-run) | nothing deleted, preview only | | | |
| **S3** `oxo-flow clean --force` (Mode B) | keep | delete | delete | keep |
| **S4** `scripts/clean_intermediates.sh --mode all --apply` | keep (tar backup of `deliver/ report/ ascat/ msi/` first) | delete | delete | keep (tar additionally captures the evidence-grade leftovers: `msi/{pair}_all`, Manta variant VCFs, scRNA per-cell results, STAR junctions) |
| **S5** manual `rm` of a declared file | destroys the result or marks its producer **stale** → next `run` recomputes it and every ancestor | | | no bookkeeping impact, but no safety net |

Notes:

- **Mode A is dormant**: the engine's tombstone loop silently skips every rule
  whose outputs contain `{wildcard}` placeholders ([#844](https://github.com/Traitome/oxo-flow/issues/844)).
  All seven B-class venus rules are wildcard rules, so today only S3/S4
  reclaim scratch. Never work around this with manual `rm` — deleting declared
  outputs without checkpoint bookkeeping makes producers stale and would
  trigger recomputation on the next `oxo-flow run`.
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

### Alignment (`rules/align.oxoflow`) — DNA

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `bwa_mem2` ×20 | `align/{sample}/{sample}.sorted.bam` | B | hundreds of GB; the biggest single scratch class |
| `mark_duplicates` ×20 | `align/{sample}/{sample}.markdup.bam` | C | regenerable from sorted.bam |
| `mark_duplicates` | `align/{sample}/{sample}.dup_metrics.txt` | **A** | library QC measurement; feeds the protected MultiQC report |
| `bqsr` ×20 | `bqsr/{sample}/{sample}.recal.table` / `.bqsr.bam` / `.bqsr.bai` | C | the BAMs every caller consumes; a full clean deletes them → realignment |

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
| — | `msi/{pair_id}_all` | D | per-locus detail, undeclared side product; survives clean by accident, do not rely on it |

### Signatures (`rules/signature.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `sbs96_fit` ×10 | `report/signatures/{pair_id}.sbs96.counts.tsv` + `.exposures.tsv` | **A** | |

### Bulk RNA (`rules/rna.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `star_align` ×10 | `rna/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam` | C | largest beneficial-to-drop RNA intermediate |
| `star_align` | `rna/star/{sample}/{sample}.Log.final.out` + `.ReadsPerGene.out.tab` | **A** | QC products |
| `featurecounts` ×10 | `rna/counts/{sample}.featurecounts.txt` + `.summary` | **A** | |
| `merge_counts` ×1 | `rna/counts/cohort.gene_counts.txt` | **A** | cohort matrix |
| `rna_qc_summary` ×1 | `rna/qc/rna_qc_summary.tsv` | **A** | |

### Single-cell RNA (`rules/scrna.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `scrna_fastqc` ×10 | `scrna/fastqc/{sample}_SC_R{1,2}_fastqc.html` | **A** | QC measurement; feeds the protected MultiQC report |
| `scrna_count` ×10 | `scrna/count/{sample}/outs/raw_matrix/{matrix.mtx,features.tsv,barcodes.tsv}.gz` | C | **re-running scrna_count is the most expensive recompute in the pipeline** — deleting these forces a full recount if any scRNA product is ever needed again |
| `scrna_count` | `outs/filter_feature.h5ad`, `outs/analysis/{cluster,marker}.csv`, `outs/metrics_summary.xls` | **A** | detection results |
| `scrna_qc_cluster` ×10 | `scrna/qc/{sample}.clusters.tsv` | **A** | |
| `scrna_integrate` ×1 | `scrna/cohort_clusters.tsv` | **A** | |
| — | `outs/anno_decon_sorted.bam` + `.bai`, `outs/filter_matrix/`, `{sample}_scRNA_report.html` | D | dnbc4tools side products |

### Reports & cohort tables (`rules/report.oxoflow`)

| Rule (×n) | Output | Class | Notes |
|---|---|---|---|
| `methods_from_rule_runs` ×1 | `report/methods_section.md` | **A** | generated from the commands that actually ran |
| `tmb_per_pair` ×10 | `report/tmb/{pair_id}.tmb.tsv` | **A** | detection result |
| `scrna_metrics_cohort` ×1 | `report/scrna_metrics.tsv` | **A** | |
| `msi_cohort` ×1 | `report/msi_cohort.tsv` | **A** | detection result |
| `cohort_tables` ×1 | `report/cohort_summary.tsv` | **A** | |
| `clinical_report` ×1 | `report/venus_clinical_report.md` + `.html` | **A** | final deliverable |

### Never-tracked files (all D)

| Path | Producer | Count (this cohort) |
|---|---|---|
| `align/{sample}/{sample}.markdup.bai` | Picard `--CREATE_INDEX` | 20 |
| `manta/{pair_id}/workspace/`, `results/{evidence,stats}/`, `runWorkflow.py`, `*.config.pickle`, `workflow.*.log.txt` | Manta | 10 dirs |
| `manta/{pair_id}/results/variants/{candidateSV,candidateSmallIndels,diploidSV}.vcf.gz(.tbi)` (paired) | Manta | ~5 files/pair |
| `cnv/{pair_id}/GRCh38*.bed` | CNVkit inputs | ~2/pair |
| `ascat/{pair_id}/{experiment}_normal{LogR,BAF}.txt` | ASCAT germline allele evidence (written by `run_ascat.R`) | 20 |
| `rna/star/{sample}/_STAR/SJ.out.tab` | STAR splice-junction inventory | 10 |
| `scrna/count/{sample}/outs/{anno_decon_sorted.bam(.bai), filter_matrix/, {sample}_scRNA_report.html}` | dnbc4tools | ~4/sample |
| `scrna/count/{sample}/outs/analysis/QC_Cluster.h5ad`, `outs/singlecell.csv` | dnbc4tools per-cell results (terminal, not consumed by `scrna_qc_cluster`) | 20 |
| `msi/{pair_id}_all` | msisensor-pro | 10 |
| `logs/**` | every rule | ~560 |
| `backups/cleanup-*.tar.gz` | `clean_intermediates.sh --apply` | per invocation |
| `.oxo-flow/` | engine state (checkpoint, env manifests) | — engine-internal, never delete |

## Deliberate losses (sanity-checked)

| Artifact | Why losing it is acceptable | Cost to regret it |
|---|---|---|
| trim FASTQs (B) | regenerate from `raw/` by re-running fastp | ~1 h/sample |
| sorted/markdup/bqsr BAMs (B/C) | realign from trim FASTQs | **days** (full BWA-MEM2 over 20 WGS) |
| `vcf.raw` scatter chain (B) | re-run Mutect2 | needs bqsr BAMs → realignment |
| STAR BAM (C) | re-align from trim FASTQs | hours |
| scRNA `raw_matrix` (C) | re-run dnbc4tools | **most expensive recompute in the pipeline** |
| manta `workspace/`, `evidence/`, cnvkit BEDs (D) | tool-internal scratch | none for results (`{pair_id}_all` MSI locus detail is NOT on this list — it backs the protected `report/msi_cohort.tsv` BH correction and is tarred by S4) |

## Coverage statement

Every `output` pattern in the 13 `rules/*.oxoflow` files (40 rules) appears
exactly once above, classified. The 561-instance census sums to the live
`completed_rules` count, and the declared-output total (1,543) equals the sum
of pattern × instance across the tables. D-class entries were verified against
the live workdir (`processed-wsx`, 2026-10) rather than inferred from tool
docs. If a rule is added, regenerate the census with the snippet at the top
and extend the matching module table.
