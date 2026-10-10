# Cleaning intermediates in the venus workdir

The venus cohort workdir (`processed-wsx`) holds ~2.7 TB across 1553 outputs
for 10 CRC pairs. Two cleanup modes are supported:

| Mode | Flag | Mechanism | Status |
|------|------|-----------|--------|
| **A — "useless" scratch only** | `--mode useless` | engine `temporary = true` tombstones | **blocked on engine fix** (see below) |
| **B — everything not protected** | `--mode all` | engine `oxo-flow clean` + `protected_output` | **works today** |

Both go through `scripts/clean_intermediates.sh`. The script never deletes
anything without `--apply`; the default invocation is a dry-run preview.

```bash
scripts/clean_intermediates.sh                    # preview mode B (dry-run)
scripts/clean_intermediates.sh --mode all --apply # mode B: backup + clean
scripts/clean_intermediates.sh --mode useless     # mode A status (no deletions yet)
```

## Mode A — clean completely useless scratch (`temporary = true`)

Seven rules produce pure scratch that is fully regenerable from inputs and only
feeds downstream rules:

| Rule | Scratch output | Approx. size |
|------|----------------|--------------|
| `fastp` (qc.oxoflow) | `qc/trim/{sample}_R{1,2}.trim.fastq.gz` | ~1 TB cohort-wide |
| `bwa_mem2` (align.oxoflow) | `align/{sample}/{sample}.sorted.bam` | hundreds of GB |
| `mutect2_paired` / `mutect2_tumoronly` | `vcf.raw/{pair_id}/{chr}.vcf.gz` + `.vcf.gz.stats` + `{chr}.tar.gz` | tens of GB |
| `gather_chr_vcfs` (varcall_merge.oxoflow) | `vcf.raw/{pair_id}.vcf.gz` + `.tbi` | GBs |
| `merge_mutect_stats` (varcall_merge.oxoflow) | `vcf.raw/{pair_id}.stats` | KBs |
| `learn_orientation` (varcall_merge.oxoflow) | `vcf.filtered/{pair_id}.orientation-bias.tar.gz` | KBs |

`fastp` keeps its KB-scale `.fastp.json`/`.fastp.html` protected even though
the rule is temporary — they are inputs to the protected MultiQC report. The
full filtered VCFs (`vcf.filtered/{pair_id}.vcf.gz`) are **not** in this
table on purpose: they are the only copy of the filtered-out variant records
(see Mode B below), so `filter_mutect` is protected instead of temporary.

These carry `temporary = true` in their rule definitions. Once the engine
fix lands, a fully successful run tombstones them automatically (delete after
success, lazy-regenerate on demand via the tombstone cascade). Leaf rules and
protected outputs are never touched.

**Known engine limitation (2026-10):** tombstone cleanup silently skips every
rule whose outputs contain `{wildcard}` placeholders — no deletion, no
tombstone, no diagnostic. All seven venus scratch rules are wildcard rules, so
mode A is currently a no-op. Do NOT work around this with manual `rm`: deleting
outputs without checkpoint bookkeeping marks the producers stale and would
trigger recomputation of the completed campaign on the next `oxo-flow run`.
Tracked upstream as [Traitome/oxo-flow#844](https://github.com/Traitome/oxo-flow/issues/844);
until it lands use `--mode all` to reclaim the same scratch (see below).

## Mode B — one-click clean of all intermediates (`protected_output` + `clean`)

Rules declare the artifacts worth keeping:

```text
protected_output = ["deliver/{pair_id}.somatic.pass.vcf.gz", ...]
```

`oxo-flow clean` collects declared outputs **plus** `protected_output`
patterns, expands wildcards, and deletes only what is not protected.

### What survives a full clean

- **Delivery** (`deliver/`): PASS somatic VCF+TBI, VEP-annotated VCF+TBI, MAF,
  SV VCF+TBI (paired + tumor-only).
- **Full filtered call sets** (`vcf.filtered/{pair_id}.vcf.gz` + `.tbi`):
  `select_pass` keeps only the PASS subset, so the records rejected by
  FilterMutectCalls (germline/strand-bias/PoN) exist only here — primary
  evidence for sensitivity audits. Regenerating them requires the whole
  `vcf.raw` chain, whose inputs (bqsr BAMs) a clean deletes, i.e. realignment.
  The tiny `.filteringStats.tsv` per pair is *not* protected.
- **Manta scored calls** (`manta/{pair_id}/results/variants/somaticSV.vcf.gz`
  + `.tbi` paired, `candidateSV.vcf.gz` + `.tbi` tumor-only): the full scored
  somatic set behind the PASS-only delivery VCF, and the unfiltered candidate
  superset. Manta re-runs need the deleted bqsr BAMs.
- **Report** (`report/`): clinical report MD+HTML, methods section, per-pair
  TMB, MSI cohort table, scRNA metrics, cohort summary, signature counts +
  exposures.
- **ASCAT** (`ascat/`): per-pair RDS, purity/ploidy, segments, metrics, LogR/BAF
  tracks, all plots, cohort TSV — the most expensive CNV step, report-facing.
- **CNVkit** (`cnv/`): per-pair `.cnr`/`.call.cns`/`.cnv.png` + reference.cnn
  (the `{pair_id}` files are symlinks to `{experiment}` files; protecting only
  the aliases would leave dangling links, so both sides are protected).
- **RNA**: STAR `Log.final.out` + `ReadsPerGene.out.tab` (QC products),
  featureCounts per-sample tables, cohort count matrix, RNA QC summary. The big
  STAR BAM is deliberately **not** protected — it is the single largest
  beneficial-to-drop RNA intermediate (regenerable from trim FASTQs).
- **scRNA**: `filter_feature.h5ad`, cluster/marker CSVs, `metrics_summary.xls`,
  QC tables, cohort clusters. The huge `raw_matrix/` is deliberately left
  unprotected (re-running `scrna_count` is the most expensive recompute in the
  pipeline — treat these protections as load-bearing).
- **QC**: multiqc report + general-stats table, plus the per-sample inputs it
  aggregates (`qc/trim/{sample}.fastp.json`/`.fastp.html`,
  `qc/fastqc/{sample}_R{1,2}.trim_fastqc.html`) so the protected report stays
  re-renderable after a clean. Picard `dup_metrics.txt` is *not* protected.

### What is deleted and what that costs

Everything declared and unprotected: trim FASTQs (~1 TB), alignment BAMs and
their declared indexes (`align/{sample}/{sample}.sorted.bam`,
`.markdup.bam` + metrics, `bqsr/{experiment}/*.bqsr.bam` + `.bai` +
`.recal.table`), the `vcf.raw/` scatter chain, the per-pair
`.filteringStats.tsv`, `msi/reference.list` + `deliver/{pair_id}.msi.tsv`,
CNVkit scratch beyond the protected set, scRNA `raw_matrix/`. On the current
cohort this is the bulk of the ~2 TB of scratch; note the undeclared leftovers
below stay behind, so `du` will not drop by the full estimate.

**Recompute-on-demand trade-off.** `oxo-flow clean` deletes files but does not
edit `.oxo-flow/checkpoint.json`. After a clean, `oxo-flow plan` reports every
cleaned rule as *stale* (outputs missing) — that is intended, and the
pre-clean invariant (`would run 0 | skip 561`) only holds *before* a clean.
Nothing reruns until you ask for it; a later `oxo-flow run` lazily regenerates
only the rules whose products you actually need (and their ancestors),
recomputing the rest of the DAG for free.

### MSI caveat

`rules/msi.oxoflow` is intentionally untouched, so a full clean deletes the
declared, unprotected MSI products: `deliver/{pair_id}.msi.tsv` (per-pair
result) and `msi/reference.list` (the cohort-shared msisensor scan, 1–3 h).
The per-pair locus detail `msi/{pair_id}_all` is an undeclared side product of
`msi_paired` and survives by accident — do not rely on it.
`clean_intermediates.sh --apply` backs up `deliver/ report/ ascat/ msi/` into
a timestamped `backups/cleanup-<stamp>.tar.gz` before deleting (skip with
`--no-backup`), which covers all of the above. If you edit `msi.oxoflow` later,
add `protected_output` there and the script's backup stays harmless.

### Undeclared leftovers not covered by `clean`

The engine only deletes *declared* rule outputs. These side products remain
after a clean and must be removed manually if space is critical:

- `align/{sample}/{sample}.markdup.bai` — created by Picard `--CREATE_INDEX`
  (which replaces the `.bam` suffix) rather than declared as an output.
- `manta/{pair_id}/workspace/`, run metadata (`runWorkflow.py`,
  `*.config.pickle`, `workflow.*.log.txt`) and `results/stats/`; paired runs
  also keep undeclared `results/variants/candidateSV.vcf.gz` and
  `diploidSV.vcf.gz` (the latter is the germline call set — worth keeping or
  archiving before deleting).
- `cnv/{pair_id}/GRCh38*.bed` — the target/antitarget BEDs copied in as
  cnvkit inputs.
- `scrna/count/{sample}/outs/singlecell.csv` and other dnbc4tools side
  products outside the declared list.

## Verifying before you clean

```bash
oxo-flow plan    # before any clean: must show "would run 0 | skip 561"
oxo-flow clean venus.oxoflow   # dry-run preview: exactly what would be deleted
```

If `plan` shows anything other than 0 pending rules *before* a clean, stop and
investigate — the fingerprint of a rule does not include `temporary` or
`protected_output`, so adding the markers does not invalidate any completed
run.
