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

Eight rules produce pure scratch that is fully regenerable from inputs and only
feeds downstream rules:

| Rule | Scratch output | Approx. size |
|------|----------------|--------------|
| `fastp` (qc.oxoflow) | `qc/trim/{sample}_R{1,2}.trim.fastq.gz` | ~1 TB cohort-wide |
| `bwa_mem2` (align.oxoflow) | `{experiment}.bam` (unsorted pre-markdup) | hundreds of GB |
| `mutect2_paired` / `mutect2_tumoronly` | per-chunk `chunks/{pair}/{chunk}.vcf.gz(.tbi)` | tens of GB |
| `gather_chr_vcfs` (varcall_merge.oxoflow) | `vcf/{pair_id}.somatic.vcf.gz(.tbi)` | GBs |
| `merge_mutect_stats` | `stats/{pair_id}.mutect_stats.tsv` | KBs |
| `learn_orientation` | `vcf/{pair_id}.learn.json` | KBs |
| `filter_mutect` | `vcf/{pair_id}.filtered.vcf.gz(.tbi)` | GBs |

These carry `temporary = true` in their rule definitions. Once the engine
fix lands, a fully successful run tombstones them automatically (delete after
success, lazy-regenerate on demand via the tombstone cascade). Leaf rules and
protected outputs are never touched.

**Known engine limitation (2026-10):** tombstone cleanup silently skips every
rule whose outputs contain `{wildcard}` placeholders — no deletion, no
tombstone, no diagnostic. All eight venus scratch rules are wildcard rules, so
mode A is currently a no-op. Do NOT work around this with manual `rm`: deleting
outputs without checkpoint bookkeeping marks the producers stale and would
trigger recomputation of the completed campaign on the next `oxo-flow run`.
Tracked upstream; until then use `--mode all` to reclaim the same scratch
(see below).

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
- **QC**: multiqc report + general-stats table.

### What is deleted and what that costs

Everything else: trim FASTQs, alignment BAMs and their indexes (BQSR, markdup),
per-chunk and intermediate VCFs, manta workspaces, CNVkit scratch, raw scRNA
matrices. On the current cohort this is the bulk of the ~2 TB of scratch.

**Recompute-on-demand trade-off.** `oxo-flow clean` deletes files but does not
edit `.oxo-flow/checkpoint.json`. After a clean, `oxo-flow plan` reports every
cleaned rule as *stale* (outputs missing) — that is intended. Nothing reruns
until you ask for it; a later `oxo-flow run` lazily regenerates only the rules
whose products you actually need (and their ancestors), recomputing the rest of
the DAG for free.

### MSI caveat

`rules/msi.oxoflow` is intentionally untouched, so `deliver/{pair_id}.msi.tsv`
and `msi/{pair_id}_all` are **not protected** and a full clean deletes them.
Re-running msisensor costs 1–3 h per pair. `clean_intermediates.sh --apply`
backs up `deliver/ report/ ascat/ msi/` into a timestamped
`backups/cleanup-<stamp>.tar.gz` before deleting (skip with `--no-backup`).
If you edit `msi.oxoflow` later, add `protected_output` there and the script's
backup stays harmless.

### Undeclared leftovers not covered by `clean`

The engine only deletes *declared* rule outputs. These side products remain
after a clean and must be removed manually if space is critical:

- `bwa/{experiment}.markdup.bam.bai` — created via `--CREATE_INDEX` inside the
  mark_duplicates command rather than declared as output.
- `sv/manta/{pair}/workspace/` and `results/` — manta scratch beyond the
  declared candidate VCFs.
- `cnv/{pair_id}/{experiment}*.target.tumor.cnn` style helper files beyond the
  protected per-pair set.

## Verifying before you clean

```bash
oxo-flow plan    # must show "would run 0 | skip 561" for the completed cohort
oxo-flow clean venus.oxoflow   # dry-run preview: exactly what would be deleted
```

If `plan` shows anything other than 0 pending rules, stop and investigate
before cleaning — the fingerprint of a rule does not include `temporary` or
`protected_output`, so adding the markers in this change does not invalidate
any completed run.
