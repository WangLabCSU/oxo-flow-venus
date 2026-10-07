# Venus (启明星)

**Venus** is a tumor multi-omics pipeline built on [oxo-flow](https://github.com/Traitome/oxo-flow), covering the full path from raw FASTQ to a cohort clinical report: somatic SNV/indel calling, copy-number analysis, bulk RNA-seq, single-cell RNA-seq, and an auto-generated methods section.

## What it does

| Module | Tools | Outputs |
|---|---|---|
| QC | fastp, FastQC, MultiQC | trimmed FASTQs, `qc/multiqc/multiqc_report.html` |
| DNA alignment | BWA-MEM2, Picard MarkDuplicates, GATK BQSR | sorted BAM + BAI per sample |
| Somatic calling | GATK Mutect2 (per-chromosome scatter) + FilterMutectCalls with orientation-bias priors | per-pair VCF (PASS), `stats`, MAF |
| Annotation | vep-rs (JSON cache) → VCF + MAF | `deliver/{pair}.somatic.pass.vep.vcf.gz`, `deliver/{pair}.maf` |
| CNV | CNVkit (`cnv_method`, default `wgs`; `--reference` for tumor-only) | per-pair CNV profiles |
| SV | Manta (paired somatic + tumor-only candidates; `--exome` for WES) | `deliver/{pair}.sv.vcf.gz` (PASS, paired) or `deliver/{pair}.sv.candidate.vcf.gz` (unscored, tumor-only) |
| MSI | MSIsensor-pro (paired only) | `deliver/{pair}.msi.tsv` |
| Signatures | in-house SBS96 counter + NNLS vs COSMIC v3.1 | `report/signatures/{pair}.sbs96.counts.tsv`, `report/signatures/{pair}.exposures.tsv` |
| Bulk RNA | 2-pass STAR, featureCounts | BAM, gene-count matrix, alignment metrics |
| Single-cell | dnbc4tools (`--chemistry {config.scrna_chemistry}`), scanpy | count matrix, QC/clustering metrics |
| Report | venus helper scripts | cohort TSV, TMB, clinical report with methods generated from what actually ran |

Analysis modes:
- **Paired tumor–normal** (`pairs.tsv` with a control sample) — Mutect2 tumor/normal, CNVkit `--normal`.
- **Tumor-only** (empty control column) — Mutect2 without `-normal`, CNVkit with a prebuilt normal `.cnn` reference. Gating is per-pair via `when = "wildcard.control != ''"`, so a single cohort may mix both modes.
- **WGS vs WES/panel** — set `target_bed` to the capture intervals for WES/panel, plus `mutect_target_args`, `cnv_method`, `cnv_targets_args`, `sv_exome_args`, and `tmb_target_args` per the checklist below; leave empty for WGS.

## Groups and assay gating

`config/groups.tsv` maps sample ids to assay modules via three columns: `name`
(group name), `samples` (comma-separated ids), and `assay` (a group metadata
column that fans out per sample as the `wildcard.assay` value). Rules are
gated per sample instance with `when`:

- `when = "wildcard.assay == 'wgs'"` — DNA chain: fastp/fastqc run for DNA and
  bulk RNA (any non-scRNA assay), but BWA-MEM2/MarkDuplicates/BQSR only for DNA.
- `when = "wildcard.assay == 'rnaseq'"` — STAR + featureCounts.
- `when = "wildcard.assay == 'scrna'"` — dnbc4tools counting (raw FASTQs; no trimming).
- Gates may combine config flags: `when = "config.run_bqsr == 'true' && wildcard.assay == 'wgs'"`.

The `assay` name is arbitrary — any metadata column in `groups.tsv` becomes a
`wildcard.<key>` available to rules. A rule that mentions `{sample}` fans out
over **every** group; the `when` gate then prunes non-matching instances at
plan time, so `star_align` never tries to align a DNA sample. Rules that
aggregate with `expand_inputs` should expand over per-group lists
(`config.samples_DNA`, `config.samples_bulkRNA`, `config.samples_scRNA` —
auto-injected from groups.tsv) rather than the global `config.samples_list`,
which now mixes assays.

Note: a sample id declared both in `groups.tsv` and in `pairs.tsv`
(tumor/control) is allowed — the engine warns about duplicate owners but the
rule paths here are group/pair-disjoint, so outputs never collide.

Platforms: BGI (DNBSEQ/MGISEQ) and Illumina FASTQs are interchangeable downstream of fastp. Single-cell libraries are counted by dnbc4tools 3.0; the chemistry is pinned per run via the `scrna_chemistry` config key (default `scRNAv2HT`; accepted values `auto`, `scRNAv1HT`, `scRNAv2HT`, `scRNAv3HT`, `scRNAv5P` — MGI stomics/DNBSEQ kits). 10x Genomics libraries are not covered by the current scRNA module and would need an additional counting tool.

## Library strandedness

`featureCounts` runs **unstranded** (`-s 0`, the rule sets no `-s` flag). For
the bulk RNA-seq data this was verified empirically rather than assumed — on a
200k read-pair subset of one sample aligned with STAR (93.4% unique, 147k
spliced reads), three independent read-geometry checks against GENCODE exon
structure all landed at ~50/50:

| Check | Same-strand as R1 | Opposite | Verdict |
|---|---|---|---|
| R1 orientation vs GTF exon strand | 50.4% | 49.6% | unstranded |
| R1 vs STAR `intronMotif` XS tag | 50.3% | 49.7% | unstranded |
| R1 vs splice-donor GT/CT motif from FASTA | 49.4% | 50.6% | unstranded |

Unstranded counting is strand-agnostic and remains correct regardless of the
true library type; if your own library is known stranded, set `-s 1` (forward)
or `-s 2` (reverse) in `rules/rna.oxoflow`.

## Repository layout

```
venus.oxoflow        # main workflow: config + include list + chromosome scatter
rules/*.oxoflow      # rule fragments (qc, align, varcall, varcall_merge, annotation,
                     #   cnv, sv, msi, signature, rna, scrna, report)
scripts/             # python helpers (MAF conversion, TMB, SBS96 fit, cohort tables,
                     #   clinical report, …)
envs/<mod>/pixi.toml # pixi environments, one dir per module (env specs in rules point here)
resources/           # vendored data tables (COSMIC v3.1 SBS96 GRCh38 signature matrix)
config/pairs.tsv     # tumor/control pairing (control may be empty for tumor-only)
config/groups.tsv    # sample groups + assay metadata (DNA, bulkRNA, scRNA)
```

## Extending venus with a new module

A module is one rule fragment plus one environment; nothing else changes:

1. **`envs/<mod>/pixi.toml`** — the tool's conda-forge/bioconda deps. Keep it
   minimal; one env per module avoids the version conflicts that a shared
   env accumulates (e.g. Manta needs python 2.7 while delivery filtering
   needs modern bcftools — hence the separate `envs/sv` and `envs/varcall`).
2. **`rules/<mod>.oxoflow`** — a `[workflow]` header + `[[rules]]`. Conventions
   the rest of the pipeline relies on:
   - inputs reference upstream deliverables by path (`bqsr/{experiment}/…`,
     `deliver/{pair_id}.maf`) — the engine resolves producers automatically;
   - per-pair deliverables land in `deliver/{pair_id}.<mod>.<ext>`, per-sample
     analyses in `report/<mod>/{pair_id}.<ext>`; run logs go to `logs/`;
   - gate paired-only rules with `when = "wildcard.control != ''"` so a
     tumor-only pair simply skips the instance instead of failing;
   - tunables live as `[config]` keys in `venus.oxoflow` (overridable at run
     time with `--arg KEY=VALUE`), not as literals inside rules.
3. **Wire it up**: `[[include]] path = "rules/<mod>.oxoflow"` in `venus.oxoflow`
   plus any new `[config]` keys. One-off data tables (signature matrices,
   panels) go in `resources/` and are declared as rule inputs so checksums
   cover them.

### WES/panel compatibility checklist

WGS is the default; a WES/panel run needs these config keys set:

| Key | Purpose |
|---|---|
| `target_bed` | capture intervals — the documented anchor for the interval set; tools receive it via the `*_args` fragments below |
| `mutect_target_args` | e.g. `-L /path/regions.bed` — restricts Mutect2 to the capture targets |
| `cnv_method` | `hybrid` for WES/panel capture (default `wgs`) |
| `cnv_targets_args` | e.g. `--targets /path/regions.bed` for CNVkit |
| `sv_exome_args` | set `--exome` so Manta's off-target depth filters keep real WES events |
| `tmb_target_args` | e.g. `--target-bed /path/regions.bed` so `all_somatic_mut_per_mb` uses merged capture Mb (the clinical `tmb_mut_per_mb` uses coding Mb for both WGS and WES) |
| `msi_coverage` | `20` for WES depth floors — capture-restricted MSI is NOT supported (msisensor-pro has no BED flag) |

### Operating on a run workdir (vendored workflow copy)

A run workdir (`processed-wsx/`) contains its own copy of `venus.oxoflow`,
`rules/`, `scripts/`, `envs/`, `resources/`. oxo-flow binds environment
manifests and checkpoint fingerprints to **workdir-resolved paths**, so:

- run and dry-run from the workdir against its own copy (`oxo-flow dry-run
  venus.oxoflow`), never against the checkout elsewhere — otherwise every
  env manifest resolves to the checkout and everything looks stale;
- after editing the repo, sync changed files into the workdir first
  (`rsync -a rules/ scripts/ envs/ venus.oxoflow <workdir>/`);
- a config edit only invalidates rules that *reference* the changed keys;
  `oxo-flow dry-run venus.oxoflow` reports the exact stale/skip/completed
  split before anything runs (the `Summary: N rules` line, by contrast,
  counts total graph instances, not stale ones).

## Requirements

- Linux x86-64, oxo-flow CLI ≥ 0.23
- [pixi](https://pixi.sh) on PATH (workflow environments are resolved per rule into `.oxo-flow/pixi/`)
- Reference data (see `venus.oxoflow [config]` for the exact keys): GRCh38 FASTA + BWA-MEM2 index, STAR index + GENCODE GTF, GATK resource bundle (dbsnp/mills/gnomad/PoN), vep-rs binary + JSON cache, single-cell reference for dnbc4tools, CNVkit normal reference for tumor-only CNV.

## Quick start

```bash
git clone https://github.com/WangLabCSU/oxo-flow-venus.git
cd oxo-flow-venus

# 1. Stage inputs (relative to your run workdir)
#    raw/{sample}_R1.fastq.gz / raw/{sample}_R2.fastq.gz  (WGS + bulk RNA)
#    single-cell FASTQs keep their library naming — see rules/scrna.oxoflow
# 2. Edit config/pairs.tsv and config/groups.tsv for your cohort
#    (groups.tsv assigns each sample id an assay: wgs / rnaseq / scrna)
# 3. Point the absolute paths in venus.oxoflow [config] at your references

# Validate, then run
oxo-flow validate
oxo-flow run --background

# Reports land in report/ ; delivery files in deliver/
```

Every config key can be overridden per run without editing the file:

```bash
oxo-flow run --arg run_bqsr=false --arg target_bed=/data/exome.bed \
  --arg mutect_target_args="-L /data/exome.bed" --arg cnv_method=hybrid \
  --arg cnv_targets_args="--targets /data/exome.bed" --arg sv_exome_args="--exome" \
  --arg tmb_target_args="--target-bed /data/exome.bed" --arg msi_coverage=20
```

## Outputs of note

- `deliver/{pair_id}.somatic.pass.vcf.gz` — PASS somatic variants (unannotated)
- `deliver/{pair_id}.somatic.pass.vep.vcf.gz` — VEP-annotated VCF
- `deliver/{pair_id}.maf` — MAF (2.4.1 subset) for downstream TMB/mutational-signature tools
- `deliver/{pair_id}.sv.vcf.gz` — Manta somatic SVs filtered to PASS (paired); `deliver/{pair_id}.sv.candidate.vcf.gz` — unscored candidates (tumor-only)
- `deliver/{pair_id}.msi.tsv` — MSIsensor-pro MSI score (paired runs only; tumor-only MSI needs a panel-of-normals baseline)
- `report/signatures/{pair_id}.sbs96.counts.tsv` / `.exposures.tsv` — SBS96 channel counts and NNLS-fitted COSMIC v3.1 exposures (with cosine similarity in the run log)
- `report/tmb/{pair_id}.tmb.tsv` — TMB under two conventions: `tmb_mut_per_mb` (coding non-silent per coding Mb — clinical, ≥10 mut/Mb cutoff) and `all_somatic_mut_per_mb` (all PASS per target Mb — WGS literature convention)
- `report/venus_clinical_report.html` — cohort report; its methods section is generated from `rule_runs` in `.oxo-flow/checkpoint.json`, so it always describes the commands that actually executed
- `qc/multiqc/multiqc_report.html`, `rna/qc/rna_qc_summary.tsv`, `report/scrna_metrics.tsv`

## License

MIT — see [LICENSE](LICENSE).
