# Venus (启明星)

**Venus** is a tumor multi-omics pipeline built on [oxo-flow](https://github.com/Traitome/oxo-flow), covering the full path from raw FASTQ to a cohort clinical report: somatic SNV/indel calling, copy-number analysis, bulk RNA-seq, single-cell RNA-seq, and an auto-generated methods section.

## What it does

| Module | Tools | Outputs |
|---|---|---|
| QC | fastp, FastQC, MultiQC | trimmed FASTQs, `qc/multiqc/multiqc_report.html` |
| DNA alignment | BWA-MEM2, Picard MarkDuplicates, GATK BQSR | sorted BAM + BAI per sample |
| Somatic calling | GATK Mutect2 (per-chromosome scatter) + FilterMutectCalls with orientation-bias priors | per-pair VCF (PASS), `stats`, MAF |
| Annotation | vep-rs (JSON cache) → VCF + MAF | `deliver/{pair}.somatic.pass.vep.vcf.gz`, `deliver/{pair}.maf` |
| CNV | CNVkit (`wgs` method; `--reference` for tumor-only) | per-pair CNV profiles |
| Bulk RNA | 2-pass STAR, featureCounts | BAM, gene-count matrix, alignment metrics |
| Single-cell | dnbc4tools (chemistry auto), scanpy | count matrix, QC/clustering metrics |
| Report | venus helper scripts | cohort TSV, TMB, clinical report with methods generated from what actually ran |

Analysis modes:
- **Paired tumor–normal** (`pairs.tsv` with a control sample) — Mutect2 tumor/normal, CNVkit `--normal`.
- **Tumor-only** (empty control column) — Mutect2 without `-normal`, CNVkit with a prebuilt normal `.cnn` reference. Gating is per-pair via `when = "wildcard.control != ''"`, so a single cohort may mix both modes.
- **WGS vs WES/panel** — set `target_bed` to the capture intervals for WES/panel (Mutect2 `-L`, CNVkit `--targets` via `cnv_targets_args`, TMB denominator `target_mb`); leave empty for WGS.

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

Platforms: BGI (DNBSEQ/MGISEQ) and Illumina FASTQs are interchangeable downstream of fastp; single-cell libraries are handled by dnbc4tools `--chemistry auto` (BGI DNBSEQ and 10x-style both supported).

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
rules/*.oxoflow      # rule fragments (qc, align, varcall, varcall_merge, annotation, cnv, rna, scrna, report)
scripts/             # python helpers (MAF conversion, TMB, cohort tables, clinical report, …)
envs/<mod>/pixi.toml # pixi environments, one dir per module (env specs in rules point here)
config/pairs.tsv     # tumor/control pairing (control may be empty for tumor-only)
config/groups.tsv    # sample groups + assay metadata (DNA, bulkRNA, scRNA)
```

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
oxo-flow run --arg run_bqsr=false --arg target_bed=/data/exome.bed --arg target_mb=34
```

## Outputs of note

- `deliver/{pair_id}.somatic.pass.vcf.gz` — PASS somatic variants (unannotated)
- `deliver/{pair_id}.somatic.pass.vep.vcf.gz` — VEP-annotated VCF
- `deliver/{pair_id}.maf` — MAF (2.4.1 subset) for downstream TMB/mutational-signature tools
- `report/tmb/{pair_id}.tmb.tsv` — TMB (mut/Mb, non-silent coding)
- `report/clinical_report.html` — cohort report; its methods section is generated from `rule_runs` in `.oxo-flow/checkpoint.json`, so it always describes the commands that actually executed
- `qc/multiqc/multiqc_report.html`, `rna/qc/rna_qc_summary.tsv`, `report/scrna_metrics.tsv`

## License

MIT — see [LICENSE](LICENSE).
