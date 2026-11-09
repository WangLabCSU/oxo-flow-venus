# Method rationale — venus multi-omics pipeline

Per-module justification: what each step does, why this specific tool was
chosen, the evidence behind it, and its known limitations. The rendered
clinical report (`report/venus_clinical_report.md`) links here; the *Methods*
section of the report is generated from the actual execution record
(`methods_from_rule_runs.py` over the oxo-flow checkpoint), so it always
describes what really ran — this document explains **why** those tools were
selected.

Versions listed are the ones verified from this run's artifacts (BAM `@PG`
records, tool logs, binary probes). All steps run under the oxo-flow engine
with per-rule pinned pixi environments, content-hashed input manifests and
stored rule fingerprints, so any change to code, parameters, references or
environment files invalidates exactly the affected downstream steps.

---

## 1. Read QC and trimming — fastp + FastQC, aggregated with MultiQC

**What.** Adapter trimming and per-base quality trimming of paired-end FASTQs
(`fastp`), independent per-sample module QC (`fastqc`), cohort aggregation of
all QC reports into one page (`multiqc`).

**Why.** fastp performs adapter detection (`--detect_adapter_for_pe`) and
sliding-window trimming (`--cut_right`, window 4, mean quality 20) in a single
pass, orders of magnitude faster than the classic cutadapt+trimmomatic chain,
and emits a machine-readable JSON report per sample. FastQC is retained as an
independent second opinion at the module level (per-base quality, GC, adapter
contamination) because it is the de-facto standard whose failure modes are
widely understood. MultiQC aggregates both into one cohort view.

**Evidence.** fastp 1.3.7 — Chen et al. 2018, *Nat Biotechnol* 36:4875;
FastQC 0.13.0 — Andrews 2010 (Babraham); MultiQC 1.35 — Ewels et al. 2016,
*Bioinformatics* 32:3671.

**Limitations.** Automatic adapter detection can miss rare/non-standard
adapters; sliding-window trimming shortens reads more aggressively than
single-end trimming, slightly affecting very short inserts.

## 2. DNA alignment — BWA-MEM2

**What.** Alignment of trimmed WGS pairs to GRCh38 (GATK bundle
`GRCh38.d1.vd1`) in MEM mode with read groups (`@RG`), sorted on the fly
(SAMtools 1.24).

**Why.** BWA-MEM is the most widely validated aligner for 100–1000 bp WGS
reads and the reference aligner of the GATK somatic toolchain. BWA-MEM2 2.2.1
reproduces BWA-MEM's alignments (identical seeds/scoring) at roughly twice
the speed with modest memory; keeping `@RG` headers preserves sample
provenance for downstream joint tools.

**Evidence.** Vasimuddin et al. 2019, *IEEE ACS* (BWA-MEM2); Li 2013,
*arXiv:1303.3997* (BWA-MEM); SAMtools 1.24 — Danecek et al. 2021,
*Gigascience* 10:giab008.

**Limitations.** BWA-MEM2's output is BWA-MEM-equivalent but not
bit-identical across versions; alt-contig handling follows the standard
decoy-aware GATK bundle. No read-marking of chimeric pairs beyond the
standard supplementary-alignment tags.

## 3. Duplicate marking — Picard MarkDuplicates

**What.** PCR/sequencing duplicates are *marked* (not removed) using the
library read-group tags.

**Why.** Marking keeps the information available for downstream tools
(Mutect2 preprocesses marked BAMs itself; CNVkit and ASCAT can ignore
duplicates) while leaving the data intact for audit. Picard is the GATK
toolchain's reference implementation.

**Evidence.** Picard 3.5.0 — Broad Institute (github.com/broadinstitute/picard).

**Limitations.** Optical duplicates are only flagged when tile information is
present; duplicate rates in WGS at this depth are low and do not drive
variant calls.

## 4. Base quality score recalibration — GATK BQSR

**What.** Machine-learning recalibration of base quality scores against
dbSNP138 and Mills/1000G gold-standard indels (`Homo_sapiens_assembly38`
resources for hg38), applied in place.

**Why.** Empirical base-quality correction is the standard GATK pre-processing
step and measurably improves SBS quality and Mutect2's error model. Known-sites
masking prevents true variants from being learned as systematic error.

**Evidence.** GATK 4.6.2.0 — McKenna et al. 2010, *Genome Res* 20:1297;
DePristo et al. 2011, *Genome Res* 21:436; Poplin et al. 2018, *Nat
Biotechnol* 36:875.

**Limitations.** Recalibration is reference-build- and known-sites-dependent;
novel systematic errors absent from the known-sites masks are partially
learned as noise. WGS at adequate depth tolerates this well.

## 5. Somatic SNV/indel calling — Mutect2 (paired + tumor-only) with PoN filtering

**What.** Tumor/normal somatic calling (`mutect2_paired`) and tumor-only
calling against the 1000 Genomes panel-of-normals for germline/artifact
filtering (`mutect2_tumoronly`), `FilterMutectCalls`, per-chromosome VCF
gathering, PASS selection, and MAF conversion for downstream reporting.

**Why.** Mutect2 is the reference somatic caller of the GATK ecosystem,
with a lesion-specific error model, orientation-bias learning
(`learn_orientation`) and germline-resource contamination estimation. The
paired mode uses the matched blood normal for subtraction; the tumor-only
mode exists for samples without a matched normal and is *only* interpreted
against the PoN. MAF format enables integration with standard clinical
reporting conventions.

**Evidence.** Cibulskis et al. 2013, *Nature* 505:183; Poplin et al. 2018,
*Nat Biotechnol* 36:875; vcf2maf — Mayakonda et al. 2018, *Leukemia* 32:1152.

**Limitations.** Sensitivity drops below ~5% VAF; PoN filtering removes
recurring artifacts but can also remove true recurrent hotspots shared with
the population; tumor-only calls are never equivalent to paired calls and are
labeled as such in delivery.

## 6. Variant annotation — vep-rs (VEP, GRCh38 cache v113)

**What.** Functional consequence annotation of the delivered VCF/MAF.

**Why.** VEP's consequence model is the community reference used by COSMIC,
ClinVar and dbNSFP integrations. The Rust re-implementation (vep-rs 0.3.1)
avoids the Perl dependency chain while using the standard Ensembl cache
(v113; cache metadata recorded in `vep_cache/info.json`).

**Evidence.** McLaren et al. 2016, *Genome Biol* 17:122 (VEP); vep-rs —
github.com/natera/vep-rs.

**Limitations.** Consequence prediction is transcript-centric (canonical
transcript reported); non-coding/regulatory annotation is limited to cache
content.

## 7. Structural variants — Manta (paired + tumor-only)

**What.** Deletion, duplication, inversion and insertion calling from paired
tumor/normal BAMs (`manta_paired`) and germline SV calling for normal tissue
(`manta_tumoronly`), with scoring of somatic candidates.

**Why.** Manta's local assembly + split-read pairing model is the most widely
validated WGS SV caller in clinical pipelines, with deterministic runtime and
built-in somatic scoring.

**Evidence.** Manta 1.6.0 — Chen et al. 2016, *Nat Methods* 13:49.

**Limitations.** Events shorter than ~50 bp are the SNV/indel caller's job;
very large (>100 kb) balanced events and complex rearrangements are detected
less completely; no mobile-element calling.

## 8. Copy-number segments — CNVkit (WGS mode)

**What.** Whole-genome binning and log2 copy-ratio segmentation
(`cnvkit_paired` with matched normal; `cnvkit_tumoronly` against a prebuilt
baseline) — `--method wgs`, no target BED in WGS mode.

**Why.** CNVkit provides robust, normalized per-segment log2 ratios and is
the community-standard hybrid/WGS CNV tool; in this pipeline it serves as the
*segment-level copy-ratio view*. Allele-specific purity/ploidy-aware
copy-number fitting is delegated to ASCAT (§13), which resolves total and
minor allele copy numbers per segment — information log2 ratios alone do not
carry.

**Evidence.** CNVkit 0.9.14 — Talevich et al. 2016, *PLoS Comput Biol*
12:e1004873.

**Limitations.** Log2 ratios are purity-blind: a 2:0 loss at 30% purity looks
like a 1:1 loss at 100%. Purity/ploidy interpretation therefore comes from
ASCAT; CNVkit segments are the sanity view. WGS bins without a matched
panel-of-normal baseline have slightly elevated noise in tumor-only mode.

## 9. Microsatellite instability — MSIsensor

**What.** Per-pair MSI scoring against the matched normal from paired WGS
(`msi_paired`), with a raw per-site table and a summary TSV; the cohort table
(`msi_cohort`) re-derives the unstable-site count with a proper
Benjamini–Hochberg FDR from the per-site output.

**Why.** MSIsensor is the original WGS-validated MSI caller (Niu et al. 2014),
using length distributions at homopolymer/microsatellite loci against a
matched normal. **Known bug handled explicitly:** MSIsensor v1.3.0 stores each
site's rank in a `uint16` while computing FDR, so above 65,535 sites the
reported unstable count saturates (observed as exactly 65,535). The pipeline
recomputes BH-FDR from the raw per-site table (`collect_msi.py`) and reports
*both* the vendor count and the corrected count side by side — for the
saturated sample the corrected count is higher (16.46% vs 10.54%), which
*strengthens* the MSI-H call.

**Evidence.** MSIsensor 1.3.0 — Niu et al. 2014, *Bioinformatics* 30:1017;
cutoffs — 3.5% unstable sites for WGS (Niu et al. 2014), ≥15% common for
WES/panel assays (msisensor-pro convention).

**Limitations.** Requires a matched normal; panel-of-normal-free; homopolymer-
heavy regions are noisier. Clinical MSI-H calls require orthogonal
confirmation (MMR immunohistochemistry or MLH1 promoter methylation).

## 10. Mutational signatures — SigProfiler SBS96 fitting

**What.** SBS96 mutational-matrix extraction from delivered SNVs and fitting
against reference COSMIC signature profiles.

**Why.** The SBS96 (Alexandrov-96) representation is the COSMIC standard;
fitting exposes APOBEC activity, aging signatures and therapy-related
signatures that contextualize TMB and MSI findings.

**Evidence.** Bergstrom et al. 2019, *Nature* 575:765; Alexandrov et al.
2013, *Nature* 500:415 (SBS96 scheme).

**Limitations.** Reliable fitting needs sufficiently many mutations; samples
with low TMB produce unstable signature attributions and are flagged rather
than over-interpreted.

## 11. Bulk RNA-seq — 2-pass STAR + featureCounts

**What.** Trimmed RNA-seq reads are aligned with two-pass STAR against the
same GRCh38 reference (shared index), then gene-level quantification with
featureCounts (`-t exon -g gene_id -p --countReadPairs`); QC metrics
(alignment rates, insert size, strand) are summarized per sample
(`rna_qc_summary`), and counts are merged into a cohort gene×sample matrix
with an added `gene_name` column (`merge_counts`).

**Why.** Two-pass STAR (pass 1 discovers novel junctions, pass 2 realigns
with them) is the ENCODE-standard alignment approach for human RNA-seq.
featureCounts is the fastest widely validated exon-level counter and pairs
naturally with STAR BAMs. Reporting both `gene_id` (Ensembl) and `gene_name`
(HGNC symbol) avoids the id-lookup ambiguity in downstream tools.

**Evidence.** STAR 2.7.11b — Dobin et al. 2013, *Genome Biol* 14:R95;
featureCounts 2.1.1 — Liao et al. 2014, *Bioinformatics* 30:923.

**Limitations.** No transcript-level quantification (gene-level only — use
Salmon/RSEM when isoform usage matters); multi-mapping at paralogous genes
follows featureCounts' default fractional/unique assignment rules; no
dedicated small-RNA handling.

## 12. Single-cell RNA-seq — dnbc4tools

**What.** Cell counting (`scrna_count`), per-sample QC/clustering
(`scrna_qc_cluster`) and cohort integration (`scrna_integrate`) for the
DNBelab C-series libraries, with cohort QC metrics aggregation
(`scrna_metrics_cohort`).

**Why.** dnbc4tools 3.0 is the vendor tool matching the exact library
chemistry (DNBelab C4); vendor chemistry-specific whitelist, correction and
UMI collapsing are not reliably reproducible by generic tooling.

**Evidence.** dnbc4tools 3.0 (BGI / DNBC4Tools).

**Limitations.** Closed-source vendor algorithms limit independent auditing;
standard outputs (matrix, metrics) are exported for cross-tool verification;
QC/clustering thresholds are library-complexity dependent.

## 13. Tumor purity, ploidy and allele-specific copy number — ASCAT (HTS chain)

**What.** For each of the 10 WGS tumor/normal pairs: allele counting at the
1000 Genomes hg38 SNP loci (`alleleCounter`), logR GC-content and
replication-timing correction, ASPCF segmentation, and the ASCAT
purity/ploidy grid fit producing per-segment total/minor-allele copy numbers,
aberrant cell fraction and goodness of fit.

**Why this tool, and why this exact chain.** ASCAT is the reference method for
joint tumor purity/ploidy estimation and allele-specific copy-number
profiling (Van Loo et al. 2010). The ASCAT project now ships an
NGS/WGS-optimized chain itself (Raine et al. 2023), which this pipeline
follows **exactly as documented by VanLoo-lab**: `prepareHTS` (allele
counting with alleleCounter, min base quality 20, min mapping quality 35) →
`loadData` → GC/RT `correctLogR` with the official G1000 hg38 WGS references
(Zenodo record 14008443) → `aspcf` (penalty 70) → `runAscat` with `gamma = 1`
(the HTS-tuned value; 0.55 is the SNP-array legacy). Orchestration and
per-pair output structure follow the `gcap::gcap.runASCAT` template
(Wang et al. 2024, *Nat Commun* 15:3231), re-implemented natively because
gcap itself is not OSI-licensed and no trusted public workflow (nf-core et
al.) ships ASCAT end-to-end. Autosomes 1:22 only: patient sex is unannotated
for this cohort, so X/Y are excluded from the fit (gender fixed to `"XX"`);
conclusions are unaffected for the arm-level CRC events reported.

**Evidence.** ASCAT 3.2.0 — Van Loo et al. 2010, *PNAS* 108:12313; NGS chain
— Raine et al. 2023, *NAR Genom Bioinform* 5:lqad041; alleleCounter 4.3.0 —
Raine et al. 2023; G1000 hg38 WGS references — Zenodo 14008443; orchestration
template — Wang et al. 2024, *Nat Commun* 15:3231 (gcap). Reference-file
provenance and the verified download/verification script: docs/refdata.md
§8 + scripts/download_ascat_refdata.sh (size+md5-gated; the resumed-transfer
corruption incident is documented in both).

**Limitations.** Purity resolution is ~5% and low-purity fits (<30%) have
reduced sensitivity for subclonal events; goodness-of-fit below ~0.80 flags a
less reliable solution; the grid fit assumes near-diploid ancestry; sex
chromosomes excluded here by design; homozygous-loss regions carry no
allele signal and are bounded by flanking segments.

## 14. Reporting — execution-derived methods, patient-level cohort table, clinical report

**What.** `methods_from_rule_runs` renders the Methods section **from the
actual execution checkpoint** (which rules ran, with which versions and
citations) rather than from prose that can drift; `cohort_tables` outer-joins
every module table to one patient-level summary (TMB + RNA QC + single-cell +
ASCAT + MSI); `clinical_report` renders the final Markdown+HTML report with
rule-based interpretation flags (TMB outliers under two conventions, indel
fraction → MSI hint, ASCAT purity/ploidy/GOF flags, MSI cutoffs under both
WGS and WES conventions).

**Why.** Execution-derived methods eliminate the classic failure mode where a
pipeline's written methods no longer match what ran. Patient-level joins make
cross-module interpretation (e.g. MSI-H + high TMB + APOBEC) explicit and
auditable; interpretation thresholds are stated inline with their citations
and convention-dependence.

**Limitations.** Interpretation is rule-based, not diagnostic; every clinical
call in the report carries its cutoff and requires orthogonal clinical
confirmation as noted inline.

### Tumor-only cohorts: paired-only analyses and report gating

ASCAT (§13) and MSIsensor (§9) are **paired-only analyses**: ASCAT's allelic
fit needs the matched normal to separate germline from somatic allele
ratios, and MSIsensor's baseline is built from the normal's length
distributions (a panel-of-normals baseline cannot be built from a single
tumor-only cohort). In a tumor-only cohort (declared with
`--arg tumor_only=true`) the per-pair rules (`ascat_paired`, `msi_paired`)
are never instantiated — their `when = "wildcard.control != ''"` gate skips
every instance — so the cohort aggregations would reference unproducible
inputs and fail, taking `cohort_tables` and `clinical_report` down with
them.

The pipeline handles this with a **header-only companion** pattern
(`ascat_cohort_tumoronly`, `msi_cohort_tumoronly`): the companion rules share
the output paths (`ascat/cohort_ascat.tsv`, `report/msi_cohort.tsv`) with
complementary `when` gates on `config.tumor_only` (the same
shared-output-path pattern as `cnvkit_paired` / `cnvkit_tumoronly`), and
write just the header row. Downstream rules keep their hard inputs, and the
report renders the ASCAT/MSI sections as "no data" — honest output (the
analyses genuinely did not run) instead of a broken run. Mixed cohorts
(some pairs paired, some tumor-only) are not resolved by this flag: the
aggregations fan out over `config.pairs_list` and still reference the
tumor-only pairs' unproducible inputs; engine support for a paired-only
pair list is tracked separately.

### Neoantigen module (P1): pVACfuse, TCR/HLA-RNA QC, HLA-LOH, MSI-H x frameshift

The P1 items (venus#31) extend the P0 chain (SpecHLA typing → VEP plugin
re-annotation → pVACseq/MHCflurry → burden report) with four opt-in layers
and one cohort page. Tool choices were verified against actual packages
(bioconda/PyPI/wheel inspection, not memory):

- **Fusion neoantigens (pVACfuse).** Manta's BND records are not in
  pVACfuse's input format; `scripts/manta_to_arriba.py` converts them to the
  Arriba TSV schema pVACfuse consumes (5'/3' transcript models, junction
  sequences fetched from the reference FASTA, PR/SR support). Paired cohorts
  feed the PASS-scored `deliver/{pair_id}.sv.vcf.gz`; tumor-only cohorts
  feed the candidateSV file (explicitly UNSCORED Manta calls) and the output
  filename says so (`fuse.candidate.filtered.tsv`). pVACfuse ships in the
  same pvactools env as pVACseq — no new env.
- **TCR repertoire (TRUST4).** TRUST4 (bioconda `trust4`, MIT) assembles
  TCRβ clonotypes from the bulk-RNA BAM. Opt-in (`run_tcr`): most venus
  cohorts are tumor-only organoids where bulk-RNA TCR adds little, and
  sensitivity is depth-dependent. Reference bring-up (bcrtcr.fa +
  human_IMGT+C.fa) is an ops step via `scripts/prepare_trust4_ref.sh`.
- **RNA-side HLA cross-check (arcasHLA).** arcasHLA (bioconda `arcashla`,
  GPL-3.0; repo now under RabadanLab) types HLA from the bulk-RNA BAM and
  `hla_dna_rna_qc` compares per-gene against the SpecHLA DNA typing
  (match / mismatch / absent_rna / absent_dna / not_run, order-insensitive
  2-field normalization). Opt-in (`run_hla_rna`); arcasHLA builds its `dat/`
  reference inside the env prefix, so `scripts/prepare_arcas_ref.sh` runs
  once per machine after the pixi install.
- **HLA-LOH — DASH infeasibility and the SpecHLA substitute.** DASH (Shukla
  et al. 2015) is the publication-standard HLA-LOH caller but has **no
  bioconda or PyPI package** (checked 2026-11; it is a standalone Perl/
  R code drop), so it cannot satisfy the pipeline's pixi lock-only env
  policy. The honest in-lock substitute is `spechla-loh`
  (SpecHLA 1.0.12's built-in `cal.hla.copy.pl`, shipped inside the existing
  envs/hla): per-allele copy ratios + LOH calls with a heterozygous-SNP
  floor (`-C`, default 5), purity/ploidy REQUIRED and fed from the ASCAT
  fit (which itself needs a matched normal — hence paired-cohort only and
  opt-in via `run_hla_loh`). Consequence: HLA-LOH calls here are copy-ratio
  inferences under ASCAT's purity/ploidy uncertainty, not DASH's
  germline-heterozygosity model; this disclosure travels with the rule
  (rules/hla_qc.oxoflow) and the artifacts doc.
- **scRNA corroboration.** %positive of each candidate gene in the pair's
  own filter_feature.h5ad, reusing envs/scrna (the only env with h5py —
  no new env, no new lock gap). Honest caveat (embedded in the report):
  the dnbc4tools h5ad carries no malignant-cell annotation, so %positive
  is over ALL filtered cells — a coarse "is the gene expressed at all"
  signal, never tumor-cell prevalence.
- **MSI-H x frameshift cohort page.** `msi_frameshift_page` joins the
  cohort-wide `frameshift_neoantigens` counts (frameshift-derived among
  strong binders, from `report/neoantigen_cohort.tsv`) with the MSI
  statuses from `report/msi_cohort.tsv` — `msi_status_wgs` at the 3.5%
  BH-corrected unstable-site cutoff (Niu et al. 2014) and
  `msi_status_wes` at the 15% WES/panel convention, shown side by side
  and never collapsed. Because the two sides live in opposite cohort
  modes (MSI is paired-only, the neoantigen module is tumor-only-gated),
  every joined cell is either a value or the explicit sentinel
  `not_available` — in a tumor-only cohort all MSI cells are honestly
  not-available (header-only companion table), and the markdown page
  states which side carries data instead of implying a join that never
  happened. A frameshift count is a lower bound: it covers only
  candidates that passed pVACseq binding prediction.

**Why.** Each addition is either license-clean in-lock (MHCflurry, TRUST4,
arcasHLA, SpecHLA tools) or explicitly disclosed as operator-side
(pVACview's R Shiny stack); the MSI-H x frameshift join makes the
MSI-H → frameshift-indel neoantigen expectation explicit instead of
requiring a manual cross-table read.

**Limitations.** The neoantigen module is tumor-only-gated (P0 scope), so
the frameshift x MSI join has both sides populated only if a future
paired-cohort extension lifts that gate; until then one side of every row
is `not_available` by construction.

---

## Reference list

- Alexandrov et al. 2013, *Nature* 500:415 — SBS96 scheme.
- Bergstrom et al. 2019, *Nature* 575:765 — mutational signature fitting.
- Chen et al. 2016, *Nat Methods* 13:49 — Manta.
- Chen et al. 2018, *Nat Biotechnol* 36:4875 — fastp.
- Cibulskis et al. 2013, *Nature* 505:183 — Mutect2.
- Danecek et al. 2021, *Gigascience* 10:giab008 — SAMtools.
- DePristo et al. 2011, *Genome Res* 21:436 — GATK/BQSR.
- Dobin et al. 2013, *Genome Biol* 14:R95 — STAR.
- Ewels et al. 2016, *Bioinformatics* 32:3671 — MultiQC.
- Hurt et al. 2019, *Genome Med* 11:78 — pVACtools (pVACseq/pVACfuse
  prediction and the pVACview review app).
- Liao et al. 2014, *Bioinformatics* 30:923 — featureCounts.
- Mayakonda et al. 2018, *Leukemia* 32:1152 — vcf2maf.
- McKenna et al. 2010, *Genome Res* 20:1297 — GATK.
- McLaren et al. 2016, *Genome Biol* 17:122 — VEP.
- Niu et al. 2014, *Bioinformatics* 30:1017 — MSIsensor; *Sci Rep*
  4:4731 (erratum/addendum) — the 3.5% WGS MSI-H cutoff used by
  collect_msi.py (the 15% WES/panel convention is carried alongside).
- Poplin et al. 2018, *Nat Biotechnol* 36:875 — Mutect2/DeepVariant-grade
  filtering.
- Raine et al. 2023, *NAR Genom Bioinform* 5:lqad041 — ASCAT NGS chain +
  alleleCounter.
- Shukla et al. 2015, *Genome Med* 7:63 — DASH (HLA-LOH). NOT deployed:
  no bioconda/PyPI package exists (checked 2026-11), so the pixi
  lock-only env policy excludes it; spechla-loh substitutes (see the
  neoantigen section above).
- SpecHLA — Luo et al. 2024, *Nat Biotechnol* 42:1318 (SpecHLA HLA typing
  and the built-in cal.hla.copy.pl LOH step).
- Talevich et al. 2016, *PLoS Comput Biol* 12:e1004873 — CNVkit.
- TRUST4 — Song et al. 2021, *Nat Commun* 12:5710 — TCR/BCR assembly from
  RNA-seq.
- arcasHLA — O'Donnell et al. 2020, *Bioinformatics* 36:4192 — RNA-side
  HLA typing for the DNA/RNA cross-check.
- Van Loo et al. 2010, *PNAS* 108:12313 — ASCAT.
- Vasimuddin et al. 2019, *IEEE ACS* — BWA-MEM2.
- Wang et al. 2024, *Nat Commun* 15:3231 — gcap ASCAT-NGS orchestration.
- Wells et al. 2020, *Immunity* 53:932 — TESLA consortium best practices
  for neoantigen prediction pipelines (binding-threshold strategy).
