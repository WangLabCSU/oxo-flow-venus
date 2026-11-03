# Full-repo audit — 2026-11-03

Scope: every rule file (`rules/*.oxoflow`, 16 files / 59 rules), `venus.oxoflow`,
all 25 scripts under `scripts/` (14 Python rule-consumed helpers + R driver +
6 ops shell scripts), all 17 `envs/*/pixi.toml`, tests, CI, README, docs and
CHANGELOG, audited against oxo-flow 0.24.x semantics (when-gates on shared
output paths, `expand_inputs`, fingerprint = rule content + input manifest,
env binding from repo checkout). Engine-internal behavior (#931 tombstone
regression, #893 semantics) is out of scope per the audit charter. The
clinical module (docs/20261027-clinical-production-audit.md) and the
neoantigen module (@326c1a4/9cd32dd/d4390ee/52cf913) were checked for
regressions only.

## Executive summary

The pipeline is architecturally sound: the when-gate discipline is complete
(all 15 shared output paths are covered by exactly one complementary rule
pair, verified programmatically), the config surface has no ghost keys, and
the delivered-analysis chain (Mutect2 scatter → FilterMutectCalls → VEP →
MAF → TMB/signature/clinical) is internally consistent end to end, with tool
parameters matching GATK/CNVkit/ASCAT/MSIsensor documented practice. The
audit found **no blocker** for the default paired-cohort run, **three
reliability defects that were fixed in place** (a checkpoint-reading race
that can nondeterministically drop the clinical module from the report's
methods section; a Manta retry deadlock on the paired path; an undeclared
ecDNA gene-list dependency that could let stale data flow into the cohort
table), plus one **documented-but-unrunnable config trap** (`run_bqsr=false`,
now fenced with warnings). The remaining findings are improvement items:
four environments violate the repo's own lock-only policy (missing
`pixi.lock`), QC completeness gaps vs. nf-core/sarek-class pipelines
(coverage/insert-size metrics, sex check, contamination, RNA fusions,
scRNA ambient/doublet), thin test coverage outside the clinical/neoantigen
modules, and a tumor-type inconsistency ("CRC" vs gastric/GEJ) that runs
through the cohort-narrative documents and reaches the delivered clinical
report — that one needs an owner decision, not a code change.

## Findings

Levels: **blocker** = breaks a supported run mode; **major** = wrong or
nondeterministic results / reproducibility hazard; **minor** = doc drift,
dead weight, polish. "Fixed" = changed in this audit's commit; "issue" =
recommended follow-up (not fixed here, with the reason).

### Major

1. **Methods section race on the clinical chain** — `rules/report.oxoflow`
   `methods_from_rule_runs` read `.oxo-flow/checkpoint.json` at execution
   time but declared only `cohort_summary.tsv` as ordering input. The
   clinical chain (`clinical_kb_prep` → `clinical_snv_annotate` → … →
   `clinical_cohort`) is *not* upstream of `cohort_summary`, so the methods
   rule could snapshot the checkpoint before the clinical rules finished and
   render their rows as "no checkpoint record" — nondeterministically, in the
   headline clinical deliverable. **Fixed**: `report/clinical_cohort.tsv` +
   `report/clinical_therapies_cohort.tsv` added as declared inputs (producers
   are ungated, no cycle). Neoantigen outputs were deliberately *not* added:
   their producer is gated `tumor_only == 'true'`, so an unconditional input
   would break paired-cohort plans (needs a gated companion — see issue N2).

2. **`manta_paired` retry deadlock** — `configManta.py` refuses an existing
   runDir; the paired rule did not clear it before configuring, so any retry
   after a failed `runWorkflow.py` failed forever on the leftover directory.
   The tumor-only sibling already carried the `rm -rf` guard ("deterministic
   retries" per its own comment) — the asymmetry was an omission, not a
   design choice. **Fixed**: `rm -rf manta/{pair_id} &&` added, mirroring
   `manta_tumoronly`. (Stales the 10 completed `manta_paired` instances on
   workdir sync — merge in a planned re-run window, same convention the
   repo documents for `star_align`.)

3. **`ecdna_cohort` undeclared input** — `collect_ecdna.py` reads
   `ecdna/<pair>/<pair>_gene_list.tsv` (the ecDNA-gene cross-flag source for
   `clinical_cnv_genes`) via a side path, but the rule only declared
   `_ecDNA_counts.tsv` + `_result_table.tsv`. A suite re-run that changed the
   gene list could not invalidate the cohort table. **Fixed**: third
   `expand_inputs` pattern added (producer `ampliconsuite` is ungated, so no
   mode-dependence).

4. **`run_bqsr=false` is a dead end advertised as an example** — `bqsr` is
   gated `config.run_bqsr == 'true' && wildcard.assay == 'wgs'` with no
   companion rule, while every consumer (Mutect2, Manta, CNVkit, ASCAT,
   MSIsensor, AmpliconSuite, SpecHLA) hard-wires
   `bqsr/{sample}/{sample}.bqsr.bam(.bai)`. The README's override example
   literally suggested `--arg run_bqsr=false`, which yields a DAG with no
   producer for the entire somatic chain. **Fixed (guardrail)**: warning in
   `venus.oxoflow`, README example replaced with safe keys + explicit note.
   **Issue N1**: a real fix needs a no-BQSR companion rule publishing the
   markdup BAM under the `bqsr/` path.

5. **Four environments violate the lock-only policy** —
   `envs/{hla,idx,pvac,vep_perl}` ship `pixi.toml` without `pixi.lock`,
   while README/#26 declare "pixi.lock is the single source of truth" and
   engine env fingerprints bind to the checkout. `bootstrap_envs.sh` would
   solve these at install time (host-dependent versions). **Issue N3** —
   not fixed here on purpose: locks must be generated from the server where
   the current envs are materialized (a Mac-side solve risks shipping a lock
   that diverges from the live envs and mass-invalidates fingerprints).

6. **Tumor-type inconsistency in cohort narrative docs** —
   `notes/interpretation.md` ("10 CRC patients", "aged-colon/SBS1-like"),
   `docs/artifacts.md` ("10 CRC pairs") and `docs/method-rationale.md` §13
   ("arm-level CRC events") describe the cohort as colorectal, while the
   clinical panels (`actionable_cnv_genes.tsv` "gastric/GEJ",
   `actionable_fusion_pairs.tsv` CLDN18 fusions, NMPA table zolbetuximab/
   disitamab) and `docs/refdata.md` are gastric/GEJ throughout.
   `notes/interpretation.md` is a rule input of `clinical_report`, so the
   **delivered** report inherits whatever the label should be. **Issue N4** —
   deliberately not auto-fixed: adjudicating the cohort's tumor type is an
   owner decision (both cancers can be MSI-H; the docs disagree, not the
   data).

### Minor (fixed in this audit)

7. `docs/artifacts.md` coverage statement claimed "13 rule files (40
   rules)"; the manifest is 16 files / 59 rules, and the ecDNA/neoantigen
   modules plus `ecdna_cohort` / `msi_cohort_tumoronly` were absent from the
   inventory. **Fixed** (statement rewritten; module tables added; the
   561-instance census explicitly marked as predating the newer modules).

8. README drift: module list missing ascat/ecDNA/clinical/neoantigen;
   scanpy credited as a pipeline tool (no rule or script imports it — the
   scRNA summaries are deliberately dependency-light); engine requirement
   "≥ 0.23" while the armed `[engine]` disk guard needs ≥ 0.24 + commit
   `a456c796`. **Fixed.**

9. `docs/refdata.md` §9 named rule `msi_scan` (real name
   `msisensor_scan`); the AmpliconSuite `aa_data_repo` — a ~27 GB host
   reference with no §-entry despite verify_refdata gating on it — and the
   Perl-VEP cache/plugin deployment were undocumented. **Fixed** (§12 added,
   §11 extended, asset table rows added).

10. `scripts/cohort_tables.py` docstring omitted the implemented
    `cohort_ecdna.tsv` route; `scripts/tmb_from_maf.py` docstring claimed a
    "'Precision'-qualified rows" drop that the code never does. **Fixed.**

11. Dead config key `workflow_file` (defined in `venus.oxoflow`, referenced
    by nothing — rules, scripts or ops tooling). **Removed.** (`target_bed`
    is also reference-free by rules but is the *documented anchor* of the
    WES/panel key set — kept, as the config comment itself explains.)

### Minor (recommended issues, not fixed)

12. **VEP version triplet** — repo says binary vep-rs 0.3.1 + cache v113
    (venus.oxoflow, refdata §5, method-rationale §6), but
    `methods_from_rule_runs.VERIFIED_VERSIONS` lists "vep-rs 115.2" with
    provenance "run log banner". Either the banner carries a VEP-compat
    release number (fine, but the row mislabels it as the tool version) or
    one of the v113 claims is stale. Needs a host-side check of the actual
    banner before editing any number. **Issue N5.**

13. **Env dead weight** — `envs/report` (matplotlib/jinja2/pyyaml/openpyxl/
    pysam unused by any report rule), `envs/clinical` (`requests` unused —
    clinical_oncokb uses urllib), `envs/rna` (`fastp` unused since trimming
    lives in qc), `envs/scrna` (scanpy/leidenalg/igraph/scrublet unused by
    rules; raw matrices are kept *for* SoupX/DecontX/scDblFinder, so the
    deps are a deliberate future hook — but they are undocumented as such).
    Not fixed: any pixi.toml edit requires a lock regen (see N3).
    **Issue N6.**

14. **`pvacseq_run` single-predictor assumption** —
    `cp .../MHC_Class_I/*.filtered.tsv {output}` breaks when
    `neo_predictors` carries more than one predictor (multiple .filtered.tsv
    sources into a file target). Safe at P0 (MHCflurry only); add a guard or
    a predictor-scoped filename before DTU predictors land. **Issue N7.**

15. **Neoantigen module is not integrated into the cohort report** —
    `neoantigen/{pair_id}.{summary,report.md}` are standalone; neither
    `cohort_tables` nor `clinical_report` consumes them, and the methods
    section can miss their checkpoint records (same race class as finding 1,
    unfixable by unconditional inputs because of the tumor_only gate).
    **Issue N2.**

16. **AmpliconSuite retry semantics unverified** — the rule asserts the
    suite's finish flag (good), but whether
    `AmpliconSuite-pipeline.py --output_directory <existing dir>` resumes,
    fails, or silently reuses stale stage outputs after a partial failure
    was not verified in this audit (MSI got an explicit scratch-dir pattern
    for exactly this class of problem). **Issue N8.**

17. **Test coverage gaps** — only clinical_tier/clinical_interpret/
    neoantigen scripts have tests. Untested and testable offline:
    `vcf2maf_light` (CSQ contract, indel trimming/coordinates),
    `tmb_from_maf` (BED merge, dual denominators), `sbs96_fit` (Fasta.fai
    reader, pyrimidine normalization, NNLS), `collect_msi` (BH math +
    uint16 saturation), `collect_ecdna`, `collect_ascat`, `cohort_tables`
    (routing/patient normalization), `merge_counts`, `summarize_star_logs`,
    `scrna_metrics` (alias/thousands-separator handling), `methods_from_rule_runs`
    (attribution longest-first, bash -c unwrap), `clinical_report`
    (md_to_html escaping), `clinical_cnv`, `clinical_sv`, `str_context`,
    `annotate_cnv_reference`. **Issue N9** (highest value: vcf2maf_light +
    collect_msi, both feed clinical decisions).

18. **QC completeness vs. standard production practice** (all additive,
    each its own issue; none is a correctness bug):
    - No per-sample alignment/coverage QC beyond dup metrics: flagstat,
      CollectInsertSizeMetrics, CollectWgsMetrics/mosdepth coverage+insert
      distribution are standard in sarek/GATK best-practice pipelines and
      feed MultiQC today only via dup metrics (the MultiQC shell already
      scans `align/`, and its comment mentions the samtools flagstat
      matcher — but no rule produces flagstat). **Issue N10.**
    - No sample-sex check (ASCAT excludes X/Y citing unannotated sex — the
      gap that motivates it is unmonitored) and no cross-contamination
      estimate (VerifyBamID2/Conpair). **Issue N11.**
    - No RNA fusion detection: `actionable_fusion_pairs.tsv` is screened
      only against Manta BND breakpoints; STAR chimeric junctions
      (`--chimSegmentMin`) or an RNA-seq caller would materially improve
      CLDN18/ARHGAP-class recall for this cohort's own panel. **Issue N12.**
    - No scRNA ambient/doublet correction despite protecting `raw_matrix/`
      explicitly "for SoupX/DecontX/scDblFinder". **Issue N13.**
    - mtDNA deliberately excluded (documented in venus.oxoflow) — no action.
    - `featureCounts` strandedness is a fixed unstranded default verified
      once for this library type; a `config.featurecounts_strandness`
      passthrough would prevent silent mis-counting on a stranded cohort.
      **Issue N14.**

## Verified-correct highlights (checked, no action)

- **When-gate completeness**: all 15 shared output paths have exactly one
  complementary producer pair (`mutect2_{paired,tumoronly}`,
  `cnvkit_{paired,tumoronly}`, `ascat_cohort(_tumoronly)`,
  `msi_cohort(_tumoronly)`, `clinical_oncokb(_skip)`,
  `clinical_sv_fusions_{paired,tumoronly}`,
  `clinical_biomarkers_{paired,tumoronly}`) — verified by script over the
  parsed manifests; no path has a single gated producer without a companion
  except `bqsr` (finding 4) and `hla_type_dna`/neoantigen rules (gated only
  by `tumor_only`/`assay`, consumed only under the same gates).
- **Config keys**: no ghost keys; the only rule-referenced keys absent from
  `[config]` are the engine-injected `pairs_list` / `samples_<group>`
  (documented in README). One dead key removed (finding 11).
- **expand_inputs coverage**: every directory the `multiqc` shell scans
  (`qc/trim qc/fastqc scrna/fastqc align rna/star rna/counts`) has a
  matching anchor pattern, so the protected report cannot freeze with
  missing sections; all cohort aggregators expand over producers that exist
  in both cohort modes.
- **MSI output contract**: `mv msi/{pair_id}` was confirmed against the
  msisensor-pro source (`sample.cpp` writes a bare-prefix summary
  `Total_Number_of_Sites/Number_of_Unstable_Sites/%` next to
  `_all/_unstable/_dis`) — the suspicious-looking move is correct, and
  `collect_msi.read_summary` consumes exactly that schema; the BH step-up
  recomputation is standard.
- **ASCAT chain**: `run_ascat.R`'s output filenames match the declared rule
  outputs one-for-one (including the `plots/{job}.{tumour}.…` prefix math
  and the loud `expected outputs missing` guard); parameters (penalty 70,
  gamma 1, minBQ 20/minMQ 35, G1000 hg38 plain-numbered loci, GC/RT
  correction, autosomes-only) match VanLoo-lab HTS documentation and the
  gcap orchestration the module cites.
- **Mutect2 scatter/merge chain**: GatherVcfs order is guaranteed by
  `config.chromosomes` (chr1…chrY, dict order); MergeMutectStats /
  LearnReadOrientationModel / FilterMutectCalls flag semantics and the
  `.vcf.gz.filteringStats.tsv` default path match GATK 4.x (and are
  production-verified per docs/artifacts.md); PoN is applied in tumor-only
  mode too (0.4.0 fix still in place — regression check passed).
- **CSQ contracts**: `vcf2maf_light` and `clinical_tier` both derive CSQ
  subfield order from the header regex `Format[:=]` — consistent with
  vep-rs's `Description="… Format: A|B|…"` spelling; allele trimming and
  MAF coordinates follow the mskcc convention; compound terms and
  frameshift typing resolve deterministically.
- **TMB conventions**: coding-Mb (34) vs target-Mb (3000) denominators and
  the FDA ≥13 / legacy ≥10 cutoffs are consistent across rule, config
  comment, script defaults, clinical_interpret and clinical_report.
- **SBS96**: channel enumeration, pyrimidine normalization, .fai random
  access and NNLS exposure fractions are internally consistent with the
  COSMIC v3.1 matrix format (`Type`-keyed rows, columns sum to 1).
- **Neoantigen module (regression check)**: post-audit fixes all present —
  plugin re-annotation feeds pvacseq, sample_name = `{experiment}`,
  G-group wide-table parsing, pVACtools 7.x column aliases, exploratory
  percentile strategy, VAF/depth filters off, MHCFLURRY_DATA_DIR wired,
  `gene_symbol_map` awk substr offsets verified correct, tests green.
- **Clinical module (regression check)**: no behavioral drift since
  docs/20261027-clinical-production-audit.md — M1 scalar-CLNSIG guard, M4
  G-group, M6 50% interval-overlap guard, Predictive-only tiering and the
  biomarker-keyed therapy matching are all intact; the only clinical-rule
  edits in this audit are the two report-side input declarations.

## Fixes implemented in this audit (working tree → `audit/20261103`)

| File | Change |
|---|---|
| `rules/report.oxoflow` | methods_from_rule_runs: + clinical cohort tables as inputs (finding 1); ecdna_cohort: + gene_list expand pattern (finding 3) |
| `rules/sv.oxoflow` | manta_paired: + `rm -rf` retry guard (finding 2) |
| `venus.oxoflow` | run_bqsr=false trap documented; dead key `workflow_file` removed; header scanpy mention corrected |
| `README.md` | module list completed; scanpy attribution dropped; engine ≥ 0.24 + a456c796; safe override example + run_bqsr note |
| `docs/refdata.md` | `msisensor_scan` name fix; §12 (aa_data_repo); §11 Perl-VEP cache/plugins; asset-table rows |
| `docs/artifacts.md` | coverage statement updated to 16 files / 59 rules; ecDNA + neoantigen module tables; missing report rows |
| `scripts/cohort_tables.py` | docstring: ecdna route documented |
| `scripts/tmb_from_maf.py` | docstring: unimplemented "Precision" claim removed |
| `CHANGELOG.md` | audit entry |

Tests: `tests/test_clinical_tier.py` (24), `tests/test_clinical_interpret.py`
(8), `tests/test_neoantigen.py` (8) — all green before and after (run with
`/tmp/venus-t-venv/bin/python` + `PYTHONNOUSERSITE=1`; note the venv's pysam
0.24.1 install was broken/incomplete and had to be reinstalled to run
test_clinical_tier). TOML validity re-checked for all 17 workflow files;
`py_compile` clean for all scripts.

Deployment note: findings 1–3 change rule inputs/shell templates → they
stale `manta_paired`, the report chain (`methods_from_rule_runs`,
`ecdna_cohort`, `cohort_tables`, `clinical_report`) and their downstream on
workdir sync. Sync together with a planned re-run window, per the repo's own
fingerprint convention.

## Recommended issue list

- **N1** (major): no-BQSR companion rule (`run_bqsr=false` end-state) or
  removal of the config key until it exists.
- **N2** (major): integrate neoantigen outputs into the cohort report +
  methods section (gated input or cohort collector).
- **N3** (major): generate `pixi.lock` for `envs/{hla,idx,pvac,vep_perl}`
  from the server (`pixi lock` per manifest; verify `bootstrap_envs.sh`
  then no-ops), restoring the lock-only policy.
- **N4** (major, owner decision): adjudicate cohort tumor type (CRC vs
  gastric/GEJ) and align `notes/interpretation.md` / `docs/artifacts.md` /
  `docs/method-rationale.md` §13 — the delivered report inherits the notes.
- **N5** (minor): reconcile the VEP version triplet (0.3.1 binary / v113
  cache / "115.2" banner) against the host banner; then fix whichever doc
  or VERIFIED_VERSIONS row is wrong.
- **N6** (minor): env slimming or explicit justification comments for unused
  deps (report/clinical/rna/scrna) — bundle with N3 lock regen.
- **N7** (minor): pvacseq multi-predictor output guard.
- **N8** (minor): verify AmpliconSuite resume/existing-outdir behavior;
  adopt the MSI scratch-dir pattern if it does not resume cleanly.
- **N9** (major, cheap): add offline unit tests for `vcf2maf_light`,
  `collect_msi`, `sbs96_fit`, `cohort_tables`, `tmb_from_maf` (fixture
  style of tests/test_neoantigen.py).
- **N10** (minor): alignment QC depth — flagstat + CollectInsertSizeMetrics
  + coverage (mosdepth/CollectWgsMetrics) into MultiQC's already-scanned
  `align/` root.
- **N11** (minor): sample-sex concordance + cross-contamination check
  (VerifyBamID2/Conpair) — would also unblock ASCAT sex-chromosome handling.
- **N12** (minor): RNA fusion support (STAR `--chimSegmentMin` + junction
  screen against `actionable_fusion_pairs.tsv`).
- **N13** (minor): scRNA ambient/doublet QC rules consuming the protected
  `raw_matrix/`.
- **N14** (minor): `featurecounts` strandedness passthrough config key.
