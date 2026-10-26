# Reference data — provenance, download and verification

Every external reference file the venus pipeline consumes: where it comes
from, which version, how it was built, how to re-download/rebuild it, and how
it is verified. **Rule:** if a step needs data that is not in the repository,
this page must say exactly how to obtain it — a pipeline whose reference
provenance lives only in someone's shell history is not reproducible.

Two categories exist by design:

1. **Repo-tracked** (`resources/`, committed) — small, stable, versioned with
   the pipeline. Currently: `resources/cosmic_v3.1_sbs96_grch38.txt` (COSMIC
   v3.1 SBS96 signature matrix, consumed by `signature_fit`).
2. **Host-path references** (absolute paths in `venus.oxoflow` config) — too
   large or too site-specific to commit. Each is listed below with origin and
   verification. Defaults assume the `lab-wsx-in` host layout
   (`/home/data/reference/`); every entry includes the rebuild/re-download
   command so another site can reconstruct them.

Host reference root: `/home/data/reference/`. Since 2026-10-11 every asset
below lives physically under this single shared root. Two config strings still
carry historical locations outside it (`scrna_reference`/`scrna_bin` under
`/home/data/data-resource/dnbc4tools`, `ascat_refdir` under
`/home/data/wsx/hnzl_data/refdata/ascat`); those paths are now **compat
symlinks** into the root. The strings stay byte-identical on purpose — rule
commands embed them, so changing them would stale every fingerprint that
touches the affected references. Flip the strings to the physical paths only
with a deliberate staleness review. Verify the whole set on the host with:

```bash
scripts/verify_refdata.sh    # from the repo root, on the reference host
```

| Asset | Config key(s) | Version | Origin |
|---|---|---|---|
| GRCh38 FASTA + BWA-MEM2 index | `genome_fa/fai/dict` | GATK bundle d1.vd1 | see §1 |
| GATK known-sites resources | `dbsnp/mills/gnomad/pon` | hg38 bundle | see §2 |
| STAR index | `star_index` | GRCh38 + GENCODE v47 | see §3 |
| GENCODE annotation | `gencode_gtf` | v47 (primary assembly) | see §4 |
| VEP cache | `vep_cache_dir` | Ensembl 113, GRCh38 | see §5 |
| dnbc4tools reference | `scrna_reference` | vendor bundle | see §6 |
| CNVkit tumor-only reference | `cnv_reference` | not deployed here | see §7 |
| ASCAT G1000 hg38 WGS set | `ascat_refdir` | ASCAT 3.2 / Zenodo 14008443 | see §8 |
| MSIsensor reference list | (derived) | generated per run | see §9 |
| Clinical KB (CIViC dump + ClinVar VCF) | `civic_*_tsv`, `clinvar_vcf` | CIViC monthly + ClinVar weekly | see §10 |

---

## 1. Genome FASTA + BWA-MEM2 index — `GRCh38.d1.vd1`

- **Path:** `/home/data/reference/GRCh38.d1.vd1.fa` (+ `.fai`, `.dict`)
- **Origin:** GATK resource bundle `GRCh38d1` (b37→hg38 decoy-aware build used
  by GATK and by the Sarek/nf-core ecosystem),
  `https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta`.
- **BWA-MEM2 index:** co-located with the FASTA
  (`GRCh38.d1.vd1.fa.{0123,amb,ann,bwt.2bit.64,pac}`), prefix = genome_fa.
  `bwa-mem2 index` is deterministic; rebuild with
  `bwa-mem2 index GRCh38.d1.vd1.fa` if missing.
- **Verify:** `samtools faidx` + `gatk CreateSequenceDictionary` regenerate
  the sidecars byte-comparably; index presence is asserted by `bwa-mem2 mem`
  failing fast otherwise.

## 2. GATK known-sites / germline resources (hg38)

All under `/home/data/reference/gatk-resources/hg38/`, from the
[Broad resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811),
each with `.tbi`:

| File | Used by |
|---|---|
| `Homo_sapiens_assembly38.dbsnp138.vcf.gz` | BQSR known-sites (SNV mask) |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | BQSR known-sites (indel mask) |
| `af-only-gnomad.hg38.vcf.gz` | Mutect2 germline resource |
| `1000g_pon.hg38.vcf.gz` | Mutect2 panel of normals (tumor-only filtering) |

- **Verify:** `bcftools index -t` succeeds; BQSR/Mutect2 fail fast on
  contig-name mismatch, which cross-checks the build.

## 3. STAR index — `star_index_GRCh38_d1vd1`

- **Path:** `/home/data/reference/star_index_GRCh38_d1vd1/`
- **Origin:** built locally from §1 FASTA + §4 GENCODE v47 GTF:
  ```bash
  STAR --runMode genomeGenerate --runThreadN 24 \
    --genomeDir star_index_GRCh38_d1vd1 \
    --genomeFastaFiles GRCh38.d1.vd1.fa \
    --sjdbGTFfile gencode.v47.primary_assembly.annotation.gtf \
    --sjdbOverhang 149        # readLength - 1 (150 bp runs)
  ```
  If read length changes, `--sjdbOverhang` must be `readLength − 1` and the
  index rebuilt (2-pass STAR's junction-insertion pass depends on it).
- **Verify:** `STAR --version` + `cat star_index_GRCh38_d1vd1/genomeParameters.txt`
  records the exact build command, genome and GTF — that file IS the
  provenance record.

## 4. GENCODE v47 annotation

- **Path:** `/home/data/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf`
- **Origin:** `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz`
  (gunzip in place). Used by STAR index build, `featureCounts -a` and `rna_qc --gtf`.
- **Verify:** `zcat | md5sum` against the release-page checksum; GENCODE
  publishes per-file md5 sums.

## 5. VEP cache — Ensembl 113 (GRCh38)

- **Path:** `/home/data/reference/vep_cache/` (`info.json` records
  `cache_version: 113`, `assembly: GRCh38`; transcripts under `transcripts/`)
- **Binary:** `/home/data/reference/vep-rs/vep-0.3.1-x86_64-unknown-linux-gnu/vep`
  (vep-rs 0.3.1, the Rust VEP re-implementation).
- **Origin:** standard Ensembl VEP cache installer content
  (`https://ftp.ensembl.org/pub/release-113/variation/vep/homo_sapiens_vep_113_GRCh38.tar.gz`),
  ingested into the vep-rs layout. The cache version must match what the
  binary supports; `info.json` is the local provenance record.
- **Verify:** annotating a dbSNP-known variant must return its rsID; cache
  version is echoed in the annotation log.

## 6. dnbc4tools single-cell reference — `ref-Homo_sapiens`

- **Path (physical):** `/home/data/reference/dnbc4tools/ref-Homo_sapiens/`
  (35 GB; sibling `/home/data/reference/dnbc4tools/dnbc4tools3.0/` holds the
  binary — config keys `scrna_reference` / `scrna_bin`). The config strings
  still read `/home/data/data-resource/dnbc4tools/...` and resolve through the
  compat symlink
  `/home/data/data-resource/dnbc4tools -> /home/data/reference/dnbc4tools`
  (layout unified 2026-10-11; see the note at the top of this page).
- **Origin:** BGI **vendor bundle** shipped with dnbc4tools 3.0 for DNBelab
  C-series chemistry (prebuilt genomeDir + whitelists). Vendor-matched
  chemistry data cannot be reconstructed generically; use the bundle the
  vendor distributes for the exact kit version.
- **Verify:** `dnbc4tools` validates the reference at startup and reports its
  identity in the count log; `scripts/verify_refdata.sh` checks the directory
  is non-empty and the binary executable *through the configured string*, so
  a broken symlink is caught before a run starts.

## 7. CNVkit tumor-only reference — configured, NOT deployed on this host

- **Config:** `cnv_reference = /home/data/reference/cnvkit_reference/hg38.wgs.reference.cnn`
- **Status:** the directory is **empty** on this host. Impact: **none for the
  paired run** — `cnvkit_paired` builds a per-pair reference from the matched
  normal and never reads this file. `cnvkit_tumoronly` WOULD fail; if
  tumor-only CNV is ever needed, first build the reference from a normal
  cohort:
  ```bash
  cnvkit.py batch tumor.sample.bam -n normals/*.bam -m reference \
    -f GRCh38.d1.vd1.fa --method wgs -d cnvkit_reference/
  ```
- **Verify:** `test -s "$cnv_reference"` before enabling any tumor-only CNV rule.

- **Gene annotation (`cnv_annotate_args`):** reference `.cnn` files built
  without `--annotate` carry an empty gene column (`-`), and `cnvkit batch`
  against such a reference propagates the gap to every `.cnr/.cns` — report
  tables then need coordinate fallbacks to name a locus (seen on GC01 with
  ERBB2). The fix is asymmetric (verified against cnvkit 0.9.14: `batch
  --reference` **hard-rejects** `--annotate/-n/-t/-a/...` as
  reference-construction flags):
  1. **Paired mode** — the reference is built per-run by `cnvkit_paired`, so
     set `cnv_annotate_args = --annotate <gene.bed>` and the flag flows into
     the batch call that constructs the reference; the annotated `.cnn` then
     names every `.cnr/.cns` bin it produces.
  2. **Tumor-only mode** — the reference is a prebuilt external `.cnn`, and
     batch cannot re-annotate it. Gene names come solely from that file's
     gene column (extracted regions → target bed → coverage → fix keep it).
     So annotate the reference itself: either rebuild it with
     `cnvkit.py reference ... --annotate <gene.bed>`, or re-annotate the
     existing file in place with
     `scripts/annotate_cnv_reference.py <ref.cnn> <gene.bed> [-o out.cnn]`
     (same `into_ranges` call cnvkit's own `target --annotate` uses; verified
     bin-for-bin identical output, ERBB2 bins named on the chr17 test slice).
     Do NOT pass `cnv_annotate_args` in tumor-only runs — the batch call
     would exit with a flag-conflict error.
  Adding the flag (or swapping in an annotated reference) changes rule
  fingerprints, so fold the switch into a planned rerun window.

## 8. ASCAT G1000 hg38 WGS reference set (Zenodo 14008443)

- **Path (physical):** `/home/data/reference/ascat/v3.2-hg38-wgs/` (config
  key `ascat_refdir`; the configured string
  `/home/data/wsx/hnzl_data/refdata/ascat/v3.2-hg38-wgs` resolves through the
  compat symlink
  `/home/data/wsx/hnzl_data/refdata/ascat -> /home/data/reference/ascat`,
  layout unified 2026-10-11 — see the note at the top of this page).
- **Origin:** official VanLoo-lab ASCAT NGS references,
  `https://zenodo.org/records/14008443` — allele/loci files for
  alleleCounter, GC-content and replication-timing tables for logR
  correction (see docs/method-rationale.md §13).
- **Download + verify (automated):**
  ```bash
  scripts/download_ascat_refdata.sh -o /home/data/reference/ascat/v3.2-hg38-wgs
  # flaky link: add -r (resume loop); verification is identical either way
  ```
  The script size-checks and md5-checks all four zips (`refdata/ascat_hg38_wgs.md5`)
  before extraction and validates the unpacked tables (22 alleles + 22 loci
  files; GC/RT are tab-separated with a header, every row 18 fields; RT spans
  chr1–X so chr22 presence is asserted). Note the zips ship inner files named
  `GC_G1000_hg38.txt` / `RT_G1000_hg38.txt` — the script copies them to the
  `*_WGS_hg38.txt` names that `rules/ascat.oxoflow` passes via `--gc`/`--rt`.
- **Two naming generations coexist on disk** (verified 2026-10-11): the
  pipeline consumes the **plain-numbered** set
  (`G1000_alleles_hg38_1..22.txt`, `G1000_loci_hg38_1..22.txt`) —
  `run_ascat.R` probes `..._1.txt` first and only falls back to the
  chr-prefixed names if plain is absent. A chr-prefixed copy
  (`G1000_alleles_hg38_chr1..22,X.txt`) is also present and harmless. Filenames
  stay plain because ASCAT's `readAlleleFiles` recovers the chromosome from
  the filename, while file **content** (and the BAM contigs it is matched
  against) is chr-prefixed. `scripts/verify_refdata.sh` gates on the plain
  set and reports the fallback set informationally.
- **Verify after any move/re-download:** `scripts/verify_refdata.sh` re-runs
  the table validation (22+22 plain files, 18-field GC/RT rows, RT chr22
  sentinel) without re-downloading.
- **Known failure mode (do not repeat):** a *resumed* download that appended
  ~12 MB of garbage once passed a naive size glance and killed ASCAT with
  `Stopped early on line 1813274. Expected 18 fields but found 51`. Never
  trust a resumed transfer without the md5 gate.
- **Cache note:** RT/GC tables are **not** rule inputs (only the `ascat_refdir`
  string is), so replacing them never invalidates other cached steps; failed
  `ascat_paired` outputs re-run automatically on the next invocation.

## 9. MSIsensor reference list — derived, not downloaded

`msi_paired` does not consume an external MSI database: the reference list of
microsatellite loci is **generated per cohort** from the genome FASTA by
`msisensor-pro scan -d $genome_fa` (rule `msi_scan`, shared across pairs).
Its provenance is therefore exactly §1's FASTA — nothing to download.

---

## 10. Clinical knowledgebases — CIViC monthly dump + ClinVar weekly VCF

Consumed by `rules/clinical.oxoflow` (`clinical_snv_annotate` and friends).
Deployed at `/home/data/reference/clinical_kb/` (override with
`CLINICAL_KB_DIR` for the prep script):

| File | Origin | Refresh |
|---|---|---|
| `<TAG>-VariantSummaries.tsv` | CIViC monthly release (CC0), `civicdb.org/releases` | monthly |
| `<TAG>-AcceptedClinicalEvidenceSummaries.tsv` | same release | monthly |
| `<TAG>-AcceptedAssertionSummaries.tsv` | same release (carries `amp_category` — the AMP tier backbone) | monthly |
| `clinvar.vcf.gz` (+`.tbi`) | NCBI `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/` | weekly |
| **`clinvar.chr.vcf.gz` (+`.tbi`)** | **derived locally**: contigs renamed `1`→`chr1` etc. (`bcftools annotate --rename-chrs`) — **this is the file `config.clinvar_vcf` points at**; NCBI's unprefixed contigs would make the join silently empty | rebuilt by the prep script when `clinvar.vcf.gz` updates |
| `.civic_tag` | marker written by the prep script (which CIViC tag this directory is complete for; fast-path key) | per refresh |

Download/refresh is one command — `scripts/prepare_clinical_kb.sh` discovers
the latest CIViC monthly release via the public GraphQL API, downloads the
TSVs + ClinVar VCF idempotently, and prints the config lines to paste into
`venus.oxoflow` (bump `clinical_kb_tag` so the report records the KB version).

**Coordinate-build caveat:** CIViC variant coordinates are GRCh37 (checked
2026-10: 1,292 GRCh37 vs 1 GRCh38 records in `VariantSummaries`), while the
pipeline is GRCh38. `clinical_snv_annotate` therefore only joins **ClinVar and
gnomAD** by GRCh38 coordinates via bcftools; CIViC matching happens in
`scripts/clinical_tier.py` build-independently — gene + normalized protein
change (primary), ClinVar VariationID bridge (secondary), GRCh38 coordinates
(tertiary) — the same strategy as the OncoKB annotator. When CIViC ships real
GRCh38 coordinates the tertiary path tightens automatically. Matching is
substitution-centric: in-frame indels/delins/exon skips (EGFR exon-19 del,
MET exon-14 class) have no portable protein key and match only via the
ClinVar VariationID bridge — a known recall gap (venus issue #30).

Verification is integrated into `scripts/verify_refdata.sh` (column-shape
checks on the TSVs; WARN-only until first download, like §7).

OncoKB (optional enrichment, `oncokb_token` in config) has **no local dump**:
the rule queries the REST API at run time with a free academic token from
<https://www.oncokb.org/api-access> (institutional email; Chakravarty et al.
2017, JCO PO).

---

## Adding a new host reference — checklist

1. Find the canonical source (vendor bundle, consortium FTP, Zenodo DOI);
   record the exact URL/version here.
2. Write/extend a download script under `scripts/` with size+md5 gates
   (pattern: `scripts/download_ascat_refdata.sh`); verify before extraction.
3. Deploy under the canonical root `/home/data/reference/` (or add a compat
   symlink if the config string must stay put) and point the `venus.oxoflow`
   config key at the deployed path.
4. Confirm the rule environment has the reader tool (e.g. `unzip`, `bcftools`).
5. Run `scripts/verify_refdata.sh` on the host — extend its checks if the new
   asset has its own structure to validate (file counts, table shapes,
   md5-known artifacts).
6. Update the table at the top of this page.
