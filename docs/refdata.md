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

Host reference root: `/home/data/reference/`

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

- **Path:** `/home/data/reference/dnbc4tools/ref-Homo_sapiens/`
- **Origin:** BGI **vendor bundle** shipped with dnbc4tools 3.0 for DNBelab
  C-series chemistry (prebuilt genomeDir + whitelists). Vendor-matched
  chemistry data cannot be reconstructed generically; use the bundle the
  vendor distributes for the exact kit version.
- **Verify:** `dnbc4tools` validates the reference at startup and reports its
  identity in the count log.

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

## 8. ASCAT G1000 hg38 WGS reference set (Zenodo 14008443)

- **Path:** `/home/data/wsx/hnzl_data/refdata/ascat/v3.2-hg38-wgs/` (config
  key `ascat_refdir`)
- **Origin:** official VanLoo-lab ASCAT NGS references,
  `https://zenodo.org/records/14008443` — allele/loci files for
  alleleCounter, GC-content and replication-timing tables for logR
  correction (see docs/method-rationale.md §13).
- **Download + verify (automated):**
  ```bash
  scripts/download_ascat_refdata.sh -o /home/data/wsx/hnzl_data/refdata/ascat/v3.2-hg38-wgs
  # flaky link: add -r (resume loop); verification is identical either way
  ```
  The script size-checks and md5-checks all four zips (`refdata/ascat_hg38_wgs.md5`)
  before extraction and validates the unpacked tables (22 alleles + 22 loci
  files; GC/RT are tab-separated with a header, every row 18 fields; RT spans
  chr1–X so chr22 presence is asserted). Note the zips ship inner files named
  `GC_G1000_hg38.txt` / `RT_G1000_hg38.txt` — the script copies them to the
  `*_WGS_hg38.txt` names that `rules/ascat.oxoflow` passes via `--gc`/`--rt`.
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

## Adding a new host reference — checklist

1. Find the canonical source (vendor bundle, consortium FTP, Zenodo DOI);
   record the exact URL/version here.
2. Write/extend a download script under `scripts/` with size+md5 gates
   (pattern: `scripts/download_ascat_refdata.sh`); verify before extraction.
3. Point the `venus.oxoflow` config key at the deployed path.
4. Confirm the rule environment has the reader tool (e.g. `unzip`, `bcftools`).
5. Update the table at the top of this page.
