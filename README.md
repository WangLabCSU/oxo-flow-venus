# Venus (启明星)

**Venus** (启明星 — "Morning Star") is a clinical-grade tumor variant detection pipeline built on [oxo-flow](https://github.com/Traitome/oxo-flow).

## Features

- **Multiple Analysis Modes**
  - `experiment_only`: Tumor-only somatic calling
  - `control_only`: Germline variant calling
  - `experiment_control`: Paired tumor-normal analysis

- **Sequencing Support**
  - Whole Genome Sequencing (WGS)
  - Whole Exome Sequencing (WES)
  - Targeted Panel Sequencing

- **Core Tools**
  - **fastp/FastQC**: Quality control and trimming
  - **BWA-MEM2**: Read alignment
  - **GATK 4**: MarkDuplicates, BQSR, Mutect2, HaplotypeCaller
  - **Strelka2**: Alternative somatic caller (paired mode)
  - **VEP**: Variant annotation

## Installation

### From Source

```bash
git clone https://github.com/Traitome/oxo-flow-venus.git
cd oxo-flow-venus
cargo install --path .
```

### Prerequisites

- Rust 1.75+ (for building)
- conda/mamba (for environments)
- oxo-flow CLI (`cargo install oxo-flow`)

## Quick Start

1. **Create a configuration file** (`config.toml`):

```toml
[venus]
mode = "experiment_control"
seq_type = "wes"
genome_build = "GRCh38"
reference_fasta = "/path/to/genome.fa"
target_bed = "/path/to/exome.bed"

[[samples]]
name = "TUMOR"
type = "tumor"
r1 = "TUMOR_R1.fq.gz"
r2 = "TUMOR_R2.fq.gz"
pair_id = "NORMAL"

[[samples]]
name = "NORMAL"
type = "normal"
r1 = "NORMAL_R1.fq.gz"
r2 = "NORMAL_R2.fq.gz"
```

2. **Generate workflow**:

```bash
venus generate config.toml -o venus.oxoflow
```

3. **Run pipeline**:

```bash
oxo-flow run venus.oxoflow -j 8
```

## CLI Commands

| Command | Description |
|---------|-------------|
| `venus generate <config>` | Generate .oxoflow workflow file |
| `venus validate <config>` | Validate configuration |
| `venus list-steps` | List pipeline steps |

## Pipeline Steps

1. **fastp** — FASTQ QC and trimming
2. **BWA-MEM2** — Read alignment
3. **MarkDuplicates** — PCR duplicate marking
4. **BQSR** — Base quality recalibration
5. **Mutect2** — Somatic variant calling
6. **FilterMutectCalls** — Variant filtering
7. **Strelka2** — Alternative somatic caller (paired mode)
8. **VEP** — Variant annotation
9. **Clinical Report** — Report generation

## Output Structure

```
results/
├── trimmed/          # Trimmed FASTQ
├── aligned/          # Aligned BAM
├── dedup/            # Deduplicated BAM
├── recal/            # Recalibrated BAM
├── variants/         # VCF files
├── annotated/        # Annotated VCF
├── tmb/              # TMB calculations
└── reports/          # Clinical reports
```

## License

Apache-2.0

## Links

- [oxo-flow](https://github.com/Traitome/oxo-flow)
- [GATK](https://gatk.broadinstitute.org/)
- [VEP](https://www.ensembl.org/vep)
