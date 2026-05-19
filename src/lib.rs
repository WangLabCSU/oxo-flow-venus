//! Venus configuration library for clinical-grade tumor variant detection.
//!
//! This library provides configuration types and validation for Venus workflows,
//! generating `.oxoflow` workflow files for the oxo-flow pipeline engine.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use thiserror::Error;

/// Errors that can occur during Venus configuration processing.
#[derive(Debug, Error)]
pub enum VenusError {
    /// A required field is missing from the configuration.
    #[error("missing required field: {0}")]
    MissingRequiredField(String),

    /// Invalid analysis mode specified.
    #[error("invalid analysis mode: {0}. Expected one of: experiment_only, control_only, experiment_control")]
    InvalidMode(String),

    /// Invalid sequencing type specified.
    #[error("invalid sequencing type: {0}. Expected one of: wgs, wes, panel")]
    InvalidSeqType(String),

    /// Target BED file required but not provided for WES/Panel modes.
    #[error("target_bed is required for {mode} sequencing mode")]
    MissingTargetBed { mode: String },

    /// A sample is missing required fields.
    #[error("sample '{name}' is missing required field: {field}")]
    SampleMissingField { name: String, field: String },

    /// A pair_id references a sample that doesn't exist.
    #[error("pair_id '{pair_id}' references non-existent sample: {sample}")]
    PairIdNotFound { pair_id: String, sample: String },

    /// A referenced file does not exist.
    #[error("file not found: {0}")]
    FileNotFound(String),

    /// Workflow generation failed.
    #[error("failed to generate oxoflow workflow: {0}")]
    GenerationFailed(String),

    /// Inconsistent analysis mode with sample types.
    #[error("inconsistent mode: analysis mode '{mode}' does not match sample types provided")]
    InconsistentModeSamples { mode: String },

    /// I/O error wrapper.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// TOML deserialization error.
    #[error("TOML parse error: {0}")]
    Toml(#[from] toml::de::Error),

    /// JSON serialization/deserialization error.
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
}

/// Generate a .oxoflow TOML workflow from a Venus configuration.
pub fn generate_oxoflow(config: &VenusConfig) -> Result<String, VenusError> {
    config.generate_oxoflow()
}

/// Analysis mode determining what types of samples are processed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum AnalysisMode {
    /// Process experiment (tumor) samples only without matched normal.
    #[default]
    ExperimentOnly,

    /// Process control (normal) samples only for germline calling.
    ControlOnly,

    /// Process tumor-normal pairs for somatic calling.
    ExperimentControl,
}

impl std::fmt::Display for AnalysisMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AnalysisMode::ExperimentOnly => write!(f, "experiment_only"),
            AnalysisMode::ControlOnly => write!(f, "control_only"),
            AnalysisMode::ExperimentControl => write!(f, "experiment_control"),
        }
    }
}

/// Sequencing type determining the capture scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum SeqType {
    /// Whole genome sequencing.
    WGS,

    /// Whole exome sequencing.
    #[default]
    WES,

    /// Targeted panel sequencing.
    Panel,
}

impl std::fmt::Display for SeqType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqType::WGS => write!(f, "wgs"),
            SeqType::WES => write!(f, "wes"),
            SeqType::Panel => write!(f, "panel"),
        }
    }
}

/// Reference genome build version.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum GenomeBuild {
    /// GRCh37 (hg19) reference build.
    GRCh37,

    /// GRCh38 (hg38) reference build.
    #[default]
    GRCh38,
}

impl std::fmt::Display for GenomeBuild {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GenomeBuild::GRCh37 => write!(f, "GRCh37"),
            GenomeBuild::GRCh38 => write!(f, "GRCh38"),
        }
    }
}

/// Sample type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum SampleType {
    /// Tumor (experiment) sample.
    Tumor,

    /// Normal (control) sample.
    Normal,
}

impl std::fmt::Display for SampleType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SampleType::Tumor => write!(f, "tumor"),
            SampleType::Normal => write!(f, "normal"),
        }
    }
}

/// A single sample with associated metadata and file paths.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sample {
    /// Unique sample identifier.
    pub name: String,

    /// Sample type (tumor or normal).
    #[serde(rename = "type")]
    pub sample_type: SampleType,

    /// Path to R1 (forward) FASTQ file.
    pub r1: String,

    /// Optional path to R2 (reverse) FASTQ file for paired-end sequencing.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r2: Option<String>,

    /// Optional pair ID linking tumor to matched normal in ExperimentControl mode.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pair_id: Option<String>,

    /// Additional sample-level metadata.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub metadata: HashMap<String, String>,
}

impl Sample {
    /// Create a new sample with required fields.
    pub fn new(name: impl Into<String>, sample_type: SampleType, r1: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            sample_type,
            r1: r1.into(),
            r2: None,
            pair_id: None,
            metadata: HashMap::new(),
        }
    }

    /// Add R2 (reverse) FASTQ path.
    pub fn with_r2(mut self, r2: impl Into<String>) -> Self {
        self.r2 = Some(r2.into());
        self
    }

    /// Add pair ID for tumor-normal pairing.
    pub fn with_pair_id(mut self, pair_id: impl Into<String>) -> Self {
        self.pair_id = Some(pair_id.into());
        self
    }

    /// Add metadata key-value pair.
    pub fn with_metadata(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }
}

/// Environment configuration for conda environments.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EnvGroup {
    /// Path to conda environment file.
    pub conda: String,
}

impl Default for EnvGroup {
    fn default() -> Self {
        Self {
            conda: "envs/default.yaml".to_string(),
        }
    }
}

/// Default resource settings applied to all rules.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Defaults {
    /// Number of threads per rule.
    #[serde(default = "default_threads")]
    pub threads: u32,

    /// Memory allocation per rule.
    #[serde(default = "default_memory")]
    pub memory: String,

    /// Environment group for tool execution.
    #[serde(default = "default_env_group")]
    pub env_group: String,
}

fn default_threads() -> u32 {
    8
}

fn default_memory() -> String {
    "16G".to_string()
}

fn default_env_group() -> String {
    "default".to_string()
}

impl Default for Defaults {
    fn default() -> Self {
        Self {
            threads: default_threads(),
            memory: default_memory(),
            env_group: default_env_group(),
        }
    }
}

/// Complete Venus workflow configuration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenusConfig {
    /// Workflow name.
    pub name: String,

    /// Workflow version.
    #[serde(default = "default_version")]
    pub version: String,

    /// Workflow description.
    #[serde(default = "default_description")]
    pub description: String,

    /// Analysis mode (experiment_only, control_only, experiment_control).
    #[serde(default)]
    pub mode: AnalysisMode,

    /// Sequencing type (wgs, wes, panel).
    #[serde(default)]
    pub seq_type: SeqType,

    /// Reference genome build.
    #[serde(default)]
    pub genome_build: GenomeBuild,

    /// Path to reference genome FASTA.
    pub reference_fasta: String,

    /// Path to known sites VCF for BQSR.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub known_sites: Option<String>,

    /// Path to target BED file (required for WES/Panel).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_bed: Option<String>,

    /// List of samples to process.
    pub samples: Vec<Sample>,

    /// Environment groups mapping names to environment configs.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub env_groups: HashMap<String, EnvGroup>,

    /// Default resource settings.
    #[serde(default)]
    pub defaults: Defaults,

    /// Output directory for results.
    #[serde(default = "default_output_dir")]
    pub output_dir: String,
}

fn default_version() -> String {
    "0.5.4".to_string()
}

fn default_description() -> String {
    "Clinical-grade tumor variant detection pipeline".to_string()
}

fn default_output_dir() -> String {
    "results".to_string()
}

impl VenusConfig {
    /// Load configuration from a TOML file.
    pub fn from_file(path: &Path) -> Result<Self, VenusError> {
        let content = std::fs::read_to_string(path)?;
        let config: Self = toml::from_str(&content)?;
        config.validate()?;
        Ok(config)
    }

    /// Validate the configuration for consistency.
    pub fn validate(&self) -> Result<(), VenusError> {
        // Validate required fields
        if self.reference_fasta.is_empty() {
            return Err(VenusError::MissingRequiredField("reference_fasta".to_string()));
        }

        if self.samples.is_empty() {
            return Err(VenusError::MissingRequiredField("samples".to_string()));
        }

        // Validate sequencing type constraints
        match self.seq_type {
            SeqType::WES | SeqType::Panel => {
                if self.target_bed.is_none() || self.target_bed.as_ref().is_none_or(|s| s.is_empty()) {
                    return Err(VenusError::MissingTargetBed {
                        mode: self.seq_type.to_string(),
                    });
                }
            }
            SeqType::WGS => {}
        }

        // Validate samples based on mode
        for sample in &self.samples {
            if sample.name.is_empty() {
                return Err(VenusError::SampleMissingField {
                    name: "<unnamed>".to_string(),
                    field: "name".to_string(),
                });
            }
            if sample.r1.is_empty() {
                return Err(VenusError::SampleMissingField {
                    name: sample.name.clone(),
                    field: "r1".to_string(),
                });
            }
        }

        // Validate mode consistency with sample types
        match self.mode {
            AnalysisMode::ExperimentOnly => {
                for sample in &self.samples {
                    if sample.sample_type == SampleType::Normal {
                        return Err(VenusError::InconsistentModeSamples {
                            mode: self.mode.to_string(),
                        });
                    }
                }
            }
            AnalysisMode::ControlOnly => {
                for sample in &self.samples {
                    if sample.sample_type == SampleType::Tumor {
                        return Err(VenusError::InconsistentModeSamples {
                            mode: self.mode.to_string(),
                        });
                    }
                }
            }
            AnalysisMode::ExperimentControl => {
                // Verify that pairs are properly defined
                let sample_names: std::collections::HashSet<_> =
                    self.samples.iter().map(|s| &s.name).collect();

                for sample in &self.samples {
                    if sample.sample_type == SampleType::Tumor
                        && let Some(pair_id) = &sample.pair_id
                        && !sample_names.contains(pair_id)
                    {
                        return Err(VenusError::PairIdNotFound {
                            pair_id: pair_id.clone(),
                            sample: sample.name.clone(),
                        });
                    }
                }
            }
        }

        Ok(())
    }

    /// Generate a TOML workflow configuration for oxo-flow.
    pub fn generate_oxoflow(&self) -> Result<String, VenusError> {
        let mut output = String::new();

        // Workflow header
        output.push_str("# Generated by Venus CLI\n");
        output.push_str("# Clinical-grade tumor variant detection pipeline\n\n");

        // Workflow metadata
        output.push_str("[workflow]\n");
        output.push_str(&format!("name = \"{}\"\n", self.name));
        output.push_str(&format!("version = \"{}\"\n", self.version));
        output.push_str(&format!("description = \"{}\"\n", self.description));
        output.push_str("author = \"Venus\"\n\n");

        // Config section
        output.push_str("[config]\n");
        output.push_str(&format!("reference_fasta = \"{}\"\n", self.reference_fasta));
        if let Some(known_sites) = &self.known_sites {
            output.push_str(&format!("known_sites = \"{}\"\n", known_sites));
        }
        output.push_str(&format!("genome_build = \"{}\"\n", self.genome_build));
        if let Some(target_bed) = &self.target_bed {
            output.push_str(&format!("target_bed = \"{}\"\n", target_bed));
        }
        output.push_str(&format!("output_dir = \"{}\"\n", self.output_dir));
        output.push_str(&format!("mode = \"{}\"\n", self.mode));
        output.push_str(&format!("seq_type = \"{}\"\n", self.seq_type));
        output.push('\n');

        // Sample list
        output.push_str("# Sample definitions\n");
        output.push_str("[config.samples]\n");
        for sample in &self.samples {
            output.push_str(&format!("{} = {{ type = \"{}\", r1 = \"{}\"",
                sample.name, sample.sample_type, sample.r1));
            if let Some(r2) = &sample.r2 {
                output.push_str(&format!(", r2 = \"{}\"", r2));
            }
            if let Some(pair_id) = &sample.pair_id {
                output.push_str(&format!(", pair_id = \"{}\"", pair_id));
            }
            output.push_str(" }\n");
        }
        output.push('\n');

        // Defaults section
        output.push_str("[defaults]\n");
        output.push_str(&format!("threads = {}\n", self.defaults.threads));
        output.push_str(&format!("memory = \"{}\"\n", self.defaults.memory));
        output.push_str(&format!("env_group = \"{}\"\n", self.defaults.env_group));
        output.push('\n');

        // Environment groups
        if !self.env_groups.is_empty() {
            output.push_str("# Environment definitions\n");
            for (name, env) in &self.env_groups {
                output.push_str(&format!("[env_groups.{}]\n", name));
                output.push_str(&format!("conda = \"{}\"\n", env.conda));
            }
            output.push('\n');
        }

        // Generate rules based on mode
        self.generate_rules(&mut output)?;

        Ok(output)
    }

    /// Generate workflow rules based on analysis mode.
    fn generate_rules(&self, output: &mut String) -> Result<(), VenusError> {
        output.push_str("# Pipeline rules\n");

        // Generate QC rules for all samples
        for sample in &self.samples {
            self.generate_sample_qc_rule(sample, output)?;
        }

        // Generate alignment rules for all samples
        for sample in &self.samples {
            self.generate_alignment_rule(sample, output)?;
        }

        // Generate variant calling rules based on mode
        match self.mode {
            AnalysisMode::ExperimentOnly => {
                // Tumor-only somatic calling
                for sample in &self.samples {
                    self.generate_tumor_only_calling(sample, output)?;
                }
            }
            AnalysisMode::ControlOnly => {
                // Germline calling only
                for sample in &self.samples {
                    self.generate_germline_calling(sample, output)?;
                }
            }
            AnalysisMode::ExperimentControl => {
                // Paired somatic calling
                for sample in &self.samples {
                    if sample.sample_type == SampleType::Tumor
                        && let Some(normal_name) = &sample.pair_id
                    {
                        self.generate_paired_calling(&sample.name, normal_name, output)?;
                    }
                }
            }
        }

        Ok(())
    }

    /// Generate QC rule for a sample.
    fn generate_sample_qc_rule(&self, sample: &Sample, output: &mut String) -> Result<(), VenusError> {
        output.push_str("\n[[rules]]\n");
        output.push_str(&format!("name = \"fastp_{}\"\n", sample.name));
        output.push_str("input = [");
        output.push_str(&format!("\"{}\"", sample.r1));
        if let Some(r2) = &sample.r2 {
            output.push_str(&format!(", \"{}\"", r2));
        }
        output.push_str("]\n");
        output.push_str(&format!("output = [\"{}/trimmed/{}_R1.fq.gz\"",
            self.output_dir, sample.name));
        if sample.r2.is_some() {
            output.push_str(&format!(", \"{}/trimmed/{}_R2.fq.gz\"", self.output_dir, sample.name));
        }
        output.push_str("]\n");
        if sample.r2.is_some() {
            output.push_str(&format!("shell = \"fastp -i {{input[0]}} -I {{input[1]}} -o {{output[0]}} -O {{output[1]}} --json {}/qc/{}.fastp.json --thread {{threads}}\"\n",
                self.output_dir, sample.name));
        } else {
            output.push_str(&format!("shell = \"fastp -i {{input[0]}} -o {{output[0]}} --json {}/qc/{}.fastp.json --thread {{threads}}\"\n",
                self.output_dir, sample.name));
        }
        output.push_str(&format!("threads = {}\n", self.defaults.threads));
        output.push_str("[rules.environment]\n");
        output.push_str("conda = \"envs/fastp.yaml\"\n");
        Ok(())
    }

    /// Generate alignment rule for a sample.
    fn generate_alignment_rule(&self, sample: &Sample, output: &mut String) -> Result<(), VenusError> {
        output.push_str("\n[[rules]]\n");
        output.push_str(&format!("name = \"align_{}\"\n", sample.name));
        if sample.r2.is_some() {
            output.push_str(&format!("input = [\"{}/trimmed/{}_R1.fq.gz\", \"{}/trimmed/{}_R2.fq.gz\"]\n",
                self.output_dir, sample.name, self.output_dir, sample.name));
            output.push_str(&format!("shell = \"bwa-mem2 mem -t {{threads}} {} {{input[0]}} {{input[1]}} | samtools sort -@ {{threads}} -o {{output[0]}}\"\n",
                self.reference_fasta));
        } else {
            output.push_str(&format!("input = [\"{}/trimmed/{}_R1.fq.gz\"]\n",
                self.output_dir, sample.name));
            output.push_str(&format!("shell = \"bwa-mem2 mem -t {{threads}} {} {{input[0]}} | samtools sort -@ {{threads}} -o {{output[0]}}\"\n",
                self.reference_fasta));
        }
        output.push_str(&format!("output = [\"{}/aligned/{}.sorted.bam\"]\n",
            self.output_dir, sample.name));
        output.push_str(&format!("threads = {}\n", self.defaults.threads));
        output.push_str(&format!("memory = \"{}\"\n", self.defaults.memory));
        output.push_str("[rules.environment]\n");
        output.push_str("conda = \"envs/bwa_mem2.yaml\"\n");
        Ok(())
    }

    /// Generate tumor-only somatic calling rule.
    fn generate_tumor_only_calling(&self, sample: &Sample, output: &mut String) -> Result<(), VenusError> {
        output.push_str("\n[[rules]]\n");
        output.push_str(&format!("name = \"mutect2_{}\"\n", sample.name));
        output.push_str(&format!("input = [\"{}/aligned/{}.sorted.bam\"]\n",
            self.output_dir, sample.name));
        output.push_str(&format!("output = [\"{}/variants/{}.mutect2.vcf.gz\"]\n",
            self.output_dir, sample.name));
        output.push_str(&format!("shell = \"gatk Mutect2 -I {{input[0]}} -R {} -O {{output[0]}}\"\n",
            self.reference_fasta));
        output.push_str(&format!("threads = {}\n", self.defaults.threads));
        output.push_str("[rules.environment]\n");
        output.push_str("conda = \"envs/gatk.yaml\"\n");
        Ok(())
    }

    /// Generate germline calling rule.
    fn generate_germline_calling(&self, sample: &Sample, output: &mut String) -> Result<(), VenusError> {
        output.push_str("\n[[rules]]\n");
        output.push_str(&format!("name = \"haplotype_{}\"\n", sample.name));
        output.push_str(&format!("input = [\"{}/aligned/{}.sorted.bam\"]\n",
            self.output_dir, sample.name));
        output.push_str(&format!("output = [\"{}/variants/{}.g.vcf.gz\"]\n",
            self.output_dir, sample.name));
        output.push_str(&format!("shell = \"gatk HaplotypeCaller -I {{input[0]}} -R {} -O {{output[0]}} -ERC GVCF\"\n",
            self.reference_fasta));
        output.push_str(&format!("threads = {}\n", self.defaults.threads));
        output.push_str("[rules.environment]\n");
        output.push_str("conda = \"envs/gatk.yaml\"\n");
        Ok(())
    }

    /// Generate paired tumor-normal calling rule.
    fn generate_paired_calling(&self, tumor: &str, normal: &str, output: &mut String) -> Result<(), VenusError> {
        output.push_str("\n[[rules]]\n");
        output.push_str(&format!("name = \"mutect2_{}_paired\"\n", tumor));
        output.push_str(&format!("input = [\"{}/aligned/{}.sorted.bam\", \"{}/aligned/{}.sorted.bam\"]\n",
            self.output_dir, tumor, self.output_dir, normal));
        output.push_str(&format!("output = [\"{}/variants/{}.mutect2.vcf.gz\"]\n",
            self.output_dir, tumor));
        output.push_str(&format!("shell = \"gatk Mutect2 -I {{input[0]}} -I {{input[1]}} -normal {} -R {} -O {{output[0]}}\"\n",
            normal, self.reference_fasta));
        output.push_str(&format!("threads = {}\n", self.defaults.threads));
        output.push_str("[rules.environment]\n");
        output.push_str("conda = \"envs/gatk.yaml\"\n");
        Ok(())
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analysis_mode_display() {
        assert_eq!(AnalysisMode::ExperimentOnly.to_string(), "experiment_only");
        assert_eq!(AnalysisMode::ControlOnly.to_string(), "control_only");
        assert_eq!(AnalysisMode::ExperimentControl.to_string(), "experiment_control");
    }

    #[test]
    fn seq_type_display() {
        assert_eq!(SeqType::WGS.to_string(), "wgs");
        assert_eq!(SeqType::WES.to_string(), "wes");
        assert_eq!(SeqType::Panel.to_string(), "panel");
    }

    #[test]
    fn genome_build_display() {
        assert_eq!(GenomeBuild::GRCh37.to_string(), "GRCh37");
        assert_eq!(GenomeBuild::GRCh38.to_string(), "GRCh38");
    }

    #[test]
    fn sample_type_display() {
        assert_eq!(SampleType::Tumor.to_string(), "tumor");
        assert_eq!(SampleType::Normal.to_string(), "normal");
    }

    #[test]
    fn validate_wes_needs_bed() {
        let config = VenusConfig {
            name: "test".to_string(),
            version: "0.5.4".to_string(),
            description: "test".to_string(),
            mode: AnalysisMode::ExperimentOnly,
            seq_type: SeqType::WES,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: None,
            target_bed: None,
            samples: vec![Sample::new("S1", SampleType::Tumor, "/data/S1_R1.fq.gz")],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(matches!(result, Err(VenusError::MissingTargetBed { .. })));
    }

    #[test]
    fn validate_produces_toml() {
        let config = VenusConfig {
            name: "test_pipeline".to_string(),
            version: "0.5.4".to_string(),
            description: "Test pipeline".to_string(),
            mode: AnalysisMode::ExperimentOnly,
            seq_type: SeqType::WGS,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: Some("/data/dbsnp.vcf.gz".to_string()),
            target_bed: None,
            samples: vec![
                Sample::new("Tumor1", SampleType::Tumor, "/data/T1_R1.fq.gz")
                    .with_r2("/data/T1_R2.fq.gz"),
            ],
            env_groups: {
                let mut map = HashMap::new();
                map.insert("gatk".to_string(), EnvGroup {
                    conda: "envs/gatk.yaml".to_string(),
                });
                map
            },
            defaults: Defaults {
                threads: 16,
                memory: "32G".to_string(),
                env_group: "gatk".to_string(),
            },
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(result.is_ok());

        let toml = config.generate_oxoflow().unwrap();
        assert!(toml.contains("[workflow]"));
        assert!(toml.contains("name = \"test_pipeline\""));
        assert!(toml.contains("[config]"));
        assert!(toml.contains("reference_fasta = \"/data/genome.fa\""));
        assert!(toml.contains("[defaults]"));
        assert!(toml.contains("threads = 16"));
        assert!(toml.contains("memory = \"32G\""));
    }

    #[test]
    fn validate_experiment_only_with_normal_fails() {
        let config = VenusConfig {
            name: "test".to_string(),
            version: "0.5.4".to_string(),
            description: "test".to_string(),
            mode: AnalysisMode::ExperimentOnly,
            seq_type: SeqType::WGS,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: None,
            target_bed: None,
            samples: vec![Sample::new("N1", SampleType::Normal, "/data/N1_R1.fq.gz")],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(matches!(result, Err(VenusError::InconsistentModeSamples { .. })));
    }

    #[test]
    fn validate_control_only_with_tumor_fails() {
        let config = VenusConfig {
            name: "test".to_string(),
            version: "0.5.4".to_string(),
            description: "test".to_string(),
            mode: AnalysisMode::ControlOnly,
            seq_type: SeqType::WGS,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: None,
            target_bed: None,
            samples: vec![Sample::new("T1", SampleType::Tumor, "/data/T1_R1.fq.gz")],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(matches!(result, Err(VenusError::InconsistentModeSamples { .. })));
    }

    #[test]
    fn validate_pair_id_not_found() {
        let config = VenusConfig {
            name: "test".to_string(),
            version: "0.5.4".to_string(),
            description: "test".to_string(),
            mode: AnalysisMode::ExperimentControl,
            seq_type: SeqType::WGS,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: None,
            target_bed: None,
            samples: vec![
                Sample::new("T1", SampleType::Tumor, "/data/T1_R1.fq.gz")
                    .with_pair_id("N_NONEXISTENT"),
                Sample::new("N1", SampleType::Normal, "/data/N1_R1.fq.gz"),
            ],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(matches!(result, Err(VenusError::PairIdNotFound { .. })));
    }

    #[test]
    fn validate_experiment_control_success() {
        let config = VenusConfig {
            name: "test".to_string(),
            version: "0.5.4".to_string(),
            description: "test".to_string(),
            mode: AnalysisMode::ExperimentControl,
            seq_type: SeqType::WGS,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: None,
            target_bed: None,
            samples: vec![
                Sample::new("T1", SampleType::Tumor, "/data/T1_R1.fq.gz")
                    .with_pair_id("N1"),
                Sample::new("N1", SampleType::Normal, "/data/N1_R1.fq.gz"),
            ],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(result.is_ok());
    }

    #[test]
    fn validate_empty_samples_fails() {
        let config = VenusConfig {
            name: "test".to_string(),
            version: "0.5.4".to_string(),
            description: "test".to_string(),
            mode: AnalysisMode::ExperimentOnly,
            seq_type: SeqType::WGS,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/data/genome.fa".to_string(),
            known_sites: None,
            target_bed: None,
            samples: vec![],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "results".to_string(),
        };

        let result = config.validate();
        assert!(matches!(result, Err(VenusError::MissingRequiredField(_))));
    }

    #[test]
    fn generate_oxoflow_includes_all_sections() {
        let config = VenusConfig {
            name: "complete_test".to_string(),
            version: "0.5.4".to_string(),
            description: "Complete test pipeline".to_string(),
            mode: AnalysisMode::ExperimentControl,
            seq_type: SeqType::Panel,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/ref/genome.fa".to_string(),
            known_sites: Some("/ref/dbsnp.vcf.gz".to_string()),
            target_bed: Some("/ref/targets.bed".to_string()),
            samples: vec![
                Sample::new("Tumor1", SampleType::Tumor, "/data/T1_R1.fq.gz")
                    .with_r2("/data/T1_R2.fq.gz")
                    .with_pair_id("Normal1")
                    .with_metadata("patient_id", "P001"),
                Sample::new("Normal1", SampleType::Normal, "/data/N1_R1.fq.gz")
                    .with_r2("/data/N1_R2.fq.gz"),
            ],
            env_groups: {
                let mut map = HashMap::new();
                map.insert("align".to_string(), EnvGroup {
                    conda: "envs/bwa.yaml".to_string(),
                });
                map.insert("call".to_string(), EnvGroup {
                    conda: "envs/gatk.yaml".to_string(),
                });
                map
            },
            defaults: Defaults {
                threads: 24,
                memory: "64G".to_string(),
                env_group: "align".to_string(),
            },
            output_dir: "output".to_string(),
        };

        config.validate().unwrap();
        let toml = config.generate_oxoflow().unwrap();

        // Verify all major sections are present
        assert!(toml.contains("[workflow]"));
        assert!(toml.contains("[config]"));
        assert!(toml.contains("[defaults]"));
        assert!(toml.contains("[env_groups"));
        assert!(toml.contains("[[rules]]"));

        // Verify key values
        assert!(toml.contains("name = \"complete_test\""));
        assert!(toml.contains("seq_type = \"panel\""));
        assert!(toml.contains("target_bed = \"/ref/targets.bed\""));
        assert!(toml.contains("threads = 24"));
        assert!(toml.contains("memory = \"64G\""));
    }

    #[test]
    fn sample_builder_pattern() {
        let sample = Sample::new("TestSample", SampleType::Tumor, "/data/R1.fq.gz")
            .with_r2("/data/R2.fq.gz")
            .with_pair_id("MatchedNormal")
            .with_metadata("patient", "P001")
            .with_metadata("batch", "B001");

        assert_eq!(sample.name, "TestSample");
        assert_eq!(sample.sample_type, SampleType::Tumor);
        assert_eq!(sample.r1, "/data/R1.fq.gz");
        assert_eq!(sample.r2, Some("/data/R2.fq.gz".to_string()));
        assert_eq!(sample.pair_id, Some("MatchedNormal".to_string()));
        assert_eq!(sample.metadata.get("patient"), Some(&"P001".to_string()));
        assert_eq!(sample.metadata.get("batch"), Some(&"B001".to_string()));
    }

    #[test]
    fn defaults_default_values() {
        let defaults = Defaults::default();
        assert_eq!(defaults.threads, 8);
        assert_eq!(defaults.memory, "16G");
        assert_eq!(defaults.env_group, "default");
    }

    #[test]
    fn enum_default_values() {
        assert!(matches!(AnalysisMode::default(), AnalysisMode::ExperimentOnly));
        assert!(matches!(SeqType::default(), SeqType::WES));
        assert!(matches!(GenomeBuild::default(), GenomeBuild::GRCh38));
    }

    #[test]
    fn serialization_roundtrip() {
        let config = VenusConfig {
            name: "roundtrip".to_string(),
            version: "0.5.4".to_string(),
            description: "Test".to_string(),
            mode: AnalysisMode::ExperimentOnly,
            seq_type: SeqType::WES,
            genome_build: GenomeBuild::GRCh38,
            reference_fasta: "/ref.fa".to_string(),
            known_sites: None,
            target_bed: Some("/targets.bed".to_string()),
            samples: vec![Sample::new("S1", SampleType::Tumor, "/data/S1.fq.gz")],
            env_groups: HashMap::new(),
            defaults: Defaults::default(),
            output_dir: "out".to_string(),
        };

        let json = serde_json::to_string(&config).unwrap();
        let parsed: VenusConfig = serde_json::from_str(&json).unwrap();

        assert_eq!(config.name, parsed.name);
        assert_eq!(config.mode, parsed.mode);
        assert_eq!(config.seq_type, parsed.seq_type);
        assert_eq!(config.samples.len(), parsed.samples.len());
    }
}
