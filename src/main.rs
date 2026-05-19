#![forbid(unsafe_code)]
//! Venus CLI — Clinical-grade tumor variant detection pipeline generator.
//!
//! This binary provides CLI commands for generating, validating, and inspecting
//! Venus pipeline configurations that produce `.oxoflow` workflow files for the
//! oxo-flow pipeline engine.

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

/// Venus — Clinical-grade tumor variant detection pipeline.
#[derive(Parser)]
#[command(name = "venus")]
#[command(author = "Traitome")]
#[command(version)]
#[command(about = "Venus (启明星) — Clinical-grade tumor variant detection pipeline", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate a .oxoflow workflow file from configuration.
    Generate {
        /// Path to Venus configuration TOML file.
        config: PathBuf,
        /// Output workflow file path.
        #[arg(short, long, default_value = "venus.oxoflow")]
        output: PathBuf,
    },
    /// Validate a Venus configuration file.
    Validate {
        /// Path to Venus configuration TOML file.
        config: PathBuf,
        /// Skip checking if files exist.
        #[arg(long)]
        skip_file_check: bool,
    },
    /// List all available pipeline steps.
    ListSteps,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Generate { config, output } => {
            let cfg = oxo_flow_venus::VenusConfig::from_file(&config)?;
            let oxoflow = oxo_flow_venus::generate_oxoflow(&cfg)?;
            std::fs::write(&output, oxoflow)?;
            println!("Generated: {}", output.display());
        }
        Commands::Validate {
            config,
            skip_file_check: _,
        } => {
            let cfg = oxo_flow_venus::VenusConfig::from_file(&config)?;
            println!("✓ Configuration valid");
            println!("  Mode: {}", cfg.mode);
            println!("  Sequencing: {}", cfg.seq_type);
            println!("  Genome: {}", cfg.genome_build);
            println!("  Samples: {}", cfg.samples.len());
        }
        Commands::ListSteps => {
            println!("Venus Pipeline Steps:");
            println!();
            println!("  1. fastp           — FASTQ quality control and trimming");
            println!("  2. bwa_mem2        — Read alignment to reference genome");
            println!("  3. mark_duplicates — PCR duplicate marking");
            println!("  4. bqsr            — Base quality score recalibration");
            println!("  5. mutect2         — Somatic variant calling (GATK)");
            println!("  6. filter_mutect   — Somatic variant filtering");
            println!("  7. strelka2        — Somatic variant calling (paired mode)");
            println!("  8. vep             — Variant effect prediction/annotation");
            println!("  9. clinical_report — Clinical report generation");
        }
    }
    Ok(())
}
