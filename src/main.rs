mod dedup;
mod network;

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(name = "umi-tools-rs")]
#[command(about = "Rust implementation of selected UMI-tools commands")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Deduplicate reads using UMI and mapping coordinates.
    Dedup(dedup::DedupArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Dedup(args) => dedup::run(args),
    }
}
