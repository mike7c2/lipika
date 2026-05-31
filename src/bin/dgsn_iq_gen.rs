//! CLI entry point. Reads a TOML scene description, generates per-observer
//! IQ, writes `.npy` + JSON sidecar files into the output directory.

use std::path::PathBuf;

use clap::Parser;
use dgsn_iq_gen::{config::SceneConfig, pipeline};

#[derive(Parser, Debug)]
#[command(
    name = "dgsn-iq-gen",
    version,
    about = "Multi-source, multi-observer IQ scene simulator",
)]
struct Cli {
    /// Scene description file (TOML).
    #[arg(short, long)]
    config: PathBuf,
    /// Output directory; per-observer .npy files written here.
    #[arg(short, long)]
    output: PathBuf,
    /// Seed for the pipeline's per-observer RNG derivation.
    #[arg(long, default_value_t = 0u64)]
    seed: u64,
    /// Backend: cpu (default) or gpu (requires --features gpu at build time).
    #[arg(long, default_value = "cpu")]
    backend: String,
}

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp_millis()
        .init();

    let cli = Cli::parse();
    if cli.backend != "cpu" {
        anyhow::bail!("backend '{}' not yet implemented (only 'cpu' for now)", cli.backend);
    }

    log::info!("loading scene from {}", cli.config.display());
    let scene = SceneConfig::from_toml_file(&cli.config)?;
    let total_bands: usize = scene.receivers.iter().map(|r| r.bands.len()).sum();
    log::info!(
        "scene: epoch={} duration={}s, {} transmitters, {} receivers ({} bands total)",
        scene.scene.epoch_unix, scene.scene.duration_s,
        scene.transmitters.len(), scene.receivers.len(), total_bands,
    );

    let start = std::time::Instant::now();
    let outputs = pipeline::run_scene(&scene, cli.seed)?;
    let elapsed = start.elapsed();
    log::info!("pipeline finished in {:.2?}", elapsed);

    log::info!("writing {} observer outputs to {}", outputs.len(), cli.output.display());
    for output in &outputs {
        output.write_raw_with_sidecar(&cli.output)?;
        log::info!("  wrote {} ({} samples)", output.name, output.iq.len());
    }
    Ok(())
}
