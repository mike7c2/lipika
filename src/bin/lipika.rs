//! CLI entry point. Reads a TOML scene description and streams per-observer
//! IQ (raw `complex64` + JSON sidecar) into the output directory.

use std::path::PathBuf;

use clap::Parser;
use lipika::{config::SceneConfig, pipeline, profile};

#[derive(Parser, Debug)]
#[command(
    name = "lipika",
    version,
    about = "Multi-source, multi-observer IQ scene simulator",
)]
struct Cli {
    /// Scene description file (TOML).
    #[arg(short, long)]
    config: PathBuf,
    /// Output directory; per-observer-band .iq.c64 files written here.
    #[arg(short, long)]
    output: PathBuf,
    /// Seed for the pipeline's per-observer RNG derivation.
    // r[impl fn.7]
    #[arg(long, default_value_t = 0u64)]
    seed: u64,
    /// Backend: cpu (default) or gpu (requires --features gpu at build time).
    #[arg(long, default_value = "cpu")]
    backend: String,
    /// Streaming chunk size in samples. Peak memory ≈ chunk × (#in-band
    /// transmitters + few) × 8 B, flat in scene duration.
    // r[impl nf.2]
    #[arg(long, default_value_t = lipika::pipeline::DEFAULT_CHUNK_SAMPLES)]
    chunk_samples: usize,
    /// Worker threads for the per-chunk transmitter fan-out (0 = all cores).
    // r[impl nf.3]
    #[arg(long, default_value_t = 0usize)]
    threads: usize,
    /// If > 0, also write per-receiver geometric truth (delay, range,
    /// range-rate for every transmitter) sampled at this rate, as
    /// `<output>/<receiver>__truth.json`. 0 disables it.
    #[arg(long, default_value_t = 0.0f64)]
    truth_hz: f64,
    /// Print a per-section CPU profile at the end of the run.
    #[arg(long, default_value_t = false)]
    profile: bool,
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

    log::info!(
        "streaming to {} (chunk={} samples, threads={})",
        cli.output.display(), cli.chunk_samples,
        if cli.threads == 0 { "all".to_string() } else { cli.threads.to_string() },
    );
    let start = std::time::Instant::now();
    let manifests = pipeline::run_scene_to_dir(
        &scene, cli.seed, &cli.output, cli.chunk_samples, cli.threads,
    )?;
    let elapsed = start.elapsed();
    log::info!("pipeline finished in {:.2?}", elapsed);

    if cli.truth_hz > 0.0 {
        log::info!("writing geometric truth at {} Hz", cli.truth_hz);
        pipeline::export_truth(&scene, &cli.output, cli.truth_hz)?;
    }

    for m in &manifests {
        log::info!("  wrote {} ({} samples) -> {}", m.name, m.n_samples, m.path.display());
    }

    if cli.profile {
        println!();
        profile::print_report();
    }
    Ok(())
}
