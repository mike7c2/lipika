//! Throughput of `lipika::pipeline::render_window`: how many seconds of IQ one
//! emulated receiver produces per wall-clock second.
//!
//!   cargo run --release --example window_throughput
//!
//! One receiver (Glasgow, 290 K, horizon cutoff) on the Hiber beacon sky
//! (`examples/hiber_beacon_sky.toml`), tuned to 401.15 MHz at 48 kHz, over a
//! 600 s window covering the pass that culminates at 35° elevation. Timed
//! with 1 worker thread and with one per core. The samples are discarded
//! (only counted and power-summed), so this is generation cost alone.

use std::ops::ControlFlow;
use std::path::Path;
use std::time::Instant;

use lipika::config::{AntennaConfig, FrontendConfig, ReceiverBandConfig, ReceiverConfig, SkyConfig};
use lipika::pipeline::{render_window, WindowOptions};

/// Window start: 300 s before the culmination at Unix 1717260990.
const START_UNIX: f64 = 1_717_260_690.0;
const DURATION_S: f64 = 600.0;
const SAMPLE_RATE_HZ: f64 = 48_000.0;

fn main() -> anyhow::Result<()> {
    let sky = SkyConfig::from_toml_file(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/hiber_beacon_sky.toml"))?;
    let receiver = ReceiverConfig::stationary_geodetic("glasgow", 55.86, -4.25, 10.0);
    // Isotropic antenna that sees nothing below the horizon.
    let band = ReceiverBandConfig::new("uhf", 401_150_000.0, SAMPLE_RATE_HZ)
        .with_noise_temperature_k(290.0)
        .with_frontend(FrontendConfig {
            antenna: Some(AntennaConfig {
                peak_dbi: 0.0, horizon_dbi: 0.0, fb_db: 0.0, cutoff_deg: Some(0.0),
            }),
            gain_db: None,
            nonlinearity: None,
            adc: None,
        });

    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    println!(
        "render_window: {DURATION_S} s at {SAMPLE_RATE_HZ} Hz, {} transmitter(s), {cores} cores",
        sky.transmitters.len(),
    );
    for (label, n_threads) in [("1 thread", 1usize), ("all cores", 0usize)] {
        let options = WindowOptions { seed: 1, n_threads, ..WindowOptions::default() };
        let mut best = f64::INFINITY;
        let mut power = 0.0f64;
        let mut samples = 0usize;
        for _ in 0..3 {
            power = 0.0;
            let t = Instant::now();
            let summary = render_window(
                &sky, &receiver, &band, START_UNIX, DURATION_S, &options,
                |chunk| {
                    power += chunk.iter().map(|v| v.norm_sqr() as f64).sum::<f64>();
                    ControlFlow::Continue(())
                },
            )?;
            best = best.min(t.elapsed().as_secs_f64());
            samples = summary.n_samples_emitted;
        }
        println!(
            "  {label:>9}: {best:.3} s wall (best of 3) -> {:.0} emulated s per wall s \
             ({:.2} Msample/s); mean power {:.3e} W over {samples} samples",
            DURATION_S / best,
            samples as f64 / best / 1e6,
            power / samples as f64,
        );
    }
    Ok(())
}
