//! Throughput of `lipika::pipeline::render_window` on a sky of many
//! transmitters, with and without the horizon gate.
//!
//!   cargo run --release --example sky_throughput -- <sky.toml> [threads] [duration_s] [--no-gate] [--chunk N] [--profile]
//!
//! One receiver (Glasgow: 55.86 N, 4.25 W, 20 m; 290 K) with an isotropic
//! antenna that sees nothing below the horizon (`cutoff_deg = 0`), tuned to
//! 401.15 MHz at 48 kHz, rendering `duration_s` (default 120) from Unix
//! 1717260690 with `threads` workers (default 1; 0 = one per core). The
//! samples are discarded (only hashed), so this is generation cost alone, and
//! it includes the per-transmitter orbit tabulation done at the start of every
//! render.
//!
//! `--no-gate` renders through the reference path that evaluates every in-band
//! transmitter for every sample. The two print the same hash: the gate does
//! not change the samples.
//!
//! `--profile` adds lipika's per-section time table (summed over threads).
//!
//! Prints the emulated seconds produced per wall-clock second and how many of
//! the sky's transmitters were in band and above the horizon at some point of
//! the window.

use std::ops::ControlFlow;
use std::time::Instant;

use lipika::config::{AntennaConfig, FrontendConfig, ReceiverBandConfig, ReceiverConfig, SkyConfig};
use lipika::pipeline::{
    render_window, render_window_without_horizon_gate, WindowOptions, DEFAULT_WINDOW_CHUNK_SAMPLES,
};

const START_UNIX: f64 = 1_717_260_690.0;
const SAMPLE_RATE_HZ: f64 = 48_000.0;
const CARRIER_HZ: f64 = 401_150_000.0;

fn main() -> anyhow::Result<()> {
    let mut positional: Vec<String> = Vec::new();
    let mut gate = true;
    let mut profile = false;
    let mut chunk_samples = DEFAULT_WINDOW_CHUNK_SAMPLES;
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--no-gate" => gate = false,
            "--profile" => profile = true,
            "--chunk" => {
                chunk_samples = args.next()
                    .ok_or_else(|| anyhow::anyhow!("--chunk needs a sample count"))?.parse()?;
            }
            _ => positional.push(a),
        }
    }
    let Some(sky_path) = positional.first() else {
        anyhow::bail!("usage: sky_throughput <sky.toml> [threads] [duration_s] [--no-gate] [--chunk N] [--profile]");
    };
    let n_threads: usize = positional.get(1).map(|s| s.parse()).transpose()?.unwrap_or(1);
    let duration_s: f64 = positional.get(2).map(|s| s.parse()).transpose()?.unwrap_or(120.0);

    let sky = SkyConfig::from_toml_file(sky_path.as_ref())?;
    let receiver = ReceiverConfig::stationary_geodetic("glasgow", 55.86, -4.25, 20.0);
    let band = ReceiverBandConfig::new("uhf", CARRIER_HZ, SAMPLE_RATE_HZ)
        .with_noise_temperature_k(290.0)
        .with_frontend(FrontendConfig {
            antenna: Some(AntennaConfig {
                peak_dbi: 0.0, horizon_dbi: 0.0, fb_db: 0.0, cutoff_deg: Some(0.0),
            }),
            gain_db: None,
            nonlinearity: None,
            adc: None,
        });
    let options = WindowOptions { seed: 1, chunk_samples, n_threads };

    // FNV-1a over the IQ bytes: equal hashes ⇒ (to overwhelming odds) equal IQ.
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    let sink = |chunk: &[num_complex::Complex32]| {
        for v in chunk {
            for b in v.re.to_le_bytes().into_iter().chain(v.im.to_le_bytes()) {
                hash = (hash ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        ControlFlow::Continue(())
    };
    let t = Instant::now();
    let summary = if gate {
        render_window(&sky, &receiver, &band, START_UNIX, duration_s, &options, sink)?
    } else {
        render_window_without_horizon_gate(
            &sky, &receiver, &band, START_UNIX, duration_s, &options, sink)?
    };
    let wall = t.elapsed().as_secs_f64();

    println!(
        "{} transmitter(s) in the sky, {} in band, {} above the horizon during the window",
        sky.transmitters.len(),
        summary.contributing_transmitters.len(),
        summary.visible_transmitters.len(),
    );
    println!(
        "horizon gate {}, {} thread(s), chunk {chunk_samples}: {duration_s} s at {SAMPLE_RATE_HZ} Hz \
         in {wall:.3} s wall -> {:.1} emulated s per wall s ({:.2} Msample/s); iq fnv1a {hash:016x}",
        if gate { "on" } else { "off" },
        if n_threads == 0 { "all".to_string() } else { n_threads.to_string() },
        duration_s / wall,
        summary.n_samples_emitted as f64 / wall / 1e6,
    );
    if profile {
        lipika::profile::print_report();
    }
    Ok(())
}
