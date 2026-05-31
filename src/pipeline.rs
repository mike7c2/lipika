//! Scene orchestration: turn a `SceneConfig` into per-(receiver, band) IQ.
//!
//! The high-level flow:
//!
//! 1. Build all transmitters (signal model + trajectory + tx power).
//! 2. Build all receivers (trajectory + clock + 1+ bands).
//! 3. For each (receiver, band) pair, run in parallel via rayon:
//!    - Realise the receiver's clock track Δt(t) (shared across this
//!      receiver's bands).
//!    - For each transmitter, evaluate its contribution at this band's
//!      output samples via `channel::transmitter_contribution_at_receiver`.
//!    - Sum across transmitters, add AWGN, package the output.
//! 4. Caller writes each band's IQ to a separate `.npy` file.

use std::sync::Arc;

use anyhow::Result;
use ndarray::Array1;
use num_complex::Complex32;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use crate::channel::{add_awgn_at_power, thermal_noise_power_w, transmitter_contribution_at_receiver};
use crate::clock::{NodeClock, NodeClockSpec};
use crate::config::{
    ReceiverBandConfig, ReceiverConfig, SceneConfig, SignalKind, TrajectoryConfig,
    TransmitterConfig,
};
use crate::geometry::{
    kepler::CircularKeplerian, precomputed::PrecomputedTrajectory, stationary::Stationary,
    tle::TleTrajectory, Trajectory,
};
use crate::output::ObserverOutput;
use crate::signal::{bpsk::BpskFrame, ca_code::CaCodeSource, cw::ContinuousWave, SignalModel};

const NEWTON_ITERATIONS: usize = 2;

/// Spacing of the precomputed-trajectory grid for transmitters. 10 ms is
/// well inside the orbit-arc-vs-chord error budget (sub-µm for LEO/GPS
/// over 250 s passes) and keeps the per-source ephemeris table small
/// (250 s × 100 Hz × 3 × 8 bytes ≈ 60 kB per source).
const TX_TRAJECTORY_GRID_DT_S: f64 = 0.01;

fn build_trajectory(cfg: &TrajectoryConfig) -> Result<Arc<dyn Trajectory>> {
    Ok(match cfg {
        TrajectoryConfig::StationaryGeodetic { lat_deg, lon_deg, h_m } => {
            Arc::new(Stationary::from_geodetic(*lat_deg, *lon_deg, *h_m))
        }
        TrajectoryConfig::StationaryEcef { x_m, y_m, z_m } => {
            Arc::new(Stationary { position_ecef: [*x_m, *y_m, *z_m] })
        }
        TrajectoryConfig::Tle { line1, line2 } => {
            Arc::new(TleTrajectory::from_tle(line1, line2)?)
        }
        TrajectoryConfig::CircularKeplerian {
            semi_major_axis_m, inclination_deg, raan_deg,
            true_anomaly_deg_at_epoch, epoch_unix,
        } => Arc::new(CircularKeplerian {
            semi_major_axis_m: *semi_major_axis_m,
            inclination_deg: *inclination_deg,
            raan_deg: *raan_deg,
            true_anomaly_deg_at_epoch: *true_anomaly_deg_at_epoch,
            epoch_unix: *epoch_unix,
        }),
    })
}

fn build_signal(cfg: &TransmitterConfig) -> Result<Arc<dyn SignalModel>> {
    Ok(match &cfg.kind {
        SignalKind::Cw { baseband_offset_hz } => Arc::new(ContinuousWave {
            carrier_hz: cfg.carrier_hz,
            baseband_offset_hz: *baseband_offset_hz,
        }),
        SignalKind::Bpsk { bit_rate_bps, payload_bits, payload_seed, sync_word } => {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(*payload_seed);
            Arc::new(BpskFrame::random(
                sync_word, *payload_bits, *bit_rate_bps,
                0.0, cfg.carrier_hz, &mut rng,
            ))
        }
        SignalKind::CaCode { prn, epoch_unix } => {
            Arc::new(CaCodeSource::new(cfg.carrier_hz, *prn, *epoch_unix))
        }
    })
}

struct Transmitter {
    #[allow(dead_code)] // for logging / future use
    name: String,
    signal: Arc<dyn SignalModel>,
    trajectory: Arc<dyn Trajectory>,
    tx_power_w: f64,
}

impl Transmitter {
    fn from_config(cfg: &TransmitterConfig, scene: &SceneConfig) -> Result<Self> {
        // Wrap the raw trajectory in a precomputed-and-interpolated decorator
        // so the per-output-sample Newton iterations only pay an array
        // lookup + interpolation per query, not a fresh SGP4 / Kepler eval.
        let raw = build_trajectory(&cfg.trajectory)?;
        let precomputed: Arc<dyn Trajectory> = match &*cfg.trajectory.kind_name() {
            // Stationary: precompute is overhead for no gain; pass through.
            "stationary_geodetic" | "stationary_ecef" => raw,
            _ => Arc::new(PrecomputedTrajectory::sample(
                &*raw,
                scene.scene.epoch_unix,
                scene.scene.duration_s,
                TX_TRAJECTORY_GRID_DT_S,
            )),
        };
        Ok(Self {
            name: cfg.name.clone(),
            signal: build_signal(cfg)?,
            trajectory: precomputed,
            tx_power_w: 10f64.powf(cfg.tx_power_dbw / 10.0),
        })
    }
}

struct Receiver {
    name: String,
    trajectory: Arc<dyn Trajectory>,
    clock_spec: NodeClockSpec,
    bands: Vec<ReceiverBandConfig>,
}

impl Receiver {
    fn from_config(cfg: &ReceiverConfig) -> Result<Self> {
        if cfg.bands.is_empty() {
            anyhow::bail!("receiver '{}' has no bands; at least one [[receiver.band]] required", cfg.name);
        }
        Ok(Self {
            name: cfg.name.clone(),
            trajectory: build_trajectory(&cfg.trajectory)?,
            clock_spec: cfg.clock,
            bands: cfg.bands.clone(),
        })
    }
}

/// Generate per-(receiver, band) IQ for an entire scene. Output naming:
/// `<receiver_name>__<band_id>` so multi-band receivers get distinct files.
pub fn run_scene(scene: &SceneConfig, rng_seed_base: u64) -> Result<Vec<ObserverOutput>> {
    let transmitters: Vec<Transmitter> = scene
        .transmitters
        .iter()
        .map(|cfg| Transmitter::from_config(cfg, scene))
        .collect::<Result<Vec<_>>>()?;
    let receivers: Vec<Receiver> = scene
        .receivers
        .iter()
        .map(Receiver::from_config)
        .collect::<Result<Vec<_>>>()?;

    // Enumerate all (receiver, band) work items.
    let mut work: Vec<(usize, usize)> = Vec::new();
    for (r_idx, rx) in receivers.iter().enumerate() {
        for (b_idx, _) in rx.bands.iter().enumerate() {
            work.push((r_idx, b_idx));
        }
    }

    let outputs: Result<Vec<ObserverOutput>> = work
        .par_iter()
        .enumerate()
        .map(|(work_idx, &(r_idx, b_idx))| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                rng_seed_base.wrapping_add(work_idx as u64)
            );
            generate_one_band(&receivers[r_idx], b_idx, &transmitters, scene, &mut rng)
        })
        .collect();
    outputs
}

fn generate_one_band(
    receiver: &Receiver,
    band_idx: usize,
    transmitters: &[Transmitter],
    scene: &SceneConfig,
    rng: &mut Xoshiro256PlusPlus,
) -> Result<ObserverOutput> {
    let band = &receiver.bands[band_idx];
    let n = (band.sample_rate_hz * scene.scene.duration_s).round() as usize;
    let t_local = Array1::from_iter(
        (0..n).map(|i| scene.scene.epoch_unix + i as f64 / band.sample_rate_hz)
    );

    // The receiver's clock realisation. Master rate = this band's sample rate.
    // (For tighter inter-band coherence in real DGSN hardware, all bands
    // should pull from the same underlying clock track; we approximate by
    // re-realising at each band's rate. The deterministic terms — c, m·t —
    // are identical between bands; only the white-FM ε(t) realisation
    // differs in its sample grid. For Allan deviations < 1e-9 this gap is
    // numerically negligible.)
    let clock = NodeClock::new(
        receiver.clock_spec, scene.scene.epoch_unix,
        band.sample_rate_hz, scene.scene.duration_s, rng,
    );
    let delta_t = clock.time_offset_at(&t_local);
    let t_utc = &t_local - &delta_t;
    let receiver_position = receiver.trajectory.position_ecef_batch(&t_utc);

    let mut accumulator = Array1::<Complex32>::zeros(n);
    for tx in transmitters {
        let contribution = transmitter_contribution_at_receiver(
            &*tx.signal, &*tx.trajectory,
            &receiver_position, &t_local, &t_utc, &delta_t,
            band.lo_carrier_hz, tx.tx_power_w, NEWTON_ITERATIONS,
        );
        accumulator = accumulator + contribution;
    }

    let noise_power_w = band.noise_power_w.unwrap_or_else(||
        thermal_noise_power_w(band.noise_temperature_k, band.sample_rate_hz)
    );
    add_awgn_at_power(&mut accumulator, noise_power_w, rng);

    let composite_name = format!("{}__{}", receiver.name, band.band_id);
    Ok(ObserverOutput {
        name: composite_name,
        epoch_unix: scene.scene.epoch_unix,
        sample_rate_hz: band.sample_rate_hz,
        lo_carrier_hz: band.lo_carrier_hz,
        iq: accumulator,
    })
}
