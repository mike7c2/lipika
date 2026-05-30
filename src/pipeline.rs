//! Scene orchestration: turn a `SceneConfig` into per-observer IQ.

use std::sync::Arc;

use anyhow::Result;
use ndarray::Array1;
use num_complex::Complex32;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use crate::channel::{add_awgn, source_contribution_at_observer};
use crate::clock::{NodeClock, NodeClockSpec};
use crate::config::{ObserverConfig, SceneConfig, SourceConfig, SourceKind, TrajectoryConfig};
use crate::geometry::{
    kepler::CircularKeplerian, stationary::Stationary, tle::TleTrajectory, Trajectory,
};
use crate::output::ObserverOutput;
use crate::signal::{bpsk::BpskFrame, ca_code::CaCodeSource, cw::ContinuousWave, SignalSource};

const NEWTON_ITERATIONS: usize = 2;

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

fn build_source(cfg: &SourceConfig) -> Result<Arc<dyn SignalSource>> {
    Ok(match &cfg.kind {
        SourceKind::Cw { baseband_offset_hz } => Arc::new(ContinuousWave {
            carrier_hz: cfg.carrier_hz,
            baseband_offset_hz: *baseband_offset_hz,
        }),
        SourceKind::Bpsk { bit_rate_bps, payload_bits, payload_seed, sync_word } => {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(*payload_seed);
            Arc::new(BpskFrame::random(
                sync_word, *payload_bits, *bit_rate_bps,
                0.0, cfg.carrier_hz, &mut rng,
            ))
        }
        SourceKind::CaCode { prn, epoch_unix } => {
            Arc::new(CaCodeSource::new(cfg.carrier_hz, *prn, *epoch_unix))
        }
    })
}

struct Observer {
    name: String,
    lo_carrier_hz: f64,
    sample_rate_hz: f64,
    snr_db_per_source: f64,
    trajectory: Arc<dyn Trajectory>,
    clock_spec: NodeClockSpec,
}

impl Observer {
    fn from_config(cfg: &ObserverConfig) -> Result<Self> {
        Ok(Self {
            name: cfg.name.clone(),
            lo_carrier_hz: cfg.lo_carrier_hz,
            sample_rate_hz: cfg.sample_rate_hz,
            snr_db_per_source: cfg.snr_db_per_source,
            trajectory: build_trajectory(&cfg.trajectory)?,
            clock_spec: cfg.clock,
        })
    }
}

struct SourceBundle {
    source: Arc<dyn SignalSource>,
    trajectory: Arc<dyn Trajectory>,
}

/// Generate per-observer IQ for an entire scene.
pub fn run_scene(scene: &SceneConfig, rng_seed_base: u64) -> Result<Vec<ObserverOutput>> {
    let observers: Vec<Observer> = scene
        .observers
        .iter()
        .map(Observer::from_config)
        .collect::<Result<Vec<_>>>()?;
    let sources: Vec<SourceBundle> = scene
        .sources
        .iter()
        .map(|sc| {
            Ok(SourceBundle {
                source: build_source(sc)?,
                trajectory: build_trajectory(&sc.trajectory)?,
            })
        })
        .collect::<Result<Vec<_>>>()?;

    // Run observers in parallel across rayon threads. Each observer is
    // independent (its own RNG seed, its own clock realisation).
    let outputs: Result<Vec<ObserverOutput>> = observers
        .par_iter()
        .enumerate()
        .map(|(idx, observer)| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(
                rng_seed_base.wrapping_add(idx as u64)
            );
            generate_one_observer(observer, &sources, scene, &mut rng)
        })
        .collect();
    outputs
}

fn generate_one_observer(
    observer: &Observer,
    sources: &[SourceBundle],
    scene: &SceneConfig,
    rng: &mut Xoshiro256PlusPlus,
) -> Result<ObserverOutput> {
    let n = (observer.sample_rate_hz * scene.scene.duration_s).round() as usize;
    // t_local — observer's nominal sample epochs in its own time frame.
    let t_local = Array1::from_iter(
        (0..n).map(|i| scene.scene.epoch_unix + i as f64 / observer.sample_rate_hz)
    );

    // Clock realisation. Master rate matches sample rate for now; could
    // be raised when the observer is multi-band.
    let clock = NodeClock::new(
        observer.clock_spec, scene.scene.epoch_unix,
        observer.sample_rate_hz, scene.scene.duration_s, rng,
    );
    let delta_t = clock.time_offset_at(&t_local);
    let t_utc = &t_local - &delta_t;
    let observer_position = observer.trajectory.position_ecef_batch(&t_utc);

    let mut accumulator = Array1::<Complex32>::zeros(n);
    for bundle in sources {
        let contribution = source_contribution_at_observer(
            &*bundle.source, &*bundle.trajectory,
            &observer_position, &t_local, &t_utc, &delta_t,
            observer.lo_carrier_hz, NEWTON_ITERATIONS,
        );
        accumulator = accumulator + contribution;
    }

    // AWGN — variance picked so each source sees `snr_db_per_source`.
    // (Composite SNR across N sources is therefore +10 log10(N) higher.)
    add_awgn(&mut accumulator, observer.snr_db_per_source, 1.0, rng);

    Ok(ObserverOutput {
        name: observer.name.clone(),
        epoch_unix: scene.scene.epoch_unix,
        sample_rate_hz: observer.sample_rate_hz,
        lo_carrier_hz: observer.lo_carrier_hz,
        iq: accumulator,
    })
}
