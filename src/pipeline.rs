//! Scene orchestration: turn a `SceneConfig` into per-(receiver, band) IQ.
//!
//! The high-level flow:
//!
//! 1. Build all transmitters (signal model + trajectory + tx power) and
//!    receivers (trajectory + clock + 1+ bands).
//! 2. For each (receiver, band), generate the IQ in sequential time-chunks:
//!    - A forward-only [`ClockStream`] (one per receiver, shared across its
//!      bands) yields Δt(t) for the chunk with O(1) memory.
//!    - The overlapping transmitters are evaluated **in parallel** within the
//!      chunk, then reduced in fixed order and AWGN added.
//!    - The chunk is pushed to a sink (in-memory or stream-to-disk).
//!
//! [`run_scene`] keeps everything in RAM; [`run_scene_to_dir`] streams each
//! band to a `.iq.c64` file with memory independent of duration.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::Result;
use ndarray::Array1;
use num_complex::Complex32;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

use crate::channel::{
    add_awgn_at_power, solve_emission_times, thermal_noise_power_w,
    transmitter_contribution_at_receiver,
};
use crate::SPEED_OF_LIGHT_M_S;
use crate::profile;
use crate::clock::{ClockStream, NodeClockSpec};
use crate::config::{
    NavMode, ReceiverBandConfig, ReceiverConfig, SceneConfig, SignalKind, TrajectoryConfig,
    TransmitterConfig,
};
use crate::geometry::{
    kepler::CircularKeplerian, precomputed::PrecomputedTrajectory, stationary::Stationary,
    tle::TleTrajectory, AnyTrajectory, Trajectory,
};
use crate::output::{ObserverOutput, StreamingC64Writer};
use crate::signal::{
    bpsk::BpskFrame, ca_code::CaCodeSource, cw::ContinuousWave,
    gps_nav::{lnav_bits, parse_rinex3_gps_nav, Ephemeris, GpsBroadcastTrajectory, GpsTime},
    SignalModel,
};

const NEWTON_ITERATIONS: usize = 2;

/// Spacing of the precomputed-trajectory grid for transmitters. 10 ms is
/// well inside the orbit-arc-vs-chord error budget (sub-µm for LEO/GPS
/// over 250 s passes) and keeps the per-source ephemeris table small
/// (250 s × 100 Hz × 3 × 8 bytes ≈ 60 kB per source).
const TX_TRAJECTORY_GRID_DT_S: f64 = 0.01;

/// Does a transmitter's occupied spectrum intersect a receiver band's
/// acceptance window? Both are expressed as a centre frequency and a total
/// bandwidth; they overlap iff the gap between centres is less than the
/// sum of the half-bandwidths. A zero-bandwidth (CW) transmitter still
/// overlaps a band that straddles its tone.
fn bands_overlap(
    tx_carrier_hz: f64, tx_bandwidth_hz: f64,
    band_lo_hz: f64, band_bandwidth_hz: f64,
) -> bool {
    (tx_carrier_hz - band_lo_hz).abs() <= 0.5 * (tx_bandwidth_hz + band_bandwidth_hz)
}

fn build_trajectory(
    cfg: &TrajectoryConfig,
    scene_epoch_unix: f64,
    gps_utc_offset_s: f64,
) -> Result<Arc<dyn Trajectory>> {
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
        TrajectoryConfig::RinexNav { file, prn } => {
            let eph = load_rinex_ephemeris(file, *prn, scene_epoch_unix, gps_utc_offset_s)?;
            Arc::new(GpsBroadcastTrajectory::new(eph, gps_utc_offset_s))
        }
    })
}

/// Read a RINEX nav file and select the GPS ephemeris for `prn` nearest the
/// scene epoch. Used for both the trajectory and (for full NAV) the broadcast
/// message, so the two are the same ephemeris.
fn load_rinex_ephemeris(
    file: &str,
    prn: usize,
    scene_epoch_unix: f64,
    gps_utc_offset_s: f64,
) -> Result<Ephemeris> {
    let text = std::fs::read_to_string(file)
        .map_err(|e| anyhow::anyhow!("reading RINEX nav '{file}': {e}"))?;
    parse_rinex3_gps_nav(&text, prn as u8, scene_epoch_unix, gps_utc_offset_s)
        .ok_or_else(|| anyhow::anyhow!("no GPS PRN {prn} ephemeris found in '{file}'"))
}

fn build_signal(
    cfg: &TransmitterConfig,
    raw_traj: &dyn Trajectory,
    scene: &SceneConfig,
) -> Result<Arc<dyn SignalModel>> {
    Ok(match &cfg.kind {
        SignalKind::Cw { baseband_offset_hz } => Arc::new(ContinuousWave {
            carrier_hz: cfg.carrier_hz,
            baseband_offset_hz: *baseband_offset_hz,
        }),
        SignalKind::Bpsk { bit_rate_bps, payload_bits, payload_seed, sync_word } => {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(*payload_seed);
            // Bit 0 is emitted at the scene epoch ⇒ epoch-relative tag 0.0.
            Arc::new(BpskFrame::random(
                sync_word, *payload_bits, *bit_rate_bps,
                0.0, cfg.carrier_hz, &mut rng,
            ))
        }
        SignalKind::CaCode { prn, epoch_unix, nav, nav_seed } => {
            // Signals are evaluated at epoch-relative emission times, so the
            // chip-0 tag must be relative to the scene epoch too.
            let epoch_rel = *epoch_unix - scene.scene.epoch_unix;
            match nav {
                NavMode::Off =>
                    Arc::new(CaCodeSource::new(cfg.carrier_hz, *prn, epoch_rel)),
                NavMode::Dummy =>
                    Arc::new(CaCodeSource::with_nav_data(
                        cfg.carrier_hz, *prn, epoch_rel, *nav_seed)),
                NavMode::Full => {
                    let bits = build_lnav_bits(cfg, raw_traj, scene, *prn)?;
                    Arc::new(CaCodeSource::with_full_nav(
                        cfg.carrier_hz, *prn, epoch_rel, bits))
                }
            }
        }
    })
}

/// Build the LNAV bit stream for the whole scene from the satellite's
/// ephemeris. The ephemeris is the broadcast record itself for `rinex_nav`
/// (exact — same object drives the orbit), exact closed-form for
/// `circular_keplerian`, and osculating-at-`toe` for any other trajectory.
fn build_lnav_bits(
    cfg: &TransmitterConfig,
    raw_traj: &dyn Trajectory,
    scene: &SceneConfig,
    prn: usize,
) -> Result<Vec<i8>> {
    let epoch = scene.scene.epoch_unix;
    let gps = GpsTime::from_unix(epoch, scene.scene.gps_utc_offset_s);
    if (gps.sow / 6.0).fract() != 0.0 {
        log::warn!(
            "scene epoch is not on a 6 s GPS boundary (sow={}); full-NAV \
             subframe edges won't align to the TOW Z-count and SV-position \
             time will be biased. Use an epoch_unix divisible by 6.",
            gps.sow,
        );
    }
    let eph = if let TrajectoryConfig::CircularKeplerian {
        semi_major_axis_m, inclination_deg, raan_deg,
        true_anomaly_deg_at_epoch, epoch_unix,
    } = &cfg.trajectory {
        let orbit = CircularKeplerian {
            semi_major_axis_m: *semi_major_axis_m,
            inclination_deg: *inclination_deg,
            raan_deg: *raan_deg,
            true_anomaly_deg_at_epoch: *true_anomaly_deg_at_epoch,
            epoch_unix: *epoch_unix,
        };
        Ephemeris::from_circular_keplerian(&orbit, prn as u8, epoch, gps)
    } else if let TrajectoryConfig::RinexNav { file, prn: rprn } = &cfg.trajectory {
        // Use the broadcast record directly — no derivation, so the
        // transmitted ephemeris is identical to the orbit we fly.
        load_rinex_ephemeris(file, *rprn, epoch, scene.scene.gps_utc_offset_s)?
    } else {
        Ephemeris::from_trajectory(raw_traj, prn as u8, epoch, gps)
    };
    Ok(lnav_bits(&eph, gps, scene.scene.duration_s))
}

struct Transmitter {
    #[allow(dead_code)] // used for the output descriptor (contributing-tx list)
    name: String,
    signal: Arc<dyn SignalModel>,
    /// Wrapped into an `AnyTrajectory` enum-dispatch instead of
    /// `Arc<dyn Trajectory>` so the per-sample Newton-loop call to
    /// `position_ecef` can be inlined.
    trajectory: AnyTrajectory,
    tx_power_w: f64,
    /// Occupied RF bandwidth (Hz), centred on `signal.carrier_hz()`. Used
    /// to decide whether this transmitter falls within a receiver band.
    bandwidth_hz: f64,
}

impl Transmitter {
    fn from_config(cfg: &TransmitterConfig, scene: &SceneConfig) -> Result<Self> {
        // Build raw trajectory, then wrap into the right AnyTrajectory
        // variant. Non-stationary ones go through PrecomputedTrajectory
        // so Newton-loop queries are just two array reads + a lerp.
        let raw = build_trajectory(
            &cfg.trajectory, scene.scene.epoch_unix, scene.scene.gps_utc_offset_s)?;
        let trajectory = match cfg.trajectory.kind_name() {
            "stationary_geodetic" | "stationary_ecef" => {
                // The Stationary impl is sub-µs to clone, no point
                // wrapping in PrecomputedTrajectory.
                let pos = raw.position_ecef(scene.scene.epoch_unix);
                AnyTrajectory::Stationary(Stationary { position_ecef: pos })
            }
            _ => profile::time("trajectory:precompute", || {
                AnyTrajectory::Precomputed(PrecomputedTrajectory::sample(
                    &*raw,
                    scene.scene.epoch_unix,
                    scene.scene.duration_s,
                    TX_TRAJECTORY_GRID_DT_S,
                ))
            }),
        };
        Ok(Self {
            name: cfg.name.clone(),
            signal: build_signal(cfg, &*raw, scene)?,
            trajectory,
            tx_power_w: 10f64.powf(cfg.tx_power_dbw / 10.0),
            bandwidth_hz: cfg.effective_bandwidth_hz(),
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
    fn from_config(cfg: &ReceiverConfig, scene: &SceneConfig) -> Result<Self> {
        if cfg.bands.is_empty() {
            anyhow::bail!("receiver '{}' has no bands; at least one [[receiver.band]] required", cfg.name);
        }
        Ok(Self {
            name: cfg.name.clone(),
            trajectory: build_trajectory(
                &cfg.trajectory, scene.scene.epoch_unix, scene.scene.gps_utc_offset_s)?,
            clock_spec: cfg.clock,
            bands: cfg.bands.clone(),
        })
    }
}

/// Default streaming chunk size (samples). ~1 Mi-sample ≈ 16 MB of c64.
pub const DEFAULT_CHUNK_SAMPLES: usize = 1 << 20;

/// A consumer of generated IQ chunks. The chunked generator drives one of
/// these so the in-memory ([`MemSink`]) and stream-to-disk ([`FileSink`])
/// paths share a single code path (and therefore produce identical bytes).
trait IqSink {
    fn push(&mut self, chunk: &Array1<Complex32>) -> Result<()>;
}

/// Collects chunks into one contiguous in-memory buffer (for `run_scene`).
struct MemSink {
    buf: Vec<Complex32>,
}
impl MemSink {
    fn with_capacity(n: usize) -> Self { Self { buf: Vec::with_capacity(n) } }
}
impl IqSink for MemSink {
    fn push(&mut self, chunk: &Array1<Complex32>) -> Result<()> {
        self.buf.extend_from_slice(chunk.as_slice().expect("contiguous chunk"));
        Ok(())
    }
}

/// Streams chunks straight to a `.iq.c64` file via [`StreamingC64Writer`].
struct FileSink {
    writer: StreamingC64Writer,
}
impl IqSink for FileSink {
    fn push(&mut self, chunk: &Array1<Complex32>) -> Result<()> {
        self.writer.push_chunk(chunk)
    }
}

/// Metadata describing one generated band, returned by the chunked generator.
struct BandMeta {
    name: String,
    epoch_unix: f64,
    sample_rate_hz: f64,
    lo_carrier_hz: f64,
    n_samples: usize,
    descriptor: serde_json::Value,
}

/// Lightweight result of a streamed band — no IQ retained.
#[derive(Debug, Clone)]
pub struct OutputManifest {
    pub name: String,
    pub path: PathBuf,
    pub n_samples: usize,
    pub descriptor: serde_json::Value,
}

/// Generate per-(receiver, band) IQ for an entire scene **in memory**. Output
/// naming: `<receiver_name>__<band_id>`. Convenience for tests/benches and
/// callers that want the samples in RAM; for large scenes prefer
/// [`run_scene_to_dir`], which streams to disk with bounded memory.
pub fn run_scene(scene: &SceneConfig, rng_seed_base: u64) -> Result<Vec<ObserverOutput>> {
    let (transmitters, receivers, work) = prepare_scene(scene)?;
    let mut outputs = Vec::with_capacity(work.len());
    for (work_idx, (r_idx, b_idx)) in work.into_iter().enumerate() {
        let receiver = &receivers[r_idx];
        let n = (receiver.bands[b_idx].sample_rate_hz * scene.scene.duration_s).round() as usize;
        let mut sink = MemSink::with_capacity(n);
        let meta = generate_band(
            receiver, r_idx, b_idx, &transmitters, scene,
            rng_seed_base, work_idx, DEFAULT_CHUNK_SAMPLES, &mut sink,
        )?;
        outputs.push(ObserverOutput {
            name: meta.name,
            epoch_unix: meta.epoch_unix,
            sample_rate_hz: meta.sample_rate_hz,
            lo_carrier_hz: meta.lo_carrier_hz,
            iq: Array1::from(sink.buf),
            descriptor: meta.descriptor,
        });
    }
    Ok(outputs)
}

/// Generate a scene **streaming to disk**: each band's IQ is written
/// incrementally to `<dir>/<name>.iq.c64` (+ JSON sidecar) in `chunk_samples`
/// chunks, so peak memory is independent of scene duration. Within each chunk
/// the (expensive) per-transmitter work is parallelised across `n_threads`
/// (0 ⇒ all cores); chunks stay sequential to keep the clock walk and noise
/// stream — and hence the output — bit-identical regardless of chunk size.
pub fn run_scene_to_dir(
    scene: &SceneConfig,
    rng_seed_base: u64,
    dir: &Path,
    chunk_samples: usize,
    n_threads: usize,
) -> Result<Vec<OutputManifest>> {
    let chunk_samples = chunk_samples.max(1);
    let (transmitters, receivers, work) = prepare_scene(scene)?;

    let pool = {
        let mut b = rayon::ThreadPoolBuilder::new();
        if n_threads > 0 { b = b.num_threads(n_threads); }
        b.build()?
    };

    // Parallelise across the independent (receiver, band) work units — each owns
    // its file and a `work_idx`-seeded RNG, so its bytes depend only on the unit,
    // not on order or concurrency (bit-exactness preserved; manifest order kept by
    // the indexed `collect`). This is the dominant parallel axis: there are many
    // units (receivers × bands) but only a handful of transmitters per band, so
    // the inner transmitter `par_iter` alone leaves most cores idle. Nested rayon
    // shares this pool, so the two levels load-balance via work-stealing.
    pool.install(|| {
        work.into_par_iter()
            .enumerate()
            .map(|(work_idx, (r_idx, b_idx))| -> Result<OutputManifest> {
                let receiver = &receivers[r_idx];
                let name = format!("{}__{}", receiver.name, receiver.bands[b_idx].band_id);
                let mut sink = FileSink { writer: StreamingC64Writer::create(dir, &name)? };
                let meta = generate_band(
                    receiver, r_idx, b_idx, &transmitters, scene,
                    rng_seed_base, work_idx, chunk_samples, &mut sink,
                )?;
                sink.writer.finalize(
                    meta.epoch_unix, meta.sample_rate_hz, meta.lo_carrier_hz, &meta.descriptor,
                )?;
                Ok(OutputManifest {
                    name: meta.name,
                    path: dir.join(format!("{name}.iq.c64")),
                    n_samples: meta.n_samples,
                    descriptor: meta.descriptor,
                })
            })
            .collect::<Result<Vec<_>>>()
    })
}

/// Export geometric truth alongside the IQ. For every (receiver, transmitter)
/// pair this writes the true light-time delay — from the *same* Sagnac-correct
/// emission solve the signal uses — plus slant range and range-rate, sampled on
/// a coarse grid at `truth_hz`. One `<dir>/<receiver>__truth.json` per receiver.
///
/// The grid is treated as ideal UTC (no receiver clock), so this is the
/// clock-free geometry baked into the signal. An analysis can therefore predict
/// each signal's arrival and Doppler without re-deriving orbits, the Sagnac
/// term or leap seconds — the residual it then measures is purely the receiver
/// clock (time transfer) or the estimator error (Doppler / TDOA).
pub fn export_truth(scene: &SceneConfig, dir: &Path, truth_hz: f64) -> Result<()> {
    if truth_hz <= 0.0 {
        return Ok(());
    }
    std::fs::create_dir_all(dir)?;
    let (transmitters, receivers, _work) = prepare_scene(scene)?;
    let epoch = scene.scene.epoch_unix;
    let dur = scene.scene.duration_s;
    let ng = ((dur * truth_hz).round() as usize).max(1) + 1;
    let dt = 1.0 / truth_hz;
    // Epoch-relative grid for the geometry solve; absolute only for the JSON.
    let t_grid_rel = Array1::from_iter((0..ng).map(|k| k as f64 * dt));
    let t_grid_abs = t_grid_rel.mapv(|x| epoch + x);

    for receiver in &receivers {
        let rx_pos = receiver.trajectory.position_ecef_batch(&t_grid_abs);
        let mut txs = Vec::new();
        for (i, tx) in transmitters.iter().enumerate() {
            // Extra Newton iterations vs the per-sample path: this is a coarse
            // one-off grid, so the cost is negligible and convergence certain.
            let (_t_emit, ranges) = solve_emission_times(
                &tx.trajectory, &rx_pos, &t_grid_rel, epoch, NEWTON_ITERATIONS.max(4),
            );
            let delay: Vec<f64> = ranges.iter().map(|&r| r / SPEED_OF_LIGHT_M_S).collect();
            let n = ranges.len();
            let mut rr = vec![0.0_f64; n];
            if n >= 2 {
                for k in 0..n {
                    rr[k] = if k == 0 {
                        (ranges[1] - ranges[0]) / dt
                    } else if k == n - 1 {
                        (ranges[n - 1] - ranges[n - 2]) / dt
                    } else {
                        (ranges[k + 1] - ranges[k - 1]) / (2.0 * dt)
                    };
                }
            }
            // True source ECEF on the grid (absolute time) — an exact trajectory
            // reference for moving emitters, where range-only reconstruction from
            // a tight receiver array is ill-conditioned in the radial direction.
            let sp = tx.trajectory.position_ecef_batch(&t_grid_abs);
            let mut d = transmitter_descriptor(&scene.transmitters[i]);
            d["range_m"] = serde_json::json!(ranges.to_vec());
            d["delay_s"] = serde_json::json!(delay);
            d["range_rate_m_s"] = serde_json::json!(rr);
            d["source_x_m"] = serde_json::json!(sp.row(0).to_vec());
            d["source_y_m"] = serde_json::json!(sp.row(1).to_vec());
            d["source_z_m"] = serde_json::json!(sp.row(2).to_vec());
            txs.push(d);
        }
        let truth = serde_json::json!({
            "receiver": receiver.name,
            "epoch_unix": epoch,
            "duration_s": dur,
            "grid_hz": truth_hz,
            "t_utc": t_grid_abs.to_vec(),
            "transmitters": txs,
        });
        std::fs::write(
            dir.join(format!("{}__truth.json", receiver.name)),
            serde_json::to_string(&truth)?,
        )?;
    }
    Ok(())
}

/// Build transmitters + receivers and enumerate (receiver, band) work items.
type PreparedScene = (Vec<Transmitter>, Vec<Receiver>, Vec<(usize, usize)>);
fn prepare_scene(scene: &SceneConfig) -> Result<PreparedScene> {
    let transmitters: Vec<Transmitter> = scene
        .transmitters
        .iter()
        .map(|cfg| Transmitter::from_config(cfg, scene))
        .collect::<Result<Vec<_>>>()?;
    let receivers: Vec<Receiver> = scene
        .receivers
        .iter()
        .map(|cfg| Receiver::from_config(cfg, scene))
        .collect::<Result<Vec<_>>>()?;
    let mut work: Vec<(usize, usize)> = Vec::new();
    for (r_idx, rx) in receivers.iter().enumerate() {
        for (b_idx, _) in rx.bands.iter().enumerate() {
            work.push((r_idx, b_idx));
        }
    }
    Ok((transmitters, receivers, work))
}

/// Clock RNG is seeded per *receiver* (shared by all its bands, so they are
/// mutually coherent); noise RNG is seeded per *(receiver, band)*. Keeping them
/// in separate sub-streams makes the result independent of chunking and of the
/// transmitter count.
fn clock_seed(base: u64, r_idx: usize) -> u64 {
    base.wrapping_add(0xC10C_0000_0000_0000).wrapping_add(r_idx as u64)
}
fn noise_seed(base: u64, work_idx: usize) -> u64 {
    base.wrapping_add(work_idx as u64)
}

/// Chunked, transmitter-parallel generation of one band, driving `sink`.
/// **`SIGNAL_MODEL.md` E5 (superposition) + E10 (full forward: Σ + noise + the
/// `RxChain` post-sum operators + ADC).**
fn generate_band(
    receiver: &Receiver,
    r_idx: usize,
    band_idx: usize,
    transmitters: &[Transmitter],
    scene: &SceneConfig,
    rng_seed_base: u64,
    work_idx: usize,
    chunk_samples: usize,
    sink: &mut dyn IqSink,
) -> Result<BandMeta> {
    let band = &receiver.bands[band_idx];
    let epoch = scene.scene.epoch_unix;
    let fs = band.sample_rate_hz;
    let n = (fs * scene.scene.duration_s).round() as usize;
    let lo = band.lo_carrier_hz;
    let band_bandwidth_hz = band.effective_bandwidth_hz();

    // One forward-only clock per receiver, realised on the receiver's highest
    // band rate so every band references the identical ε(t) walk (the high-rate
    // band lands on the grid exactly; lower-rate bands interpolate the same
    // walk). Seeded per receiver ⇒ coherent across bands.
    let master_rate = receiver.bands.iter()
        .map(|b| b.sample_rate_hz)
        .fold(0.0_f64, f64::max);
    // The clock runs in epoch-relative time (start = 0.0): all per-sample
    // timing is kept at small magnitudes so sub-ns offsets survive f64.
    let mut clock = ClockStream::new(
        receiver.clock_spec, 0.0, master_rate,
        Xoshiro256PlusPlus::seed_from_u64(clock_seed(rng_seed_base, r_idx)),
    );
    let mut noise_rng = Xoshiro256PlusPlus::seed_from_u64(noise_seed(rng_seed_base, work_idx));
    let noise_power_w = band.noise_power_w.unwrap_or_else(||
        thermal_noise_power_w(band.noise_temperature_k, fs)
    );

    // Band-overlap gate, evaluated once: which transmitters land in this band.
    // Skipping the rest stops a far-out-of-band emitter (e.g. a 437 MHz beacon)
    // aliasing into the band. Partial overlap is admitted unfiltered (no
    // anti-alias) — documented limitation.
    let mut active: Vec<usize> = Vec::new();
    let mut contributing: Vec<serde_json::Value> = Vec::new();
    for (i, tx) in transmitters.iter().enumerate() {
        if bands_overlap(tx.signal.carrier_hz(), tx.bandwidth_hz, lo, band_bandwidth_hz) {
            active.push(i);
            contributing.push(transmitter_descriptor(&scene.transmitters[i]));
        }
    }

    // Receiver front-end operator chain (built once; identity if unconfigured).
    let rx_chain = crate::rx::RxChain::from_config(band.frontend.as_ref());

    // Sequential chunk loop (clock walk + noise stream + file order); the
    // expensive per-transmitter work is parallelised within each chunk.
    let mut i0 = 0usize;
    while i0 < n {
        let i1 = (i0 + chunk_samples).min(n);
        let clen = i1 - i0;

        // Epoch-relative sample times from the GLOBAL index (never accumulate
        // time, and keep magnitudes small so sub-ns terms survive f64). The
        // global index keeps these identical regardless of chunking.
        let t_local_rel = Array1::from_iter((i0..i1).map(|i| i as f64 / fs));
        let delta_t = profile::time("clock:stream", || {
            let mut d = Array1::<f64>::zeros(clen);
            for (j, &t) in t_local_rel.iter().enumerate() {
                d[j] = clock.delta_t_at(t);
            }
            d
        });
        let t_utc_rel = &t_local_rel - &delta_t;
        // Trajectory queries need absolute time; a 238 ns rounding here is
        // sub-mm of orbit motion, so reconstruct the absolute UTC grid only here.
        let receiver_position = profile::time("receiver:position_batch",
            || receiver.trajectory.position_ecef_batch(&t_utc_rel.mapv(|x| epoch + x)));

        // Parallel over the overlapping transmitters → one contribution buffer
        // each; reduce in fixed transmitter order (deterministic f32 sum).
        let mut acc = profile::time("tx_contribution:total", || {
            let contributions: Vec<Array1<Complex32>> = active
                .par_iter()
                .map(|&i| {
                    let tx = &transmitters[i];
                    transmitter_contribution_at_receiver(
                        &*tx.signal, &tx.trajectory,
                        &receiver_position, &t_local_rel, &t_utc_rel, &delta_t,
                        lo, tx.tx_power_w, epoch, NEWTON_ITERATIONS,
                        rx_chain.antenna.as_ref(),
                    )
                })
                .collect();
            let mut acc = Array1::<Complex32>::zeros(clen);
            for c in &contributions {
                acc += c;
            }
            acc
        });

        profile::time("noise:awgn", || add_awgn_at_power(&mut acc, noise_power_w, &mut noise_rng));
        // Receiver front-end: post-summation operator chain (gain, ADC collapse …).
        // Identity when no `[receiver.band.frontend]` is configured. Deterministic
        // per sample ⇒ chunk/thread invariance preserved.
        profile::time("rx:frontend", || rx_chain.apply_post_sum(&mut acc));
        sink.push(&acc)?;
        i0 = i1;
    }

    let descriptor = serde_json::json!({
        "receiver": receiver.name,
        "band_id": band.band_id,
        "duration_s": scene.scene.duration_s,
        "band": {
            "lo_carrier_hz": lo,
            "sample_rate_hz": fs,
            "bandwidth_hz": band_bandwidth_hz,
            "noise_temperature_k": band.noise_temperature_k,
            "noise_power_w": noise_power_w,
        },
        "contributing_transmitters": contributing,
        "scene_config": serde_json::to_value(scene).unwrap_or(serde_json::Value::Null),
    });

    Ok(BandMeta {
        name: format!("{}__{}", receiver.name, band.band_id),
        epoch_unix: epoch,
        sample_rate_hz: fs,
        lo_carrier_hz: lo,
        n_samples: n,
        descriptor,
    })
}

/// Per-transmitter detail captured in the output descriptor: enough to know
/// exactly which emitters (and, for GNSS, which PRNs) landed in a band.
fn transmitter_descriptor(cfg: &TransmitterConfig) -> serde_json::Value {
    let mut v = serde_json::json!({
        "name": cfg.name,
        "kind": cfg.kind.kind_name(),
        "carrier_hz": cfg.carrier_hz,
        "tx_power_dbw": cfg.tx_power_dbw,
        "bandwidth_hz": cfg.effective_bandwidth_hz(),
        "trajectory_kind": cfg.trajectory.kind_name(),
    });
    if let SignalKind::CaCode { prn, nav, .. } = &cfg.kind {
        v["prn"] = serde_json::json!(prn);
        v["nav"] = serde_json::json!(format!("{nav:?}").to_lowercase());
    }
    v
}

#[cfg(test)]
mod tests {
    use super::bands_overlap;

    // GNSS L1 band: 1575.42 MHz LO, 2.048 MHz acceptance window.
    const GNSS_LO: f64 = 1_575_420_000.0;
    const GNSS_BW: f64 = 2_048_000.0;

    #[test]
    fn far_out_of_band_beacon_is_excluded() {
        // The 437 MHz ISS UHF beacon must NOT contribute to the GNSS band —
        // this is the aliasing bug the gate fixes.
        let beacon_bw = 2.0 * 9600.0; // BPSK @ 9.6 kbps
        assert!(!bands_overlap(437_000_000.0, beacon_bw, GNSS_LO, GNSS_BW));
    }

    #[test]
    fn co_carrier_ca_code_is_included() {
        // A C/A transmitter at the L1 carrier overlaps the GNSS band.
        let ca_bw = 2.0 * 1_023_000.0;
        assert!(bands_overlap(GNSS_LO, ca_bw, GNSS_LO, GNSS_BW));
    }

    #[test]
    fn edge_overlap_is_inclusive() {
        // A tone sitting exactly on the band edge counts as overlapping.
        let edge = GNSS_LO + GNSS_BW / 2.0;
        assert!(bands_overlap(edge, 0.0, GNSS_LO, GNSS_BW));
        // Just beyond the edge does not.
        assert!(!bands_overlap(edge + 1.0, 0.0, GNSS_LO, GNSS_BW));
    }

    #[test]
    fn wideband_tx_can_reach_a_nearby_band() {
        // A wide emission whose skirts reach into the band overlaps even
        // if its centre is outside the band's own bandwidth.
        let centre = GNSS_LO + 5_000_000.0; // 5 MHz away
        assert!(bands_overlap(centre, 20_000_000.0, GNSS_LO, GNSS_BW));
    }

    #[test]
    fn truth_export_delay_matches_geometry() {
        // For a stationary tx/rx pair the exported delay must equal the
        // straight-line range / c, and the range-rate must be ~0.
        use crate::config::SceneConfig;
        let scene: SceneConfig = toml::from_str(r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.1

[[transmitter]]
name = "src"
carrier_hz = 1_000_000_000.0
kind = "cw"
baseband_offset_hz = 0.0
trajectory = { kind = "stationary_ecef", x_m = 7000000.0, y_m = 0.0, z_m = 0.0 }

[[receiver]]
name = "rx"
trajectory = { kind = "stationary_ecef", x_m = 7000000.0, y_m = 0.0, z_m = 300000.0 }
[[receiver.band]]
band_id = "b"
lo_carrier_hz = 1_000_000_000.0
sample_rate_hz = 1_000_000.0
"#).unwrap();
        let dir = std::env::temp_dir().join("lipika_truth_export_test");
        super::export_truth(&scene, &dir, 100.0).unwrap();
        let v: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("rx__truth.json")).unwrap()).unwrap();
        let tx0 = &v["transmitters"][0];
        let delay0 = tx0["delay_s"][0].as_f64().unwrap();
        let expected = 300_000.0 / crate::SPEED_OF_LIGHT_M_S;
        assert!((delay0 - expected).abs() < 1e-9,
                "exported delay {delay0} != range/c {expected}");
        let rr0 = tx0["range_rate_m_s"][0].as_f64().unwrap();
        assert!(rr0.abs() < 1e-3, "stationary range-rate {rr0} should be ~0");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
