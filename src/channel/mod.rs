//! Propagation channel: emission-time lookup + delay + (optional) AWGN.
//!
//! Each observer sample is computed by solving for the source emission time
//! `t_emit` whose signal arrives at the sample's UTC instant, then evaluating
//! the source baseband there. Doppler (including code-rate Doppler and
//! sample-rate skew) is the per-sample variation of the propagation delay.

use ndarray::{Array1, Array2};
use num_complex::Complex32;
use rand::Rng;
use rand_distr::{Distribution, Normal};

use crate::geometry::{AnyTrajectory, Trajectory};
use crate::profile;
use crate::signal::SignalModel;
use crate::SPEED_OF_LIGHT_M_S;

/// Solve for emission times `t_emit` such that
/// `t_emit + |P_src(t_emit) − P_obs(t_utc)| / c = t_utc`
/// for each entry of `t_utc`. Returns also the per-sample slant range (m).
///
/// Uses fixed-point iteration from an initial guess of `t_utc − 70 ms`
/// (median LEO/GPS one-way delay); two iterations reach ns precision across
/// the LEO–GPS altitude range.
/// Generic monomorphised solver — each concrete `T: Trajectory` gets its
/// own compiled copy with `T::position_ecef` inlined into the Newton
/// inner loop. Callers usually go through the [`AnyTrajectory`]-dispatch
/// helper below, but this one is exposed for hand-tuned codepaths.
pub fn solve_emission_times_for<T: Trajectory + ?Sized>(
    source_traj: &T,
    observer_position_utc: &Array2<f64>,
    t_utc: &Array1<f64>,
    iterations: usize,
) -> (Array1<f64>, Array1<f64>) {
    let n = t_utc.len();
    assert_eq!(observer_position_utc.dim(), (3, n),
               "observer_position_utc must be (3, n) matching t_utc length");

    // Per-sample Newton loop: no (3, N) intermediate allocations, all the
    // hot data stays in registers / L1 cache. Trajectory queries go via
    // `position_ecef` (returns `[f64; 3]`, no heap), which for the
    // PrecomputedTrajectory decorator is just two array reads + a lerp.
    let mut t_emit = Array1::<f64>::zeros(n);
    let mut ranges = Array1::<f64>::zeros(n);
    let inv_c = 1.0 / SPEED_OF_LIGHT_M_S;
    let t_utc_s = t_utc.as_slice().expect("contiguous t_utc");
    let t_emit_s = t_emit.as_slice_mut().expect("contiguous t_emit");
    let ranges_s = ranges.as_slice_mut().expect("contiguous ranges");
    for i in 0..n {
        let tu = t_utc_s[i];
        let ox = observer_position_utc[[0, i]];
        let oy = observer_position_utc[[1, i]];
        let oz = observer_position_utc[[2, i]];
        let mut t = tu - 0.07;
        let mut r = 0.0;
        for _ in 0..iterations {
            let [sx, sy, sz] = source_traj.position_ecef(t);
            let dx = sx - ox;
            let dy = sy - oy;
            let dz = sz - oz;
            r = (dx * dx + dy * dy + dz * dz).sqrt();
            t = tu - r * inv_c;
        }
        t_emit_s[i] = t;
        ranges_s[i] = r;
    }
    (t_emit, ranges)
}

/// Thin AnyTrajectory-dispatching wrapper. Each match arm calls the
/// monomorphised generic, so the inner Newton loop inlines the concrete
/// `position_ecef` body for whatever trajectory type is in use.
pub fn solve_emission_times(
    source_traj: &AnyTrajectory,
    observer_position_utc: &Array2<f64>,
    t_utc: &Array1<f64>,
    iterations: usize,
) -> (Array1<f64>, Array1<f64>) {
    match source_traj {
        AnyTrajectory::Stationary(s) =>
            solve_emission_times_for(s, observer_position_utc, t_utc, iterations),
        AnyTrajectory::Precomputed(p) =>
            solve_emission_times_for(p, observer_position_utc, t_utc, iterations),
        AnyTrajectory::Dyn(d) =>
            solve_emission_times_for(&**d, observer_position_utc, t_utc, iterations),
    }
}

/// Per-transmitter contribution at one receiver band, **accumulated**
/// directly into a caller-supplied buffer. Avoids allocating an
/// intermediate per-source contribution array and the second loop pass
/// that summing into a band accumulator would otherwise need.
///
/// `t_local` / `t_utc` are the receiver's nominal sample times in
/// receiver-local / UTC frames. `receiver_lo_carrier_hz` is the receiver-
/// band LO. `tx_power_w` is the transmitter's EIRP in watts (FSPL).
/// `accumulator` is `n`-long and is added to in place.
pub fn accumulate_transmitter_contribution(
    accumulator: &mut Array1<Complex32>,
    signal: &dyn SignalModel,
    transmitter_traj: &AnyTrajectory,
    receiver_position_utc: &Array2<f64>,
    t_local: &Array1<f64>,
    t_utc: &Array1<f64>,
    delta_t_receiver: &Array1<f64>,
    receiver_lo_carrier_hz: f64,
    tx_power_w: f64,
    newton_iters: usize,
) {
    let n = t_utc.len();
    let (t_emit, ranges) = profile::time("tx_contribution:newton",
        || solve_emission_times(transmitter_traj, receiver_position_utc, t_utc, newton_iters));
    let baseband = profile::time("tx_contribution:baseband",
        || signal.evaluate_baseband_batch(&t_emit));
    let src_carrier = signal.carrier_hz();
    let wavelength_m = SPEED_OF_LIGHT_M_S / src_carrier;
    let tx_amplitude_v = tx_power_w.sqrt();

    debug_assert_eq!(accumulator.len(), n,
        "accumulator length must match per-sample arrays");
    profile::time("tx_contribution:mix", move || {

    // Phase at the observer's antenna mixer at sample n:
    //   φ_obs(n) = φ_src_emit(t_emit) - 2π · f_src · (t_local - t_emit)
    //                                  + 2π · (f_src - f_lo_obs) · t_local
    //                                  - 2π · f_lo_obs · Δt_o(t_local)
    //
    // The first term is the source's baseband (already in `baseband`).
    // The (t_local - t_emit) term = delay + Δt_o (since t_utc = t_local - Δt_o
    // and t_emit = t_utc - delay). Carrier-phase Doppler is the time
    // derivative of this and falls out automatically as we walk n.
    //
    // For an observer tuned to the source carrier (f_src == f_lo_obs) this
    // collapses to:
    //   φ_obs(n) = baseband_phase(t_emit) - 2π · f_src · delay - 2π · f_src · Δt_o(t_local)
    //
    // We allow f_lo_obs ≠ f_src to keep the IF-tuning case open.
    let tau = std::f64::consts::TAU;
    // Loop invariants, hoisted.
    let neg_tau_src = -tau * src_carrier;
    let tau_if = tau * (src_carrier - receiver_lo_carrier_hz);
    // λ/(4π) · √P_tx — constant factor of the FSPL amplitude term;
    // per-sample we only divide by range.
    let amp_factor = (wavelength_m / (2.0 * tau)) * tx_amplitude_v;
    let inv_c = 1.0 / SPEED_OF_LIGHT_M_S;
    let acc = accumulator.as_slice_mut().expect("contiguous accumulator");
    let bb = baseband.as_slice().expect("contiguous baseband");
    let ranges_s = ranges.as_slice().expect("contiguous ranges");
    let delta_t_s = delta_t_receiver.as_slice().expect("contiguous delta_t");
    let t_local_s = t_local.as_slice().expect("contiguous t_local");
    for i in 0..n {
        let r = ranges_s[i];
        let delay = r * inv_c;
        let phase = neg_tau_src * (delay + delta_t_s[i])
                  + tau_if * t_local_s[i];
        let (sp, cp) = phase.sin_cos();
        let amplitude = (amp_factor / r) as f32;
        // Complex multiply baseband × (cos + j sin), scale by amplitude,
        // accumulate into the band buffer in one fused expression.
        let b = bb[i];
        acc[i] = Complex32::new(
            acc[i].re + amplitude * (b.re * (cp as f32) - b.im * (sp as f32)),
            acc[i].im + amplitude * (b.re * (sp as f32) + b.im * (cp as f32)),
        );
    }
    })  // end of "tx_contribution:mix" profile scope
}

/// Allocating wrapper: produces a fresh `n`-long contribution buffer.
/// Used by tests / benches; the pipeline uses
/// [`accumulate_transmitter_contribution`] directly to avoid the
/// allocation + second-pass.
pub fn transmitter_contribution_at_receiver(
    signal: &dyn SignalModel,
    transmitter_traj: &AnyTrajectory,
    receiver_position_utc: &Array2<f64>,
    t_local: &Array1<f64>,
    t_utc: &Array1<f64>,
    delta_t_receiver: &Array1<f64>,
    receiver_lo_carrier_hz: f64,
    tx_power_w: f64,
    newton_iters: usize,
) -> Array1<Complex32> {
    let mut out = Array1::<Complex32>::zeros(t_utc.len());
    accumulate_transmitter_contribution(
        &mut out, signal, transmitter_traj, receiver_position_utc, t_local, t_utc,
        delta_t_receiver, receiver_lo_carrier_hz, tx_power_w, newton_iters,
    );
    out
}

/// Boltzmann's constant in J/K.
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// Thermal noise power in a bandwidth `bandwidth_hz` at temperature
/// `temperature_k` (K). Returns W.
///
/// `P_noise = k · T · B`. The bandwidth is normally the receiver's
/// sample rate (Nyquist), but you might pass a narrower value if the
/// observation bandwidth is filtered further before sampling.
pub fn thermal_noise_power_w(temperature_k: f64, bandwidth_hz: f64) -> f64 {
    BOLTZMANN_J_PER_K * temperature_k * bandwidth_hz
}

/// Add complex AWGN at a fixed receiver-side noise power (W).
///
/// The noise power is set by the receiver's thermal noise floor (kTB)
/// — independent of how much signal happens to be arriving. SNR
/// emerges from the physics: it's the ratio of received signal power
/// (after FSPL) to this fixed noise floor.
///
/// Mutates `iq` in place. Caller provides the RNG.
pub fn add_awgn_at_power(
    iq: &mut Array1<Complex32>,
    noise_power_w: f64,
    rng: &mut impl Rng,
) {
    let noise_std = (noise_power_w / 2.0).sqrt() as f32;
    let normal = Normal::new(0.0f32, noise_std).expect("valid σ");
    for v in iq.iter_mut() {
        v.re += normal.sample(rng);
        v.im += normal.sample(rng);
    }
}
