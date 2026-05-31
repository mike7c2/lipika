//! Propagation channel: emission-time lookup + delay + (optional) AWGN.
//!
//! The core innovation vs frequency-shift simulators: we compute each
//! observer sample by solving for the source's emission time `t_emit` such
//! that the signal reaches the observer's antenna exactly at the UTC
//! instant the observer's sample is taken, and then look up the source
//! baseband at `t_emit`. The propagation delay's per-sample variation
//! *is* the Doppler — it falls out automatically, including code-rate
//! Doppler and sample-rate skew effects.

use ndarray::{Array1, Array2};
use num_complex::Complex32;
use rand::Rng;
use rand_distr::{Distribution, Normal};

use crate::geometry::Trajectory;
use crate::signal::SignalModel;
use crate::SPEED_OF_LIGHT_M_S;

pub mod free_space;

/// Solve for emission times `t_emit` such that
/// `t_emit + |P_src(t_emit) − P_obs(t_utc)| / c = t_utc`
/// for each entry of `t_utc`. Returns also the per-sample slant range (m).
///
/// Uses fixed-point iteration from an initial guess of `t_utc − 70 ms`
/// (median LEO/GPS one-way delay). Two iterations bring us to ns precision
/// across the LEO–GPS altitude range; three for paranoia.
pub fn solve_emission_times(
    source_traj: &dyn Trajectory,
    observer_position_utc: &Array2<f64>,
    t_utc: &Array1<f64>,
    iterations: usize,
) -> (Array1<f64>, Array1<f64>) {
    let n = t_utc.len();
    assert_eq!(observer_position_utc.dim(), (3, n),
               "observer_position_utc must be (3, n) matching t_utc length");

    let mut t_emit = t_utc - 0.07;
    let mut ranges = Array1::<f64>::zeros(n);

    for _ in 0..iterations {
        let src_pos = source_traj.position_ecef_batch(&t_emit);
        for i in 0..n {
            let dx = src_pos[[0, i]] - observer_position_utc[[0, i]];
            let dy = src_pos[[1, i]] - observer_position_utc[[1, i]];
            let dz = src_pos[[2, i]] - observer_position_utc[[2, i]];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            ranges[i] = r;
            t_emit[i] = t_utc[i] - r / SPEED_OF_LIGHT_M_S;
        }
    }
    (t_emit, ranges)
}

/// Per-transmitter contribution at one receiver band.
///
/// `t_local` and `t_utc` are the receiver's nominal sample times in
/// receiver-local and UTC frames respectively. `receiver_lo_carrier_hz`
/// is the receiver-band LO (used to mix to baseband). `tx_power_w`
/// is the transmitter's effective isotropic radiated power in watts
/// (used for free-space-path-loss amplitude scaling).
pub fn transmitter_contribution_at_receiver(
    signal: &dyn SignalModel,
    transmitter_traj: &dyn Trajectory,
    receiver_position_utc: &Array2<f64>,
    t_local: &Array1<f64>,
    t_utc: &Array1<f64>,
    delta_t_receiver: &Array1<f64>,
    receiver_lo_carrier_hz: f64,
    tx_power_w: f64,
    newton_iters: usize,
) -> Array1<Complex32> {
    let n = t_utc.len();
    let (t_emit, ranges) = solve_emission_times(transmitter_traj, receiver_position_utc, t_utc, newton_iters);
    let baseband = signal.evaluate_baseband_batch(&t_emit);
    let src_carrier = signal.carrier_hz();
    let wavelength_m = SPEED_OF_LIGHT_M_S / src_carrier;
    let tx_amplitude_v = tx_power_w.sqrt();

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
    let mut out = Array1::<Complex32>::zeros(n);
    let tau = std::f64::consts::TAU;
    for i in 0..n {
        let r = ranges[i];
        let delay = r / SPEED_OF_LIGHT_M_S;
        let if_offset = src_carrier - receiver_lo_carrier_hz;
        // Phase of the received carrier at sample i, after baseband mix:
        let phase = -tau * src_carrier * (delay + delta_t_receiver[i])
                  + tau * if_offset * t_local[i];
        let (sp, cp) = phase.sin_cos();
        // Free-space path loss: |E| ∝ λ / (4π·r). Combined with transmitter
        // EIRP √(P_tx) we get receiver-side voltage amplitude.
        // Isotropic Rx antenna assumed; per-band gain patterns are a future
        // extension (`channel::antenna_gain` placeholder).
        let amplitude = (wavelength_m / (2.0 * tau * r)) * tx_amplitude_v;
        let mix = Complex32::new((amplitude * cp) as f32, (amplitude * sp) as f32);
        out[i] = baseband[i] * mix;
    }
    out
}

/// Keep the legacy name as a thin shim during the refactor.
#[deprecated(note = "renamed to transmitter_contribution_at_receiver; use that instead")]
pub fn source_contribution_at_observer(
    source: &dyn SignalModel,
    source_traj: &dyn Trajectory,
    observer_position_utc: &Array2<f64>,
    t_local: &Array1<f64>,
    t_utc: &Array1<f64>,
    delta_t_observer: &Array1<f64>,
    observer_lo_carrier_hz: f64,
    newton_iters: usize,
) -> Array1<Complex32> {
    transmitter_contribution_at_receiver(
        source, source_traj, observer_position_utc, t_local, t_utc,
        delta_t_observer, observer_lo_carrier_hz,
        1.0, // 1 W default for backward-compat
        newton_iters,
    )
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

/// Compatibility wrapper for callers that want "noise calibrated to a
/// declared per-source SNR against a 1 W reference". Kept so existing
/// tests don't break, but new scenes should specify a noise temperature
/// or noise floor explicitly.
#[deprecated(note = "specify receiver noise via thermal_noise_power_w + add_awgn_at_power; absolute noise model is physically meaningful")]
pub fn add_awgn(
    iq: &mut Array1<Complex32>,
    snr_db: f64,
    _expected_signal_power: f64,
    rng: &mut impl Rng,
) {
    let noise_power = 1.0 / 10f64.powf(snr_db / 10.0);
    add_awgn_at_power(iq, noise_power, rng);
}
