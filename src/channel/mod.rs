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
use crate::signal::SignalSource;
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

/// Per-source contribution at an observer.
///
/// `t_local` and `t_utc` are the observer's nominal sample times in
/// observer-local and UTC frames respectively. `observer_lo_carrier_hz`
/// is the observer's LO frequency (used to mix the carrier-phase
/// contribution to baseband).
pub fn source_contribution_at_observer(
    source: &dyn SignalSource,
    source_traj: &dyn Trajectory,
    observer_position_utc: &Array2<f64>,
    t_local: &Array1<f64>,
    t_utc: &Array1<f64>,
    delta_t_observer: &Array1<f64>,
    observer_lo_carrier_hz: f64,
    newton_iters: usize,
) -> Array1<Complex32> {
    let n = t_utc.len();
    let (t_emit, ranges) = solve_emission_times(source_traj, observer_position_utc, t_utc, newton_iters);
    let baseband = source.evaluate_baseband_batch(&t_emit);
    let src_carrier = source.carrier_hz();

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
        let delay = ranges[i] / SPEED_OF_LIGHT_M_S;
        let if_offset = src_carrier - observer_lo_carrier_hz;
        // Phase of the received carrier at sample i, after baseband mix:
        let phase = -tau * src_carrier * (delay + delta_t_observer[i])
                  + tau * if_offset * t_local[i];
        let (sp, cp) = phase.sin_cos();
        let mix = Complex32::new(cp as f32, sp as f32);
        out[i] = baseband[i] * mix;
    }
    out
}

/// Add complex AWGN such that the per-source SNR is `snr_db`.
/// Each call mutates `iq` in place and uses `rng` for noise samples.
pub fn add_awgn(
    iq: &mut Array1<Complex32>,
    snr_db: f64,
    expected_signal_power: f64,
    rng: &mut impl Rng,
) {
    let noise_power = expected_signal_power / 10f64.powf(snr_db / 10.0);
    let noise_std = (noise_power / 2.0).sqrt() as f32;
    let normal = Normal::new(0.0f32, noise_std).expect("valid σ");
    for v in iq.iter_mut() {
        v.re += normal.sample(rng);
        v.im += normal.sample(rng);
    }
}
