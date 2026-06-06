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

/// Earth rotation rate (rad/s), shared with the trajectory and GPS NAV code so
/// the Sagnac-corrected geometry is self-consistent across the stack.
const EARTH_ROTATION_RATE_RAD_S: f64 = 7.292_115_146_7e-5;

/// Solve for emission times `t_emit` such that
/// `t_emit + |P_src(t_emit) − P_obs(t_utc)| / c = t_utc`
/// for each entry of `t_rel_utc`. Returns also the per-sample slant range (m).
///
/// **`SIGNAL_MODEL.md` E1 (light-time solve) + E2 (Sagnac).**
///
/// **Times are epoch-relative** (seconds since the scene epoch). Keeping the
/// per-sample arithmetic in small magnitudes preserves sub-ns precision — at
/// absolute Unix time (~1.7e9 s) the f64 ULP is ~238 ns, which would quantise
/// every delay and clock offset. The absolute epoch is added back only to query
/// the trajectory (`epoch + t`), where a 238 ns error moves a satellite ~0.2 mm.
///
/// Uses fixed-point iteration from an initial guess of `t_rel_utc − 70 ms`
/// (median LEO/GPS one-way delay); two iterations reach ns precision across
/// the LEO–GPS altitude range.
/// Generic monomorphised solver — each concrete `T: Trajectory` gets its
/// own compiled copy with `T::position_ecef` inlined into the Newton
/// inner loop. Callers usually go through the [`AnyTrajectory`]-dispatch
/// helper below, but this one is exposed for hand-tuned codepaths.
pub fn solve_emission_times_for<T: Trajectory + ?Sized>(
    source_traj: &T,
    observer_position_utc: &Array2<f64>,
    t_rel_utc: &Array1<f64>,
    epoch: f64,
    iterations: usize,
) -> (Array1<f64>, Array1<f64>) {
    let n = t_rel_utc.len();
    assert_eq!(observer_position_utc.dim(), (3, n),
               "observer_position_utc must be (3, n) matching t_rel_utc length");

    // Per-sample Newton loop: no (3, N) intermediate allocations, all the
    // hot data stays in registers / L1 cache. Trajectory queries go via
    // `position_ecef` (returns `[f64; 3]`, no heap), which for the
    // PrecomputedTrajectory decorator is just two array reads + a lerp.
    let mut t_emit = Array1::<f64>::zeros(n);
    let mut ranges = Array1::<f64>::zeros(n);
    let inv_c = 1.0 / SPEED_OF_LIGHT_M_S;
    let t_utc_s = t_rel_utc.as_slice().expect("contiguous t_rel_utc");
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
            // Absolute time only for the orbit query (epoch-tolerant).
            let [sx, sy, sz] = source_traj.position_ecef(epoch + t);
            // Earth-rotation (Sagnac) correction: the source's ECEF coordinates
            // are taken at emission time `t`, but the geometric light path must
            // be measured in one frame. Rotate the source into the ECEF frame
            // at reception time `tu` by the Earth rotation over the transit
            // `tu − t`. The angle is ≤ ~6 µrad, so the linearised rotation is
            // exact to sub-mm and avoids a per-sample sin_cos. Without this the
            // range is wrong by ~tens of metres (~tens of ns) — significant for
            // a timing reference.
            let wtau = EARTH_ROTATION_RATE_RAD_S * (tu - t);
            let sx_rot = sx + wtau * sy;
            let sy_rot = sy - wtau * sx;
            let dx = sx_rot - ox;
            let dy = sy_rot - oy;
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
    t_rel_utc: &Array1<f64>,
    epoch: f64,
    iterations: usize,
) -> (Array1<f64>, Array1<f64>) {
    match source_traj {
        AnyTrajectory::Stationary(s) =>
            solve_emission_times_for(s, observer_position_utc, t_rel_utc, epoch, iterations),
        AnyTrajectory::Precomputed(p) =>
            solve_emission_times_for(p, observer_position_utc, t_rel_utc, epoch, iterations),
        AnyTrajectory::Dyn(d) =>
            solve_emission_times_for(&**d, observer_position_utc, t_rel_utc, epoch, iterations),
    }
}

/// Per-transmitter contribution at one receiver band, **accumulated**
/// directly into a caller-supplied buffer. Avoids allocating an
/// intermediate per-source contribution array and the second loop pass
/// that summing into a band accumulator would otherwise need.
///
/// **`SIGNAL_MODEL.md` E3 (amplitude), E4 (received phase); one term of E5.**
///
/// `t_local_rel` / `t_utc_rel` are the receiver's nominal sample times in the
/// receiver-local / UTC frames, **relative to the scene `epoch`** (seconds).
/// Working in relative time keeps sub-ns precision (see
/// [`solve_emission_times_for`]); `epoch` is forwarded so the solver can query
/// trajectories at absolute time. `receiver_lo_carrier_hz` is the receiver-band
/// LO. `tx_power_w` is the transmitter's EIRP in watts (FSPL). `accumulator` is
/// `n`-long and is added to in place.
pub fn accumulate_transmitter_contribution(
    accumulator: &mut Array1<Complex32>,
    signal: &dyn SignalModel,
    transmitter_traj: &AnyTrajectory,
    receiver_position_utc: &Array2<f64>,
    t_local_rel: &Array1<f64>,
    t_utc_rel: &Array1<f64>,
    delta_t_receiver: &Array1<f64>,
    receiver_lo_carrier_hz: f64,
    tx_power_w: f64,
    epoch: f64,
    newton_iters: usize,
    antenna: Option<&crate::rx::Antenna>,
) {
    let n = t_utc_rel.len();
    let (t_emit, ranges) = profile::time("tx_contribution:newton",
        || solve_emission_times(transmitter_traj, receiver_position_utc, t_utc_rel, epoch, newton_iters));
    let baseband = profile::time("tx_contribution:baseband",
        || signal.evaluate_baseband_batch(&t_emit));

    // Geometry-dependent per-sample amplitude gain: the **antenna pattern**
    // √G(elevation) from the source's arrival elevation (sin elevation = LOS·up,
    // up = O/|O|). This also carries **Earth occlusion** — an antenna with a
    // `cutoff_deg` returns 0 below it (a set satellite is blocked). `None` (no
    // antenna) ⇒ no work, byte-identical to before. Deterministic per sample, so
    // chunk/thread invariance holds. Uses the batch query at the solved emission times.
    let geom_gain: Option<Vec<f32>> = antenna.map(|a| {
        let abs_emit = t_emit.mapv(|t| epoch + t);
        let src = transmitter_traj.position_ecef_batch(&abs_emit);
        (0..n).map(|i| {
            let (ox, oy, oz) = (receiver_position_utc[[0, i]],
                                receiver_position_utc[[1, i]],
                                receiver_position_utc[[2, i]]);
            let on = (ox * ox + oy * oy + oz * oz).sqrt();
            let (dx, dy, dz) = (src[[0, i]] - ox, src[[1, i]] - oy, src[[2, i]] - oz);
            let dn = (dx * dx + dy * dy + dz * dz).sqrt();
            let sin_el = (dx * ox + dy * oy + dz * oz) / (dn * on);
            a.amplitude_at(sin_el.clamp(-1.0, 1.0).asin().to_degrees())
        }).collect()
    });
    let src_carrier = signal.carrier_hz();
    let wavelength_m = SPEED_OF_LIGHT_M_S / src_carrier;
    let tx_amplitude_v = tx_power_w.sqrt();

    debug_assert_eq!(accumulator.len(), n,
        "accumulator length must match per-sample arrays");
    profile::time("tx_contribution:mix", move || {

    // Phase at the observer's antenna mixer at sample n (SIGNAL_MODEL.md E4).
    // Fundamentally just (arriving carrier phase) − (LO phase). The arriving
    // carrier was set at the source at emission, so it rides f_src·t_emit; the
    // LO is the node oscillator, so it rides f_lo·t_local (clock time, exactly):
    //   φ_obs(n) = arg s(t_emit) + 2π·f_src·t_emit − 2π·f_lo_obs·t_local
    //
    // Rearranged into the delay/IF terms this loop evaluates:
    //   φ_obs(n) = baseband(t_emit) − 2π·f_src·(t_local − t_emit)
    //                                + 2π·(f_src − f_lo_obs)·t_local
    //
    // The (t_local − t_emit) term = delay + Δt_o. Δt_o is in it because the
    // emission solve was run against true time t_utc = t_local − Δt_o (the
    // pipeline forms `t_utc_rel = t_local_rel − delta_t`), so t_emit =
    // t_utc − delay. Hence the ONE Δt_o here drives both the code phase (via
    // t_emit → baseband) and the carrier phase: pseudorange = r + c·Δt_o, GNSS.
    // Carrier-phase Doppler is the time derivative of this, falling out as n walks.
    //
    // On-tune (f_src == f_lo_obs) this collapses to:
    //   φ_obs(n) = baseband(t_emit) − 2π·f_src·(delay + Δt_o)
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
    let t_local_s = t_local_rel.as_slice().expect("contiguous t_local_rel");
    let gg = geom_gain.as_deref();
    for i in 0..n {
        // Antenna pattern gain (carries occlusion: 0 below the antenna cutoff).
        let g = match gg { Some(v) => v[i], None => 1.0 };
        if g == 0.0 { continue; }                 // occluded ⇒ no signal (noise stays)
        let r = ranges_s[i];
        let delay = r * inv_c;
        let phase = neg_tau_src * (delay + delta_t_s[i])
                  + tau_if * t_local_s[i];
        let (sp, cp) = phase.sin_cos();
        let amplitude = (amp_factor / r) as f32 * g;
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
    t_local_rel: &Array1<f64>,
    t_utc_rel: &Array1<f64>,
    delta_t_receiver: &Array1<f64>,
    receiver_lo_carrier_hz: f64,
    tx_power_w: f64,
    epoch: f64,
    newton_iters: usize,
    antenna: Option<&crate::rx::Antenna>,
) -> Array1<Complex32> {
    let mut out = Array1::<Complex32>::zeros(t_utc_rel.len());
    accumulate_transmitter_contribution(
        &mut out, signal, transmitter_traj, receiver_position_utc, t_local_rel, t_utc_rel,
        delta_t_receiver, receiver_lo_carrier_hz, tx_power_w, epoch, newton_iters,
        antenna,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::stationary::Stationary;

    #[test]
    fn emission_solve_applies_earth_rotation_correction() {
        // GPS-altitude source with an east-west line of sight, where the Sagnac
        // term is largest. The solved range must (a) differ from the naive
        // same-frame distance by ~tens of metres and (b) match the exact
        // rotate-into-reception-frame computation to sub-cm.
        let sat = Stationary { position_ecef: [25.0e6, 10.0e6, 0.0] };
        let obs = [6_378_137.0, 0.0, 0.0];
        let mut obs_arr = Array2::<f64>::zeros((3, 1));
        for k in 0..3 { obs_arr[[k, 0]] = obs[k]; }
        // Epoch-relative reception time (the stationary source ignores absolute
        // time; the Sagnac term depends only on the transit r/c).
        let t_rel = Array1::from_vec(vec![0.0]);
        let traj = AnyTrajectory::Stationary(sat);
        let (_te, ranges) = solve_emission_times(&traj, &obs_arr, &t_rel, 1_717_260_840.0, 4);
        let solved = ranges[0];

        let s = [25.0e6, 10.0e6, 0.0];
        let naive = ((s[0]-obs[0]).powi(2) + (s[1]-obs[1]).powi(2) + (s[2]-obs[2]).powi(2)).sqrt();
        // Exact rotation by the Earth angle over the (solved) transit.
        let th = EARTH_ROTATION_RATE_RAD_S * solved / SPEED_OF_LIGHT_M_S;
        let (sn, cs) = th.sin_cos();
        let sr = [s[0]*cs + s[1]*sn, -s[0]*sn + s[1]*cs, s[2]];
        let exact = ((sr[0]-obs[0]).powi(2) + (sr[1]-obs[1]).powi(2) + (sr[2]-obs[2]).powi(2)).sqrt();

        assert!((solved - naive).abs() > 1.0,
            "Sagnac correction not applied: solved {solved} vs naive {naive}");
        assert!((solved - exact).abs() < 0.01,
            "linearised correction off from exact rotation: {solved} vs {exact}");
    }

    #[test]
    fn antenna_cutoff_occludes_below_horizon() {
        use crate::signal::ca_code::CaCodeSource;
        use crate::rx::Antenna;
        // Observer on the equator; local up is +x.
        let obs = [6_378_137.0, 0.0, 0.0];
        let n = 8;
        let mut obs_arr = Array2::<f64>::zeros((3, n));
        for i in 0..n { for k in 0..3 { obs_arr[[k, i]] = obs[k]; } }
        let t = Array1::<f64>::zeros(n);
        let dt = Array1::<f64>::zeros(n);

        // Overhead source (elevation +90°) and a source clearly below the horizon
        // (on the far side of the Earth ⇒ LOS·up < 0).
        let above = AnyTrajectory::Stationary(Stationary { position_ecef: [6.9e6, 0.0, 0.0] });
        let below = AnyTrajectory::Stationary(Stationary { position_ecef: [-7.0e6, 0.0, 0.0] });
        let src = CaCodeSource::new(1_575_420_000.0, 1, 0.0);
        // Unit-gain antennas: one with a 0° cutoff (the occlusion), one without.
        let occ = Antenna { peak_dbi: 0.0, horizon_dbi: 0.0, fb_db: 0.0, cutoff_deg: Some(0.0) };
        let open = Antenna { peak_dbi: 0.0, horizon_dbi: 0.0, fb_db: 0.0, cutoff_deg: None };
        let power = |traj, ant: Option<&Antenna>| {
            let c = transmitter_contribution_at_receiver(
                &src, traj, &obs_arr, &t, &t, &dt, 1_575_420_000.0, 20.0, 1_717_260_840.0, 4, ant);
            c.iter().map(|v| v.norm_sqr() as f64).sum::<f64>()
        };

        // No antenna ⇒ no occlusion; an open antenna radiates the back-lobe.
        assert!(power(&below, None) > 0.0, "below-horizon source radiates with no antenna");
        assert!(power(&below, Some(&open)) > 0.0, "open antenna: below-horizon via back-lobe");
        // A 0° cutoff occludes below the horizon, leaves the overhead source.
        assert_eq!(power(&below, Some(&occ)), 0.0, "0° cutoff must drop a below-horizon source");
        assert!(power(&above, Some(&occ)) > 0.0, "overhead source survives the cutoff");
    }
}
