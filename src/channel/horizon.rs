//! Horizon gate: decide, cheaply and without ever being wrong, that a
//! transmitter is below the receive antenna's cutoff elevation for a whole
//! stretch of samples, so the per-sample work for it can be skipped.
//!
//! The per-sample path ([`super::accumulate_transmitter_contribution`]) gives
//! a sample exactly zero contribution when the arrival elevation it computes
//! is below `Antenna::cutoff_deg`. That elevation is taken between the
//! receiver position `O` at the sample's UTC instant and the source position
//! `S` at the solved emission time, against the geocentric up `O/|O|`. The gate
//! bounds that same quantity over a stretch of samples:
//!
//! - **Receiver.** Every `O` of the stretch lies in a ball: centre `O_c` (the
//!   first sample's position), radius `ρ_o` (the largest distance of any other
//!   sample's position from it, measured on the very array the per-sample path
//!   reads). `ρ_o = 0` for an Earth-fixed receiver.
//! - **Emission time.** The solve leaves `t_emit = t_utc − r/c` with `r` the
//!   last iterate's range, so `t_utc − L ≤ t_emit ≤ t_utc` for any `L ≥ r/c`.
//!   `r ≤ 1.01·R_max + |O|`, where `R_max` bounds the source's distance from
//!   the Earth's centre and 1.01 covers the Sagnac rotation of the solve for
//!   any light time up to 100 s (the gate declines beyond that). So every
//!   emission time of the stretch lies in `[min t_utc − L, max t_utc]`.
//! - **Source.** The trajectory the per-sample path queries is either fixed or
//!   a [`PrecomputedTrajectory`](crate::geometry::precomputed::PrecomputedTrajectory),
//!   whose position moves no faster than its `max_speed_m_s` (the longest grid
//!   segment over the grid step; it is constant outside the grid). So every
//!   `S` of the stretch lies in a ball: centre `S_c` (the position at the
//!   middle of that time interval), radius `ρ_s = v_max · (half the interval)`
//!   plus 1 m for float rounding.
//! - **Elevation.** With `d_c = S_c − O_c` and `ρ = ρ_s + ρ_o`, the direction
//!   `S − O` is within `asin(ρ/|d_c|)` of `d_c` and the up direction `O` is
//!   within `asin(ρ_o/|O_c|)` of `O_c`, so every elevation of the stretch is at
//!   most `el(d_c, O_c) + asin(ρ/|d_c|) + asin(ρ_o/|O_c|)`.
//!
//! The stretch is skipped only if that bound, plus [`MARGIN_DEG`] for the
//! rounding of the per-sample `asin`, is below the cutoff. Anything else
//! (including any non-finite input, a trajectory kind without a speed bound,
//! or a source that is merely close to the cutoff) goes to the per-sample
//! path, which decides sample by sample exactly as it always has. The bound
//! uses the tabulated motion itself, not an assumption about orbits, so it
//! holds for every trajectory kind and for moving receivers.
//!
//! How tight it is: for a LEO (7.6 km/s) and the default 0.1 s stretch at
//! 48 kHz, `ρ_s ≈ 7.6 km/s × (0.1 s + 0.07 s)/2 ≈ 0.65 km`; a satellite at the
//! horizon is ≥ 2000 km away, so the spread is ≈ 0.02°. A satellite is handed
//! to the per-sample path a fraction of a second before it reaches the cutoff.

use std::ops::Range;

use ndarray::{Array1, Array2};

use super::EMISSION_INITIAL_LAG_S;
use crate::geometry::AnyTrajectory;
use crate::SPEED_OF_LIGHT_M_S;

/// Nominal length of one gate stretch (s).
pub const STRETCH_S: f64 = 0.1;

/// Shortest gate stretch (samples), so low sample rates do not pay a gate
/// test every few samples.
pub const MIN_STRETCH_SAMPLES: usize = 256;

/// Extra elevation (degrees) the bound must clear the cutoff by. The
/// per-sample elevation is `asin` of an f64 ratio; near ±90° a 1-ulp error in
/// the ratio moves the angle by up to ~3e-6°, and the gate's own `asin` calls
/// have the same error. 1e-3° covers both several hundred times over and
/// costs ~35 m of range at the horizon.
pub const MARGIN_DEG: f64 = 1.0e-3;

/// Longest light time (s) for which the 1 % Sagnac allowance in the range
/// bound holds (`ω⊕ · 100 s = 0.0073 < 0.01`). Farther sources are not gated.
const MAX_LIGHT_TIME_S: f64 = 100.0;

/// Gate stretch length in samples for a band sampled at `sample_rate_hz`.
pub fn stretch_samples(sample_rate_hz: f64) -> usize {
    let n = (sample_rate_hz * STRETCH_S).round();
    if n.is_finite() && n > MIN_STRETCH_SAMPLES as f64 { n as usize } else { MIN_STRETCH_SAMPLES }
}

/// What the gate needs to know about the receiver over one stretch of a
/// chunk: a ball holding every receiver position and the span of UTC sample
/// times. Computed once per chunk and shared by all transmitters.
#[derive(Debug, Clone)]
pub struct ReceiverStretch {
    /// Sample range within the chunk.
    pub samples: Range<usize>,
    /// Receiver position at the first sample of the stretch (ECEF, m).
    centre: [f64; 3],
    centre_norm: f64,
    /// Largest distance of any sample's receiver position from `centre` (m).
    radius: f64,
    /// Smallest and largest epoch-relative UTC sample time (s).
    t_utc_min: f64,
    t_utc_max: f64,
    /// False if any input of the stretch was not finite: never skip.
    finite: bool,
}

/// Cut a chunk into stretches of at most `stretch_samples` and summarise the
/// receiver over each. `receiver_position_utc` is the `(3, n)` array and
/// `t_utc_rel` the `n` epoch-relative UTC sample times that the per-sample
/// path is given for the same chunk.
#[mantra_macros::req_satisfied("nf-5")]
pub fn receiver_stretches(
    receiver_position_utc: &Array2<f64>,
    t_utc_rel: &Array1<f64>,
    stretch_samples: usize,
) -> Vec<ReceiverStretch> {
    let n = t_utc_rel.len();
    assert_eq!(receiver_position_utc.dim(), (3, n),
               "receiver_position_utc must be (3, n) matching t_utc_rel length");
    let step = stretch_samples.max(1);
    let mut out = Vec::with_capacity(n.div_ceil(step));
    let mut a = 0usize;
    while a < n {
        let b = (a + step).min(n);
        let p = |i: usize| [receiver_position_utc[[0, i]], receiver_position_utc[[1, i]],
                            receiver_position_utc[[2, i]]];
        let centre = p(a);
        let centre_norm = norm(centre);
        let mut finite = centre_norm.is_finite();
        let mut r2_max = 0.0_f64;
        let (mut t_min, mut t_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for i in a..b {
            let q = p(i);
            let d = [q[0] - centre[0], q[1] - centre[1], q[2] - centre[2]];
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let t = t_utc_rel[i];
            // `is_finite` first: a NaN would slip through the `>`/`<` below.
            finite &= r2.is_finite() & t.is_finite();
            if r2 > r2_max { r2_max = r2; }
            if t < t_min { t_min = t; }
            if t > t_max { t_max = t; }
        }
        out.push(ReceiverStretch {
            samples: a..b,
            centre,
            centre_norm,
            radius: r2_max.sqrt(),
            t_utc_min: t_min,
            t_utc_max: t_max,
            finite,
        });
        a = b;
    }
    out
}

/// True only if the per-sample path would give **every** sample of the
/// stretch a zero antenna gain for this transmitter, i.e. the arrival
/// elevation it computes is below `cutoff_deg` throughout. False means "not
/// certain", never "above". See the module documentation for the bound.
///
/// `epoch` is the scene epoch the relative times are measured from, exactly as
/// passed to [`super::accumulate_transmitter_contribution`].
#[mantra_macros::req_satisfied("nf-5")]
pub fn certainly_below_cutoff(
    transmitter_traj: &AnyTrajectory,
    stretch: &ReceiverStretch,
    epoch: f64,
    cutoff_deg: f64,
) -> bool {
    if !stretch.finite {
        return false;
    }
    // (speed bound, bound on distance from the Earth's centre, rounding slack)
    let (max_speed, max_radius, slack_m) = match transmitter_traj {
        AnyTrajectory::Stationary(s) => (0.0, norm(s.position_ecef), 0.0),
        AnyTrajectory::Precomputed(p) => (p.max_speed_m_s(), p.max_radius_m(), 1.0),
        // No bound on how an arbitrary trajectory moves between two queries.
        AnyTrajectory::Dyn(_) => return false,
    };

    // Light-time bound L: t_utc − L ≤ t_emit ≤ t_utc for every sample. Each
    // solver iterate has r ≤ |S|·(1 + ω⊕·Δ) + |O| with Δ the previous
    // iterate's lag (70 ms at the start); by induction Δ ≤ L ≤ 100 s, so
    // ω⊕·Δ < 0.01. With zero iterations the lag is the 70 ms guess itself.
    let light_time = ((1.01 * max_radius + stretch.centre_norm + stretch.radius)
        / SPEED_OF_LIGHT_M_S).max(EMISSION_INITIAL_LAG_S);
    if !(light_time <= MAX_LIGHT_TIME_S) {
        return false; // too far, or not finite
    }
    // Rounding of `t_utc − r/c` and of the subtraction below: far inside 1 µs
    // plus a part in 1e12 of the time magnitude.
    let t_mag = stretch.t_utc_min.abs().max(stretch.t_utc_max.abs());
    let lag = light_time * (1.0 + 1e-9) + 1e-6 + 1e-12 * t_mag;

    // Absolute emission-time interval. Float addition is monotone, so
    // `epoch + t_emit` (what the per-sample path queries) lies in
    // [t_lo, t_hi] for every t_emit in [t_utc_min − lag, t_utc_max].
    let t_lo = epoch + (stretch.t_utc_min - lag);
    let t_hi = epoch + stretch.t_utc_max;
    let t_mid = t_lo + 0.5 * (t_hi - t_lo);
    let half_span = (t_mid - t_lo).max(t_hi - t_mid);

    // Ball holding every source position of the stretch.
    let s_c = transmitter_traj.position_ecef(t_mid);
    let source_radius = max_speed * half_span * (1.0 + 1e-9) + slack_m;

    let o_c = stretch.centre;
    let d_c = [s_c[0] - o_c[0], s_c[1] - o_c[1], s_c[2] - o_c[2]];
    let range = norm(d_c);
    let rho = source_radius + stretch.radius;
    // Written so that a NaN anywhere answers "not certain".
    if !(rho < range) || !(stretch.radius < stretch.centre_norm) {
        return false;
    }
    let sin_el = (d_c[0] * o_c[0] + d_c[1] * o_c[1] + d_c[2] * o_c[2])
        / (range * stretch.centre_norm);
    let el_centre_deg = sin_el.clamp(-1.0, 1.0).asin().to_degrees();
    let spread_deg = (rho / range).asin().to_degrees()
        + (stretch.radius / stretch.centre_norm).asin().to_degrees();
    el_centre_deg + spread_deg + MARGIN_DEG < cutoff_deg
}

/// The sample ranges of a chunk that must go through the per-sample path for
/// this transmitter: every stretch not [`certainly_below_cutoff`], adjacent
/// ones merged. Empty ⇒ the transmitter contributes exactly zero to the whole
/// chunk.
#[mantra_macros::req_satisfied("nf-5")]
pub fn kept_ranges(
    transmitter_traj: &AnyTrajectory,
    stretches: &[ReceiverStretch],
    epoch: f64,
    cutoff_deg: f64,
) -> Vec<Range<usize>> {
    let mut kept: Vec<Range<usize>> = Vec::new();
    for s in stretches {
        if certainly_below_cutoff(transmitter_traj, s, epoch, cutoff_deg) {
            continue;
        }
        match kept.last_mut() {
            Some(last) if last.end == s.samples.start => last.end = s.samples.end,
            _ => kept.push(s.samples.clone()),
        }
    }
    kept
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::channel::{solve_emission_times, transmitter_contribution_at_receiver};
    use crate::geometry::kepler::CircularKeplerian;
    use crate::geometry::precomputed::PrecomputedTrajectory;
    use crate::geometry::stationary::Stationary;
    use crate::geometry::Trajectory;
    use crate::rx::Antenna;
    use crate::signal::cw::ContinuousWave;

    const EPOCH: f64 = 1_717_260_690.0;

    /// Per-sample inputs for `n` samples at `fs` from the epoch, for a
    /// receiver on `rx` with a constant clock offset.
    fn sample_inputs(rx: &dyn Trajectory, n: usize, fs: f64, clock_offset_s: f64)
        -> (Array2<f64>, Array1<f64>, Array1<f64>, Array1<f64>)
    {
        let t_local = Array1::from_iter((0..n).map(|i| i as f64 / fs));
        let delta_t = Array1::from_elem(n, clock_offset_s);
        let t_utc = &t_local - &delta_t;
        let pos = rx.position_ecef_batch(&t_utc.mapv(|t| EPOCH + t));
        (pos, t_local, t_utc, delta_t)
    }

    /// The elevation (degrees) the per-sample path computes for each sample:
    /// source at the solved emission time against the receiver's geocentric up.
    fn per_sample_elevation_deg(traj: &AnyTrajectory, pos: &Array2<f64>, t_utc: &Array1<f64>)
        -> Vec<f64>
    {
        let (t_emit, _) = solve_emission_times(traj, pos, t_utc, EPOCH, 2);
        let src = traj.position_ecef_batch(&t_emit.mapv(|t| EPOCH + t));
        (0..t_utc.len()).map(|i| {
            let o = [pos[[0, i]], pos[[1, i]], pos[[2, i]]];
            let d = [src[[0, i]] - o[0], src[[1, i]] - o[1], src[[2, i]] - o[2]];
            ((d[0] * o[0] + d[1] * o[1] + d[2] * o[2]) / (norm(d) * norm(o)))
                .clamp(-1.0, 1.0).asin().to_degrees()
        }).collect()
    }

    fn unit_antenna(cutoff_deg: f64) -> Antenna {
        Antenna { peak_dbi: 0.0, horizon_dbi: 0.0, fb_db: 0.0, cutoff_deg: Some(cutoff_deg) }
    }

    /// Sweep LEOs past a receiver and check every verdict of the gate against
    /// the per-sample path. Returns (stretches skipped, stretches wholly more
    /// than `slack_deg` below the cutoff, stretches in which the source
    /// crossed the cutoff).
    fn sweep(rx: &dyn Trajectory, cutoff_deg: f64, clock_offset_s: f64, slack_deg: f64)
        -> (usize, usize, usize)
    {
        let (fs, n, stretch) = (500.0, 10_000usize, 50usize); // 20 s, 0.1 s stretches
        let (pos, t_local, t_utc, delta_t) = sample_inputs(rx, n, fs, clock_offset_s);
        let stretches = receiver_stretches(&pos, &t_utc, stretch);
        assert_eq!(stretches.len(), n / stretch);
        let antenna = unit_antenna(cutoff_deg);
        let tone = ContinuousWave { carrier_hz: 401.15e6, baseband_offset_hz: 0.0, anchor_phase_rad: None };
        let (mut skipped, mut clearly_below, mut crossing) = (0, 0, 0);
        // Geometric elevation of an orbit from the receiver at a sample.
        let duration_s = n as f64 / fs;
        let elevation_at = |orbit: &CircularKeplerian, i: usize| {
            let o = [pos[[0, i]], pos[[1, i]], pos[[2, i]]];
            let p = orbit.position_ecef(EPOCH + t_utc[i]);
            let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
            ((d[0] * o[0] + d[1] * o[1] + d[2] * o[2]) / (norm(d) * norm(o))).asin().to_degrees()
        };
        // 40 orbits taken as they come (nearly all far below the cutoff) and
        // 80 picked because they are within 1.5° of it at one end of the
        // window, so the sweep is dense in rises, sets and grazing passes.
        let (mut any, mut near, mut k) = (0, 0, 0usize);
        while any < 40 || near < 80 {
            k += 1;
            assert!(k < 200_000, "ran out of candidate orbits");
            let orbit = CircularKeplerian {
                semi_major_axis_m: 6_878_137.0 + 40_000.0 * (k % 7) as f64,
                inclination_deg: [51.6, 65.0, 87.9, 97.5][k % 4],
                raan_deg: 37.0 * k as f64,
                true_anomaly_deg_at_epoch: 53.0 * k as f64,
                epoch_unix: EPOCH,
            };
            let is_near = (elevation_at(&orbit, 0) - cutoff_deg).abs() < 1.5
                || (elevation_at(&orbit, n - 1) - cutoff_deg).abs() < 1.5;
            if is_near && near < 80 {
                near += 1;
            } else if any < 40 {
                any += 1;
            } else {
                continue;
            }
            let traj = AnyTrajectory::Precomputed(
                PrecomputedTrajectory::sample(&orbit, EPOCH, duration_s, 0.01));
            let elevation = per_sample_elevation_deg(&traj, &pos, &t_utc);
            let contribution = transmitter_contribution_at_receiver(
                &tone, &traj, &pos, &t_local, &t_utc, &delta_t, 401.15e6, 1.0, EPOCH, 2,
                Some(&antenna));
            for s in &stretches {
                let el = &elevation[s.samples.clone()];
                let max_el = el.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let min_el = el.iter().cloned().fold(f64::INFINITY, f64::min);
                if min_el < cutoff_deg && max_el >= cutoff_deg { crossing += 1; }
                let skip = certainly_below_cutoff(&traj, s, EPOCH, cutoff_deg);
                if skip {
                    skipped += 1;
                    // Sound: nothing the per-sample path would have produced.
                    assert!(max_el < cutoff_deg,
                        "orbit {k}: skipped a stretch that reaches {max_el}° (cutoff {cutoff_deg}°)");
                    assert!(contribution.slice(ndarray::s![s.samples.clone()])
                        .iter().all(|v| v.re == 0.0 && v.im == 0.0));
                }
                if max_el < cutoff_deg - slack_deg {
                    clearly_below += 1;
                    // Tight: a source this far below is always skipped.
                    assert!(skip, "orbit {k}: kept a stretch that peaks at {max_el}° \
                                   (cutoff {cutoff_deg}°, slack {slack_deg}°)");
                }
            }
        }
        (skipped, clearly_below, crossing)
    }

    #[mantra_macros::req_test("nf-5")]
    #[test]
    fn gate_skips_only_what_the_per_sample_path_zeroes() {
        let glasgow = Stationary::from_geodetic(55.86, -4.25, 20.0);
        for (cutoff_deg, clock_offset_s) in [(0.0, 0.0), (10.0, 2.5e-3), (-5.0, -40.0e-3)] {
            // 0.1 s stretches of a LEO seen from the ground: within 0.1° of
            // the cutoff is the most the gate ever gives away.
            let (skipped, clearly_below, crossing) = sweep(&glasgow, cutoff_deg, clock_offset_s, 0.1);
            assert!(crossing >= 20, "cutoff {cutoff_deg}°: only {crossing} rises/sets in the sweep");
            assert!(clearly_below > 5_000 && skipped >= clearly_below,
                "cutoff {cutoff_deg}°: skipped {skipped}, clearly below {clearly_below}");
        }
    }

    #[mantra_macros::req_test("nf-5")]
    #[test]
    fn gate_is_sound_for_a_moving_receiver() {
        // A receiver that is itself in orbit (below the sources): over a
        // stretch its positions fill a ball ~0.8 km in radius, the sources can
        // be near, and both the line of sight and the local up swing. The gate
        // widens its bound accordingly (so it gives more away, here up to 2°)
        // and stays sound.
        let rx = CircularKeplerian {
            semi_major_axis_m: 6_778_137.0,
            inclination_deg: 80.0,
            raan_deg: 200.0,
            true_anomaly_deg_at_epoch: 10.0,
            epoch_unix: EPOCH,
        };
        let (skipped, clearly_below, crossing) = sweep(&rx, 0.0, 0.0, 2.0);
        assert!(crossing >= 20 && clearly_below > 5_000 && skipped >= clearly_below,
            "skipped {skipped}, clearly below {clearly_below}, crossing {crossing}");
    }

    /// A fixed source at `elevation_deg` and 2000 km from an equatorial
    /// observer whose up is +x.
    fn fixed_source_at(elevation_deg: f64) -> AnyTrajectory {
        let (s, c) = elevation_deg.to_radians().sin_cos();
        AnyTrajectory::Stationary(Stationary {
            position_ecef: [6_378_137.0 + 2.0e6 * s, 2.0e6 * c, 0.0],
        })
    }

    #[mantra_macros::req_test("nf-5")]
    #[test]
    fn fixed_source_is_skipped_exactly_when_it_clears_the_margin() {
        let rx = Stationary { position_ecef: [6_378_137.0, 0.0, 0.0] };
        let (pos, _, t_utc, _) = sample_inputs(&rx, 300, 1000.0, 0.0);
        let stretches = receiver_stretches(&pos, &t_utc, 300);
        assert_eq!(stretches.len(), 1);
        // Both ends fixed: the bound is the elevation itself plus the margin.
        for cutoff_deg in [-20.0, 0.0, 5.0, 60.0] {
            let verdict = |offset_deg: f64| certainly_below_cutoff(
                &fixed_source_at(cutoff_deg + offset_deg), &stretches[0], EPOCH, cutoff_deg);
            assert!(verdict(-1.0));
            assert!(verdict(-2.0 * MARGIN_DEG), "just beyond the margin is skipped");
            assert!(!verdict(-0.5 * MARGIN_DEG), "inside the margin is kept");
            assert!(!verdict(0.0), "at the cutoff is kept (the antenna admits it)");
            assert!(!verdict(1.0));
        }
        // The kept range is all or nothing here.
        let below = fixed_source_at(-3.0);
        let above = fixed_source_at(3.0);
        assert!(kept_ranges(&below, &stretches, EPOCH, 0.0).is_empty());
        assert_eq!(kept_ranges(&above, &stretches, EPOCH, 0.0), [0..300]);
    }

    #[mantra_macros::req_test("nf-5")]
    #[test]
    fn gate_declines_whenever_it_cannot_bound() {
        let rx = Stationary { position_ecef: [6_378_137.0, 0.0, 0.0] };
        let (pos, _, t_utc, _) = sample_inputs(&rx, 100, 1000.0, 0.0);
        let stretch = &receiver_stretches(&pos, &t_utc, 100)[0];
        let far_below = || Stationary { position_ecef: [-7.0e6, 0.0, 0.0] }; // antipodal
        assert!(certainly_below_cutoff(
            &AnyTrajectory::Stationary(far_below()), stretch, EPOCH, 0.0));

        // A trajectory kind with no speed bound is never skipped.
        let dynamic = AnyTrajectory::Dyn(std::sync::Arc::new(far_below()));
        assert!(!certainly_below_cutoff(&dynamic, stretch, EPOCH, 0.0));
        // A cutoff that is not a number occludes nothing per sample either.
        assert!(!certainly_below_cutoff(
            &AnyTrajectory::Stationary(far_below()), stretch, EPOCH, f64::NAN));
        // A source more than 100 light-seconds away is outside the bound.
        let remote = AnyTrajectory::Stationary(Stationary { position_ecef: [-4.0e10, 0.0, 0.0] });
        assert!(!certainly_below_cutoff(&remote, stretch, EPOCH, 0.0));
        // A tabulated path holding a non-finite position has no speed bound.
        struct Broken;
        impl Trajectory for Broken {
            fn position_ecef(&self, t: f64) -> [f64; 3] {
                if t > EPOCH + 0.5 { [f64::NAN, 0.0, 0.0] } else { [-7.0e6, 0.0, 0.0] }
            }
        }
        let broken = AnyTrajectory::Precomputed(PrecomputedTrajectory::sample(&Broken, EPOCH, 1.0, 0.01));
        assert!(!certainly_below_cutoff(&broken, stretch, EPOCH, 0.0));

        // Non-finite receiver position or sample time anywhere in the stretch.
        for (row, bad) in [(0usize, f64::NAN), (2, f64::INFINITY)] {
            let mut p = pos.clone();
            p[[row, 57]] = bad;
            let s = &receiver_stretches(&p, &t_utc, 100)[0];
            assert!(!certainly_below_cutoff(&AnyTrajectory::Stationary(far_below()), s, EPOCH, 0.0));
        }
        let mut t = t_utc.clone();
        t[99] = f64::NAN;
        let s = &receiver_stretches(&pos, &t, 100)[0];
        assert!(!certainly_below_cutoff(&AnyTrajectory::Stationary(far_below()), s, EPOCH, 0.0));
        // A receiver at the Earth's centre has no up direction.
        let centre = Array2::<f64>::zeros((3, 100));
        let s = &receiver_stretches(&centre, &t_utc, 100)[0];
        assert!(!certainly_below_cutoff(&AnyTrajectory::Stationary(far_below()), s, EPOCH, 0.0));
    }

    #[mantra_macros::req_test("nf-5")]
    #[test]
    fn stretches_tile_the_chunk_and_kept_ranges_merge() {
        assert_eq!(stretch_samples(48_000.0), 4800);
        assert_eq!(stretch_samples(1_000.0), MIN_STRETCH_SAMPLES);
        assert_eq!(stretch_samples(f64::NAN), MIN_STRETCH_SAMPLES);

        // A receiver that steps 3 m per sample along z, with times running
        // backwards (a clock can do that): the summary must not assume order.
        let n = 25;
        let mut pos = Array2::<f64>::zeros((3, n));
        for i in 0..n {
            pos[[0, i]] = 6.4e6;
            pos[[2, i]] = 3.0 * i as f64;
        }
        let t_utc = Array1::from_iter((0..n).map(|i| 10.0 - i as f64));
        let st = receiver_stretches(&pos, &t_utc, 10);
        let ranges: Vec<_> = st.iter().map(|s| s.samples.clone()).collect();
        assert_eq!(ranges, [0..10, 10..20, 20..25]);
        assert_eq!(st[0].radius, 27.0);
        assert_eq!(st[2].radius, 12.0);
        assert_eq!((st[1].t_utc_min, st[1].t_utc_max), (-9.0, 0.0));
        assert!(receiver_stretches(&Array2::zeros((3, 0)), &Array1::zeros(0), 10).is_empty());

        // Three stretches for a fixed receiver; a source that is below, a
        // source that is above, and the merged result.
        let rx = Stationary { position_ecef: [6_378_137.0, 0.0, 0.0] };
        let (pos, _, t_utc, _) = sample_inputs(&rx, 25, 1000.0, 0.0);
        let st = receiver_stretches(&pos, &t_utc, 10);
        assert_eq!(kept_ranges(&fixed_source_at(20.0), &st, EPOCH, 0.0), [0..25]);
        assert!(kept_ranges(&fixed_source_at(-20.0), &st, EPOCH, 0.0).is_empty());
        // Make only the middle stretch unusable: the two kept neighbours of a
        // skipped stretch stay separate, and a kept run merges.
        let mut t_bad = t_utc.clone();
        t_bad[12] = f64::NAN;
        let st_bad = receiver_stretches(&pos, &t_bad, 10);
        assert_eq!(kept_ranges(&fixed_source_at(-20.0), &st_bad, EPOCH, 0.0), [10..20]);
        let mut t_bad = t_utc.clone();
        t_bad[3] = f64::NAN;
        t_bad[22] = f64::NAN;
        let st_bad = receiver_stretches(&pos, &t_bad, 10);
        assert_eq!(kept_ranges(&fixed_source_at(-20.0), &st_bad, EPOCH, 0.0), [0..10, 20..25]);
        t_bad[15] = f64::NAN;
        let st_bad = receiver_stretches(&pos, &t_bad, 10);
        assert_eq!(kept_ranges(&fixed_source_at(-20.0), &st_bad, EPOCH, 0.0), [0..25]);
    }
}
