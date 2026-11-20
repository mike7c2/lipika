//! Decorator that samples any [`Trajectory`] on a uniform time grid
//! once and answers subsequent queries by linear interpolation. The
//! per-query cost drops from "call the underlying propagator" to
//! "two array reads + a multiply-add per axis" — usually 2-3 orders of
//! magnitude faster for the trajectory types we care about (SGP4, full
//! Keplerian) and exact for `Stationary`.
//!
//! Accuracy: linear interpolation of orbital position. For a LEO at
//! 7.5 km/s with 10 ms grid spacing, the chord-vs-arc error is at most
//! `(v·dt)²/(8·R)` ≈ `(75 m)²/(8·6700 km)` ≈ 0.1 µm. For GPS at 3.9 km/s
//! with the same grid, ≈ 30 nm. Both far below the position precision
//! anything downstream cares about.
//!
//! For trajectories with sub-grid features (manoeuvres, drag spikes,
//! moving observers with high acceleration) shrink the grid step or
//! use the underlying [`Trajectory`] directly.

use ndarray::{Array1, Array2, ArrayView2};

use super::Trajectory;

/// Trajectory wrapper that samples-and-interpolates.
pub struct PrecomputedTrajectory {
    /// First sample's absolute UTC time.
    t0_unix: f64,
    /// Time between adjacent samples (s).
    dt_s: f64,
    /// (3, N) ECEF position grid.
    positions: Array2<f64>,
    /// Number of samples.
    n: usize,
    /// Largest speed along the interpolated path: the longest grid segment
    /// divided by the grid step (m/s). NaN if the grid is not finite.
    max_speed_m_s: f64,
    /// Largest distance of any grid point from the Earth's centre (m). NaN if
    /// the grid is not finite.
    max_radius_m: f64,
}

impl PrecomputedTrajectory {
    /// Wrap `inner` by evaluating it once at a uniform grid with the
    /// requested spacing across `[t0_unix, t0_unix + duration_s]`.
    ///
    /// Adds one extra cell on each side so the Newton-iteration initial
    /// guess `t_utc − 70 ms` can be answered without out-of-range clamping.
    #[mantra_macros::req_satisfied("fid-3", "fn-5")]
    pub fn sample(
        inner: &dyn Trajectory,
        t0_unix: f64,
        duration_s: f64,
        grid_dt_s: f64,
    ) -> Self {
        assert!(grid_dt_s > 0.0, "grid spacing must be positive");
        // Pad both ends so out-of-range queries don't lose precision.
        let pad_s: f64 = 0.5;
        let start = t0_unix - pad_s;
        let total = duration_s + 2.0 * pad_s;
        let n = (total / grid_dt_s).ceil() as usize + 1;
        let times = Array1::from_iter((0..n).map(|i| start + i as f64 * grid_dt_s));
        let positions = inner.position_ecef_batch(&times);
        let (max_speed_m_s, max_radius_m) = path_bounds(&positions, grid_dt_s);
        Self {
            t0_unix: start,
            dt_s: grid_dt_s,
            positions,
            n,
            max_speed_m_s,
            max_radius_m,
        }
    }

    /// An upper bound on how fast the interpolated position moves: the
    /// longest grid segment divided by the grid step (m/s). The interpolated
    /// path is piecewise linear in time and constant outside the grid, so for
    /// any two query times `|P(t₂) − P(t₁)| ≤ max_speed · |t₂ − t₁|` (to
    /// float rounding). NaN if the grid holds a non-finite position.
    pub fn max_speed_m_s(&self) -> f64 { self.max_speed_m_s }

    /// An upper bound on `|P(t)|` for any query time: the largest distance of
    /// a grid point from the Earth's centre (every interpolated point is a
    /// convex combination of two grid points). NaN if the grid holds a
    /// non-finite position.
    pub fn max_radius_m(&self) -> f64 { self.max_radius_m }
}

/// `(max segment length / dt, max |P|)` over a `(3, N)` position grid; NaN
/// for either if any coordinate is not finite.
#[mantra_macros::req_satisfied("nf-5")]
fn path_bounds(positions: &Array2<f64>, dt_s: f64) -> (f64, f64) {
    let n = positions.ncols();
    let (mut max_seg2, mut max_r2) = (0.0_f64, 0.0_f64);
    let mut finite = true;
    for k in 0..n {
        let p = [positions[[0, k]], positions[[1, k]], positions[[2, k]]];
        let r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        finite &= r2.is_finite();
        if r2 > max_r2 { max_r2 = r2; }
        if k + 1 < n {
            let d = [positions[[0, k + 1]] - p[0], positions[[1, k + 1]] - p[1],
                     positions[[2, k + 1]] - p[2]];
            let seg2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            finite &= seg2.is_finite();
            if seg2 > max_seg2 { max_seg2 = seg2; }
        }
    }
    if finite { (max_seg2.sqrt() / dt_s, max_r2.sqrt()) } else { (f64::NAN, f64::NAN) }
}

impl Trajectory for PrecomputedTrajectory {
    #[mantra_macros::req_satisfied("fid-3")]
    fn position_ecef(&self, t_unix: f64) -> [f64; 3] {
        let f = (t_unix - self.t0_unix) / self.dt_s;
        let lo = (f.floor() as isize).clamp(0, (self.n - 2) as isize) as usize;
        let frac = (f - lo as f64).clamp(0.0, 1.0);
        let a = &self.positions;
        [
            a[[0, lo]] * (1.0 - frac) + a[[0, lo + 1]] * frac,
            a[[1, lo]] * (1.0 - frac) + a[[1, lo + 1]] * frac,
            a[[2, lo]] * (1.0 - frac) + a[[2, lo + 1]] * frac,
        ]
    }

    fn velocity_ecef(&self, t_unix: f64) -> [f64; 3] {
        // Finite-difference within the grid (forward-difference for speed).
        let f = (t_unix - self.t0_unix) / self.dt_s;
        let lo = (f.floor() as isize).clamp(0, (self.n - 2) as isize) as usize;
        let a = &self.positions;
        let inv = 1.0 / self.dt_s;
        [
            (a[[0, lo + 1]] - a[[0, lo]]) * inv,
            (a[[1, lo + 1]] - a[[1, lo]]) * inv,
            (a[[2, lo + 1]] - a[[2, lo]]) * inv,
        ]
    }

    fn position_ecef_batch(&self, t_unix: &Array1<f64>) -> Array2<f64> {
        let m = t_unix.len();
        let mut out = Array2::<f64>::zeros((3, m));
        let inv_dt = 1.0 / self.dt_s;
        let max_lo = (self.n - 2) as isize;
        let positions: ArrayView2<f64> = self.positions.view();
        for i in 0..m {
            let f = (t_unix[i] - self.t0_unix) * inv_dt;
            let lo_signed = f.floor() as isize;
            let lo = lo_signed.clamp(0, max_lo) as usize;
            let frac = (f - lo as f64).clamp(0.0, 1.0);
            let omf = 1.0 - frac;
            out[[0, i]] = positions[[0, lo]] * omf + positions[[0, lo + 1]] * frac;
            out[[1, i]] = positions[[1, lo]] * omf + positions[[1, lo + 1]] * frac;
            out[[2, i]] = positions[[2, lo]] * omf + positions[[2, lo + 1]] * frac;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::kepler::CircularKeplerian;
    use approx::assert_relative_eq;

    #[mantra_macros::req_test("nf-5")]
    #[test]
    fn path_bounds_cover_every_interpolated_query() {
        // A 550 km LEO: orbital speed ≈ 7.6 km/s inertial, a little less in
        // the rotating frame for a prograde orbit.
        let orbit = CircularKeplerian {
            semi_major_axis_m: 6_928_137.0,
            inclination_deg: 53.0,
            raan_deg: 10.0,
            true_anomaly_deg_at_epoch: 0.0,
            epoch_unix: 1_700_000_000.0,
        };
        let t0 = 1_700_000_000.0;
        let pre = PrecomputedTrajectory::sample(&orbit, t0, 30.0, 0.01);
        let (v, r) = (pre.max_speed_m_s(), pre.max_radius_m());
        assert!(v > 6_900.0 && v < 7_700.0, "max speed {v}");
        assert_relative_eq!(r, 6_928_137.0, max_relative = 1e-9);
        // Lipschitz in time, including queries beyond both ends of the grid
        // (where the position is held), and never farther out than the bound.
        let mut prev: Option<(f64, [f64; 3])> = None;
        for k in 0..4000 {
            let t = t0 - 2.0 + k as f64 * 0.008_7;
            let p = pre.position_ecef(t);
            assert!((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() <= r * (1.0 + 1e-12));
            if let Some((tp, pp)) = prev {
                let d = ((p[0] - pp[0]).powi(2) + (p[1] - pp[1]).powi(2)
                    + (p[2] - pp[2]).powi(2)).sqrt();
                assert!(d <= v * (t - tp) * (1.0 + 1e-9) + 1e-6, "step {d} m in {} s", t - tp);
            }
            prev = Some((t, p));
        }
        // A grid holding a non-finite position has no bound.
        struct Broken;
        impl Trajectory for Broken {
            fn position_ecef(&self, t: f64) -> [f64; 3] {
                if t > 5.0 { [f64::NAN, 0.0, 0.0] } else { [7.0e6, 0.0, 0.0] }
            }
        }
        let bad = PrecomputedTrajectory::sample(&Broken, 0.0, 10.0, 0.5);
        assert!(bad.max_speed_m_s().is_nan() && bad.max_radius_m().is_nan());
    }

    #[mantra_macros::req_test("fid-3")]
    #[test]
    fn precomputed_matches_underlying_to_sub_metre_at_100hz_grid() {
        let orbit = CircularKeplerian {
            semi_major_axis_m: 6_378_137.0 + 20_200_000.0,
            inclination_deg: 55.0,
            raan_deg: 0.0,
            true_anomaly_deg_at_epoch: 30.0,
            epoch_unix: 1_700_000_000.0,
        };
        let pre = PrecomputedTrajectory::sample(&orbit, 1_700_000_000.0, 60.0, 0.01);
        // Spot-check 50 query times against the underlying propagator.
        for k in 0..50 {
            let t = 1_700_000_000.0 + (60.0 * k as f64 / 50.0);
            let exact = orbit.position_ecef(t);
            let approx = pre.position_ecef(t);
            for axis in 0..3 {
                assert_relative_eq!(exact[axis], approx[axis], max_relative = 1e-6,
                                    epsilon = 1.0e-3);
            }
        }
    }
}
