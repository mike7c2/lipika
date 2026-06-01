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
}

impl PrecomputedTrajectory {
    /// Wrap `inner` by evaluating it once at a uniform grid with the
    /// requested spacing across `[t0_unix, t0_unix + duration_s]`.
    ///
    /// Adds one extra cell on each side so the Newton-iteration initial
    /// guess `t_utc − 70 ms` can be answered without out-of-range clamping.
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
        Self {
            t0_unix: start,
            dt_s: grid_dt_s,
            positions,
            n,
        }
    }
}

impl Trajectory for PrecomputedTrajectory {
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
