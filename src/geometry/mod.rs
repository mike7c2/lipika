//! Geometry: positions, trajectories, frame conversions.
//!
//! Internal representation: ECEF metres. Coordinate-frame conversions
//! (geodetic, ECI) are at the input/output boundaries.

use ndarray::{Array1, Array2};

pub mod ecef;
pub mod kepler;
pub mod precomputed;
pub mod stationary;
pub mod tle;

/// Enum-dispatched trajectory for the perf-critical hot path.
///
/// The pipeline wraps every transmitter trajectory into one of these
/// variants (almost always `Precomputed` after the trajectory-precompute
/// step). The Newton-loop inner call goes through a `match` on this
/// enum rather than a `&dyn Trajectory` vtable, which the compiler can
/// inline and (where possible) autovectorise.
///
/// New trajectory types extend the `Trajectory` trait freely; the
/// pipeline downcasts/wraps known concrete types here for speed, and
/// falls back to a boxed-trait variant for everything else.
pub enum AnyTrajectory {
    Stationary(stationary::Stationary),
    Precomputed(precomputed::PrecomputedTrajectory),
    /// Catch-all for trajectories the pipeline hasn't been specialised
    /// for. Each call pays a vtable hop; fine for cold paths and tests.
    Dyn(std::sync::Arc<dyn Trajectory>),
}

impl AnyTrajectory {
    /// One-sample position query. The match is fully visible to the
    /// compiler so the variant body inlines into callers.
    #[mantra_macros::req_satisfied("fn-5")]
    #[inline(always)]
    pub fn position_ecef(&self, t_unix: f64) -> [f64; 3] {
        match self {
            AnyTrajectory::Stationary(s) => s.position_ecef(t_unix),
            AnyTrajectory::Precomputed(p) => p.position_ecef(t_unix),
            AnyTrajectory::Dyn(a) => a.position_ecef(t_unix),
        }
    }

    /// Batched query (kept so existing callers that take whole arrays
    /// still work).
    pub fn position_ecef_batch(&self, t_unix: &Array1<f64>) -> Array2<f64> {
        match self {
            AnyTrajectory::Stationary(s) => s.position_ecef_batch(t_unix),
            AnyTrajectory::Precomputed(p) => p.position_ecef_batch(t_unix),
            AnyTrajectory::Dyn(a) => a.position_ecef_batch(t_unix),
        }
    }
}

/// A trajectory whose ECEF position can be queried at any UTC time.
///
/// Implementations should be cheap to evaluate (precomputed tables OK).
/// Velocity is optional; if not provided, finite-differencing the position
/// is fine for diagnostics.
pub trait Trajectory: Send + Sync {
    /// ECEF position at a single time. Returns `[x, y, z]` in metres.
    fn position_ecef(&self, t_unix: f64) -> [f64; 3];

    /// ECEF velocity at a single time. Returns `[vx, vy, vz]` in m/s.
    /// Default impl uses central differencing if not overridden.
    fn velocity_ecef(&self, t_unix: f64) -> [f64; 3] {
        let h = 0.5; // 0.5 s differencing step; conservative
        let pp = self.position_ecef(t_unix + h);
        let pm = self.position_ecef(t_unix - h);
        [(pp[0] - pm[0]) / (2.0 * h),
         (pp[1] - pm[1]) / (2.0 * h),
         (pp[2] - pm[2]) / (2.0 * h)]
    }

    /// Batched query — important for performance. Default loops, but
    /// implementations like TLE/Kepler should override for vectorisation.
    fn position_ecef_batch(&self, t_unix: &Array1<f64>) -> Array2<f64> {
        let n = t_unix.len();
        let mut out = Array2::<f64>::zeros((3, n));
        for (i, &t) in t_unix.iter().enumerate() {
            let p = self.position_ecef(t);
            out[[0, i]] = p[0];
            out[[1, i]] = p[1];
            out[[2, i]] = p[2];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[mantra_macros::req_test("fn-5")]
    #[test]
    fn any_trajectory_dispatch_matches_the_wrapped_trajectory() {
        let t0 = 1_700_000_000.0;
        let times = Array1::from(vec![t0, t0 + 1.5, t0 + 9.99]);
        let orbit = kepler::CircularKeplerian {
            semi_major_axis_m: 26_578_137.0,
            inclination_deg: 55.0,
            raan_deg: 0.0,
            true_anomaly_deg_at_epoch: 30.0,
            epoch_unix: t0,
        };
        let fixed = stationary::Stationary::from_geodetic(55.86, -4.25, 10.0);
        let pre = precomputed::PrecomputedTrajectory::sample(&orbit, t0, 10.0, 0.01);

        let expect_orbit: Vec<[f64; 3]> = times.iter().map(|&t| orbit.position_ecef(t)).collect();
        let expect_pre: Vec<[f64; 3]> = times.iter().map(|&t| pre.position_ecef(t)).collect();
        let expect_fixed = vec![fixed.position_ecef(t0); 3];

        let wrapped = [
            (AnyTrajectory::Stationary(fixed), expect_fixed),
            (AnyTrajectory::Precomputed(pre), expect_pre),
            (AnyTrajectory::Dyn(std::sync::Arc::new(orbit)), expect_orbit.clone()),
        ];
        for (any, expected) in &wrapped {
            let batch = any.position_ecef_batch(&times);
            for (i, &t) in times.iter().enumerate() {
                let p = any.position_ecef(t);
                assert_eq!(p, expected[i]);
                for axis in 0..3 {
                    assert_abs_diff_eq!(batch[[axis, i]], p[axis], epsilon = 1e-6);
                }
            }
        }
        // A moving source moves; an Earth-fixed one does not.
        assert_ne!(expect_orbit[0], expect_orbit[2]);
    }
}
