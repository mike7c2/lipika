//! Geometry: positions, trajectories, frame conversions.
//!
//! Internal representation: ECEF metres. Coordinate-frame conversions
//! (geodetic, ECI) are at the input/output boundaries.

use ndarray::{Array1, Array2};

pub mod ecef;
pub mod stationary;
pub mod kepler;
pub mod tle;

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
