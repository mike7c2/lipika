//! Stationary observer / source on the rotating Earth (ECEF-fixed).

use ndarray::{Array1, Array2};

use super::ecef::geodetic_to_ecef;
use super::Trajectory;

/// A point fixed in the ECEF frame.
// r[impl fn.5.stationary]
#[derive(Debug, Clone, Copy)]
pub struct Stationary {
    pub position_ecef: [f64; 3],
}

impl Stationary {
    pub fn from_geodetic(lat_deg: f64, lon_deg: f64, h_m: f64) -> Self {
        Self { position_ecef: geodetic_to_ecef(lat_deg, lon_deg, h_m) }
    }
}

impl Trajectory for Stationary {
    fn position_ecef(&self, _t_unix: f64) -> [f64; 3] {
        self.position_ecef
    }

    fn velocity_ecef(&self, _t_unix: f64) -> [f64; 3] {
        [0.0, 0.0, 0.0]
    }

    fn position_ecef_batch(&self, t_unix: &Array1<f64>) -> Array2<f64> {
        let n = t_unix.len();
        let mut out = Array2::<f64>::zeros((3, n));
        for i in 0..n {
            out[[0, i]] = self.position_ecef[0];
            out[[1, i]] = self.position_ecef[1];
            out[[2, i]] = self.position_ecef[2];
        }
        out
    }
}
