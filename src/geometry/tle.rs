//! SGP4 propagation from a TLE pair via the `sgp4` crate.

use ndarray::{Array1, Array2};

use super::Trajectory;

/// A satellite trajectory propagated by SGP4 from a TLE.
pub struct TleTrajectory {
    sat: sgp4::Constants,
    epoch_unix: f64,
}

impl TleTrajectory {
    /// Build from the two TLE lines (the optional name line is not used).
    pub fn from_tle(line1: &str, line2: &str) -> anyhow::Result<Self> {
        let elements = sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes())
            .map_err(|e| anyhow::anyhow!("TLE parse: {:?}", e))?;
        let epoch_unix = elements.datetime.and_utc().timestamp() as f64
            + elements.datetime.and_utc().timestamp_subsec_micros() as f64 * 1e-6;
        let sat = sgp4::Constants::from_elements(&elements)
            .map_err(|e| anyhow::anyhow!("sgp4 constants from elements: {:?}", e))?;
        Ok(Self { sat, epoch_unix })
    }
}

impl Trajectory for TleTrajectory {
    fn position_ecef(&self, t_unix: f64) -> [f64; 3] {
        // sgp4 wants minutes-since-epoch; output is TEME (km, km/s). For
        // metre-precision DGSN work we'd want a proper TEME→ECEF rotation;
        // for our research-grade tooling, sgp4's TEME ≈ ECI ≈ ECEF after
        // a sidereal rotation. Implementing the proper conversion is on
        // the to-do list; for now we return the TEME position scaled to m.
        let dt_min = (t_unix - self.epoch_unix) / 60.0;
        let pred = self.sat.propagate(sgp4::MinutesSinceEpoch(dt_min))
            .expect("sgp4 propagation failed");
        [pred.position[0] * 1000.0, pred.position[1] * 1000.0, pred.position[2] * 1000.0]
    }

    fn velocity_ecef(&self, t_unix: f64) -> [f64; 3] {
        let dt_min = (t_unix - self.epoch_unix) / 60.0;
        let pred = self.sat.propagate(sgp4::MinutesSinceEpoch(dt_min))
            .expect("sgp4 propagation failed");
        [pred.velocity[0] * 1000.0, pred.velocity[1] * 1000.0, pred.velocity[2] * 1000.0]
    }

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
