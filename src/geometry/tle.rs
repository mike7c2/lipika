//! SGP4 propagation from a TLE pair via the `sgp4` crate, rotated TEME → ECEF.

use ndarray::{Array1, Array2};

use super::Trajectory;

/// Earth rotation rate (rad/s), matching the GPS NAV ephemeris algorithm so a
/// TLE trajectory and its derived broadcast ephemeris stay consistent.
const OMEGA_E: f64 = 7.292_115_146_7e-5;

/// Greenwich Mean Sidereal Time (rad) at a UTC unix time (IAU 1982; UT1≈UTC).
/// Sufficient for a self-consistent simulation — the same rotation is applied
/// to the signal geometry and to the trajectory the ephemeris is fitted from,
/// so sub-arcsecond frame approximations cancel.
fn gmst_rad(t_unix: f64) -> f64 {
    let jd = t_unix / 86400.0 + 2_440_587.5;
    let tc = (jd - 2_451_545.0) / 36_525.0;
    let gmst_sec = 67_310.548_41
        + (876_600.0 * 3600.0 + 8_640_184.812_866) * tc
        + 0.093_104 * tc * tc
        - 6.2e-6 * tc * tc * tc;
    (gmst_sec / 86_400.0).rem_euclid(1.0) * std::f64::consts::TAU
}

/// A satellite trajectory propagated by SGP4 from a TLE.
pub struct TleTrajectory {
    sat: sgp4::Constants,
    epoch_unix: f64,
}

impl TleTrajectory {
    /// Build from the two TLE lines (the optional name line is not used).
    #[mantra_macros::req_satisfied("fn-5-tle")]
    pub fn from_tle(line1: &str, line2: &str) -> anyhow::Result<Self> {
        let elements = sgp4::Elements::from_tle(None, line1.as_bytes(), line2.as_bytes())
            .map_err(|e| anyhow::anyhow!("TLE parse: {:?}", e))?;
        let epoch_unix = elements.datetime.and_utc().timestamp() as f64
            + elements.datetime.and_utc().timestamp_subsec_micros() as f64 * 1e-6;
        let sat = sgp4::Constants::from_elements(&elements)
            .map_err(|e| anyhow::anyhow!("sgp4 constants from elements: {:?}", e))?;
        Ok(Self { sat, epoch_unix })
    }

    /// UTC unix time of the TLE epoch.
    pub fn epoch_unix(&self) -> f64 { self.epoch_unix }
}

impl Trajectory for TleTrajectory {
    fn position_ecef(&self, t_unix: f64) -> [f64; 3] {
        // SGP4 wants minutes-since-epoch; output is TEME (km). Rotate about z
        // by GMST to Earth-fixed (ECEF), in metres.
        let dt_min = (t_unix - self.epoch_unix) / 60.0;
        let pred = self.sat.propagate(sgp4::MinutesSinceEpoch(dt_min))
            .expect("sgp4 propagation failed");
        let (s, c) = gmst_rad(t_unix).sin_cos();
        let (xt, yt, zt) = (pred.position[0] * 1000.0, pred.position[1] * 1000.0, pred.position[2] * 1000.0);
        [xt * c + yt * s, -xt * s + yt * c, zt]
    }

    fn velocity_ecef(&self, t_unix: f64) -> [f64; 3] {
        let dt_min = (t_unix - self.epoch_unix) / 60.0;
        let pred = self.sat.propagate(sgp4::MinutesSinceEpoch(dt_min))
            .expect("sgp4 propagation failed");
        let (s, c) = gmst_rad(t_unix).sin_cos();
        let (xt, yt) = (pred.position[0] * 1000.0, pred.position[1] * 1000.0);
        let (vxt, vyt, vzt) = (pred.velocity[0] * 1000.0, pred.velocity[1] * 1000.0, pred.velocity[2] * 1000.0);
        let (xe, ye) = (xt * c + yt * s, -xt * s + yt * c);
        // v_ecef = R3(θ)·v_teme − ω⊕ × r_ecef
        [vxt * c + vyt * s + OMEGA_E * ye,
         -vxt * s + vyt * c - OMEGA_E * xe,
         vzt]
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

#[cfg(test)]
mod tests {
    use super::*;

    // A current GPS satellite (PRN 22).
    const L1: &str = "1 26407U 00040A   26150.58242625  .00000046  00000+0  00000+0 0  9999";
    const L2: &str = "2 26407  54.8538 215.8617 0121278 302.3883 230.2243  2.00557636189594";

    #[mantra_macros::req_test("fn-5-tle")]
    #[test]
    fn gps_tle_gives_gps_altitude_and_consistent_velocity() {
        let t = TleTrajectory::from_tle(L1, L2).unwrap();
        let t0 = t.epoch_unix;
        let p = t.position_ecef(t0);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        // GPS orbital radius ≈ 26 560 km (± eccentricity).
        assert!((r - 26_560_000.0).abs() < 1_500_000.0, "radius {r} not GPS-like");

        // Finite-difference position vs. analytic ECEF velocity should agree.
        let h = 0.1;
        let pp = t.position_ecef(t0 + h);
        let pm = t.position_ecef(t0 - h);
        let v_fd = [(pp[0] - pm[0]) / (2.0 * h), (pp[1] - pm[1]) / (2.0 * h), (pp[2] - pm[2]) / (2.0 * h)];
        let v = t.velocity_ecef(t0);
        for k in 0..3 {
            assert!((v[k] - v_fd[k]).abs() < 1.0, "velocity axis {k}: {} vs fd {}", v[k], v_fd[k]);
        }
    }
}
