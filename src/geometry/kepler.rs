//! Synthetic Keplerian propagator — port of the Python `gnss_constellation`
//! module. Two-body, no J2/drag. Used for deterministic GPS-like
//! constellations in unit tests and demo scenes; for real GPS use
//! [`super::tle`] with a current TLE/SP3 instead.

use ndarray::{Array1, Array2};

use super::Trajectory;

// IS-GPS-200 WGS-84 values, shared with the GPS NAV ephemeris algorithm so a
// circular-Keplerian orbit and its broadcast ephemeris agree exactly.
const GM_EARTH: f64 = 3.986_005e14;
const EARTH_ROTATION_RATE_RAD_S: f64 = 7.292_115_146_7e-5;

/// Minimal element set for a circular Keplerian orbit (no eccentricity,
/// no argument-of-perigee; gives a GPS-like circle that the propagator
/// can evaluate analytically without iteration).
// r[impl fn.5.keplerian]
#[derive(Debug, Clone, Copy)]
pub struct CircularKeplerian {
    pub semi_major_axis_m: f64,
    pub inclination_deg: f64,
    pub raan_deg: f64,
    pub true_anomaly_deg_at_epoch: f64,
    pub epoch_unix: f64,
}

impl CircularKeplerian {
    fn mean_motion_rad_s(&self) -> f64 {
        (GM_EARTH / self.semi_major_axis_m.powi(3)).sqrt()
    }

    fn eci_at(&self, t_unix: f64) -> ([f64; 3], [f64; 3]) {
        let r = self.semi_major_axis_m;
        let n = self.mean_motion_rad_s();
        let dt = t_unix - self.epoch_unix;
        let nu = self.true_anomaly_deg_at_epoch.to_radians() + n * dt;
        let (s_nu, c_nu) = nu.sin_cos();
        let (s_i, c_i) = self.inclination_deg.to_radians().sin_cos();
        let (s_o, c_o) = self.raan_deg.to_radians().sin_cos();

        // Perifocal -> ECI via R_z(RAAN) * R_x(incl), argument-of-perigee 0.
        let xp = r * c_nu;
        let yp = r * s_nu;
        let y1 = yp * c_i;
        let z1 = yp * s_i;
        let x_eci = xp * c_o - y1 * s_o;
        let y_eci = xp * s_o + y1 * c_o;
        let z_eci = z1;

        // Velocity in perifocal frame: r·dν/dt perpendicular to position.
        let v_mag = r * n;
        let vxp = -v_mag * s_nu;
        let vyp = v_mag * c_nu;
        let vy1 = vyp * c_i;
        let vz1 = vyp * s_i;
        let vx_eci = vxp * c_o - vy1 * s_o;
        let vy_eci = vxp * s_o + vy1 * c_o;
        let vz_eci = vz1;

        ([x_eci, y_eci, z_eci], [vx_eci, vy_eci, vz_eci])
    }

    /// Rotate an ECI vector into ECEF (z-axis rotation by GMST since epoch).
    /// Velocity transform: `v_ecef = R · v_eci − ω × r_ecef`.
    fn eci_to_ecef(&self, t_unix: f64, p_eci: [f64; 3], v_eci: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let theta = EARTH_ROTATION_RATE_RAD_S * (t_unix - self.epoch_unix);
        let (s, c) = theta.sin_cos();
        let xe = p_eci[0] * c + p_eci[1] * s;
        let ye = -p_eci[0] * s + p_eci[1] * c;
        let ze = p_eci[2];
        let mut vxe = v_eci[0] * c + v_eci[1] * s;
        let mut vye = -v_eci[0] * s + v_eci[1] * c;
        let vze = v_eci[2];
        vxe += EARTH_ROTATION_RATE_RAD_S * ye;
        vye -= EARTH_ROTATION_RATE_RAD_S * xe;
        ([xe, ye, ze], [vxe, vye, vze])
    }
}

impl Trajectory for CircularKeplerian {
    fn position_ecef(&self, t_unix: f64) -> [f64; 3] {
        let (p_eci, v_eci) = self.eci_at(t_unix);
        self.eci_to_ecef(t_unix, p_eci, v_eci).0
    }

    fn velocity_ecef(&self, t_unix: f64) -> [f64; 3] {
        let (p_eci, v_eci) = self.eci_at(t_unix);
        self.eci_to_ecef(t_unix, p_eci, v_eci).1
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
    use approx::assert_relative_eq;

    // r[verify fn.5.keplerian]
    #[test]
    fn gps_altitude_orbit_distance_correct() {
        let gps_radius = 6_378_137.0 + 20_200_000.0;
        let orbit = CircularKeplerian {
            semi_major_axis_m: gps_radius,
            inclination_deg: 55.0,
            raan_deg: 0.0,
            true_anomaly_deg_at_epoch: 0.0,
            epoch_unix: 0.0,
        };
        let p = orbit.position_ecef(0.0);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert_relative_eq!(r, gps_radius, max_relative = 1e-9);
    }
}
