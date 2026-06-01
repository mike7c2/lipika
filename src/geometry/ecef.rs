//! WGS84 ↔ ECEF conversions.

const WGS84_A: f64 = 6_378_137.0;  // semi-major axis
const WGS84_F: f64 = 1.0 / 298.257_223_563;
const WGS84_E2: f64 = WGS84_F * (2.0 - WGS84_F);

/// Convert geodetic (latitude °, longitude °, elevation m above ellipsoid)
/// to ECEF (m).
pub fn geodetic_to_ecef(lat_deg: f64, lon_deg: f64, h_m: f64) -> [f64; 3] {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
    let x = (n + h_m) * cos_lat * lon.cos();
    let y = (n + h_m) * cos_lat * lon.sin();
    let z = (n * (1.0 - WGS84_E2) + h_m) * sin_lat;
    [x, y, z]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn glasgow_matches_reference() {
        // Cross-checked against a reference WGS-84 geodetic→ECEF implementation:
        //   geodetic_to_ecef(55.86, -4.25, 10.0)
        //   -> [3577894.5228608, -265883.49155448, 5255718.21468305]
        let p = geodetic_to_ecef(55.86, -4.25, 10.0);
        assert_relative_eq!(p[0],  3_577_894.5228608, max_relative = 1e-9);
        assert_relative_eq!(p[1],   -265_883.49155448, max_relative = 1e-9);
        assert_relative_eq!(p[2],  5_255_718.21468305, max_relative = 1e-9);
    }

    #[test]
    fn equator_at_zero_longitude_lies_on_x_axis() {
        let p = geodetic_to_ecef(0.0, 0.0, 0.0);
        assert_relative_eq!(p[0], WGS84_A, max_relative = 1e-12);
        assert_relative_eq!(p[1], 0.0, epsilon = 1e-6);
        assert_relative_eq!(p[2], 0.0, epsilon = 1e-6);
    }
}
