//! Regression: transmitter contributions sum coherently in the band
//! accumulator, then noise is added once at the receiver's thermal floor.
//!
//! Verifies the architectural contract:
//!   iq_band(t) = Σ_tx s_tx(t) + n(t)
//! with ONE noise realisation per band (not per-transmitter). This is
//! load-bearing for every downstream SNR claim — if noise were
//! re-realised per transmitter, the band's effective noise floor would
//! silently scale with transmitter count.

use dgsn_iq_gen::config::SceneConfig;
use dgsn_iq_gen::pipeline::run_scene;

const TOML_NO_TX: &str = r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.5

[[receiver]]
name = "lab"
trajectory = { kind = "stationary_geodetic", lat_deg = 55.86, lon_deg = -4.25, h_m = 10.0 }
clock = { allan_deviation_1s = 0.0 }
[[receiver.band]]
band_id = "test"
lo_carrier_hz = 1_000_000_000.0
sample_rate_hz = 1_000_000.0
noise_power_w = 1.0e-15
"#;

const TOML_TWO_TX_TEMPLATE: &str = r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.5

[[transmitter]]
name = "alpha"
carrier_hz = 1_000_000_000.0
tx_power_dbw = 30.0   # 1 kW just to push above the test noise floor
kind = "cw"
baseband_offset_hz = 0.0
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 1000.0 }

[[transmitter]]
name = "beta"
carrier_hz = 1_005_000_000.0   # different carrier so they don't beat
tx_power_dbw = 30.0
kind = "cw"
baseband_offset_hz = 0.0
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 1000.0, z_m = 1000.0 }

[[receiver]]
name = "lab"
trajectory = { kind = "stationary_geodetic", lat_deg = 55.86, lon_deg = -4.25, h_m = 10.0 }
clock = { allan_deviation_1s = 0.0 }
[[receiver.band]]
band_id = "test"
lo_carrier_hz = 1_000_000_000.0
sample_rate_hz = 1_000_000.0
noise_power_w = 1.0e-15
"#;

#[test]
fn noise_only_band_has_power_within_floor_tolerance() {
    let scene: SceneConfig = toml::from_str(TOML_NO_TX).unwrap();
    let outputs = run_scene(&scene, /*seed=*/ 12345).unwrap();
    assert_eq!(outputs.len(), 1);
    let iq = &outputs[0].iq;
    let mean_power: f64 =
        iq.iter().map(|v| (v.re as f64).powi(2) + (v.im as f64).powi(2)).sum::<f64>()
            / iq.len() as f64;
    assert!((mean_power / 1e-15 - 1.0).abs() < 0.02,
            "noise-only mean power {} W deviates >2% from 1e-15 W", mean_power);
}

#[test]
fn band_iq_is_sum_of_per_transmitter_runs_plus_one_noise_draw() {
    // The architectural contract:
    //   out_both = signal_alpha + signal_beta + noise
    //   out_alpha = signal_alpha + noise           (same noise realisation)
    //   out_beta  = signal_beta  + noise           (same noise realisation)
    //   out_noise = noise                          (same noise realisation)
    //
    // Therefore: out_both = out_alpha + out_beta - out_noise (sample-by-sample).
    //
    // This holds iff (a) signals add linearly into the band accumulator and
    // (b) AWGN is drawn ONCE per band, not per transmitter. The same RNG
    // seed across all four runs makes the noise realisation deterministic.

    let scene_both: SceneConfig = toml::from_str(TOML_TWO_TX_TEMPLATE).unwrap();
    let mut scene_alpha = scene_both.clone();
    scene_alpha.transmitters.retain(|t| t.name == "alpha");
    let mut scene_beta = scene_both.clone();
    scene_beta.transmitters.retain(|t| t.name == "beta");
    let scene_noise: SceneConfig = toml::from_str(TOML_NO_TX).unwrap();

    let seed = 99;
    let a    = run_scene(&scene_alpha, seed).unwrap().remove(0).iq;
    let b    = run_scene(&scene_beta,  seed).unwrap().remove(0).iq;
    let both = run_scene(&scene_both,  seed).unwrap().remove(0).iq;
    let n    = run_scene(&scene_noise, seed).unwrap().remove(0).iq;
    let len = a.len();
    assert_eq!(b.len(), len);
    assert_eq!(both.len(), len);
    assert_eq!(n.len(), len);

    // First sanity: noise-only run reproduces deterministically across
    // calls to run_scene (same seed → same RNG → same noise samples).
    // We don't verify this directly here; if the additivity test below
    // passes, determinism is implicit.

    // Sample-by-sample: out_both ≈ out_alpha + out_beta − out_noise.
    // Allow tiny absolute slop for f32 round-off.
    let mut max_abs_err: f32 = 0.0;
    for i in 0..len {
        let lhs = both[i];
        let rhs = a[i] + b[i] - n[i];
        let err = (lhs - rhs).norm();
        if err > max_abs_err { max_abs_err = err; }
    }
    let avg_signal = both.iter().map(|v| v.norm()).sum::<f32>() / len as f32;
    eprintln!("max |both - (a + b - n)|: {:.3e};  mean |both|: {:.3e}",
              max_abs_err, avg_signal);
    // The signals are ~1e-7 magnitude; round-off must be far below that.
    assert!(max_abs_err < 1e-10,
            "additivity violated: max abs err {:.3e} exceeds f32 round-off",
            max_abs_err);
}
