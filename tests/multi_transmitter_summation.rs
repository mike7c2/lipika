//! Regression: transmitter contributions sum coherently in the band
//! accumulator, then noise is added once at the receiver's thermal floor.
//!
//! Verifies the architectural contract:
//!   iq_band(t) = Σ_tx s_tx(t) + n(t)
//! with ONE noise realisation per band (not per-transmitter). This is
//! load-bearing for every downstream SNR claim — if noise were
//! re-realised per transmitter, the band's effective noise floor would
//! silently scale with transmitter count.

use lipika::config::SceneConfig;
use lipika::pipeline::run_scene;

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

// A BPSK emitter at a realistic Unix epoch. Regression guard: the bit index
// is `floor((t_emit − signal_epoch)·bit_rate)`, so if the signal epoch is not
// the scene epoch (it was once hardcoded to 0.0), every emission time lands
// astronomically out of range and the baseband is all zeros — the band then
// collapses to the bare noise floor. We assert the band carries real signal.
const TOML_BPSK: &str = r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.05

[[transmitter]]
name = "beacon"
carrier_hz = 1_000_000_000.0
tx_power_dbw = 30.0
kind = "bpsk"
bit_rate_bps = 100_000.0
payload_bits = 6000
payload_seed = 7
sync_word = [1, -1, 1, 1, -1, -1, 1, -1]
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 0.0 }

[[receiver]]
name = "lab"
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 100000.0 }
clock = { allan_deviation_1s = 0.0 }
[[receiver.band]]
band_id = "test"
lo_carrier_hz = 1_000_000_000.0
sample_rate_hz = 1_000_000.0
noise_power_w = 1.0e-18
"#;

// Two co-located receivers see one CW source; they differ only by a 130 ns
// clock offset. The carrier-phase difference between their IQ must encode that
// 130 ns. Regression guard for the f64 absolute-time precision bug: doing the
// timing maths on absolute Unix times (~1.7e9 s) quantises to the ~238 ns ULP,
// so a 130 ns offset would collapse to 0 (or snap to 238 ns). Working in
// epoch-relative time preserves it.
const TOML_CW_TWO_CLOCKS: &str = r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.02

[[transmitter]]
name = "cw"
carrier_hz = 1000000.0
tx_power_dbw = 0.0
kind = "cw"
baseband_offset_hz = 0.0
trajectory = { kind = "stationary_ecef", x_m = 6478137.0, y_m = 0.0, z_m = 0.0 }

[[receiver]]
name = "rxa"
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 0.0 }
clock = { initial_time_offset_s = 0.0 }
[[receiver.band]]
band_id = "b"
lo_carrier_hz = 1000000.0
sample_rate_hz = 100000.0
noise_power_w = 1.0e-30

[[receiver]]
name = "rxb"
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 0.0 }
clock = { initial_time_offset_s = 1.3e-7 }
[[receiver.band]]
band_id = "b"
lo_carrier_hz = 1000000.0
sample_rate_hz = 100000.0
noise_power_w = 1.0e-30
"#;

#[test]
fn clock_offset_survives_at_realistic_epoch() {
    let scene: SceneConfig = toml::from_str(TOML_CW_TWO_CLOCKS).unwrap();
    let outs = run_scene(&scene, /*seed=*/ 7).unwrap();
    let iq = |name: &str| &outs.iter().find(|o| o.name == name).unwrap().iq;
    let a = iq("rxb__b"); // Δt = 130 ns
    let b = iq("rxa__b"); // Δt = 0
    // mean(a · conj(b)); the carrier phase delay term is identical for the two
    // co-located receivers, so the residual phase is −2π·f·(Δt_a − Δt_b).
    let (mut acc_re, mut acc_im) = (0.0_f64, 0.0_f64);
    for i in 0..a.len() {
        let (ar, ai) = (a[i].re as f64, a[i].im as f64);
        let (br, bi) = (b[i].re as f64, b[i].im as f64);
        acc_re += ar * br + ai * bi;
        acc_im += ai * br - ar * bi;
    }
    let phase = acc_im.atan2(acc_re);
    let f = 1.0e6;
    let recovered = -phase / (std::f64::consts::TAU * f); // f·Δt = 0.13 < 0.5 ⇒ unambiguous
    assert!((recovered - 1.3e-7).abs() < 1e-9,
            "recovered clock offset {:.2} ns != 130 ns — absolute-time quantisation regressed?",
            recovered * 1e9);
}

#[test]
fn bpsk_radiates_at_realistic_epoch() {
    let scene: SceneConfig = toml::from_str(TOML_BPSK).unwrap();
    let iq = run_scene(&scene, /*seed=*/ 1).unwrap().remove(0).iq;
    let mean_power: f64 =
        iq.iter().map(|v| (v.re as f64).powi(2) + (v.im as f64).powi(2)).sum::<f64>()
            / iq.len() as f64;
    // Noise floor is 1e-18 W; a 1 kW emitter at 100 km gives ~1e-11 W. If the
    // baseband were zeroed the mean power would sit at ~1e-18. Demand a margin
    // of many orders of magnitude so the test fails loudly if BPSK goes silent.
    assert!(mean_power > 1e-13,
            "BPSK band power {mean_power:.3e} W is at the noise floor — \
             the baseband is silent (signal-epoch regression?)");
}

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
