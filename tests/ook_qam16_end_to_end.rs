//! OOK and 16-QAM through the whole chain (TOML → signal → channel → IQ).
//!
//! The scene is stationary and noise-free to working precision, so every
//! output sample is one transmitted symbol scaled by a common path gain and
//! rotated by a common carrier phase. The envelope power therefore takes only
//! the constellation's discrete levels, and a CW source in the same geometry
//! calibrates what the configured transmit power means for each kind.

use lipika::config::SceneConfig;
use lipika::pipeline::run_scene;

const TOML_TEMPLATE: &str = r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.05

[[transmitter]]
name = "src"
carrier_hz = 1_000_000_000.0
tx_power_dbw = 30.0
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 0.0 }
@SIGNAL@

[[receiver]]
name = "lab"
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 100000.0 }
clock = { allan_deviation_1s = 0.0 }
[[receiver.band]]
band_id = "test"
lo_carrier_hz = 1_000_000_000.0
sample_rate_hz = 1_000_000.0
noise_power_w = 1.0e-30
"#;

const CW: &str = r#"kind = "cw"
baseband_offset_hz = 0.0"#;

const OOK: &str = r#"kind = "ook"
bit_rate_bps = 100_000.0
payload_bits = 6000
payload_seed = 7
sync_word = [1, 0, 1, 1, 0, 0, 1, 0]"#;

const QAM16: &str = r#"kind = "qam16"
symbol_rate_baud = 100_000.0
payload_symbols = 6000
payload_seed = 7
sync_word = [0, 0, 0, 0, 1, 0, 1, 0]"#;

/// Per-sample envelope power once the signal has arrived (the first ~334
/// samples precede the 100 km light-time and are silent).
fn envelope_power(signal: &str) -> Vec<f64> {
    let scene: SceneConfig = toml::from_str(&TOML_TEMPLATE.replace("@SIGNAL@", signal)).unwrap();
    let iq = run_scene(&scene, /*seed=*/ 1).unwrap().remove(0).iq;
    iq.iter().skip(1000).map(|v| (v.re as f64).powi(2) + (v.im as f64).powi(2)).collect()
}

fn mean(p: &[f64]) -> f64 { p.iter().sum::<f64>() / p.len() as f64 }

#[mantra_macros::req_test("fn-4-ook")]
#[test]
fn ook_keys_the_carrier_between_off_and_cw_power() {
    let cw = mean(&envelope_power(CW));
    let p = envelope_power(OOK);
    let mut on = 0usize;
    for &v in &p {
        let r = v / cw;
        if (r - 1.0).abs() < 1e-3 {
            on += 1;
        } else {
            assert!(r < 1e-6, "OOK sample at {r:.4}× CW power is neither on nor off");
        }
    }
    // Balanced pseudo-random payload: the carrier is on about half the time.
    let on_fraction = on as f64 / p.len() as f64;
    assert!((on_fraction - 0.5).abs() < 0.05, "on fraction {on_fraction:.3}");
}

#[mantra_macros::req_test("fn-4-qam16")]
#[test]
fn qam16_envelope_sits_on_the_three_rings_at_unit_mean_power() {
    let cw = mean(&envelope_power(CW));
    let p = envelope_power(QAM16);
    // |s|² ∈ {2, 10, 18}/10 for the inner, middle and corner points.
    let rings = [0.2, 1.0, 1.8];
    let mut counts = [0usize; 3];
    for &v in &p {
        let r = v / cw;
        let k = rings.iter().position(|ring| (r - ring).abs() < 1e-3)
            .unwrap_or_else(|| panic!("16-QAM sample at {r:.4}× CW power is off-constellation"));
        counts[k] += 1;
    }
    // 4 inner, 8 middle, 4 corner points, equiprobable.
    let n = p.len() as f64;
    for (k, expect) in [0.25, 0.5, 0.25].into_iter().enumerate() {
        let frac = counts[k] as f64 / n;
        assert!((frac - expect).abs() < 0.05, "ring {} fraction {frac:.3}", rings[k]);
    }
    // Average power equals the configured transmit power (the CW reference).
    assert!((mean(&p) / cw - 1.0).abs() < 0.05, "mean power {:.4}× CW", mean(&p) / cw);
}
