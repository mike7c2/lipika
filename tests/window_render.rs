//! `pipeline::render_window`: one receiver, one band, one absolute-time window
//! of a sky, streamed to a sink.
//!
//! - it is the same render as the equivalent one-receiver scene;
//! - chunk size and thread count do not change it, and the sink can stop it;
//! - two windows that overlap in absolute time see the same emitted signal;
//! - a burst gate keys the received signal on and off at the light-delayed
//!   edges;
//! - an empty sky is the receiver's own noise.

use std::ops::ControlFlow;
use std::path::Path;

use lipika::clock::NodeClockSpec;
use lipika::config::{
    ReceiverBandConfig, ReceiverConfig, SceneConfig, SignalKind, SkyConfig, TrajectoryConfig,
};
use lipika::pipeline::{render_window, run_scene, WindowOptions};
use num_complex::Complex32;

const FS: f64 = 48_000.0;
const HIBER_HZ: f64 = 401_150_000.0;
/// A whole multiple of 15 s (so a gate period of the example sky starts here),
/// 300 s before the example pass culminates over Glasgow.
const PASS_START_UNIX: f64 = 1_717_260_690.0;

fn hiber_sky_text() -> String {
    std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/hiber_beacon_sky.toml")).unwrap()
}

fn hiber_sky() -> SkyConfig {
    SkyConfig::from_toml_str(&hiber_sky_text()).unwrap()
}

fn glasgow() -> ReceiverConfig {
    ReceiverConfig::stationary_geodetic("glasgow", 55.86, -4.25, 10.0)
}

fn options(seed: u64, chunk_samples: usize, n_threads: usize) -> WindowOptions {
    WindowOptions { seed, chunk_samples, n_threads }
}

/// Render a whole window into memory.
fn collect(
    sky: &SkyConfig,
    receiver: &ReceiverConfig,
    band: &ReceiverBandConfig,
    start_unix: f64,
    duration_s: f64,
    opts: &WindowOptions,
) -> Vec<Complex32> {
    let mut iq = Vec::new();
    let summary = render_window(sky, receiver, band, start_unix, duration_s, opts, |chunk| {
        iq.extend_from_slice(chunk);
        ControlFlow::Continue(())
    })
    .unwrap();
    assert!(!summary.stopped_early);
    assert_eq!(summary.n_samples_emitted, iq.len());
    assert_eq!(summary.n_samples_window, iq.len());
    iq
}

fn bytes(iq: &[Complex32]) -> Vec<u8> {
    iq.iter().flat_map(|v| [v.re.to_le_bytes(), v.im.to_le_bytes()]).flatten().collect()
}

fn mean_power(iq: &[Complex32]) -> f64 {
    iq.iter().map(|v| v.norm_sqr() as f64).sum::<f64>() / iq.len() as f64
}

#[mantra_macros::req_test("fn-10", "fn-9", "fn-8")]
#[test]
fn window_equals_the_equivalent_scene_render() {
    // The example sky (a gated, repeating BPSK beacon on a TLE, with an
    // explicit signal epoch) seen by a receiver with a non-trivial clock and
    // a thermal floor. The scene is the same thing written as a scene file:
    // epoch = window start, the sky's transmitters, the one receiver + band.
    // The window opens 0.3 s before a gate period, so it spans a burst edge.
    let start = PASS_START_UNIX - 0.3;
    let duration = 1.0;
    let scene_text = format!(
        r#"
[scene]
epoch_unix = {start}
duration_s = {duration}
{sky}
[[receiver]]
name = "glasgow"
trajectory = {{ kind = "stationary_geodetic", lat_deg = 55.86, lon_deg = -4.25, h_m = 10.0 }}
clock = {{ allan_deviation_1s = 1e-9, frequency_offset_y0 = 0.5e-6, initial_time_offset_s = 1e-7 }}
[[receiver.band]]
band_id = "uhf"
lo_carrier_hz = {HIBER_HZ}
sample_rate_hz = {FS}
noise_temperature_k = 290.0
"#,
        sky = hiber_sky_text(),
    );
    let scene: SceneConfig = toml::from_str(&scene_text).unwrap();
    let from_scene = run_scene(&scene, /*seed=*/ 99).unwrap().remove(0).iq.to_vec();

    let receiver = glasgow().with_clock(NodeClockSpec {
        allan_deviation_1s: 1e-9,
        flicker_floor: 0.0,
        frequency_offset_y0: 0.5e-6,
        initial_time_offset_s: 1e-7,
    });
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS).with_noise_temperature_k(290.0);
    let from_window = collect(&hiber_sky(), &receiver, &band, start, duration, &options(99, 5000, 1));

    assert_eq!(from_window.len(), 48_000);
    assert!(bytes(&from_window) == bytes(&from_scene),
        "window render differs from the equivalent scene render");
    // The burst is in there: the second half (keyed on) is well above the
    // noise-only first quarter.
    assert!(mean_power(&from_window[24_000..]) > 3.0 * mean_power(&from_window[..12_000]));
}

#[mantra_macros::req_test("fn-10", "fn-7", "nf-3")]
#[test]
fn chunk_size_and_thread_count_do_not_change_the_window() {
    let sky = hiber_sky();
    let receiver = glasgow().with_clock(NodeClockSpec {
        allan_deviation_1s: 1e-9,
        flicker_floor: 1e-10,
        frequency_offset_y0: -0.3e-6,
        initial_time_offset_s: 0.0,
    });
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS);
    let (start, duration) = (PASS_START_UNIX - 0.2, 0.8);

    let whole = collect(&sky, &receiver, &band, start, duration, &options(7, 1 << 20, 1));
    assert_eq!(whole.len(), 38_400);
    for (chunk, threads) in [(1usize, 1usize), (997, 1), (4096, 4), (38_400, 0), (0, 2)] {
        let other = collect(&sky, &receiver, &band, start, duration, &options(7, chunk, threads));
        assert!(bytes(&whole) == bytes(&other),
            "output changed with chunk_samples = {chunk}, n_threads = {threads}");
    }
    // Same arguments, same bytes; another seed, another noise realisation.
    assert!(bytes(&whole)
        == bytes(&collect(&sky, &receiver, &band, start, duration, &options(7, 1 << 20, 1))));
    assert!(bytes(&whole)
        != bytes(&collect(&sky, &receiver, &band, start, duration, &options(8, 1 << 20, 1))));
}

#[mantra_macros::req_test("fn-10")]
#[test]
fn chunks_arrive_in_order_and_the_sink_can_stop_the_render() {
    let sky = hiber_sky();
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS);
    let whole = collect(&sky, &glasgow(), &band, PASS_START_UNIX, 0.5, &options(3, 1000, 1));
    assert_eq!(whole.len(), 24_000);

    // Chunks are `chunk_samples` long (the last one shorter) and concatenate
    // to the window.
    let mut sizes = Vec::new();
    render_window(&sky, &glasgow(), &band, PASS_START_UNIX, 0.5, &options(3, 7000, 1), |c| {
        sizes.push(c.len());
        ControlFlow::Continue(())
    })
    .unwrap();
    assert_eq!(sizes, [7000, 7000, 7000, 3000]);

    // Break on the third chunk: exactly three chunks are delivered, and they
    // are the first 3000 samples of the full window.
    let mut got = Vec::new();
    let mut calls = 0;
    let summary =
        render_window(&sky, &glasgow(), &band, PASS_START_UNIX, 0.5, &options(3, 1000, 1), |c| {
            got.extend_from_slice(c);
            calls += 1;
            if calls == 3 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
        })
        .unwrap();
    assert_eq!(calls, 3);
    assert!(summary.stopped_early);
    assert_eq!(summary.n_samples_emitted, 3000);
    assert_eq!(summary.n_samples_window, 24_000);
    assert_eq!(summary.contributing_transmitters, ["hiber_beacon"]);
    assert!(bytes(&got) == bytes(&whole[..3000]));

    // Breaking on the very last chunk is a complete render, not an early stop.
    let summary =
        render_window(&sky, &glasgow(), &band, PASS_START_UNIX, 0.5, &options(3, 24_000, 1), |_| {
            ControlFlow::Break(())
        })
        .unwrap();
    assert!(!summary.stopped_early);
    assert_eq!(summary.n_samples_emitted, 24_000);
}

/// Largest per-sample difference between two captures of the same instants,
/// relative to the RMS of the keyed-on signal.
fn max_relative_difference(a: &[Complex32], b: &[Complex32]) -> f64 {
    assert_eq!(a.len(), b.len());
    let on: Vec<f64> =
        a.iter().map(|v| v.norm_sqr() as f64).filter(|&p| p > 0.0).collect();
    assert!(!on.is_empty(), "the overlap should contain signal");
    let rms = (on.iter().sum::<f64>() / on.len() as f64).sqrt();
    a.iter().zip(b).map(|(x, y)| (x - y).norm() as f64).fold(0.0, f64::max) / rms
}

#[mantra_macros::req_test("fn-9", "fn-8", "fn-10")]
#[test]
fn overlapping_windows_agree_on_a_leo_beacon() {
    // Two captures of the example beacon by the same (noise-free, ideal-clock)
    // receiver: A covers [T − 0.5, T + 2.5), B covers [T + 0.75, T + 3.0).
    // Over the shared 1.75 s — burst 1's tail is before it; burst 2
    // (T + 1.5 … T + 2.0) is inside it — the samples must agree. The only
    // difference allowed is the f64 step of absolute time in the orbit query
    // (≤ 119 ns ⇒ ≤ ~1 mm of range ⇒ ≤ ~0.01 rad of 401 MHz carrier).
    let sky = hiber_sky();
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS).with_noise_power_w(0.0);
    let a = collect(&sky, &glasgow(), &band, PASS_START_UNIX - 0.5, 3.0, &options(1, 4096, 1));
    let b = collect(&sky, &glasgow(), &band, PASS_START_UNIX + 0.75, 2.25, &options(2, 9000, 1));
    let lead = (1.25 * FS) as usize; // B starts 60 000 samples after A
    let overlap = a.len() - lead;
    assert_eq!(overlap, 84_000);
    let (a_ov, b_ov) = (&a[lead..], &b[..overlap]);

    let diff = max_relative_difference(a_ov, b_ov);
    assert!(diff < 0.03, "overlap differs by {diff} of the signal RMS");
    // Both saw the same on/off pattern, sample for sample: on for burst 2
    // (0.5 s = 24 000 samples), silent otherwise.
    let on_a: Vec<bool> = a_ov.iter().map(|v| v.norm() > 0.0).collect();
    let on_b: Vec<bool> = b_ov.iter().map(|v| v.norm() > 0.0).collect();
    assert!(on_a == on_b, "gate edges differ between the two windows");
    let n_on = on_a.iter().filter(|&&x| x).count();
    assert!((23_999..=24_001).contains(&n_on), "{n_on} keyed-on samples in the overlap");
    // And the data agrees, not just the envelope: the BPSK sign flips line up.
    let flips = |x: &[Complex32]| -> Vec<usize> {
        (1..x.len())
            .filter(|&i| x[i].norm() > 0.0 && x[i - 1].norm() > 0.0
                && (x[i] * x[i - 1].conj()).re < 0.0)
            .collect()
    };
    let (fa, fb) = (flips(a_ov), flips(b_ov));
    assert!(fa.len() > 500, "expected ~1200 bit transitions, saw {}", fa.len());
    assert_eq!(fa, fb);
}

#[mantra_macros::req_test("fn-9", "fn-10")]
#[test]
fn overlapping_windows_agree_for_every_anchored_signal_kind() {
    // Stationary emitters, so nothing but the emission itself can differ: one
    // transmitter of each kind that takes an anchor, each on its own band
    // centre, default anchoring (no signal_epoch_unix ⇒ Unix 0) except the
    // BPSK one, which names an epoch that is not a whole number of seconds.
    // Short frames, so both windows span many repeats.
    let sky = SkyConfig::from_toml_str(r#"
[[transmitter]]
name = "tone"
carrier_hz = 433_000_000.0
kind = "cw"
baseband_offset_hz = 1234.5
trajectory = { kind = "stationary_ecef", x_m = 3600000.0, y_m = -270000.0, z_m = 5300000.0 }

[[transmitter]]
name = "bpsk"
carrier_hz = 434_000_000.0
kind = "bpsk"
bit_rate_bps = 1200.0
payload_bits = 95
payload_seed = 4
sync_word = [-1, 1, 1, -1, 1]
signal_epoch_unix = 1000000007.25
trajectory = { kind = "stationary_ecef", x_m = 3600000.0, y_m = -270000.0, z_m = 5300000.0 }

[[transmitter]]
name = "ook"
carrier_hz = 435_000_000.0
kind = "ook"
bit_rate_bps = 2400.0
payload_bits = 61
payload_seed = 5
sync_word = [1, 0, 1]
trajectory = { kind = "stationary_ecef", x_m = 3600000.0, y_m = -270000.0, z_m = 5300000.0 }

[[transmitter]]
name = "qam"
carrier_hz = 436_000_000.0
kind = "qam16"
symbol_rate_baud = 4800.0
payload_symbols = 50
payload_seed = 6
sync_word = [1, 0, 1, 1]
trajectory = { kind = "stationary_ecef", x_m = 3600000.0, y_m = -270000.0, z_m = 5300000.0 }

[[transmitter]]
name = "code"
carrier_hz = 437_000_000.0
kind = "ca_code"
prn = 7
epoch_unix = 1000000000.0
nav = "off"
trajectory = { kind = "stationary_ecef", x_m = 3600000.0, y_m = -270000.0, z_m = 5300000.0 }
"#).unwrap();
    let start_a = 1_717_260_840.0;
    let lead_s = 0.75;
    for (name, lo, fs) in [
        ("tone", 433.0e6, 48_000.0), ("bpsk", 434.0e6, 48_000.0), ("ook", 435.0e6, 48_000.0),
        ("qam", 436.0e6, 48_000.0), ("code", 437.0e6, 4_000_000.0),
    ] {
        // Narrow acceptance window so only the one transmitter lands in band.
        let mut band = ReceiverBandConfig::new(name, lo, fs).with_noise_power_w(0.0);
        band.bandwidth_hz = Some(100_000.0f64.min(fs));
        let dur_a = if fs > 1.0e6 { 0.8 } else { 1.5 };
        let a = collect(&sky, &glasgow(), &band, start_a, dur_a, &options(1, 1 << 16, 1));
        let b = collect(&sky, &glasgow(), &band, start_a + lead_s, dur_a - lead_s + 0.1,
                        &options(1, 1 << 16, 1));
        let lead = (lead_s * fs) as usize;
        let overlap = a.len() - lead;
        let diff = max_relative_difference(&a[lead..], &b[..overlap]);
        assert!(diff < 1e-3, "{name}: overlap differs by {diff} of the signal RMS");
    }
    // The NAV overlays are scene-only: a window render refuses rather than
    // emit data bits that depend on where the window starts.
    let mut with_nav = sky.clone();
    for tx in &mut with_nav.transmitters {
        if let SignalKind::CaCode { nav, .. } = &mut tx.kind {
            *nav = lipika::config::NavMode::Dummy;
        }
    }
    let band = ReceiverBandConfig::new("code", 437.0e6, 4_000_000.0);
    let err = render_window(&with_nav, &glasgow(), &band, start_a, 0.001,
                            &WindowOptions::default(), |_| ControlFlow::Continue(()));
    assert!(err.is_err());
}

#[mantra_macros::req_test("fn-8", "eq-e11", "fn-10")]
#[test]
fn gate_keys_the_received_signal_at_the_light_delayed_edges() {
    // A gated tone 300 km straight above a receiver on the equator: one-way
    // delay 300 km / c = 1.000 69 ms = 48.03 samples at 48 kHz. Gate: 20 s
    // period starting at Unix multiples of 20 s plus 5 s, on for [2.0, 2.5)
    // and [7.0, 7.25). The window opens 1.9 s into a period.
    let sky = SkyConfig::from_toml_str(r#"
[[transmitter]]
name = "keyed_tone"
carrier_hz = 401_150_000.0
tx_power_dbw = 0.0
kind = "cw"
baseband_offset_hz = 0.0
gate = { period_s = 20, epoch_unix = 5.0, bursts = [
    { offset_s = 2.0, duration_s = 0.5 },
    { offset_s = 7.0, duration_s = 0.25 },
] }
trajectory = { kind = "stationary_ecef", x_m = 6678137.0, y_m = 0.0, z_m = 0.0 }
"#).unwrap();
    let receiver = ReceiverConfig {
        name: "equator".into(),
        trajectory: TrajectoryConfig::StationaryEcef { x_m: 6_378_137.0, y_m: 0.0, z_m: 0.0 },
        clock: NodeClockSpec::default(),
        bands: Vec::new(),
    };
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS).with_noise_power_w(0.0);
    let period_start = 1_717_260_840.0 + 5.0; // 1717260840 is a multiple of 20
    let iq = collect(&sky, &receiver, &band, period_start + 1.9, 6.0, &options(0, 10_000, 1));

    let delay_samples = 300_000.0 / lipika::SPEED_OF_LIGHT_M_S * FS; // 48.03
    // Sample n is keyed on iff its emission time n/fs − delay falls in a burst.
    let edge = |t_rel_s: f64| (t_rel_s * FS + delay_samples).ceil() as usize;
    let (on1, off1) = (edge(0.1), edge(0.6));
    let (on2, off2) = (edge(5.1), edge(5.35));
    assert_eq!((on1, off1, on2, off2), (4849, 28_849, 244_849, 256_849));
    let expected_amplitude =
        (lipika::SPEED_OF_LIGHT_M_S / HIBER_HZ) / (4.0 * std::f64::consts::PI * 300_000.0);
    for (n, v) in iq.iter().enumerate() {
        let on = (on1..off1).contains(&n) || (on2..off2).contains(&n);
        if on {
            let a = v.norm() as f64;
            assert!((a / expected_amplitude - 1.0).abs() < 1e-4, "sample {n}: amplitude {a}");
        } else {
            assert_eq!(*v, Complex32::new(0.0, 0.0), "sample {n} should be silent");
        }
    }

    // Without the gate the same transmitter is on throughout.
    let mut ungated = sky.clone();
    ungated.transmitters[0].gate = None;
    let iq = collect(&ungated, &receiver, &band, period_start + 1.9, 0.2, &options(0, 10_000, 1));
    assert!(iq.iter().all(|v| v.norm() > 0.0));

    // A burst that does not fit its period is a configuration error.
    let mut bad = sky.clone();
    bad.transmitters[0].gate.as_mut().unwrap().bursts[1].offset_s = 19.9;
    assert!(render_window(&bad, &receiver, &band, period_start, 0.01,
                          &WindowOptions::default(), |_| ControlFlow::Continue(())).is_err());
}

#[mantra_macros::req_test("fn-8")]
#[test]
fn example_sky_is_three_bursts_every_fifteen_seconds() {
    let sky = hiber_sky();
    assert_eq!(sky.transmitters.len(), 1);
    let tx = &sky.transmitters[0];
    assert_eq!(tx.carrier_hz, 401_150_000.0);
    assert!(matches!(tx.kind, SignalKind::Bpsk { bit_rate_bps, .. } if bit_rate_bps == 4800.0));
    assert!(matches!(tx.trajectory, TrajectoryConfig::Tle { .. }));
    let gate = tx.gate.as_ref().unwrap();
    assert_eq!((gate.period_s, gate.bursts.len()), (15.0, 3));

    // One whole period, opened 0.2 s early so the first burst's edge is in
    // view: exactly three keyed-on runs of 0.5 s (24 000 samples, give or
    // take one for the Doppler stretch), then silence to the end.
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS).with_noise_power_w(0.0);
    let iq = collect(&sky, &glasgow(), &band, PASS_START_UNIX - 0.2, 15.0, &options(0, 1 << 16, 1));
    let mut runs: Vec<(usize, usize)> = Vec::new(); // (first sample, length)
    for (n, v) in iq.iter().enumerate() {
        if v.norm() > 0.0 {
            match runs.last_mut() {
                Some((start, len)) if *start + *len == n => *len += 1,
                _ => runs.push((n, 1)),
            }
        }
    }
    assert_eq!(runs.len(), 3, "runs: {runs:?}");
    for (k, &(start, len)) in runs.iter().enumerate() {
        assert!((23_999..=24_001).contains(&len), "burst {k} is {len} samples");
        // Bursts start 1.5 s apart, 0.2 s + light time (2–8 ms) into the window.
        let expected = (0.2 + 1.5 * k as f64) * FS;
        let late = start as f64 - expected;
        assert!((0.002 * FS..0.008 * FS).contains(&late), "burst {k} starts {late} samples late");
    }
}

#[mantra_macros::req_test("fn-10", "fid-4")]
#[test]
fn empty_sky_is_receiver_noise_only() {
    // An empty document is a valid sky.
    let sky = SkyConfig::from_toml_str("").unwrap();
    assert!(sky.transmitters.is_empty());
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS).with_noise_temperature_k(290.0);
    let mut iq = Vec::new();
    let summary = render_window(&sky, &glasgow(), &band, PASS_START_UNIX, 2.0,
                                &options(5, 8192, 1), |c| {
        iq.extend_from_slice(c);
        ControlFlow::Continue(())
    }).unwrap();
    assert!(summary.contributing_transmitters.is_empty());
    assert_eq!(iq.len(), 96_000);

    // White, zero-mean, at the thermal floor k·T·B.
    let floor = 1.380_649e-23 * 290.0 * FS;
    let p = mean_power(&iq);
    assert!((p / floor - 1.0).abs() < 0.02, "noise power {p:.3e} W vs k·T·B {floor:.3e} W");
    let mean = iq.iter().fold(Complex32::new(0.0, 0.0), |a, v| a + v) / iq.len() as f32;
    assert!((mean.norm() as f64) < 0.02 * floor.sqrt(), "noise mean {mean}");
    let lag1: f64 = iq.windows(2).map(|w| (w[0] * w[1].conj()).re as f64).sum::<f64>()
        / (iq.len() - 1) as f64;
    assert!(lag1.abs() < 0.02 * floor, "lag-1 correlation {lag1:.3e}");

    // A sky whose only transmitter is far out of band is the same recording:
    // it is skipped, and the noise stream does not depend on the sky.
    let far = SkyConfig::from_toml_str(r#"
[[transmitter]]
name = "elsewhere"
carrier_hz = 1_575_420_000.0
kind = "cw"
baseband_offset_hz = 0.0
trajectory = { kind = "stationary_ecef", x_m = 7000000.0, y_m = 0.0, z_m = 0.0 }
"#).unwrap();
    assert!(bytes(&iq)
        == bytes(&collect(&far, &glasgow(), &band, PASS_START_UNIX, 2.0, &options(5, 8192, 1))));
    // With the noise switched off, an empty sky is silence.
    let quiet = band.clone().with_noise_power_w(0.0);
    let iq = collect(&sky, &glasgow(), &quiet, PASS_START_UNIX, 0.1, &options(5, 8192, 1));
    assert!(iq.iter().all(|v| *v == Complex32::new(0.0, 0.0)));
    // A zero-length window is valid and empty.
    assert!(collect(&sky, &glasgow(), &band, PASS_START_UNIX, 0.0, &options(5, 8192, 1)).is_empty());
}

#[test]
fn a_sky_file_holds_transmitters_only() {
    // Scene-only tables are rejected, so a scene file is not silently
    // accepted as a sky with its receivers and epoch dropped.
    assert!(SkyConfig::from_toml_str("[scene]\nepoch_unix = 0.0\nduration_s = 1.0\n").is_err());
    assert!(SkyConfig::from_toml_str(
        "[[receiver]]\nname = \"rx\"\ntrajectory = { kind = \"stationary_ecef\", x_m = 1.0, y_m = 0.0, z_m = 0.0 }\nband = []\n"
    ).is_err());
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/hiber_beacon_sky.toml");
    assert_eq!(SkyConfig::from_toml_file(&path).unwrap().transmitters[0].name, "hiber_beacon");
    assert!(SkyConfig::from_toml_file(Path::new("/nonexistent/sky.toml")).is_err());
    // Bad window arguments are errors, not panics.
    let band = ReceiverBandConfig::new("uhf", HIBER_HZ, FS);
    for (start, dur) in [(f64::NAN, 1.0), (0.0, -1.0), (0.0, f64::INFINITY)] {
        assert!(render_window(&SkyConfig::empty(), &glasgow(), &band, start, dur,
                              &WindowOptions::default(), |_| ControlFlow::Continue(())).is_err());
    }
}
