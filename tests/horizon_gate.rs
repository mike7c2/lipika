//! The horizon gate (`channel::horizon`, nf-5) must never change a sample.
//!
//! Every test renders the same thing twice, through the gate and through the
//! reference path that evaluates every in-band transmitter for every sample
//! (`render_window_without_horizon_gate`, `run_scene_without_horizon_gate`),
//! and compares the bytes:
//!
//! - a sky of LEOs that rise, set, graze the horizon from either side, stay
//!   up and stay down during the window, across chunk sizes and thread counts;
//! - other cutoffs, and no cutoff at all (nothing may be skipped);
//! - the scene path, in memory and streamed, with a noiseless band so the
//!   skipped contributions are compared as exact zeros;
//! - a receiver that is itself in orbit;
//! - an SGP4 (TLE) transmitter rendered across its rise;
//! - Earth-fixed transmitters.

use std::ops::ControlFlow;

use lipika::clock::NodeClockSpec;
use lipika::config::{
    AdcConfig, AntennaConfig, BurstConfig, FrontendConfig, GateConfig, ReceiverBandConfig,
    ReceiverConfig, SceneConfig, SceneMeta, SignalKind, SkyConfig, TrajectoryConfig,
    TransmitterConfig,
};
use lipika::geometry::ecef::geodetic_to_ecef;
use lipika::geometry::kepler::CircularKeplerian;
use lipika::geometry::tle::TleTrajectory;
use lipika::geometry::Trajectory;
use lipika::pipeline::{
    render_window, render_window_without_horizon_gate, run_scene, run_scene_to_dir,
    run_scene_without_horizon_gate, WindowOptions, WindowSummary,
};
use num_complex::Complex32;

const CARRIER_HZ: f64 = 401_150_000.0;
const FS: f64 = 1_000.0;
const START_UNIX: f64 = 1_717_260_690.0;
const DURATION_S: f64 = 120.0;
const GLASGOW: (f64, f64, f64) = (55.86, -4.25, 20.0);

fn glasgow() -> ReceiverConfig {
    ReceiverConfig::stationary_geodetic("glasgow", GLASGOW.0, GLASGOW.1, GLASGOW.2)
}

/// A band at the test carrier with an isotropic antenna and this cutoff
/// (`None` ⇒ an antenna that occludes nothing).
fn band(cutoff_deg: Option<f64>) -> ReceiverBandConfig {
    ReceiverBandConfig::new("uhf", CARRIER_HZ, FS)
        .with_noise_temperature_k(290.0)
        .with_frontend(FrontendConfig {
            antenna: Some(AntennaConfig { peak_dbi: 2.0, horizon_dbi: -1.0, fb_db: 6.0, cutoff_deg }),
            gain_db: None,
            nonlinearity: None,
            adc: None,
        })
}

fn bytes(iq: &[Complex32]) -> Vec<u8> {
    iq.iter().flat_map(|v| [v.re.to_le_bytes(), v.im.to_le_bytes()]).flatten().collect()
}

/// Geometric elevation (degrees) of `p` seen from `o`, against `o`'s
/// geocentric up: the angle the antenna cutoff is applied to.
fn elevation_deg(o: [f64; 3], p: [f64; 3]) -> f64 {
    let d = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
    let n = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    ((d[0] * o[0] + d[1] * o[1] + d[2] * o[2]) / (n(d) * n(o))).asin().to_degrees()
}

/// What a satellite does over Glasgow during the test window.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Pass {
    /// Below −2° at the start, above +2° at the end.
    Rises,
    /// Above +2° at the start, below −2° at the end.
    Sets,
    /// Never below 5°.
    StaysUp,
    /// Peaks between 0.05° and 1° above the horizon.
    GrazesAbove,
    /// Peaks between 0.05° and 1° below the horizon: never visible.
    GrazesBelow,
    /// Never within 10° of the horizon.
    StaysDown,
}

fn classify(orbit: &CircularKeplerian) -> Option<Pass> {
    let o = geodetic_to_ecef(GLASGOW.0, GLASGOW.1, GLASGOW.2);
    let el = |t: f64| elevation_deg(o, orbit.position_ecef(START_UNIX + t));
    let (first, last) = (el(0.0), el(DURATION_S));
    if first < -25.0 && last < -25.0 {
        return (first < -40.0).then_some(Pass::StaysDown);
    }
    let track: Vec<f64> = (0..=(DURATION_S as usize * 4)).map(|k| el(k as f64 * 0.25)).collect();
    let max = track.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = track.iter().cloned().fold(f64::INFINITY, f64::min);
    if first < -2.0 && last > 2.0 { Some(Pass::Rises) }
    else if first > 2.0 && last < -2.0 { Some(Pass::Sets) }
    else if min > 5.0 { Some(Pass::StaysUp) }
    else if max > 0.05 && max < 1.0 { Some(Pass::GrazesAbove) }
    else if max < -0.05 && max > -1.0 { Some(Pass::GrazesBelow) }
    else { None }
}

/// A sky of circular LEOs chosen by what they do during the window: three of
/// each kind of [`Pass`] (four that stay down). Every third one is a
/// burst-gated BPSK beacon, the rest CW tones at distinct offsets.
fn leo_sky() -> (SkyConfig, Vec<(String, Pass)>) {
    let want = |p: Pass| if p == Pass::StaysDown { 4 } else { 3 };
    let mut have: std::collections::HashMap<Pass, usize> = Default::default();
    let mut transmitters = Vec::new();
    let mut passes = Vec::new();
    let mut k = 0usize;
    while transmitters.len() < 19 {
        k += 1;
        assert!(k < 2_000_000, "ran out of candidate orbits; found {have:?}");
        let (a, inc, raan, nu) = (
            6_878_137.0 + 30_000.0 * (k % 5) as f64,
            [53.0, 65.0, 87.9, 97.5][k % 4],
            (k as f64 * 137.508) % 360.0,
            (k as f64 * 61.803) % 360.0,
        );
        let orbit = CircularKeplerian {
            semi_major_axis_m: a,
            inclination_deg: inc,
            raan_deg: raan,
            true_anomaly_deg_at_epoch: nu,
            epoch_unix: START_UNIX,
        };
        let Some(pass) = classify(&orbit) else { continue };
        let n = have.entry(pass).or_insert(0);
        if *n >= want(pass) { continue; }
        *n += 1;
        let i = transmitters.len();
        let name = format!("leo-{i:02}");
        let (kind, gate) = if i % 3 == 2 {
            (
                SignalKind::Bpsk {
                    bit_rate_bps: 200.0,
                    payload_bits: 268,
                    payload_seed: 40 + i as u64,
                    sync_word: vec![-1, 1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1,
                                    -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1],
                },
                Some(GateConfig {
                    period_s: 5.0,
                    epoch_unix: 0.0,
                    bursts: vec![
                        BurstConfig { offset_s: 0.3 * (i % 5) as f64, duration_s: 1.5 },
                        BurstConfig { offset_s: 3.0, duration_s: 0.7 },
                    ],
                }),
            )
        } else {
            (SignalKind::Cw { baseband_offset_hz: -400.0 + 43.0 * i as f64 }, None)
        };
        transmitters.push(TransmitterConfig {
            name: name.clone(),
            carrier_hz: CARRIER_HZ,
            tx_power_dbw: 0.0,
            bandwidth_hz: None,
            trajectory: TrajectoryConfig::CircularKeplerian {
                semi_major_axis_m: a,
                inclination_deg: inc,
                raan_deg: raan,
                true_anomaly_deg_at_epoch: nu,
                epoch_unix: START_UNIX,
            },
            signal_epoch_unix: Some(0.0),
            gate,
            kind,
        });
        passes.push((name, pass));
    }
    (SkyConfig::new(transmitters), passes)
}

fn collect(
    gate: bool,
    sky: &SkyConfig,
    receiver: &ReceiverConfig,
    band: &ReceiverBandConfig,
    start_unix: f64,
    duration_s: f64,
    opts: &WindowOptions,
) -> (Vec<Complex32>, WindowSummary) {
    let mut iq = Vec::new();
    let sink = |chunk: &[Complex32]| {
        iq.extend_from_slice(chunk);
        ControlFlow::Continue(())
    };
    let summary = if gate {
        render_window(sky, receiver, band, start_unix, duration_s, opts, sink)
    } else {
        render_window_without_horizon_gate(sky, receiver, band, start_unix, duration_s, opts, sink)
    }
    .unwrap();
    assert_eq!(summary.n_samples_emitted, iq.len());
    (iq, summary)
}

fn options(chunk_samples: usize, n_threads: usize) -> WindowOptions {
    WindowOptions { seed: 11, chunk_samples, n_threads }
}

#[mantra_macros::req_test("nf-5", "fn-7", "nf-3", "fn-10")]
#[test]
fn gate_does_not_change_a_sky_of_rising_and_setting_satellites() {
    let (sky, passes) = leo_sky();
    // A receiver clock that is off and drifting: the UTC sample instants the
    // gate bounds are not the nominal ones.
    let receiver = glasgow().with_clock(NodeClockSpec {
        allan_deviation_1s: 1e-9,
        flicker_floor: 1e-10,
        frequency_offset_y0: 2.0e-6,
        initial_time_offset_s: 3.0e-3,
    });
    let band = band(Some(0.0));
    let n = (FS * DURATION_S) as usize;

    let (reference, ref_summary) =
        collect(false, &sky, &receiver, &band, START_UNIX, DURATION_S, &options(1 << 20, 1));
    assert_eq!(reference.len(), n);
    let reference = bytes(&reference);

    // Chunks shorter than a gate stretch (256 samples here), not a multiple
    // of it, longer, and the whole window; one, several and all threads. The
    // satellites cross the horizon wherever they do, so the crossings land
    // mid-chunk and mid-stretch.
    for (chunk, threads) in
        [(1usize << 20, 1usize), (100, 1), (257, 3), (1000, 1), (4096, 4), (65_536, 0), (n, 2)]
    {
        let (iq, summary) =
            collect(true, &sky, &receiver, &band, START_UNIX, DURATION_S, &options(chunk, threads));
        assert!(bytes(&iq) == reference,
            "the horizon gate changed the output at chunk_samples = {chunk}, n_threads = {threads}");
        assert_eq!(summary, ref_summary, "summary differs at chunk {chunk}, threads {threads}");
    }

    // The summary: every transmitter is in band; the visible ones are exactly
    // those whose pass takes them above the horizon.
    let names = |keep: &dyn Fn(Pass) -> bool| -> Vec<String> {
        passes.iter().filter(|(_, p)| keep(*p)).map(|(n, _)| n.clone()).collect()
    };
    assert_eq!(ref_summary.contributing_transmitters, names(&|_| true));
    assert_eq!(
        ref_summary.visible_transmitters,
        names(&|p| !matches!(p, Pass::StaysDown | Pass::GrazesBelow)),
        "passes: {passes:?}",
    );
    assert_eq!(ref_summary.visible_transmitters.len(), 12);
}

#[mantra_macros::req_test("nf-5")]
#[test]
fn gate_respects_the_configured_cutoff_and_its_absence() {
    let (sky, _) = leo_sky();
    let receiver = glasgow();
    let mut visible_counts = Vec::new();
    for cutoff_deg in [Some(-3.0), Some(0.0), Some(10.0), None] {
        let band = band(cutoff_deg);
        let (reference, ref_summary) =
            collect(false, &sky, &receiver, &band, START_UNIX, DURATION_S, &options(1 << 20, 1));
        let (gated, summary) =
            collect(true, &sky, &receiver, &band, START_UNIX, DURATION_S, &options(777, 2));
        assert!(bytes(&gated) == bytes(&reference), "cutoff {cutoff_deg:?}: output changed");
        assert_eq!(summary, ref_summary, "cutoff {cutoff_deg:?}");
        visible_counts.push(summary.visible_transmitters.len());
        if cutoff_deg.is_none() {
            // No cutoff: nothing is occluded (the back-lobe hears through the
            // Earth), nothing may be skipped.
            assert_eq!(summary.visible_transmitters, summary.contributing_transmitters);
        }
    }
    // A higher cutoff hides more: −3° sees the grazers that 0° does not, 10°
    // loses some that 0° sees, no cutoff sees all 19.
    assert!(visible_counts[0] > visible_counts[1] && visible_counts[1] > visible_counts[2],
        "visible per cutoff: {visible_counts:?}");
    assert_eq!(visible_counts[3], 19);

    // No antenna at all: likewise nothing is occluded or skipped.
    let bare = ReceiverBandConfig::new("uhf", CARRIER_HZ, FS);
    let (reference, _) = collect(false, &sky, &receiver, &bare, START_UNIX, 20.0, &options(1 << 20, 1));
    let (gated, summary) = collect(true, &sky, &receiver, &bare, START_UNIX, 20.0, &options(999, 2));
    assert!(bytes(&gated) == bytes(&reference));
    assert_eq!(summary.visible_transmitters.len(), 19);
}

#[mantra_macros::req_test("nf-5", "fn-7")]
#[test]
fn gate_does_not_change_a_scene_in_memory_or_streamed() {
    let (sky, _) = leo_sky();
    let duration_s = 60.0;
    // Two receivers. Glasgow has a noiseless band behind an ADC (so a skipped
    // contribution is compared as an exact zero, sign included) and a second
    // band without a cutoff; the other site has a 5° cutoff and thermal noise.
    let noiseless = ReceiverBandConfig {
        frontend: Some(FrontendConfig {
            adc: Some(AdcConfig { bits: 12, full_scale_v: 1.0e-6 }),
            ..band(Some(0.0)).frontend.unwrap()
        }),
        ..band(Some(0.0)).with_noise_power_w(0.0)
    };
    let open = ReceiverBandConfig { band_id: "open".into(), ..band(None) };
    let scene = SceneConfig {
        scene: SceneMeta { epoch_unix: START_UNIX, duration_s, gps_utc_offset_s: 18.0 },
        transmitters: sky.transmitters.clone(),
        receivers: vec![
            ReceiverConfig { bands: vec![noiseless, open], ..glasgow() },
            ReceiverConfig {
                bands: vec![band(Some(5.0))],
                ..ReceiverConfig::stationary_geodetic("svalbard", 78.23, 15.39, 450.0)
            },
        ],
    };

    let reference = run_scene_without_horizon_gate(&scene, 5).unwrap();
    let gated = run_scene(&scene, 5).unwrap();
    assert_eq!(reference.len(), 3);
    for (r, g) in reference.iter().zip(&gated) {
        assert_eq!(r.name, g.name);
        assert!(bytes(g.iq.as_slice().unwrap()) == bytes(r.iq.as_slice().unwrap()),
            "the horizon gate changed band {}", r.name);
        assert_eq!(r.descriptor, g.descriptor, "sidecar descriptor of {}", r.name);
    }
    // The noiseless band really is exercised: it holds signal.
    assert!(gated[0].iq.iter().any(|v| v.norm_sqr() > 0.0));

    // Streamed to disk in small chunks on several threads: the same bytes.
    let dir = std::env::temp_dir().join(format!("lipika_horizon_gate_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let manifests = run_scene_to_dir(&scene, 5, &dir, 1234, 3).unwrap();
    assert_eq!(manifests.len(), 3);
    for (m, r) in manifests.iter().zip(&reference) {
        assert_eq!(m.name, r.name);
        assert!(std::fs::read(&m.path).unwrap() == bytes(r.iq.as_slice().unwrap()),
            "streamed band {} differs from the ungated reference", r.name);
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[mantra_macros::req_test("nf-5")]
#[test]
fn gate_does_not_change_what_an_orbiting_receiver_records() {
    // The receiver is a satellite 400 km up, below the transmitters; its
    // antenna ignores everything more than 15° under its local horizontal
    // (roughly: behind the Earth's limb). Both ends move.
    let (sky, _) = leo_sky();
    let receiver = ReceiverConfig {
        name: "orbiter".into(),
        trajectory: TrajectoryConfig::CircularKeplerian {
            semi_major_axis_m: 6_778_137.0,
            inclination_deg: 70.0,
            raan_deg: 320.0,
            true_anomaly_deg_at_epoch: 50.0,
            epoch_unix: START_UNIX,
        },
        clock: NodeClockSpec::default(),
        bands: Vec::new(),
    };
    let band = band(Some(-15.0));
    let (reference, ref_summary) =
        collect(false, &sky, &receiver, &band, START_UNIX, DURATION_S, &options(1 << 20, 1));
    for (chunk, threads) in [(1usize << 20, 1usize), (300, 2), (5000, 0)] {
        let (gated, summary) =
            collect(true, &sky, &receiver, &band, START_UNIX, DURATION_S, &options(chunk, threads));
        assert!(bytes(&gated) == bytes(&reference),
            "orbiting receiver: output changed at chunk {chunk}, threads {threads}");
        assert_eq!(summary, ref_summary);
    }
    // The scenario has both: transmitters it sees and transmitters it does not.
    let visible = ref_summary.visible_transmitters.len();
    assert!(visible > 0 && visible < 19, "{visible} of 19 visible");
}

#[mantra_macros::req_test("nf-5", "fn-10")]
#[test]
fn gate_does_not_change_a_tle_transmitter_across_its_rise() {
    // The example sky's SGP4 orbit, carrying a CW tone and the example's
    // burst-gated beacon. Find the instant it rises over Glasgow before the
    // pass that culminates at Unix 1717260990, and render 6 s around it at
    // 48 kHz with no noise: silence, then signal.
    let (l1, l2) = (
        "1 25544U 98067A   24001.50000000  .00010000  00000+0  18000-3 0  9994",
        "2 25544  51.6400 100.0000 0001000  90.0000 270.0000 15.50000000000009",
    );
    let orbit = TleTrajectory::from_tle(l1, l2).unwrap();
    let o = geodetic_to_ecef(GLASGOW.0, GLASGOW.1, GLASGOW.2);
    let el = |t: f64| elevation_deg(o, orbit.position_ecef(t));
    let culmination = 1_717_260_990.0;
    let mut rise = culmination;
    while el(rise) > 0.0 {
        rise -= 1.0;
        assert!(culmination - rise < 900.0, "no rise found");
    }
    let start = rise - 3.0; // whole seconds: `rise` is within 1 s after the crossing
    let tle = || TrajectoryConfig::Tle { line1: l1.into(), line2: l2.into() };
    let example: SkyConfig = SkyConfig::from_toml_file(
        &std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/hiber_beacon_sky.toml"))
        .unwrap();
    let mut transmitters = example.transmitters.clone();
    transmitters.push(TransmitterConfig {
        name: "tone".into(),
        carrier_hz: CARRIER_HZ,
        tx_power_dbw: 0.0,
        bandwidth_hz: None,
        trajectory: tle(),
        signal_epoch_unix: None,
        gate: None,
        kind: SignalKind::Cw { baseband_offset_hz: 1_000.0 },
    });
    let sky = SkyConfig::new(transmitters);
    let fs = 48_000.0;
    let band = ReceiverBandConfig::new("uhf", CARRIER_HZ, fs)
        .with_noise_power_w(0.0)
        .with_frontend(FrontendConfig {
            antenna: Some(AntennaConfig { peak_dbi: 0.0, horizon_dbi: 0.0, fb_db: 0.0, cutoff_deg: Some(0.0) }),
            gain_db: None,
            nonlinearity: None,
            adc: None,
        });

    let (reference, ref_summary) =
        collect(false, &sky, &glasgow(), &band, start, 6.0, &options(1 << 20, 1));
    assert_eq!(reference.len(), 288_000);
    // Silent until the satellite clears the horizon, then never silent again
    // (the tone is not burst-gated); the rise is inside the window.
    let first_signal = reference.iter().position(|v| v.norm_sqr() > 0.0).expect("it rises");
    assert!(first_signal > 48_000 && first_signal < 240_000, "rise at sample {first_signal}");
    assert!(reference[first_signal..].iter().all(|v| v.norm_sqr() > 0.0));
    assert_eq!(ref_summary.visible_transmitters, ["hiber_beacon", "tone"]);

    // Default chunks (the rise falls inside one), chunks of exactly one gate
    // stretch, and odd ones.
    for (chunk, threads) in [(1usize << 16, 1usize), (4800, 2), (4799, 1), (100_000, 0)] {
        let (gated, summary) =
            collect(true, &sky, &glasgow(), &band, start, 6.0, &options(chunk, threads));
        assert!(bytes(&gated) == bytes(&reference),
            "TLE rise: output changed at chunk {chunk}, threads {threads}");
        assert_eq!(summary, ref_summary);
    }

    // A window that ends before the rise: in band, never visible, all zeros.
    let (before, summary) = collect(true, &sky, &glasgow(), &band, start - 30.0, 2.0, &options(1 << 16, 1));
    assert!(before.iter().all(|v| v.re.to_bits() == 0 && v.im.to_bits() == 0));
    assert_eq!(summary.contributing_transmitters, ["hiber_beacon", "tone"]);
    assert!(summary.visible_transmitters.is_empty());
}

#[mantra_macros::req_test("nf-5")]
#[test]
fn gate_handles_earth_fixed_transmitters() {
    // A mast 30 km from the receiver (above its horizon), a balloon overhead,
    // a site on the far side of the Earth, and one just beyond the horizon.
    let fixed = |name: &str, lat: f64, lon: f64, h: f64, offset_hz: f64| TransmitterConfig {
        name: name.into(),
        carrier_hz: CARRIER_HZ,
        tx_power_dbw: -60.0,
        bandwidth_hz: None,
        trajectory: TrajectoryConfig::StationaryGeodetic { lat_deg: lat, lon_deg: lon, h_m: h },
        signal_epoch_unix: Some(0.0),
        gate: None,
        kind: SignalKind::Cw { baseband_offset_hz: offset_hz },
    };
    let sky = SkyConfig::new(vec![
        fixed("mast", 55.86, -3.77, 3_000.0, 100.0),
        fixed("balloon", 55.86, -4.25, 30_000.0, -200.0),
        fixed("antipode", -55.86, 175.75, 100.0, 300.0),
        fixed("over-the-horizon", 50.0, -4.25, 100.0, -350.0),
    ]);
    let band = band(Some(0.0));
    let (reference, ref_summary) =
        collect(false, &sky, &glasgow(), &band, START_UNIX, 5.0, &options(1 << 20, 1));
    let (gated, summary) = collect(true, &sky, &glasgow(), &band, START_UNIX, 5.0, &options(700, 2));
    assert!(bytes(&gated) == bytes(&reference));
    assert_eq!(summary, ref_summary);
    assert_eq!(summary.contributing_transmitters.len(), 4);
    assert_eq!(summary.visible_transmitters, ["mast", "balloon"]);
}
