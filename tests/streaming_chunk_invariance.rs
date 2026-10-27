//! The load-bearing streaming guarantee: **chunking must not change the
//! result.** Generating a scene in one giant chunk vs. many tiny chunks must
//! produce byte-for-byte identical `.iq.c64` files, and the streamed bytes
//! must equal the in-memory `run_scene` output. This exercises every
//! cross-sample stateful path at once: the clock Wiener walk (shared across
//! two bands at different rates), the per-band noise stream, the absolute-time
//! geometry/signal models, and the fixed-order parallel transmitter reduction.

use std::path::{Path, PathBuf};

use lipika::config::SceneConfig;
use lipika::pipeline::{run_scene, run_scene_to_dir};

// One receiver, a non-trivial clock (allan + drift + offset), TWO bands at
// different sample rates off the same clock, and three transmitters — two
// landing in band A (exercises the multi-tx ordered reduction) and one in
// band B. Finite noise floors on both bands exercise the noise stream.
const SCENE: &str = r#"
[scene]
epoch_unix = 1717260840.0
duration_s = 0.05

[[transmitter]]
name = "a0"
carrier_hz = 1_000_000_000.0
tx_power_dbw = 30.0
kind = "cw"
baseband_offset_hz = 1234.0
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 0.0, z_m = 1000.0 }

[[transmitter]]
name = "a1"
carrier_hz = 1_000_050_000.0
tx_power_dbw = 28.0
kind = "cw"
baseband_offset_hz = -777.0
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 2000.0, z_m = 500.0 }

[[transmitter]]
name = "b0"
carrier_hz = 1_500_000_000.0
tx_power_dbw = 30.0
kind = "cw"
baseband_offset_hz = 333.0
trajectory = { kind = "stationary_ecef", x_m = 6378137.0, y_m = 1000.0, z_m = 1000.0 }

[[receiver]]
name = "rx"
trajectory = { kind = "stationary_geodetic", lat_deg = 55.86, lon_deg = -4.25, h_m = 10.0 }
clock = { allan_deviation_1s = 1e-9, frequency_offset_y0 = 0.5e-6, initial_time_offset_s = 1e-7 }
[[receiver.band]]
band_id = "a"
lo_carrier_hz = 1_000_000_000.0
sample_rate_hz = 1_000_000.0
noise_power_w = 1.0e-15
[[receiver.band]]
band_id = "b"
lo_carrier_hz = 1_500_000_000.0
sample_rate_hz = 200_000.0
noise_power_w = 1.0e-15
"#;

fn unique_dir(tag: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir()
        .join(format!("lipika_stream_{tag}_{}_{nanos}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn read_c64(dir: &Path, name: &str) -> Vec<u8> {
    std::fs::read(dir.join(format!("{name}.iq.c64"))).unwrap()
}

#[mantra_macros::req_test("fn-7", "nf-2")]
#[test]
fn chunked_output_is_byte_identical_regardless_of_chunk_size() {
    let scene: SceneConfig = toml::from_str(SCENE).unwrap();
    let seed = 4242;

    let dir_one = unique_dir("one");
    let dir_many = unique_dir("many");

    // One chunk covering the whole band vs. tiny 4096-sample chunks.
    let big = scene.scene.duration_s as usize * 4_000_000; // >> any band's n
    run_scene_to_dir(&scene, seed, &dir_one, big, /*threads=*/ 0).unwrap();
    run_scene_to_dir(&scene, seed, &dir_many, 4096, /*threads=*/ 0).unwrap();

    for name in ["rx__a", "rx__b"] {
        let a = read_c64(&dir_one, name);
        let b = read_c64(&dir_many, name);
        assert_eq!(a.len(), b.len(), "{name}: length differs between chunk sizes");
        assert!(a == b, "{name}: bytes differ between single-chunk and many-chunk runs");
        assert!(!a.is_empty(), "{name}: empty output");
    }

    let _ = std::fs::remove_dir_all(&dir_one);
    let _ = std::fs::remove_dir_all(&dir_many);
}

#[mantra_macros::req_test("fn-7", "eq-e10")]
#[test]
fn streamed_bytes_equal_in_memory_run_scene() {
    let scene: SceneConfig = toml::from_str(SCENE).unwrap();
    let seed = 4242;

    let dir = unique_dir("mem");
    run_scene_to_dir(&scene, seed, &dir, 4096, 0).unwrap();
    let outputs = run_scene(&scene, seed).unwrap();

    for out in &outputs {
        // Encode the in-memory IQ the same way the writer does (LE f32 re/im).
        let mut mem_bytes = Vec::with_capacity(out.iq.len() * 8);
        for v in out.iq.iter() {
            mem_bytes.extend_from_slice(&v.re.to_le_bytes());
            mem_bytes.extend_from_slice(&v.im.to_le_bytes());
        }
        let file_bytes = read_c64(&dir, &out.name);
        assert_eq!(mem_bytes.len(), file_bytes.len(), "{}: length mismatch", out.name);
        assert!(mem_bytes == file_bytes,
            "{}: in-memory run_scene bytes differ from streamed file", out.name);
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[mantra_macros::req_test("fn-7", "nf-3")]
#[test]
fn thread_count_does_not_change_output() {
    // Fixed-order transmitter reduction ⇒ result independent of worker count.
    let scene: SceneConfig = toml::from_str(SCENE).unwrap();
    let seed = 7;

    let dir1 = unique_dir("t1");
    let dir4 = unique_dir("t4");
    run_scene_to_dir(&scene, seed, &dir1, 8192, /*threads=*/ 1).unwrap();
    run_scene_to_dir(&scene, seed, &dir4, 8192, /*threads=*/ 4).unwrap();

    for name in ["rx__a", "rx__b"] {
        assert!(read_c64(&dir1, name) == read_c64(&dir4, name),
            "{name}: output depends on thread count");
    }

    let _ = std::fs::remove_dir_all(&dir1);
    let _ = std::fs::remove_dir_all(&dir4);
}
