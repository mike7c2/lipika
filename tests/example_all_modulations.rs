//! `examples/all_modulations.toml` must keep its promise: one transmitter of
//! every signal kind, each landing in the band and radiating.

use std::path::Path;

use lipika::config::SceneConfig;
use lipika::pipeline::run_scene;

#[test]
fn every_signal_kind_is_present_and_in_band() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/all_modulations.toml");
    let mut scene = SceneConfig::from_toml_file(&path).unwrap();

    let mut kinds: Vec<&str> = scene.transmitters.iter().map(|t| t.kind.kind_name()).collect();
    kinds.sort_unstable();
    assert_eq!(kinds, ["bpsk", "ca_code", "cw", "ook", "qam16"]);

    // Run each transmitter alone for 10 ms against a negligible noise floor:
    // if it were gated out of the band (or silent) the output would be empty.
    scene.scene.duration_s = 0.01;
    scene.receivers[0].bands[0].noise_power_w = Some(1.0e-30);
    for tx in &scene.transmitters {
        let mut solo = scene.clone();
        solo.transmitters.retain(|t| t.name == tx.name);
        let iq = run_scene(&solo, /*seed=*/ 1).unwrap().remove(0).iq;
        let mean_power: f64 =
            iq.iter().map(|v| (v.re as f64).powi(2) + (v.im as f64).powi(2)).sum::<f64>()
                / iq.len() as f64;
        assert!(mean_power > 1e-15, "{} contributes {mean_power:.3e} W to the band", tx.name);
    }
}
