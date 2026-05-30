//! IQ file writers.
//!
//! v1 supports raw NumPy `.npy` complex64 dumps + a `.json` sidecar with
//! per-observer metadata (epoch, sample rate, etc.). SigMF support is a
//! later addition once the schema settles.

use std::path::Path;

use ndarray::Array1;
use num_complex::Complex32;

pub mod npy;
pub mod sigmf;

/// One observer's output: complex-baseband IQ plus the metadata needed
/// to interpret it.
pub struct ObserverOutput {
    pub name: String,
    pub epoch_unix: f64,
    pub sample_rate_hz: f64,
    pub lo_carrier_hz: f64,
    pub iq: Array1<Complex32>,
}

impl ObserverOutput {
    pub fn write_npy_with_sidecar(&self, dir: &Path) -> Result<(), anyhow::Error> {
        std::fs::create_dir_all(dir)?;
        let iq_path = dir.join(format!("{}.iq.npy", self.name));
        let meta_path = dir.join(format!("{}.iq.json", self.name));
        npy::write_complex64_npy(&iq_path, &self.iq)?;
        let meta = serde_json::json!({
            "name": self.name,
            "epoch_unix": self.epoch_unix,
            "sample_rate_hz": self.sample_rate_hz,
            "lo_carrier_hz": self.lo_carrier_hz,
            "n_samples": self.iq.len(),
            "dtype": "complex64",
        });
        std::fs::write(&meta_path, serde_json::to_string_pretty(&meta)?)?;
        Ok(())
    }
}
