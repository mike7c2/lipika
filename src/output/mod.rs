//! IQ file writers.
//!
//! Output is raw little-endian complex64 (two f32 each: real, imag
//! interleaved) plus a `.json` sidecar with the metadata a reader needs
//! to interpret it (epoch, sample rate, LO frequency, sample count).
//!
//! Why raw rather than `.npy`: NumPy reads it just as easily via
//! `np.fromfile(path, dtype=np.complex64)` — no header parsing — and
//! every other SDR tool (GNU Radio file_source, gnss-sdr, inspectrum,
//! etc.) reads raw complex64 natively.

use std::path::Path;

use ndarray::Array1;
use num_complex::Complex32;

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
    /// Write `<dir>/<name>.iq.c64` (raw complex64 little-endian) plus
    /// `<dir>/<name>.iq.json` (metadata sidecar). Creates the directory
    /// if it doesn't exist.
    pub fn write_raw_with_sidecar(&self, dir: &Path) -> Result<(), anyhow::Error> {
        std::fs::create_dir_all(dir)?;
        let iq_path = dir.join(format!("{}.iq.c64", self.name));
        let meta_path = dir.join(format!("{}.iq.json", self.name));
        write_raw_complex64(&iq_path, &self.iq)?;
        let meta = serde_json::json!({
            "name": self.name,
            "epoch_unix": self.epoch_unix,
            "sample_rate_hz": self.sample_rate_hz,
            "lo_carrier_hz": self.lo_carrier_hz,
            "n_samples": self.iq.len(),
            "format": "raw complex64 little-endian (interleaved f32 re, f32 im)",
            "numpy_load": "np.fromfile(path, dtype=np.complex64)",
        });
        std::fs::write(&meta_path, serde_json::to_string_pretty(&meta)?)?;
        Ok(())
    }
}

/// Write a 1-D complex64 array to disk as raw little-endian bytes.
pub fn write_raw_complex64(path: &Path, iq: &Array1<Complex32>) -> Result<(), anyhow::Error> {
    use std::fs::File;
    use std::io::{BufWriter, Write};
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let mut buf = [0u8; 8];
    for v in iq.iter() {
        buf[..4].copy_from_slice(&v.re.to_le_bytes());
        buf[4..].copy_from_slice(&v.im.to_le_bytes());
        writer.write_all(&buf)?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_complex64_roundtrip_via_byte_check() {
        let path = std::env::temp_dir().join("dgsn_iq_gen_raw_test.c64");
        let iq = Array1::from(vec![
            Complex32::new(1.0, 2.0),
            Complex32::new(3.0, -4.0),
            Complex32::new(0.5, 0.25),
        ]);
        write_raw_complex64(&path, &iq).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(bytes.len(), 3 * 8);
        // Confirm first sample's real part decodes to 1.0
        let re_bytes: [u8; 4] = bytes[..4].try_into().unwrap();
        assert_eq!(f32::from_le_bytes(re_bytes), 1.0);
        let im_bytes: [u8; 4] = bytes[4..8].try_into().unwrap();
        assert_eq!(f32::from_le_bytes(im_bytes), 2.0);
        let _ = std::fs::remove_file(&path);
    }
}
