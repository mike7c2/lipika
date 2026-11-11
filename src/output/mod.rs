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

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

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
    /// Extended, "greedy" descriptor: which receiver/band this is, the list
    /// of transmitters (with PRNs) that actually contributed to this band,
    /// and the full scene config used to generate it. Written under the
    /// `descriptor` key of the JSON sidecar; the core fields above stay at
    /// the top level so simple readers can ignore it.
    pub descriptor: serde_json::Value,
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
        let meta = sidecar_json(
            &self.name, self.epoch_unix, self.sample_rate_hz, self.lo_carrier_hz,
            self.iq.len(), &self.descriptor,
        );
        std::fs::write(&meta_path, serde_json::to_string_pretty(&meta)?)?;
        Ok(())
    }
}

/// Build the JSON sidecar value. Shared by the in-memory and streaming writers
/// so the on-disk metadata is identical regardless of how the IQ was produced.
#[mantra_macros::req_satisfied("fn-6")]
fn sidecar_json(
    name: &str,
    epoch_unix: f64,
    sample_rate_hz: f64,
    lo_carrier_hz: f64,
    n_samples: usize,
    descriptor: &serde_json::Value,
) -> serde_json::Value {
    serde_json::json!({
        "name": name,
        "epoch_unix": epoch_unix,
        "sample_rate_hz": sample_rate_hz,
        "lo_carrier_hz": lo_carrier_hz,
        "n_samples": n_samples,
        "format": "raw complex64 little-endian (interleaved f32 re, f32 im)",
        "numpy_load": "np.fromfile(path, dtype=np.complex64)",
        "descriptor": descriptor,
    })
}

/// Encode a complex64 slice as raw little-endian bytes into `writer`.
#[mantra_macros::req_satisfied("fn-6", "nf-4")]
fn write_complex64_into<W: Write>(writer: &mut W, iq: &Array1<Complex32>) -> Result<(), anyhow::Error> {
    let mut buf = [0u8; 8];
    for v in iq.iter() {
        buf[..4].copy_from_slice(&v.re.to_le_bytes());
        buf[4..].copy_from_slice(&v.im.to_le_bytes());
        writer.write_all(&buf)?;
    }
    Ok(())
}

/// Write a 1-D complex64 array to disk as raw little-endian bytes.
pub fn write_raw_complex64(path: &Path, iq: &Array1<Complex32>) -> Result<(), anyhow::Error> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_complex64_into(&mut writer, iq)?;
    writer.flush()?;
    Ok(())
}

/// Incremental writer for streaming generation: opens `<dir>/<name>.iq.c64`
/// (truncating any prior file), appends complex64 chunks as they are produced,
/// then `finalize`s by flushing and writing the `<name>.iq.json` sidecar. Raw
/// c64 is headerless, so there is no header to patch as samples accrue.
pub struct StreamingC64Writer {
    writer: BufWriter<File>,
    dir: PathBuf,
    name: String,
    samples_written: usize,
}

impl StreamingC64Writer {
    pub fn create(dir: &Path, name: &str) -> Result<Self, anyhow::Error> {
        std::fs::create_dir_all(dir)?;
        let iq_path = dir.join(format!("{name}.iq.c64"));
        Ok(Self {
            writer: BufWriter::new(File::create(&iq_path)?),
            dir: dir.to_path_buf(),
            name: name.to_string(),
            samples_written: 0,
        })
    }

    /// Append one chunk of IQ to the raw file.
    #[mantra_macros::req_satisfied("nf-2")]
    pub fn push_chunk(&mut self, chunk: &Array1<Complex32>) -> Result<(), anyhow::Error> {
        write_complex64_into(&mut self.writer, chunk)?;
        self.samples_written += chunk.len();
        Ok(())
    }

    pub fn samples_written(&self) -> usize { self.samples_written }

    /// Flush the IQ file and write the JSON sidecar (now that the total sample
    /// count is known).
    pub fn finalize(
        mut self,
        epoch_unix: f64,
        sample_rate_hz: f64,
        lo_carrier_hz: f64,
        descriptor: &serde_json::Value,
    ) -> Result<(), anyhow::Error> {
        self.writer.flush()?;
        let meta = sidecar_json(
            &self.name, epoch_unix, sample_rate_hz, lo_carrier_hz,
            self.samples_written, descriptor,
        );
        let meta_path = self.dir.join(format!("{}.iq.json", self.name));
        std::fs::write(&meta_path, serde_json::to_string_pretty(&meta)?)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[mantra_macros::req_test("fn-6", "nf-4")]
    #[test]
    fn raw_complex64_roundtrip_via_byte_check() {
        let path = std::env::temp_dir().join("lipika_raw_test.c64");
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
