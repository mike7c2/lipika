//! Continuous-wave (single-tone) source. Useful for sanity checks.

use ndarray::Array1;
use num_complex::Complex32;

use super::SignalSource;

pub struct ContinuousWave {
    pub carrier_hz: f64,
    /// Baseband-frequency offset relative to the nominal carrier (Hz).
    /// Zero by default — a clean carrier at `carrier_hz` produces a pure
    /// DC tone in the source's baseband frame.
    pub baseband_offset_hz: f64,
}

impl SignalSource for ContinuousWave {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        let omega = std::f64::consts::TAU * self.baseband_offset_hz;
        let mut out = Array1::<Complex32>::zeros(t_emit.len());
        for (i, &t) in t_emit.iter().enumerate() {
            let phase = (omega * t) as f32;
            out[i] = Complex32::new(phase.cos(), phase.sin());
        }
        out
    }
}
