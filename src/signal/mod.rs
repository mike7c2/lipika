//! Signal sources — anything that can be evaluated at an emission time
//! `t_emit` to produce a complex baseband sample.
//!
//! The contract is **time-domain**, not frequency-domain. A source knows
//! its nominal carrier frequency (used downstream for the mix to baseband
//! at the observer) and produces baseband IQ at whatever emission times
//! the channel asks for. This makes Doppler emerge naturally from delay-
//! varying queries rather than being a separate frequency-shift step.

use ndarray::Array1;
use num_complex::Complex32;

pub mod bpsk;
pub mod cw;
pub mod ca_code;

pub trait SignalSource: Send + Sync {
    /// Nominal centre frequency of the source's RF emission (Hz).
    /// The observer mixes against its own LO at (this − offset) before
    /// recording IQ, so the Doppler emerges from delay-varying queries.
    fn carrier_hz(&self) -> f64;

    /// Evaluate the source's complex baseband signal at the requested
    /// emission times (one per output sample at the observer). Returns
    /// the same length as the input.
    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32>;
}
