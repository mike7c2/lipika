//! Signal models — purely algebraic descriptions of what a transmitter is
//! emitting, evaluable at any emission time `t_emit`.
//!
//! The contract is **time-domain and stateless**: the model is a function
//! `s(t)`, not a pre-sampled buffer. This makes the receiver-side clock
//! interpolation trivial (just call `s(t_corrected)`) and avoids any
//! resampling artefacts that pre-sampled-then-interpolated signals would
//! introduce.
//!
//! A model exposes its nominal RF centre frequency via `carrier_hz()`;
//! it does NOT itself apply the carrier — that's done downstream when
//! mixing to the receiver's baseband against the receiver's LO. The
//! model only produces complex baseband samples in the transmitter's
//! own frame.

use ndarray::Array1;
use num_complex::Complex32;

pub mod bpsk;
pub mod cw;
pub mod ca_code;

pub trait SignalModel: Send + Sync {
    /// Nominal centre frequency of the RF emission (Hz). The receiver
    /// will mix against its own LO at `(this − IF)` before recording IQ,
    /// so the Doppler shift falls out of delay-varying emission-time
    /// queries rather than being a separate frequency-shift step.
    fn carrier_hz(&self) -> f64;

    /// Evaluate the complex baseband signal at the requested emission
    /// times. Returns the same length as the input.
    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32>;
}
