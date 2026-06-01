//! BPSK source: deterministic bit pattern (sync word + payload) modulated
//! at a fixed bit rate.

use ndarray::Array1;
use num_complex::Complex32;
use rand::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::SignalModel;

pub struct BpskFrame {
    /// Total bit stream: sync word followed by payload, +1 / −1 valued.
    pub bits: Vec<i8>,
    pub bit_rate_bps: f64,
    pub epoch_unix: f64,
    pub carrier_hz: f64,
}

impl BpskFrame {
    pub fn random(
        sync_word: &[i8],
        n_payload_bits: usize,
        bit_rate_bps: f64,
        epoch_unix: f64,
        carrier_hz: f64,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Self {
        let mut bits = Vec::with_capacity(sync_word.len() + n_payload_bits);
        bits.extend_from_slice(sync_word);
        for _ in 0..n_payload_bits {
            bits.push(if rng.gen_bool(0.5) { 1 } else { -1 });
        }
        Self { bits, bit_rate_bps, epoch_unix, carrier_hz }
    }
}

impl SignalModel for BpskFrame {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        // Bit index = floor((t_emit - epoch) * bit_rate).
        let n_bits = self.bits.len() as i64;
        let mut out = Array1::<Complex32>::zeros(t_emit.len());
        for (i, &t) in t_emit.iter().enumerate() {
            let bit_idx_f = (t - self.epoch_unix) * self.bit_rate_bps;
            let bit_idx = bit_idx_f.floor() as i64;
            if bit_idx < 0 || bit_idx >= n_bits {
                out[i] = Complex32::new(0.0, 0.0);
            } else {
                out[i] = Complex32::new(self.bits[bit_idx as usize] as f32, 0.0);
            }
        }
        out
    }
}
