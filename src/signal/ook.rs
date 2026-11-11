//! OOK (on-off keying) source: deterministic bit pattern (sync word + payload)
//! keying the carrier on and off at a fixed bit rate.

use anyhow::{bail, Result};
use ndarray::Array1;
use num_complex::Complex32;
use rand::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::SignalModel;

pub struct OokFrame {
    /// Total bit stream: sync word followed by payload, 1 (carrier on) / 0
    /// (carrier off) valued.
    pub bits: Vec<u8>,
    pub bit_rate_bps: f64,
    /// Bit-0 time tag in **seconds relative to the scene epoch** (see
    /// [`super::bpsk::BpskFrame::epoch_rel_s`]).
    pub epoch_rel_s: f64,
    pub carrier_hz: f64,
}

impl OokFrame {
    pub fn random(
        sync_word: &[u8],
        n_payload_bits: usize,
        bit_rate_bps: f64,
        epoch_rel_s: f64,
        carrier_hz: f64,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Result<Self> {
        if let Some(bad) = sync_word.iter().find(|&&b| b > 1) {
            bail!("ook sync_word entries must be 0 or 1, got {bad}");
        }
        let mut bits = Vec::with_capacity(sync_word.len() + n_payload_bits);
        bits.extend_from_slice(sync_word);
        for _ in 0..n_payload_bits {
            bits.push(rng.gen_bool(0.5) as u8);
        }
        Ok(Self { bits, bit_rate_bps, epoch_rel_s, carrier_hz })
    }
}

impl SignalModel for OokFrame {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    #[mantra_macros::req_satisfied("fn-4", "fn-4-ook")]
    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        // Bit index = floor((t_emit - epoch) * bit_rate). The on state has unit
        // amplitude, so the transmitter's EIRP is the keyed-on carrier power.
        let n_bits = self.bits.len() as i64;
        let mut out = Array1::<Complex32>::zeros(t_emit.len());
        for (i, &t) in t_emit.iter().enumerate() {
            let bit_idx_f = (t - self.epoch_rel_s) * self.bit_rate_bps;
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn frame(sync: &[u8], n_payload: usize, seed: u64) -> OokFrame {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
        OokFrame::random(sync, n_payload, 1_000.0, 0.0, 433.92e6, &mut rng).unwrap()
    }

    #[mantra_macros::req_test("fn-4-ook")]
    #[test]
    fn carrier_is_keyed_by_sync_word_then_payload() {
        let sync = [1, 0, 1, 1, 0, 0, 1, 0];
        let f = frame(&sync, 64, 3);
        assert_eq!(f.carrier_hz(), 433.92e6);
        assert_eq!(f.bits.len(), sync.len() + 64);
        assert_eq!(&f.bits[..sync.len()], &sync);
        // Sample mid-bit: the baseband is exactly the bit, on = 1 + 0j, off = 0.
        let t = Array1::from((0..f.bits.len()).map(|k| (k as f64 + 0.5) * 1e-3).collect::<Vec<_>>());
        let s = f.evaluate_baseband_batch(&t);
        for (k, &b) in f.bits.iter().enumerate() {
            assert_eq!(s[k], Complex32::new(b as f32, 0.0), "bit {k}");
        }
        // The seeded payload is not degenerate: both states occur.
        let ones = f.bits[sync.len()..].iter().filter(|&&b| b == 1).count();
        assert!(ones > 16 && ones < 48, "payload has {ones}/64 ones");
    }

    #[mantra_macros::req_test("fn-4-ook")]
    #[test]
    fn silent_outside_the_frame_and_seed_deterministic() {
        let f = frame(&[1, 1], 6, 9);
        // Before bit 0 and after the last bit the carrier is off.
        let t = Array1::from(vec![-1e-6, 0.0, 7.999e-3, 8.0e-3, 1.0]);
        let s = f.evaluate_baseband_batch(&t);
        assert_eq!(s[0], Complex32::new(0.0, 0.0));
        assert_eq!(s[1], Complex32::new(1.0, 0.0));
        assert_eq!(s[2], Complex32::new(f.bits[7] as f32, 0.0));
        assert_eq!(s[3], Complex32::new(0.0, 0.0));
        assert_eq!(s[4], Complex32::new(0.0, 0.0));
        assert_eq!(f.bits, frame(&[1, 1], 6, 9).bits);
    }

    #[test]
    fn non_binary_sync_word_is_rejected() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        assert!(OokFrame::random(&[1, 0, 2], 8, 1_000.0, 0.0, 1.0e9, &mut rng).is_err());
    }
}
