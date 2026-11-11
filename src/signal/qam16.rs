//! 16-QAM source: deterministic symbol pattern (sync word + payload) on a
//! square, Gray-coded constellation at a fixed symbol rate.

use anyhow::{bail, Result};
use ndarray::Array1;
use num_complex::Complex32;
use rand::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;

use super::SignalModel;

/// Bits carried by one 16-QAM symbol.
pub const BITS_PER_SYMBOL: usize = 4;

/// Scale that gives the `{±1, ±3}²` constellation unit mean power
/// (`E|s|² = 10` unscaled), so a transmitter's EIRP is its average power.
const UNIT_POWER_SCALE: f32 = 0.316_227_77; // 1/√10

/// Gray-coded 4-level axis: adjacent levels differ in one bit.
fn gray_level(b0: u8, b1: u8) -> f32 {
    match (b0, b1) {
        (0, 0) => -3.0,
        (0, 1) => -1.0,
        (1, 1) => 1.0,
        _ => 3.0,
    }
}

/// Map four bits `[b0, b1, b2, b3]` to a constellation point: `b0 b1` select
/// the I level and `b2 b3` the Q level, each Gray-coded, scaled to unit mean
/// power.
pub fn map_symbol(bits: [u8; BITS_PER_SYMBOL]) -> Complex32 {
    Complex32::new(gray_level(bits[0], bits[1]), gray_level(bits[2], bits[3])) * UNIT_POWER_SCALE
}

pub struct Qam16Frame {
    /// Total symbol stream: sync word followed by payload, already mapped to
    /// unit-mean-power constellation points.
    pub symbols: Vec<Complex32>,
    pub symbol_rate_baud: f64,
    /// Symbol-0 time tag in **seconds relative to the scene epoch** (see
    /// [`super::bpsk::BpskFrame::epoch_rel_s`]).
    pub epoch_rel_s: f64,
    pub carrier_hz: f64,
}

impl Qam16Frame {
    /// `sync_word` is a 0/1 bit sequence, four bits per symbol.
    pub fn random(
        sync_word: &[u8],
        n_payload_symbols: usize,
        symbol_rate_baud: f64,
        epoch_rel_s: f64,
        carrier_hz: f64,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Result<Self> {
        if let Some(bad) = sync_word.iter().find(|&&b| b > 1) {
            bail!("qam16 sync_word entries must be 0 or 1, got {bad}");
        }
        if sync_word.len() % BITS_PER_SYMBOL != 0 {
            bail!(
                "qam16 sync_word must hold a whole number of symbols \
                 ({BITS_PER_SYMBOL} bits each), got {} bits",
                sync_word.len()
            );
        }
        let mut symbols =
            Vec::with_capacity(sync_word.len() / BITS_PER_SYMBOL + n_payload_symbols);
        for b in sync_word.chunks_exact(BITS_PER_SYMBOL) {
            symbols.push(map_symbol([b[0], b[1], b[2], b[3]]));
        }
        for _ in 0..n_payload_symbols {
            let v: u8 = rng.gen_range(0..16);
            symbols.push(map_symbol([(v >> 3) & 1, (v >> 2) & 1, (v >> 1) & 1, v & 1]));
        }
        Ok(Self { symbols, symbol_rate_baud, epoch_rel_s, carrier_hz })
    }
}

impl SignalModel for Qam16Frame {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    #[mantra_macros::req_satisfied("fn-4", "fn-4-qam16")]
    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        // Symbol index = floor((t_emit - epoch) * symbol_rate).
        let n_symbols = self.symbols.len() as i64;
        let mut out = Array1::<Complex32>::zeros(t_emit.len());
        for (i, &t) in t_emit.iter().enumerate() {
            let sym_idx_f = (t - self.epoch_rel_s) * self.symbol_rate_baud;
            let sym_idx = sym_idx_f.floor() as i64;
            if sym_idx < 0 || sym_idx >= n_symbols {
                out[i] = Complex32::new(0.0, 0.0);
            } else {
                out[i] = self.symbols[sym_idx as usize];
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use rand::SeedableRng;

    fn all_points() -> Vec<([u8; 4], Complex32)> {
        (0..16u8)
            .map(|v| {
                let bits = [(v >> 3) & 1, (v >> 2) & 1, (v >> 1) & 1, v & 1];
                (bits, map_symbol(bits))
            })
            .collect()
    }

    #[mantra_macros::req_test("fn-4-qam16")]
    #[test]
    fn constellation_is_square_unit_power_and_gray_coded() {
        let pts = all_points();
        // 16 distinct points on the {±1, ±3}/√10 grid.
        let levels = [-3.0_f32, -1.0, 1.0, 3.0].map(|l| l / 10.0_f32.sqrt());
        for (i, (_, p)) in pts.iter().enumerate() {
            assert!(levels.iter().any(|l| (l - p.re).abs() < 1e-6), "I level {}", p.re);
            assert!(levels.iter().any(|l| (l - p.im).abs() < 1e-6), "Q level {}", p.im);
            for (_, q) in &pts[..i] {
                assert!((p - q).norm() > 0.1, "duplicate constellation point");
            }
        }
        let mean_power: f32 = pts.iter().map(|(_, p)| p.norm_sqr()).sum::<f32>() / 16.0;
        assert_abs_diff_eq!(mean_power, 1.0, epsilon = 1e-6);
        // Gray: nearest neighbours (one grid step apart) differ in exactly one bit.
        let step = 2.0 / 10.0_f32.sqrt();
        for (ba, pa) in &pts {
            for (bb, pb) in &pts {
                if ((pa - pb).norm() - step).abs() < 1e-5 {
                    let differing = ba.iter().zip(bb).filter(|(x, y)| x != y).count();
                    assert_eq!(differing, 1, "{ba:?} vs {bb:?}");
                }
            }
        }
        // The frame emits exactly these points: a sync word walking all 16.
        let sync: Vec<u8> = pts.iter().flat_map(|(b, _)| *b).collect();
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        let f = Qam16Frame::random(&sync, 0, 1_000.0, 0.0, 2.4e9, &mut rng).unwrap();
        let t = Array1::from((0..16).map(|k| (k as f64 + 0.5) * 1e-3).collect::<Vec<_>>());
        let s = f.evaluate_baseband_batch(&t);
        for (k, (_, p)) in pts.iter().enumerate() {
            assert_eq!(s[k], *p);
        }
    }

    #[mantra_macros::req_test("fn-4-qam16")]
    #[test]
    fn sync_word_then_seeded_payload_at_symbol_rate() {
        let sync = [0, 0, 0, 0, 1, 0, 1, 0]; // (−3,−3), (+3,+3)
        let mk = || {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
            Qam16Frame::random(&sync, 4000, 2_000.0, 0.0, 2.4e9, &mut rng).unwrap()
        };
        let f = mk();
        assert_eq!(f.carrier_hz(), 2.4e9);
        assert_eq!(f.symbols.len(), 2 + 4000);
        let a = 3.0 / 10.0_f32.sqrt();
        // Symbols last 1/2000 s; silent before symbol 0 and after the last.
        let t = Array1::from(vec![-1e-9, 0.0, 0.499e-3, 0.5e-3, 0.999e-3, 2001.0, 2001.5]);
        let s = f.evaluate_baseband_batch(&t);
        let zero = Complex32::new(0.0, 0.0);
        assert_eq!(s[0], zero);
        for i in [1, 2] {
            assert_abs_diff_eq!(s[i].re, -a, epsilon = 1e-6);
            assert_abs_diff_eq!(s[i].im, -a, epsilon = 1e-6);
        }
        for i in [3, 4] {
            assert_abs_diff_eq!(s[i].re, a, epsilon = 1e-6);
            assert_abs_diff_eq!(s[i].im, a, epsilon = 1e-6);
        }
        assert_eq!(s[5], zero);
        assert_eq!(s[6], zero);
        // Seeded payload: reproducible, uses the whole constellation, ~unit power.
        assert_eq!(f.symbols, mk().symbols);
        let payload = &f.symbols[2..];
        let mean_power: f32 =
            payload.iter().map(|p| p.norm_sqr()).sum::<f32>() / payload.len() as f32;
        assert_abs_diff_eq!(mean_power, 1.0, epsilon = 0.05);
        for (_, p) in all_points() {
            assert!(payload.iter().any(|q| (q - p).norm() < 1e-6), "point {p} never drawn");
        }
    }

    #[test]
    fn malformed_sync_word_is_rejected() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);
        assert!(Qam16Frame::random(&[1, 0, 1], 8, 1_000.0, 0.0, 1.0e9, &mut rng).is_err());
        assert!(Qam16Frame::random(&[1, 0, 2, 0], 8, 1_000.0, 0.0, 1.0e9, &mut rng).is_err());
    }
}
