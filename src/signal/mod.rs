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
pub mod gate;
pub mod gps_nav;
pub mod ook;
pub mod qam16;

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

/// Split the number of symbol periods elapsed over `elapsed_s` at `rate_hz`
/// into a whole count and a fraction in `[0, 1)`.
///
/// This is how a signal is anchored to an **absolute** epoch without carrying
/// Unix-magnitude times into the per-sample arithmetic (`fid-2`): the render
/// epoch sits `elapsed_s` after the signal's own epoch, the whole count becomes
/// an integer symbol offset and the fraction a sub-symbol time tag. The product
/// is formed exactly (an FMA recovers its rounding error), so for a given
/// absolute instant the result does not depend on which render epoch it was
/// reached from.
#[mantra_macros::req_satisfied("fn-9")]
pub fn elapsed_periods(elapsed_s: f64, rate_hz: f64) -> anyhow::Result<(i64, f64)> {
    let p = elapsed_s * rate_hz;
    if !p.is_finite() || p.abs() >= 9.0e18 {
        anyhow::bail!(
            "cannot anchor a signal {elapsed_s} s from its epoch at {rate_hz} Hz \
             (symbol count is not representable)"
        );
    }
    let err = elapsed_s.mul_add(rate_hz, -p);
    let mut whole = p.floor();
    let mut frac = (p - whole) + err;
    if frac < 0.0 {
        whole -= 1.0;
        frac += 1.0;
    } else if frac >= 1.0 {
        whole += 1.0;
        frac -= 1.0;
    }
    Ok((whole as i64, frac))
}

/// Where a framed signal (`bpsk`, `ook`, `qam16`) sits relative to the render
/// epoch when it is anchored to an absolute epoch and repeated for all time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameAnchor {
    /// Time tag (render-epoch-relative seconds, in `(−1/rate, 0]`) of the
    /// symbol boundary at or just before the render epoch.
    pub epoch_rel_s: f64,
    /// Index within the frame of the symbol that starts at `epoch_rel_s`.
    pub symbol_offset: i64,
}

impl FrameAnchor {
    /// Anchor for a frame of `n_symbols` at `rate_hz` whose symbol 0 was
    /// emitted at `signal_epoch_unix` (and every `n_symbols / rate_hz` seconds
    /// before and after), as seen from a render starting at `render_epoch_unix`.
    #[mantra_macros::req_satisfied("fn-9")]
    pub fn at(
        render_epoch_unix: f64,
        signal_epoch_unix: f64,
        rate_hz: f64,
        n_symbols: usize,
    ) -> anyhow::Result<Self> {
        if !(rate_hz > 0.0) {
            anyhow::bail!("symbol rate must be positive, got {rate_hz}");
        }
        let (whole, frac) = elapsed_periods(render_epoch_unix - signal_epoch_unix, rate_hz)?;
        Ok(Self {
            epoch_rel_s: -frac / rate_hz,
            symbol_offset: whole.rem_euclid((n_symbols as i64).max(1)),
        })
    }
}

/// Resolve a raw symbol count `⌊(t − epoch)·rate⌋` to an index into a frame of
/// `n` symbols. `repeat_from = None` is a one-shot frame (silent before symbol
/// 0 and after the last); `Some(k)` is the periodic extension in which raw
/// count 0 is frame symbol `k`.
#[inline]
pub(crate) fn frame_index(raw: i64, n: i64, repeat_from: Option<i64>) -> Option<usize> {
    match repeat_from {
        None => (raw >= 0 && raw < n).then_some(raw as usize),
        Some(_) if n <= 0 => None,
        Some(k) => Some(raw.wrapping_add(k).rem_euclid(n) as usize),
    }
}

#[cfg(test)]
mod tests {
    use super::bpsk::BpskFrame;
    use super::cw::ContinuousWave;
    use super::ook::OokFrame;
    use super::qam16::Qam16Frame;
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[mantra_macros::req_test("fn-9")]
    #[test]
    fn elapsed_periods_splits_exactly() {
        // Whole seconds at an integer rate: an exact whole count.
        assert_eq!(elapsed_periods(1_717_260_840.0, 4800.0).unwrap(), (8_242_852_032_000, 0.0));
        // A fractional part survives at Unix magnitude: 2⁻²⁰ s is one f64 step
        // or more there, and at 4800 Hz is 4800/2²⁰ of a period.
        let (whole, frac) = elapsed_periods(1_717_260_840.0 + 1.0 / 1_048_576.0, 4800.0).unwrap();
        assert_eq!(whole, 8_242_852_032_000);
        assert!((frac - 4800.0 / 1_048_576.0).abs() < 1e-15, "frac = {frac}");
        // Negative elapsed time (signal epoch after the render epoch) floors down.
        let (whole, frac) = elapsed_periods(-0.25, 2.0).unwrap();
        assert_eq!((whole, frac), (-1, 0.5));
        // A rate that is not a binary fraction: the FMA keeps the fraction
        // consistent with the exact product, so stepping the elapsed time by
        // whole periods changes only the whole count.
        let rate = 9600.0 / 7.0;
        let (w0, f0) = elapsed_periods(1.0e9, rate).unwrap();
        let (w1, f1) = elapsed_periods(1.0e9 + 7.0, rate).unwrap();
        assert_eq!(w1 - w0, 9600);
        assert!((f1 - f0).abs() < 1e-9, "{f0} vs {f1}");
        assert!(elapsed_periods(f64::INFINITY, 1.0).is_err());
        assert!(elapsed_periods(1.0e12, 1.0e12).is_err());
    }

    /// Evaluate `signal_at(render_epoch)` at the same absolute instants from
    /// two render epochs and require the same values.
    fn assert_same_from_two_epochs<S: SignalModel>(signal_at: impl Fn(f64) -> S, rate_hz: f64) {
        let epoch_a = 1_717_260_840.0;
        let shift = 1234.5; // exact in f64, not a whole number of frames
        let (a, b) = (signal_at(epoch_a), signal_at(epoch_a + shift));
        // Mid-symbol instants spanning several frame repeats, on both sides of
        // render epoch B (negative times for B are before its window start).
        let t_a = Array1::from(
            (0..4000).map(|k| shift - 1.0 + (k as f64 + 0.5) / rate_hz).collect::<Vec<_>>());
        let t_b = t_a.mapv(|t| t - shift);
        let (sa, sb) = (a.evaluate_baseband_batch(&t_a), b.evaluate_baseband_batch(&t_b));
        for k in 0..t_a.len() {
            assert!((sa[k] - sb[k]).norm() < 1e-4, "sample {k}: {} vs {}", sa[k], sb[k]);
        }
        assert!(sa.iter().any(|v| v.norm() > 0.5), "signal should be radiating");
    }

    #[mantra_macros::req_test("fn-9", "fn-4-bpsk")]
    #[test]
    fn anchored_bpsk_repeats_for_all_time_and_ignores_the_render_epoch() {
        let sync = [-1i8, 1, 1, -1, 1];
        let make = |render_epoch: f64, signal_epoch: f64| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(11);
            let f = BpskFrame::random(&sync, 59, 1_000.0, 0.0, 401.15e6, &mut rng);
            let n = f.bits.len();
            f.repeating(FrameAnchor::at(render_epoch, signal_epoch, 1_000.0, n).unwrap())
        };
        // Render epoch == signal epoch: s(t) = bits[⌊t·R⌋ mod n], for negative
        // times and far past the end of the frame alike (never silent).
        let f = make(1.0e9, 1.0e9);
        let n = f.bits.len() as i64; // 64
        for k in [-200i64, -65, -64, -1, 0, 1, 63, 64, 65, 1_000_003] {
            let s = f.evaluate_baseband_batch(&Array1::from(vec![(k as f64 + 0.5) * 1e-3]));
            assert_eq!(s[0].re, f.bits[k.rem_euclid(n) as usize] as f32, "bit {k}");
        }
        // Signal epoch 2.5 bits + 10 frames before the render epoch: the bit
        // under way at t = 0 is frame bit 2, and bit 3 starts half a bit later.
        let g = make(1.0e9, 1.0e9 - 0.64 - 0.0025);
        assert_eq!(g.repeat_from, Some(2));
        // (The epoch itself is only as fine as an f64 Unix time: ~119 ns here.)
        assert!((g.epoch_rel_s + 0.0005).abs() < 1.2e-7, "{}", g.epoch_rel_s);
        let s = g.evaluate_baseband_batch(&Array1::from(vec![0.0, 0.0003, 0.0007]));
        assert_eq!([s[0].re, s[1].re, s[2].re],
                   [g.bits[2] as f32, g.bits[2] as f32, g.bits[3] as f32]);
        assert_same_from_two_epochs(|e| make(e, 1_000_000_007.25), 1_000.0);
    }

    #[mantra_macros::req_test("fn-9", "fn-4-ook", "fn-4-qam16", "fn-4-cw")]
    #[test]
    fn anchored_ook_qam16_and_cw_ignore_the_render_epoch() {
        let signal_epoch = 1_000_000_007.25;
        assert_same_from_two_epochs(|e| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(5);
            let f = OokFrame::random(&[1, 0, 1, 1], 93, 2_400.0, 0.0, 433.92e6, &mut rng).unwrap();
            let n = f.bits.len();
            f.repeating(FrameAnchor::at(e, signal_epoch, 2_400.0, n).unwrap())
        }, 2_400.0);
        assert_same_from_two_epochs(|e| {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(6);
            let f = Qam16Frame::random(&[1, 0, 1, 1], 30, 1_200.0, 0.0, 915.0e6, &mut rng).unwrap();
            let n = f.symbols.len();
            f.repeating(FrameAnchor::at(e, signal_epoch, 1_200.0, n).unwrap())
        }, 1_200.0);
        assert_same_from_two_epochs(|e| ContinuousWave {
            carrier_hz: 437.0e6,
            baseband_offset_hz: 1_234.5,
            anchor_phase_rad: Some(ContinuousWave::anchor_phase(e, signal_epoch, 1_234.5).unwrap()),
        }, 20_000.0);
        // At its own epoch an anchored tone has zero phase.
        let p0 = ContinuousWave::anchor_phase(signal_epoch, signal_epoch, 1_234.5).unwrap();
        assert_eq!(p0, 0.0);
    }

    #[test]
    fn empty_repeating_frame_is_silent() {
        assert_eq!(frame_index(5, 0, Some(0)), None);
        assert_eq!(frame_index(5, 0, None), None);
    }
}
