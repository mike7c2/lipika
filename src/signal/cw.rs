//! Continuous-wave (single-tone) source. Useful for sanity checks.

use ndarray::Array1;
use num_complex::Complex32;

use super::SignalModel;

pub struct ContinuousWave {
    pub carrier_hz: f64,
    /// Baseband-frequency offset relative to the nominal carrier (Hz).
    /// Zero by default — a clean carrier at `carrier_hz` produces a pure
    /// DC tone in the source's baseband frame.
    pub baseband_offset_hz: f64,
}

impl SignalModel for ContinuousWave {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    #[mantra_macros::req_satisfied("fn-4", "fn-4-cw")]
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[mantra_macros::req_test("fn-4-cw")]
    #[test]
    fn tone_is_exp_j_two_pi_f_offset_t() {
        let cw = ContinuousWave { carrier_hz: 437.0e6, baseband_offset_hz: 1_000.0 };
        assert_eq!(cw.carrier_hz(), 437.0e6);
        let t = Array1::from(vec![0.0, 0.25e-3, 0.5e-3, 0.7e-3, 1.0e-3]);
        let s = cw.evaluate_baseband_batch(&t);
        for (i, &ti) in t.iter().enumerate() {
            let ph = std::f64::consts::TAU * 1_000.0 * ti;
            assert_abs_diff_eq!(s[i].re, ph.cos() as f32, epsilon = 1e-5);
            assert_abs_diff_eq!(s[i].im, ph.sin() as f32, epsilon = 1e-5);
            assert_abs_diff_eq!(s[i].norm(), 1.0, epsilon = 1e-5);
        }
    }

    #[mantra_macros::req_test("fn-4")]
    #[test]
    fn evaluation_is_stateless_and_batch_independent() {
        // Algebraic s(t): the value at a time depends on that time alone, not on
        // which other times share the batch or on the order they are asked in.
        let cw = ContinuousWave { carrier_hz: 1.0e9, baseband_offset_hz: 123.0 };
        let all = Array1::from((0..100).map(|i| i as f64 * 1e-4).collect::<Vec<_>>());
        let whole = cw.evaluate_baseband_batch(&all);
        let (a, b) = all.view().split_at(ndarray::Axis(0), 37);
        let part_a = cw.evaluate_baseband_batch(&a.to_owned());
        let part_b = cw.evaluate_baseband_batch(&b.to_owned());
        for i in 0..100 {
            let split = if i < 37 { part_a[i] } else { part_b[i - 37] };
            assert_eq!(whole[i], split);
        }
        let reversed = Array1::from(all.iter().rev().copied().collect::<Vec<_>>());
        let back = cw.evaluate_baseband_batch(&reversed);
        for i in 0..100 {
            assert_eq!(whole[i], back[99 - i]);
        }
    }
}
