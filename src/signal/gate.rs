//! Burst gate: a periodic on/off envelope that can wrap **any** signal model.
//!
//! **`SIGNAL_MODEL.md` E11.** The gate is a function of emission time, like
//! the signal it wraps: `s_g(t) = g(t)·s(t)` with `g(t) = 1` while
//! `(t − t_g) mod T` lies in one of the configured on-intervals
//! `[o_k, o_k + d_k)` and `0` otherwise. `t_g` is an absolute epoch, so the
//! burst pattern is fixed in absolute time; because the gate acts at emission,
//! a receiver sees every burst edge delayed by the light time, as it would a
//! real keyed transmitter.

use std::sync::Arc;

use anyhow::{bail, Result};
use ndarray::Array1;
use num_complex::Complex32;

use super::SignalModel;

/// The on/off pattern `g(t)`: a period and the on-intervals inside it.
#[derive(Debug, Clone, PartialEq)]
pub struct BurstGate {
    period_s: f64,
    /// On-intervals as `(start, end)` seconds within the period, `end` exclusive.
    bursts: Vec<(f64, f64)>,
    /// Time since the last period start at the render epoch (`t = 0`), in
    /// `[0, period)`.
    phase_s: f64,
}

impl BurstGate {
    /// `bursts` are `(offset_s, duration_s)` pairs within a period of
    /// `period_s`; `phase_s` is where in the period the render epoch falls.
    pub fn new(period_s: f64, bursts: &[(f64, f64)], phase_s: f64) -> Result<Self> {
        if !(period_s > 0.0) || !period_s.is_finite() {
            bail!("gate period_s must be positive and finite, got {period_s}");
        }
        let mut on = Vec::with_capacity(bursts.len());
        for &(offset_s, duration_s) in bursts {
            if !(offset_s >= 0.0) || !(duration_s > 0.0) || !(offset_s + duration_s <= period_s) {
                bail!(
                    "gate burst (offset_s = {offset_s}, duration_s = {duration_s}) must lie \
                     within one period: 0 <= offset_s, 0 < duration_s, \
                     offset_s + duration_s <= period_s = {period_s}"
                );
            }
            on.push((offset_s, offset_s + duration_s));
        }
        if !phase_s.is_finite() {
            bail!("gate phase_s must be finite, got {phase_s}");
        }
        Ok(Self { period_s, bursts: on, phase_s: phase_s.rem_euclid(period_s) })
    }

    /// A gate whose periods start at `gate_epoch_unix` (and every `period_s`
    /// before and after), as seen from a render starting at `render_epoch_unix`.
    pub fn anchored(
        period_s: f64,
        bursts: &[(f64, f64)],
        gate_epoch_unix: f64,
        render_epoch_unix: f64,
    ) -> Result<Self> {
        // `%` on f64 is exact, so the phase carries no rounding beyond the one
        // subtraction (itself exact for epochs within a factor of two).
        Self::new(period_s, bursts, render_epoch_unix - gate_epoch_unix)
    }

    /// `g(t)`: is the transmitter keyed on at render-epoch-relative emission
    /// time `t`?
    #[mantra_macros::req_satisfied("fn-8", "eq-e11")]
    #[inline]
    pub fn is_on(&self, t: f64) -> bool {
        let mut u = (t + self.phase_s).rem_euclid(self.period_s);
        if u >= self.period_s {
            u = 0.0; // rem_euclid of a tiny negative can round up to the period
        }
        self.bursts.iter().any(|&(start, end)| u >= start && u < end)
    }
}

/// A signal model seen through a [`BurstGate`]: the inner `s(t)` while the gate
/// is on, zero while it is off. The carrier is the inner signal's.
pub struct Gated {
    pub inner: Arc<dyn SignalModel>,
    pub gate: BurstGate,
}

impl SignalModel for Gated {
    fn carrier_hz(&self) -> f64 { self.inner.carrier_hz() }

    #[mantra_macros::req_satisfied("fn-4", "fn-8", "eq-e11")]
    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        let on: Vec<bool> = t_emit.iter().map(|&t| self.gate.is_on(t)).collect();
        if !on.iter().any(|&b| b) {
            // Whole batch inside a gap: nothing to evaluate.
            return Array1::<Complex32>::zeros(t_emit.len());
        }
        let mut out = self.inner.evaluate_baseband_batch(t_emit);
        for (v, &keep) in out.iter_mut().zip(&on) {
            if !keep {
                *v = Complex32::new(0.0, 0.0);
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signal::cw::ContinuousWave;

    fn three_bursts(phase_s: f64) -> BurstGate {
        BurstGate::new(15.0, &[(0.0, 0.5), (1.5, 0.5), (3.0, 0.5)], phase_s).unwrap()
    }

    #[mantra_macros::req_test("fn-8", "eq-e11")]
    #[test]
    fn on_exactly_inside_the_configured_intervals() {
        let g = three_bursts(0.0);
        // Interval starts are inclusive, ends exclusive.
        for (t, on) in [
            (0.0, true), (0.25, true), (0.499_999, true), (0.5, false), (1.0, false),
            (1.499_999, false), (1.5, true), (1.999_999, true), (2.0, false),
            (3.0, true), (3.25, true), (3.5, false), (10.0, false), (14.999_999, false),
        ] {
            assert_eq!(g.is_on(t), on, "t = {t}");
        }
        // Periodic for all time, both directions.
        for k in [-3.0, -1.0, 1.0, 2.0, 1.0e6] {
            assert!(g.is_on(15.0 * k + 0.25), "period {k}");
            assert!(g.is_on(15.0 * k + 1.75), "period {k}");
            assert!(!g.is_on(15.0 * k + 0.75), "period {k}");
            assert!(!g.is_on(15.0 * k - 0.25), "period {k}");
        }
    }

    #[mantra_macros::req_test("fn-8")]
    #[test]
    fn pattern_is_fixed_in_absolute_time() {
        // Period starts at Unix multiples of 15 s. A render starting 4 s into a
        // period sees the tail of that period's silence, then the next bursts.
        let epoch = 1_717_260_840.0; // a multiple of 15
        let g = BurstGate::anchored(15.0, &[(0.0, 0.5), (1.5, 0.5), (3.0, 0.5)], 0.0, epoch + 4.0)
            .unwrap();
        assert!(!g.is_on(0.0));
        assert!(!g.is_on(10.9));
        assert!(g.is_on(11.0));
        assert!(g.is_on(11.4));
        assert!(!g.is_on(11.5));
        assert!(g.is_on(12.5));
        // The same absolute instant, reached from two render epochs, agrees.
        let a = BurstGate::anchored(15.0, &[(2.0, 1.0)], 7.0, epoch).unwrap();
        let b = BurstGate::anchored(15.0, &[(2.0, 1.0)], 7.0, epoch + 123.0).unwrap();
        for i in 0..3000 {
            let t_abs_rel = 123.0 + i as f64 * 0.01;
            assert_eq!(a.is_on(t_abs_rel), b.is_on(t_abs_rel - 123.0), "t = {t_abs_rel}");
        }
    }

    #[mantra_macros::req_test("fn-8", "eq-e11", "fn-4")]
    #[test]
    fn gated_signal_is_inner_when_on_and_zero_when_off() {
        let inner = Arc::new(ContinuousWave {
            carrier_hz: 401.15e6, baseband_offset_hz: 50.0, anchor_phase_rad: None,
        });
        let gated = Gated { inner: inner.clone(), gate: three_bursts(0.0) };
        assert_eq!(gated.carrier_hz(), 401.15e6);
        let t = Array1::from((0..2000).map(|i| i as f64 * 0.01).collect::<Vec<_>>());
        let s = gated.evaluate_baseband_batch(&t);
        let raw = inner.evaluate_baseband_batch(&t);
        let mut n_on = 0;
        for i in 0..t.len() {
            if gated.gate.is_on(t[i]) {
                assert_eq!(s[i], raw[i], "on at t = {}", t[i]);
                n_on += 1;
            } else {
                assert_eq!(s[i], Complex32::new(0.0, 0.0), "off at t = {}", t[i]);
            }
        }
        // 20 s = one full period (150 on-samples) + the first 5 s of the next (150).
        assert_eq!(n_on, 300);
        // A batch wholly inside a gap is all zero.
        let gap = Array1::from(vec![5.0, 6.0, 14.0]);
        assert!(gated.evaluate_baseband_batch(&gap).iter().all(|v| v.norm() == 0.0));
    }

    #[test]
    fn bursts_outside_the_period_are_rejected() {
        assert!(BurstGate::new(15.0, &[(14.8, 0.5)], 0.0).is_err());
        assert!(BurstGate::new(15.0, &[(-0.1, 0.5)], 0.0).is_err());
        assert!(BurstGate::new(15.0, &[(1.0, 0.0)], 0.0).is_err());
        assert!(BurstGate::new(0.0, &[], 0.0).is_err());
        // No bursts at all is valid: a transmitter that never keys on.
        assert!(!BurstGate::new(15.0, &[], 0.0).unwrap().is_on(1.0));
    }
}
