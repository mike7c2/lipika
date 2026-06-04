//! Receiver front-end as an **algebraic operator chain** — the dual of
//! [`crate::signal::SignalModel`].
//!
//! lipika models a transmitter as a continuous function `s(t)` and lets Doppler /
//! delay fall out of *where* it is evaluated; the only collapse to discrete samples
//! is the output grid. The receiver is the mirror image: a composition of
//! continuous-time operators on the incident field, with the **ADC sample+quantise
//! as the single discrete collapse** and everything upstream exact.
//!
//! Operators split by where they act:
//!   - **per-source** (woven into the channel contribution): antenna gain, the
//!     oscillator time-warp `δ(t)` and LO phase `φ(t)` — *added in later steps*;
//!   - **post-summation** (on the digitised band buffer): gain/AGC, nonlinearity,
//!     IQ imbalance, and the ADC collapse — the [`PostOp`] chain applied here.
//!
//! `RxChain` is a **value** built from config; this module is the shared definition
//! that lipika interprets forward and (Phase B) gnss-rx will interpret inverse, so
//! the simulator and the receiver can never drift.

use crate::config::FrontendConfig;
use ndarray::Array1;
use num_complex::Complex32;

/// A post-summation operator: acts on the summed band buffer (volts).
#[derive(Debug, Clone, Copy)]
pub enum PostOp {
    /// Linear RF/IF gain + AGC in dB — brings the (tiny, FSPL-attenuated) signal
    /// into the ADC's input range. The matching inverse for an estimator is trivial.
    Gain { db: f64 },
    /// The **ADC collapse**: clip to ±`full_scale_v`, then quantise to `bits`
    /// (symmetric mid-tread). This is the one discrete, lossy, non-invertible step
    /// in the entire receiver chain. Deterministic per sample ⇒ chunk/thread
    /// invariance is preserved.
    Adc { bits: u32, full_scale_v: f64 },
}

impl PostOp {
    fn apply(&self, iq: &mut Array1<Complex32>) {
        match *self {
            PostOp::Gain { db } => {
                let g = 10f32.powf((db / 20.0) as f32);
                iq.iter_mut().for_each(|c| *c *= g);
            }
            PostOp::Adc { bits, full_scale_v } => {
                let fs = full_scale_v as f32;
                let levels = ((1u64 << (bits.max(2) - 1)) - 1) as f32;
                let lsb = fs / levels;
                let q = |x: f32| ((x.clamp(-fs, fs)) / lsb).round() * lsb;
                iq.iter_mut().for_each(|c| *c = Complex32::new(q(c.re), q(c.im)));
            }
        }
    }
}

/// The receiver chain. For now the post-summation operator list; per-source
/// operators (antenna, oscillator) attach in the channel contribution in a later
/// step. Empty ⇒ identity (output byte-identical to no front-end).
#[derive(Debug, Clone, Default)]
pub struct RxChain {
    pub post_sum: Vec<PostOp>,
}

impl RxChain {
    /// Build the chain from a band's optional `[receiver.band.frontend]` config.
    pub fn from_config(fe: Option<&FrontendConfig>) -> Self {
        let mut post_sum = Vec::new();
        if let Some(fe) = fe {
            if let Some(db) = fe.gain_db {
                post_sum.push(PostOp::Gain { db });
            }
            if let Some(adc) = &fe.adc {
                post_sum.push(PostOp::Adc { bits: adc.bits, full_scale_v: adc.full_scale_v });
            }
        }
        RxChain { post_sum }
    }

    /// Forward-interpret the post-summation operators onto the band buffer, in
    /// chain order. No-op when the chain is empty.
    pub fn apply_post_sum(&self, iq: &mut Array1<Complex32>) {
        for op in &self.post_sum {
            op.apply(iq);
        }
    }

    pub fn is_identity(&self) -> bool {
        self.post_sum.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{AdcConfig, FrontendConfig};

    #[test]
    fn empty_chain_is_identity() {
        let chain = RxChain::from_config(None);
        assert!(chain.is_identity());
        let mut iq = Array1::from_vec(vec![Complex32::new(0.3, -0.7), Complex32::new(1.1, 0.2)]);
        let before = iq.clone();
        chain.apply_post_sum(&mut iq);
        assert_eq!(iq, before);
    }

    #[test]
    fn adc_quantises_to_levels_and_clips() {
        // 4-bit ADC, full scale ±1 V: 7 positive levels, LSB = 1/7.
        let chain = RxChain::from_config(Some(&FrontendConfig {
            gain_db: None,
            adc: Some(AdcConfig { bits: 4, full_scale_v: 1.0 }),
        }));
        let lsb = 1.0f32 / 7.0;
        let mut iq = Array1::from_vec(vec![
            Complex32::new(0.5, -0.5),  // mid range -> nearest level
            Complex32::new(2.0, -3.0),  // beyond full scale -> clips to ±1
        ]);
        chain.apply_post_sum(&mut iq);
        // every output component is an integer multiple of the LSB, within ±1 V
        for c in iq.iter() {
            for v in [c.re, c.im] {
                assert!(v.abs() <= 1.0 + 1e-6);
                assert!(((v / lsb).round() * lsb - v).abs() < 1e-6, "v={v} not on LSB grid");
            }
        }
        assert_eq!(iq[1], Complex32::new(1.0, -1.0)); // clipped
    }

    #[test]
    fn gain_scales_linearly() {
        let chain = RxChain::from_config(Some(&FrontendConfig {
            gain_db: Some(20.0),  // ×10
            adc: None,
        }));
        let mut iq = Array1::from_vec(vec![Complex32::new(0.1, -0.2)]);
        chain.apply_post_sum(&mut iq);
        assert!((iq[0].re - 1.0).abs() < 1e-5 && (iq[0].im + 2.0).abs() < 1e-5);
    }
}
