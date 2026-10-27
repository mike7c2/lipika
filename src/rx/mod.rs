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
//!   - **per-source** (woven into the channel contribution): antenna gain (E6) and
//!     the oscillator `δ(t)`/LO phase (E7, in [`crate::clock`]);
//!   - **post-summation** (on the digitised band buffer): gain/AGC, the IP3
//!     nonlinearity (E8) and the ADC collapse (E9) — the [`PostOp`] chain here.
//!
//! `RxChain` is a **value** built from config; written so an estimator can invert
//! the same operators (E-tags refer to `SIGNAL_MODEL.md`).

use crate::config::FrontendConfig;
use ndarray::Array1;
use num_complex::Complex32;

/// Parameterised antenna gain pattern (`SIGNAL_MODEL.md` E6; azimuthally
/// symmetric — gain vs elevation),
/// modelling the chosen self-phased QFH: peak at zenith, roll-off to the horizon,
/// a front/back ratio below it. A **per-source** operator: applied in the channel
/// contribution using each source's arrival elevation (the same angle the horizon
/// mask uses), so it scales the link budget per source, per band.
#[derive(Debug, Clone, Copy)]
pub struct Antenna {
    /// Boresight (zenith) gain, dBi.
    pub peak_dbi: f64,
    /// Gain at the local horizon (elevation 0°), dBi.
    pub horizon_dbi: f64,
    /// Front/back ratio (dB): below-horizon gain = `horizon_dbi − fb_db`.
    pub fb_db: f64,
    /// Optional lower cutoff elevation (degrees): below it the gain is 0. This is
    /// how **Earth line-of-sight occlusion** is modelled — as a property of the
    /// antenna, not a separate channel pass. `None` = no cutoff (back-lobe to −90°).
    pub cutoff_deg: Option<f64>,
}

impl Antenna {
    /// Voltage gain factor √G at arrival elevation `elev_deg` ∈ [−90, 90].
    /// Linear interpolation in dBi from horizon→zenith above the horizon, a flat
    /// back-lobe below, and **0 below `cutoff_deg`** (occlusion). (√G = 10^(dBi/20).)
    // r[impl eq.e6] antenna gain pattern G_ant(θ)
    pub fn amplitude_at(&self, elev_deg: f64) -> f32 {
        if let Some(c) = self.cutoff_deg {
            if elev_deg < c { return 0.0; }
        }
        let dbi = if elev_deg >= 0.0 {
            self.horizon_dbi + (self.peak_dbi - self.horizon_dbi) * (elev_deg / 90.0)
        } else {
            self.horizon_dbi - self.fb_db
        };
        10f64.powf(dbi / 20.0) as f32
    }
}

/// A post-summation operator: acts on the summed band buffer (volts).
#[derive(Debug, Clone, Copy)]
pub enum PostOp {
    /// Linear RF/IF gain + AGC in dB — brings the (tiny, FSPL-attenuated) signal
    /// into the ADC's input range. The matching inverse for an estimator is trivial.
    Gain { db: f64 },
    /// Memoryless third-order nonlinearity (`SIGNAL_MODEL.md` E8; LNA/mixer):
    /// `y = x·(1 − (4/3)|x|²/Pᵢₚ₃)`,
    /// `iip3_w` = the input third-order intercept *power* (watts; lipika amplitude²
    /// is power, so this is direct from the dBm spec). Small-signal linear; models
    /// gain compression + IM3 at strong input. The binding spec for terrestrial
    /// rejection (antenna doc §R2). Valid up to ~compression; the ADC bounds extremes.
    // r[impl eq.e8]
    Nonlinearity { iip3_w: f64 },
    /// The **ADC collapse** (`SIGNAL_MODEL.md` E9): clip to ±`full_scale_v`, then
    /// quantise to `bits` (symmetric mid-tread). The one discrete, lossy, non-invertible step
    /// in the entire receiver chain. Deterministic per sample ⇒ chunk/thread
    /// invariance is preserved.
    // r[impl eq.e9]
    Adc { bits: u32, full_scale_v: f64 },
}

impl PostOp {
    fn apply(&self, iq: &mut Array1<Complex32>) {
        match *self {
            PostOp::Gain { db } => {
                let g = 10f32.powf((db / 20.0) as f32);
                iq.iter_mut().for_each(|c| *c *= g);
            }
            PostOp::Nonlinearity { iip3_w } => {
                let k = (4.0 / 3.0 / iip3_w) as f32;
                iq.iter_mut().for_each(|c| *c *= 1.0 - k * c.norm_sqr());
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
    /// Per-source operator: antenna gain (applied in the channel contribution).
    pub antenna: Option<Antenna>,
    /// Post-summation operators (gain, ADC collapse, …) on the band buffer.
    pub post_sum: Vec<PostOp>,
}

impl RxChain {
    /// Build the chain from a band's optional `[receiver.band.frontend]` config.
    pub fn from_config(fe: Option<&FrontendConfig>) -> Self {
        let mut post_sum = Vec::new();
        let mut antenna = None;
        if let Some(fe) = fe {
            if let Some(a) = &fe.antenna {
                antenna = Some(Antenna { peak_dbi: a.peak_dbi, horizon_dbi: a.horizon_dbi, fb_db: a.fb_db, cutoff_deg: a.cutoff_deg });
            }
            if let Some(db) = fe.gain_db {
                post_sum.push(PostOp::Gain { db });
            }
            if let Some(nl) = &fe.nonlinearity {
                // dBm → watts (lipika amplitude² is power).
                let iip3_w = 10f64.powf((nl.iip3_dbm - 30.0) / 10.0);
                post_sum.push(PostOp::Nonlinearity { iip3_w });
            }
            if let Some(adc) = &fe.adc {
                post_sum.push(PostOp::Adc { bits: adc.bits, full_scale_v: adc.full_scale_v });
            }
        }
        RxChain { antenna, post_sum }
    }

    /// Forward-interpret the post-summation operators onto the band buffer, in
    /// chain order. No-op when the chain is empty.
    // r[impl eq.e10] post-sum operator chain, in order
    pub fn apply_post_sum(&self, iq: &mut Array1<Complex32>) {
        for op in &self.post_sum {
            op.apply(iq);
        }
    }

    pub fn is_identity(&self) -> bool {
        self.antenna.is_none() && self.post_sum.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{AdcConfig, FrontendConfig, NonlinearityConfig};

    // r[verify eq.e10]
    #[test]
    fn empty_chain_is_identity() {
        let chain = RxChain::from_config(None);
        assert!(chain.is_identity());
        let mut iq = Array1::from_vec(vec![Complex32::new(0.3, -0.7), Complex32::new(1.1, 0.2)]);
        let before = iq.clone();
        chain.apply_post_sum(&mut iq);
        assert_eq!(iq, before);
    }

    // r[verify eq.e9]
    #[test]
    fn adc_quantises_to_levels_and_clips() {
        // 4-bit ADC, full scale ±1 V: 7 positive levels, LSB = 1/7.
        let chain = RxChain::from_config(Some(&FrontendConfig {
            antenna: None,
            gain_db: None,
            nonlinearity: None,
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

    // r[verify eq.e6]
    #[test]
    fn antenna_pattern_peaks_at_zenith_and_rejects_back() {
        let a = Antenna { peak_dbi: 3.9, horizon_dbi: -3.0, fb_db: 8.5, cutoff_deg: None };
        let g = |e: f64| 20.0 * (a.amplitude_at(e) as f64).log10(); // back to dBi
        assert!((g(90.0) - 3.9).abs() < 1e-3, "zenith {}", g(90.0));   // peak
        assert!((g(0.0) + 3.0).abs() < 1e-3, "horizon {}", g(0.0));    // horizon
        assert!((g(-10.0) - (-3.0 - 8.5)).abs() < 1e-3);               // back-lobe
        assert!(g(45.0) > g(0.0) && g(45.0) < g(90.0));                // monotone up
        // cutoff occludes below it (Earth LOS), leaves above untouched.
        let ac = Antenna { cutoff_deg: Some(0.0), ..a };
        assert_eq!(ac.amplitude_at(-0.1), 0.0, "below cutoff occluded");
        assert!(ac.amplitude_at(0.0) > 0.0, "at cutoff present");
        assert_eq!(ac.amplitude_at(45.0), a.amplitude_at(45.0), "above cutoff unchanged");
    }

    // r[verify eq.e10]
    #[test]
    fn gain_scales_linearly() {
        let chain = RxChain::from_config(Some(&FrontendConfig {
            antenna: None,
            gain_db: Some(20.0),  // ×10
            nonlinearity: None,
            adc: None,
        }));
        let mut iq = Array1::from_vec(vec![Complex32::new(0.1, -0.2)]);
        chain.apply_post_sum(&mut iq);
        assert!((iq[0].re - 1.0).abs() < 1e-5 && (iq[0].im + 2.0).abs() < 1e-5);
    }

    // r[verify eq.e8]
    #[test]
    fn nonlinearity_is_linear_small_and_compresses_large() {
        // IIP3 = 0 dBm = 1 mW. Small signal (|x|² ≪ Pᵢₚ₃) passes ~unchanged;
        // a signal near the intercept is compressed (|y| < |x|).
        let chain = RxChain::from_config(Some(&FrontendConfig {
            antenna: None,
            gain_db: None,
            nonlinearity: Some(NonlinearityConfig { iip3_dbm: 0.0 }),
            adc: None,
        }));
        let small = Complex32::new(1e-4, 0.0);            // |x|² = 1e-8 W ≪ 1e-3
        let large = Complex32::new(0.02, 0.0);            // |x|² = 4e-4 W ~ 0.4·Pᵢₚ₃
        let mut iq = Array1::from_vec(vec![small, large]);
        chain.apply_post_sum(&mut iq);
        assert!(((iq[0] - small).norm() / small.norm()) < 1e-3, "small not linear");
        assert!(iq[1].norm() < large.norm(), "large not compressed");
        // exact coefficient: gain factor = 1 − (4/3)|x|²/Pᵢₚ₃
        let expect = large.re * (1.0 - (4.0 / 3.0) * 4e-4 / 1e-3) as f32;
        assert!((iq[1].re - expect).abs() < 1e-6);
    }
}
