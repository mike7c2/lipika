//! Per-observer oscillator model.
//!
//! `Δt(t) = c + m·t + ε(t)` where
//!   - `c` is the constant offset (residual PPS-edge jitter etc.)
//!   - `m` is the static fractional frequency offset (ppm-scale)
//!   - `ε(t)` is the white-FM Wiener walk, calibrated to a target σ_y(τ=1s)
//!
//! Real OCXO behaviour also has flicker and random-walk FM at longer τ;
//! those aren't modelled here. For seconds-to-minutes timescales (DGSN
//! pass durations), white-FM + drift captures the dominant behaviour.

use ndarray::Array1;
use rand_distr::{Distribution, Normal};
use rand_xoshiro::Xoshiro256PlusPlus;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct NodeClockSpec {
    /// White-FM Allan deviation at τ = 1 s (dimensionless). Zero → no walk.
    #[serde(default)]
    pub allan_deviation_1s: f64,
    /// Static fractional frequency offset (dimensionless; 1e-6 = 1 ppm).
    #[serde(default)]
    pub frequency_offset_y0: f64,
    /// Static initial clock bias Δt(0) in seconds.
    #[serde(default)]
    pub initial_time_offset_s: f64,
}

impl Default for NodeClockSpec {
    fn default() -> Self {
        Self {
            allan_deviation_1s: 0.0,
            frequency_offset_y0: 0.0,
            initial_time_offset_s: 0.0,
        }
    }
}

/// Materialised clock track. Build once per observer per pass, query
/// at the master sample grid or interpolate to other grids.
pub struct NodeClock {
    pub start_epoch_unix: f64,
    pub master_sample_rate_hz: f64,
    pub time_offset_track: Array1<f64>,
}

impl NodeClock {
    pub fn new(
        spec: NodeClockSpec,
        start_epoch_unix: f64,
        master_sample_rate_hz: f64,
        duration_s: f64,
        rng: &mut Xoshiro256PlusPlus,
    ) -> Self {
        let n = (master_sample_rate_hz * duration_s).round() as usize;
        let mut walk = Array1::<f64>::zeros(n);
        if spec.allan_deviation_1s > 0.0 {
            let step_std = spec.allan_deviation_1s / master_sample_rate_hz.sqrt();
            let normal = Normal::new(0.0, step_std).expect("valid σ");
            let mut acc = 0.0;
            for v in walk.iter_mut() {
                acc += normal.sample(rng);
                *v = acc;
            }
        }
        // Add deterministic offset + linear drift.
        let mut out = walk;
        for (i, v) in out.iter_mut().enumerate() {
            let elapsed = (i as f64) / master_sample_rate_hz;
            *v += spec.initial_time_offset_s + spec.frequency_offset_y0 * elapsed;
        }
        Self { start_epoch_unix, master_sample_rate_hz, time_offset_track: out }
    }

    /// Linearly interpolate Δt(t) at arbitrary UTC times.
    pub fn time_offset_at(&self, times_unix: &Array1<f64>) -> Array1<f64> {
        let fs = self.master_sample_rate_hz;
        let t0 = self.start_epoch_unix;
        let track = &self.time_offset_track;
        let n_max = track.len() as f64;
        let mut out = Array1::<f64>::zeros(times_unix.len());
        for (i, &t) in times_unix.iter().enumerate() {
            let idx_f = (t - t0) * fs;
            if idx_f <= 0.0 {
                out[i] = track[0];
            } else if idx_f >= n_max - 1.0 {
                out[i] = track[track.len() - 1];
            } else {
                let lo = idx_f.floor() as usize;
                let frac = idx_f - (lo as f64);
                out[i] = track[lo] * (1.0 - frac) + track[lo + 1] * frac;
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;

    #[test]
    fn frequency_offset_produces_exact_linear_drift() {
        let spec = NodeClockSpec {
            allan_deviation_1s: 0.0,
            frequency_offset_y0: 5e-6,
            initial_time_offset_s: 0.0,
        };
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(42);
        let clock = NodeClock::new(spec, 0.0, 1e6, 1.0, &mut rng);
        let expected_at_end = 5e-6 * (1.0 - 1.0 / 1e6);
        assert_relative_eq!(
            clock.time_offset_track[clock.time_offset_track.len() - 1],
            expected_at_end,
            max_relative = 1e-12
        );
    }

    #[test]
    fn allan_deviation_recovers_to_target() {
        let spec = NodeClockSpec {
            allan_deviation_1s: 1e-11,
            frequency_offset_y0: 0.0,
            initial_time_offset_s: 0.0,
        };
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(7);
        let fs = 1e6;
        let duration = 60.0;
        let clock = NodeClock::new(spec, 0.0, fs, duration, &mut rng);
        // Compute σ_y(1s) as std of (Δt_track[t+1s] - Δt_track[t]).
        let step = fs as usize;
        let mut deltas: Vec<f64> = Vec::new();
        let mut i = 0;
        while i + step < clock.time_offset_track.len() {
            deltas.push(clock.time_offset_track[i + step] - clock.time_offset_track[i]);
            i += step;
        }
        let mean: f64 = deltas.iter().sum::<f64>() / deltas.len() as f64;
        let var: f64 = deltas.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / deltas.len() as f64;
        let measured = var.sqrt();
        assert!(
            (measured / spec.allan_deviation_1s - 1.0).abs() < 0.4,
            "measured σ_y(1s) {} not within 40% of target {}", measured, spec.allan_deviation_1s
        );
    }
}
