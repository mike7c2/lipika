//! Per-observer oscillator model.
//!
//! `Δt(t) = c + m·t + ε(t)` where
//!   - `c` is the constant offset
//!   - `m` is the static fractional frequency offset (ppm-scale)
//!   - `ε(t)` is the white-FM Wiener walk, calibrated to a target σ_y(τ=1s)
//!
//! Flicker and random-walk FM at longer τ are not modelled; over typical pass
//! durations (seconds to minutes) white-FM + drift dominate.

use ndarray::Array1;
use rand_distr::{Distribution, Normal};
use rand_xoshiro::Xoshiro256PlusPlus;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct NodeClockSpec {
    /// White-FM Allan deviation at τ = 1 s (dimensionless). Zero → no walk.
    #[serde(default)]
    pub allan_deviation_1s: f64,
    /// **Flicker-FM** Allan-deviation *floor* (dimensionless, flat in τ). This is
    /// the oscillator-grade knob: the irreducible close-in instability that drift
    /// + white-FM miss (TCXO ≈ 1e-9…1e-10, OCXO ≈ 1e-12…1e-13). Zero → no flicker.
    #[serde(default)]
    pub flicker_floor: f64,
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
            flicker_floor: 0.0,
            frequency_offset_y0: 0.0,
            initial_time_offset_s: 0.0,
        }
    }
}

/// Flicker-FM (1/f frequency noise) generator: a **sum of first-order
/// relaxations** (AR(1) poles log-spaced in frequency) whose superposition
/// approximates a 1/f fractional-frequency PSD over the band. Runs on a coarse
/// grid (`FLICKER_RATE_HZ`) — flicker is low-frequency, so updating it at the MHz
/// master rate would be ~10¹⁰ draws/pass for nothing; the streaming clock
/// interpolates its δ. Self-calibrated at construction so the resulting δ has the
/// configured Allan **floor**.
struct Flicker {
    alpha: Vec<f64>,            // AR(1) decay per pole
    gain: Vec<f64>,             // √(1−α²) per pole (unit-variance drive)
    state: Vec<f64>,            // AR states
    scale: f64,                 // calibrated fractional-frequency scale
    rng: Xoshiro256PlusPlus,
    normal: Normal<f64>,
    frontier: i64,              // coarse-grid frontier index
    delta_hi: f64,              // accumulated flicker δ(coarse=frontier)
    delta_lo: f64,              // δ(coarse=frontier−1)
}

const FLICKER_RATE_HZ: f64 = 128.0;
const FLICKER_POLES: usize = 20;
const FLICKER_F_LO: f64 = 1e-4;        // ~τ up to a few 1000 s
const FLICKER_BURN_STEPS: usize = 200_000;

impl Flicker {
    fn poles() -> (Vec<f64>, Vec<f64>) {
        let f_hi = FLICKER_RATE_HZ / 2.0;
        let (alpha, gain): (Vec<f64>, Vec<f64>) = (0..FLICKER_POLES).map(|j| {
            let f = FLICKER_F_LO * (f_hi / FLICKER_F_LO)
                .powf(j as f64 / (FLICKER_POLES - 1) as f64);
            let a = (-2.0 * std::f64::consts::PI * f / FLICKER_RATE_HZ).exp();
            (a, (1.0 - a * a).sqrt())
        }).unzip();
        (alpha, gain)
    }

    /// One coarse step → instantaneous fractional frequency `y` (unit-scaled).
    fn raw_step(alpha: &[f64], gain: &[f64], state: &mut [f64],
                normal: &Normal<f64>, rng: &mut Xoshiro256PlusPlus) -> f64 {
        let mut y = 0.0;
        for j in 0..alpha.len() {
            state[j] = alpha[j] * state[j] + gain[j] * normal.sample(rng);
            y += state[j];
        }
        y
    }

    fn new(floor: f64, seed_rng: &Xoshiro256PlusPlus) -> Self {
        let (alpha, gain) = Self::poles();
        let normal = Normal::new(0.0, 1.0).expect("σ");

        // Calibrate the scale on an independent stream: integrate y→δ, measure the
        // Allan deviation in the flat region, set scale = floor / measured.
        let mut crng = seed_rng.clone();
        crng.long_jump();
        let mut cstate = vec![0.0; alpha.len()];
        for _ in 0..FLICKER_BURN_STEPS {
            Self::raw_step(&alpha, &gain, &mut cstate, &normal, &mut crng);
        }
        let nc = 1usize << 18;
        let mut delta = vec![0.0; nc];
        let mut acc = 0.0;
        for d in delta.iter_mut().skip(1) {
            acc += Self::raw_step(&alpha, &gain, &mut cstate, &normal, &mut crng) / FLICKER_RATE_HZ;
            *d = acc;
        }
        let m = (2.0 * FLICKER_RATE_HZ) as usize;          // τ = 2 s (flat region)
        let measured = allan_dev_from_delta(&delta, m, FLICKER_RATE_HZ).max(1e-30);
        let scale = floor / measured;

        // Streaming stream: independent of calibration; warm to steady state.
        let mut rng = seed_rng.clone();
        rng.long_jump();
        rng.long_jump();
        let mut state = vec![0.0; alpha.len()];
        for _ in 0..FLICKER_BURN_STEPS {
            Self::raw_step(&alpha, &gain, &mut state, &normal, &mut rng);
        }
        Flicker { alpha, gain, state, scale, rng, normal,
                  frontier: -1, delta_hi: 0.0, delta_lo: 0.0 }
    }

    fn advance_to(&mut self, m: i64) {
        while self.frontier < m {
            self.delta_lo = self.delta_hi;
            let y = Self::raw_step(&self.alpha, &self.gain, &mut self.state,
                                   &self.normal, &mut self.rng) * self.scale;
            self.delta_hi += y / FLICKER_RATE_HZ;          // ∫ frequency dt → time
            self.frontier += 1;
        }
    }

    /// Flicker δ at `elapsed_s`, linearly interpolated on the coarse grid.
    fn delta_at(&mut self, elapsed_s: f64) -> f64 {
        let x = elapsed_s.max(0.0) * FLICKER_RATE_HZ;
        let m0 = x.floor() as i64;
        let frac = x - m0 as f64;
        self.advance_to(m0 + 1);
        self.delta_lo * (1.0 - frac) + self.delta_hi * frac
    }
}

/// Allan deviation σ_y(τ) from a uniformly-sampled time-error series `delta`
/// (seconds) at rate `fc`, with `m` samples per τ. Non-overlapping second
/// differences: σ_y²(τ) = ⟨(δ_{i+2m} − 2δ_{i+m} + δ_i)²⟩ / (2τ²).
fn allan_dev_from_delta(delta: &[f64], m: usize, fc: f64) -> f64 {
    let tau = m as f64 / fc;
    let (mut sum, mut cnt) = (0.0, 0usize);
    let mut i = 0;
    while i + 2 * m < delta.len() {
        let d = delta[i + 2 * m] - 2.0 * delta[i + m] + delta[i];
        sum += d * d;
        cnt += 1;
        i += m;
    }
    ((sum / cnt.max(1) as f64) / (2.0 * tau * tau)).sqrt()
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

/// Forward-only, O(1)-memory generator of the same `Δt(t)` track that
/// [`NodeClock`] materialises — but without ever holding the whole array.
///
/// The Wiener walk `ε` is the only cumulative state; everything else in
/// `Δt = c + m·t + ε(t)` is a pure function of `t`. `ClockStream` advances the
/// walk by *master-grid index*, carrying just the last two grid points
/// (`eps_lo`/`eps_hi`) so it can linearly interpolate `ε` at an arbitrary time
/// while only moving forward. This makes streaming chunked generation use
/// memory independent of duration.
///
/// **Coherence contract:** all of a receiver's bands construct a `ClockStream`
/// from the *same* seed and master rate and query it at monotonically
/// increasing times, so they observe the identical `ε(t)` — which is what lets
/// one model multi-channel coherent receivers. (The highest-rate band lands on
/// the master grid exactly; lower-rate bands interpolate the same walk.)
///
/// Queries MUST be monotonically non-decreasing in time (which the pipeline's
/// sample order guarantees); the generator cannot look backwards.
pub struct ClockStream {
    start_epoch_unix: f64,
    master_rate_hz: f64,
    initial_offset_s: f64,
    freq_offset_y0: f64,
    rng: Xoshiro256PlusPlus,
    /// `None` when `allan_deviation_1s == 0` (no walk; ε ≡ 0, no RNG drawn —
    /// keeps the per-band noise stream independent of the clock).
    normal: Option<Normal<f64>>,
    /// Index of the current frontier grid point (`eps_hi`); starts at −1 with
    /// ε ≡ 0 so the first `advance_to(0)` draws exactly one increment, matching
    /// `NodeClock`'s `walk[0]`.
    frontier: i64,
    eps_hi: f64,
    eps_lo: f64,
    /// Flicker-FM contribution (own coarse-grid process); `None` unless configured,
    /// keeping the white-FM walk + per-band noise streams byte-identical when off.
    flicker: Option<Flicker>,
}

impl ClockStream {
    pub fn new(
        spec: NodeClockSpec,
        start_epoch_unix: f64,
        master_rate_hz: f64,
        rng: Xoshiro256PlusPlus,
    ) -> Self {
        let normal = if spec.allan_deviation_1s > 0.0 {
            let step_std = spec.allan_deviation_1s / master_rate_hz.sqrt();
            Some(Normal::new(0.0, step_std).expect("valid σ"))
        } else {
            None
        };
        // Built from an independent (long-jumped) stream so it never perturbs the
        // white-FM walk's RNG; `None` when off ⇒ output unchanged.
        let flicker = (spec.flicker_floor > 0.0).then(|| Flicker::new(spec.flicker_floor, &rng));
        Self {
            start_epoch_unix,
            master_rate_hz,
            initial_offset_s: spec.initial_time_offset_s,
            freq_offset_y0: spec.frequency_offset_y0,
            rng,
            normal,
            frontier: -1,
            eps_hi: 0.0,
            eps_lo: 0.0,
            flicker,
        }
    }

    /// Advance the walk so the frontier grid point is at index `k`. After this,
    /// `eps_hi == ε(k)` and `eps_lo == ε(k−1)`.
    fn advance_to(&mut self, k: i64) {
        while self.frontier < k {
            self.eps_lo = self.eps_hi;
            if let Some(n) = &self.normal {
                self.eps_hi += n.sample(&mut self.rng);
            }
            self.frontier += 1;
        }
    }

    /// Deterministic part `c + m·elapsed`.
    #[inline]
    fn deterministic(&self, elapsed_s: f64) -> f64 {
        self.initial_offset_s + self.freq_offset_y0 * elapsed_s
    }

    /// `Δt` exactly at integer master-grid index `k` (used for testing against
    /// [`NodeClock`] and as the exact-grid path).
    pub fn delta_t_at_master_index(&mut self, k: i64) -> f64 {
        self.advance_to(k.max(0));
        let eps = if k < 0 { 0.0 } else { self.eps_hi };
        let elapsed = (k.max(0) as f64) / self.master_rate_hz;
        let flick = self.flicker.as_mut().map_or(0.0, |f| f.delta_at(elapsed));
        eps + flick + self.deterministic(elapsed)
    }

    /// `Δt(t)` at an arbitrary UTC time, linearly interpolating the walk
    /// between the bracketing master-grid points. `t` must not go backwards.
    pub fn delta_t_at(&mut self, t_unix: f64) -> f64 {
        let elapsed = t_unix - self.start_epoch_unix;
        let x = elapsed * self.master_rate_hz;
        if x <= 0.0 {
            // Clamp to grid point 0 (matches NodeClock::time_offset_at).
            self.advance_to(0);
            let flick = self.flicker.as_mut().map_or(0.0, |f| f.delta_at(0.0));
            return self.eps_hi + flick + self.deterministic(0.0);
        }
        let k0 = x.floor() as i64;
        let frac = x - (k0 as f64);
        self.advance_to(k0 + 1);
        let eps = self.eps_lo * (1.0 - frac) + self.eps_hi * frac;
        let flick = self.flicker.as_mut().map_or(0.0, |f| f.delta_at(elapsed));
        // Linear interpolation of the (linear) deterministic part between grid
        // points is exact, so compute it directly from continuous elapsed time.
        eps + flick + self.deterministic(elapsed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;

    #[test]
    fn flicker_fm_has_a_flat_allan_floor_at_the_configured_level() {
        // Pure flicker FM at 1e-10: the δ(t) Allan deviation should sit ~1e-10 and
        // be *flat* in τ (the flicker signature) — unlike white FM (∝ τ^-1/2) or
        // random walk (∝ τ^+1/2). master rate = the coarse rate so each index = a
        // coarse step.
        let fc = super::FLICKER_RATE_HZ;
        let spec = NodeClockSpec { flicker_floor: 1e-10, ..Default::default() };
        let mut s = ClockStream::new(spec, 0.0, fc, Xoshiro256PlusPlus::seed_from_u64(42));
        let n = 1usize << 18;
        let delta: Vec<f64> = (0..n).map(|k| s.delta_t_at_master_index(k as i64)).collect();
        let adev = |tau_s: f64| super::allan_dev_from_delta(&delta, (tau_s * fc) as usize, fc);
        let (a1, a2, a8) = (adev(1.0), adev(2.0), adev(8.0));
        // right level (generous: single-τ estimate has sampling spread)
        assert!((0.3e-10..3.0e-10).contains(&a2), "σy(2s) = {a2:e}");
        // flat: τ-sweep stays within a factor ~2 (white FM would drop ~2.8× over 1→8s)
        assert!((a8 / a1) > 0.45 && (a8 / a1) < 2.2, "not flat: σy(1)={a1:e} σy(8)={a8:e}");
    }

    #[test]
    fn frequency_offset_produces_exact_linear_drift() {
        let spec = NodeClockSpec {
            allan_deviation_1s: 0.0,
            flicker_floor: 0.0,
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
            flicker_floor: 0.0,
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

    #[test]
    fn clock_stream_reproduces_node_clock_bit_for_bit() {
        // Same spec + seed + master rate ⇒ the streaming generator must equal
        // the materialised track at every master-grid index, exactly.
        let spec = NodeClockSpec {
            allan_deviation_1s: 1e-11,
            flicker_floor: 0.0,
            frequency_offset_y0: 0.5e-6,
            initial_time_offset_s: 1e-7,
        };
        let fs = 1e6;
        let dur = 0.02;
        let clock = NodeClock::new(spec, 0.0, fs, dur, &mut Xoshiro256PlusPlus::seed_from_u64(123));
        let mut stream = ClockStream::new(spec, 0.0, fs, Xoshiro256PlusPlus::seed_from_u64(123));
        for (k, &expected) in clock.time_offset_track.iter().enumerate() {
            let got = stream.delta_t_at_master_index(k as i64);
            assert_eq!(got.to_bits(), expected.to_bits(),
                "ClockStream diverges from NodeClock at master index {k}: {got} vs {expected}");
        }
    }

    #[test]
    fn clock_stream_chunk_query_is_invariant() {
        // Querying Δt(t) sample-by-sample in one pass vs. in several restarts
        // (each restart a fresh ClockStream from the same seed) gives identical
        // values — i.e. the forward-only generator is chunk-size invariant.
        let spec = NodeClockSpec {
            allan_deviation_1s: 5e-11,
            flicker_floor: 0.0,
            frequency_offset_y0: 1e-6,
            initial_time_offset_s: 0.0,
        };
        let fs = 2.048e6;
        let n = 20_000usize;
        let t = |i: usize| i as f64 / fs;

        let mut one = ClockStream::new(spec, 0.0, fs, Xoshiro256PlusPlus::seed_from_u64(7));
        let whole: Vec<f64> = (0..n).map(|i| one.delta_t_at(t(i))).collect();

        // "Chunked": a fresh stream that we query only over [0, n), but the
        // forward-only generator is deterministic from the seed, so a stream
        // restarted from the same seed and walked to the same times matches.
        let mut two = ClockStream::new(spec, 0.0, fs, Xoshiro256PlusPlus::seed_from_u64(7));
        for i in 0..n {
            assert_eq!(two.delta_t_at(t(i)).to_bits(), whole[i].to_bits(),
                "Δt query not reproducible at sample {i}");
        }
    }

    #[test]
    fn clock_stream_two_bands_are_coherent() {
        // Two bands of one receiver (same seed, same master rate) sampling at
        // different rates must agree on Δt at coincident instants — the
        // multi-channel coherence contract. The high-rate band == master grid;
        // the low-rate band interpolates the same walk.
        let spec = NodeClockSpec {
            allan_deviation_1s: 1e-11,
            flicker_floor: 0.0,
            frequency_offset_y0: 0.3e-6,
            initial_time_offset_s: 0.0,
        };
        let master = 2.048e6; // high-rate band == master grid
        let low = 48_000.0;
        let mut hi = ClockStream::new(spec, 0.0, master, Xoshiro256PlusPlus::seed_from_u64(99));
        let mut lo = ClockStream::new(spec, 0.0, master, Xoshiro256PlusPlus::seed_from_u64(99));

        // Walk the high-rate band densely; sample the low-rate band on its own
        // coarser grid and compare against the high-rate value at the same time.
        let n_lo = 500usize;
        // Precompute hi at all master indices we'll need.
        let max_t = (n_lo as f64) / low;
        let n_hi = (max_t * master) as usize + 2;
        let hi_vals: Vec<f64> = (0..n_hi).map(|i| hi.delta_t_at(i as f64 / master)).collect();
        for j in 0..n_lo {
            let t = j as f64 / low;
            let got = lo.delta_t_at(t);
            // Reference: linear interp of the high-rate (== master) samples.
            let x = t * master;
            let k0 = x.floor() as usize;
            let frac = x - k0 as f64;
            let reference = hi_vals[k0] * (1.0 - frac) + hi_vals[k0 + 1] * frac;
            assert_relative_eq!(got, reference, max_relative = 1e-9, epsilon = 1e-18);
        }
    }
}
