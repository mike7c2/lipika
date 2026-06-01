//! `lipika` — multi-source, multi-observer IQ scene simulator.
//!
//! The library exposes the building blocks the CLI binary orchestrates:
//!
//! - [`geometry`] — ECEF math, observer locations, source trajectories
//!   (stationary, TLE/SGP4, synthetic Keplerian).
//! - [`signal`] — source signal generators (BPSK, CW, GPS C/A code, file-IQ).
//! - [`clock`] — per-observer oscillator model: constant offset, frequency
//!   offset (ppm scale), and white-FM phase walk.
//! - [`channel`] — free-space propagation (delay + Doppler via
//!   time-domain resampling, not frequency shift), AWGN.
//! - [`backend`] — execution backend: CPU (rayon-parallel across the
//!   transmitters within each time chunk). A wgpu/WGSL GPU backend was
//!   prototyped and abandoned; see `backend::mod` for the precision analysis.
//! - [`config`] — TOML scene description (serde structs).
//! - [`output`] — IQ file writers (raw little-endian `complex64` + JSON sidecar).
//!
//! ## Architectural principle: time-domain resampling for Doppler
//!
//! Many SDR simulators apply Doppler as a frequency shift (multiply IQ by
//! `exp(j·2π·f_d·t)`). That's a narrowband approximation and leaks: the
//! chip rate of a wideband signal also Doppler-shifts, and a sample-rate
//! skew between observer clocks would be a separate effect to bolt on.
//!
//! Instead, every output sample `n` at observer `o` is built by:
//! 1. Computing the observer's *local* time `t_local = n / fs_o + epoch_o`.
//! 2. Mapping it to UTC via `t_utc = t_local − Δt_o(t_local)` where
//!    `Δt_o = c + m·t_local + ε(t_local)` is the observer's clock state.
//! 3. For each source `s`, solving for the emission time `t_emit` that
//!    satisfies `t_emit + |P_s(t_emit) − P_o(t_utc)| / c = t_utc`
//!    (one or two Newton iterations from a `t_utc − 70 ms` initial guess).
//! 4. Looking up the source signal at `t_emit`.
//! 5. Mixing to baseband against the observer's local oscillator.
//!
//! Doppler, chip-rate Doppler, sample-rate skew, and code-phase drift all
//! emerge from this one mechanism, by construction self-consistent.

pub mod backend;
pub mod channel;
pub mod clock;
pub mod config;
pub mod geometry;
pub mod output;
pub mod pipeline;
pub mod profile;
pub mod signal;

/// Speed of light in m/s (IAU 2009 / WGS84 conventional value).
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;
