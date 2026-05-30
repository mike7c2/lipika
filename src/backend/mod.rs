//! Execution backend abstraction.
//!
//! v1 ships with a CPU backend (rayon-parallel across observers + sources).
//! GPU backend (wgpu / WGSL) lives behind the `gpu` feature flag and ports
//! the per-output-sample inner loop into a compute shader.

pub mod cpu;

#[cfg(feature = "gpu")]
pub mod wgpu_backend;

/// What kind of execution to use for a scene.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Cpu,
    #[cfg(feature = "gpu")]
    Gpu,
}
