//! Execution backend abstraction.
//!
//! v0.1 ships CPU-only (rayon-parallel across observers). A GPU backend
//! via wgpu/WGSL was prototyped and abandoned: WGSL has no native f64,
//! and f32 carrier-phase precision at L-band (1.5 GHz × 70 ms delay
//! ≈ 1e8 cycles) loses sub-cycle accuracy before any shader code runs.
//! The f64-via-Vulkan/CUDA alternative gives only ~3-5× over CPU on
//! consumer GPUs (1/64 f64 throttle), not worth the cross-platform loss.
//!
//! Left here for future revisit if we ever target datacenter GPUs.

pub mod cpu;

/// What kind of execution to use for a scene.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    Cpu,
}
