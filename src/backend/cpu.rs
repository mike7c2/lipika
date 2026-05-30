//! CPU backend — rayon-parallel across (observer, source) pairs.
//!
//! Today this is the only backend that exists. Once the wgpu backend is
//! online they should produce bit-comparable IQ at f32 precision; the
//! CPU path keeps a reference implementation tests can validate against.
