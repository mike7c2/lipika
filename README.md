# dgsn-iq-gen

**A multi-source, multi-observer IQ generator.** That's all it does:
takes a TOML scene description (sources with trajectories + signal models;
observers with trajectories + clock models + RF frontends; channel) and
writes per-observer baseband IQ as NumPy `.npy` files. Receiver-side work
— acquisition, tracking, clock recovery, positioning — is intentionally
out of scope and belongs in downstream analysis (Python, GNU Radio, real
GNSS receivers, etc.).

Designed for DGSN (Distributed Ground Station Network) but useful for any
multi-receiver SDR work where you want Doppler, code-rate Doppler, and
inter-observer sample-rate skew handled coherently in one mechanism.

## Why another scene simulator?

GPS-SDR-SIM and friends solve part of this for GPS L1 only. We need:

- **Multiple bands per observer** with one coherent local oscillator
- **Realistic per-observer clock model** (constant offset + ppm-scale
  frequency offset + white-FM jitter, plus future flicker/RW)
- **Arbitrary source signals** (BPSK, CW, C/A code, file-IQ)
- **TLE trajectories alongside synthetic Keplerian** for both sources and
  observers
- **Time-domain Doppler** (resampling the source at delay-varying instants)
  so chip-rate effects and inter-observer sample-rate skew fall out
  automatically rather than being separate frequency shifts to bolt on

## Status

v0.1 — CPU-only library + CLI. Targets the exp/009-equivalent scene as the
first thing to reproduce. A wgpu/WGSL GPU backend was prototyped and
abandoned because WGSL has no f64 and carrier-phase math at L-band
(1.5 GHz × 70 ms = ~1e8 cycles) loses sub-cycle accuracy in f32 before
any shader code runs. See `src/backend/mod.rs` for the analysis. The
CPU path is already ~4× faster than the Python pipeline with rayon
across observers.

## Build

```bash
cd dgsn-iq-gen
cargo build --release
```

## Use

```bash
cargo run --release --bin dgsn-iq-gen -- \
    --config examples/iss_uhf_pass.toml \
    --output ./out/ \
    --seed 2026
```

This writes `<observer>.iq.npy` (complex64) + `<observer>.iq.json`
(metadata) per observer into `./out/`.

## Architectural principle: time-domain resampling for Doppler

Every output sample `n` at observer `o` is built by:

1. `t_local = epoch_o + n / fs_o` — observer's nominal sample time.
2. `t_utc = t_local − Δt_o(t_local)` — clock-correct to UTC.
3. For each source `s`, solve `t_emit + |P_s(t_emit) − P_o(t_utc)|/c = t_utc`
   (Newton iteration, 2–3 steps from a `t_utc − 70 ms` initial guess).
4. Look up source baseband at `t_emit`.
5. Mix to observer baseband against `f_lo_obs`.

Doppler, chip-rate Doppler, code-phase walk, sample-rate skew, and clock-
noise carrier wander all emerge from this one mechanism, by construction
self-consistent.

## License

MIT.
