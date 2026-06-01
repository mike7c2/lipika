# lipika experiments

Analysis experiments that consume lipika-generated IQ. Each experiment builds a
scene, runs the release `lipika` binary, then processes the raw `complex64`
output. Because lipika's geometry is exact, **ground truth is known**, so what
these measure is the *estimator's* error, not the data's.

Large IQ is written to the mounted drive
(`/media/.../lipika-test/exp_*`), not the repo.

## Running

```bash
cargo build --release           # from the repo root, once
cd experiments
python3 exp1_tdoa.py             # generates data + prints results
```

`common.py` holds the shared kit: WGS84 ⇄ ECEF, raw-IQ + sidecar loading,
FFT cross-correlation TDOA (sub-sample peak), Levenberg–Marquardt TDOA
multilateration, and an ENU error decomposition. Stationary scenes use
`stationary_ecef` nodes whose coordinates are computed in `common.py` and fed
to lipika verbatim, so the truth positions and the positions baked into the
signal are identical by construction.

## Exp 1 — TDOA multilateration of a fixed emitter

Six receivers spread across central Scotland (37–85 km baselines) observe one
wideband BPSK emitter (137 MHz, 2 Mcps, 1 s, sampled at 8 MHz). Pairwise FFT
cross-correlation of the baseband IQ gives time-difference-of-arrival; LM
multilateration solves the emitter position. Two scenes:

| scene | horizontal error | vertical error | fit residual |
|-------|-----------------:|---------------:|-------------:|
| ideal clocks            | **32.6 m** | 19 km (GDOP) | 21 m |
| unsynchronised clocks   | **78.9 m** | 20 km (GDOP) | 23 m |

Findings:

- **Horizontal geolocation works to ~30 m** with perfect clocks. The floor is
  the correlation resolution: at 8 MHz one sample is 37.5 m, and sub-sample
  interpolation reaches roughly one sample. Wider bandwidth / higher sample
  rate would tighten it — a natural follow-on sweep.
- **Vertical is GDOP-limited (~19 km).** A regional ground array has almost no
  vertical leverage on a 10 km-altitude emitter; height is essentially
  unobservable. This is geometry, not a data defect — the per-pair TDOAs match
  truth to 15–40 m. A wider-baseline or altitude-diverse array, or an
  altitude-constrained solve, is the fix.
- **Unsynchronised clocks degrade the fix to ~79 m.** Per-node offsets of tens
  to hundreds of ns map to ±40–100 m of range-difference error that propagates
  straight into the position. This is the distributed-network synchronisation
  problem, quantified — and it motivates GNSS-disciplined timing.

### Bug found

This experiment surfaced a latent lipika bug: the BPSK signal epoch was
hardcoded to `0.0` in `build_signal`, so for any real Unix-epoch scene every
emission time fell out of the bit-stream range and BPSK emitted **pure zeros**.
Fixed to use the scene epoch; covered by the
`bpsk_radiates_at_realistic_epoch` regression test.
