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

## Exp 3 — GNSS-disciplined timing rescues TDOA

The keystone result. Each node from Exp 1 additionally carries a **GNSS L1 band
sharing its clock** (1575.42 MHz, 8.192 Msps, 6 GPS SVs with full LNAV from
`gps_synth.rnx`), alongside the VHF band. The clock offset measured on the GNSS
band is — by construction in lipika — the same offset corrupting the VHF TDOA.
We recover each node's clock by code-phase time transfer against the
truth-export geometry (`gps.py`), then de-bias the VHF TDOA with it.

| | horizontal | notes |
|-|-----------:|-------|
| raw (unsynchronised)   | **76 m** | clocks uncorrected |
| GNSS-disciplined       | **17 m** | clock recovered from GPS, removed |

Per-node clock recovery lands within **~4 ns** of the configured offsets
(e.g. ayr 342.4 vs 341.0 ns, stirling −214.6 vs −210.2 ns). This is the full
distributed-timing loop: GNSS gives each node a common time reference, and that
reference makes coherent cross-node geolocation work.

### The bug this surfaced — 238 ns timing quantisation

Exp 3 only produced clean numbers *after* fixing a fundamental precision bug it
exposed. lipika had been doing per-sample timing arithmetic on **absolute Unix
times** (~1.717×10⁹ s); the f64 ULP there is exactly 2⁻²² s = **238.4 ns**, so
every sub-238 ns quantity — clock offsets, the Allan walk, the Sagnac term,
sub-sample delay — was rounded away. Measured directly: configured offsets of
25/50 ns read back as 0; 130/250 ns as 238 ns; 500 ns as 477 ns — snapping to
the ULP. (GPS still validated because its delay is range-derived and the
sample-time jitter averages out in tracking.)

Fixed by moving the precise per-sample maths to an **epoch-relative timebase**
(seconds since the scene epoch → femtosecond ULP), reconstructing absolute time
only to query trajectories (where 238 ns ≈ 0.2 mm of orbit motion). The lever is
magnitude, not bit width — i128/f128 wouldn't help, since precision dies at the
f64 boundary unless the value crossing it is small. Guarded by
`clock_offset_survives_at_realistic_epoch` (a 130 ns offset must survive at a
1.7×10⁹ epoch), and the chunk/thread invariance tests still pass because the
relative time keys on the global sample index.

## Exp 4 — VHF carrier-phase interferometry

Moving from envelope timing to the **carrier phase** as the cross-node
observable (−2π·Δr/λ; at VHF λ≈2.2 m, ~1000× finer than Exp 1's TDOA). The
target emits a co-located CW tone (140 MHz) for phase plus a BPSK (137 MHz) for
a coarse prior; each node's clock is removed using the GNSS time transfer from
Exp 3. A position-domain coherent search (distributed beamformer)
`J(x)=|Σᵢ exp(j(ψᵢ + 2π f·rᵢ(x)/c))|` then locates the emitter — with the
Sagnac term included to match lipika's geometry.

Two findings:

- **lipika's carrier phase is faithful.** Clock-corrected node phases match the
  truth-export geometry to **±0.3° (0.01 ns)**, the coherence at the true
  position is **1.0000**, and a search within ±1 λ of truth geolocates to
  **2 mm**. This is the headline simulator validation — the interferometric
  observable is correct to the millimetre.
- **Static single-frequency geolocation is ambiguity-limited.** A metre-class
  envelope prior is ~8 wavelengths, and the 6-node coherence has near-unity
  sidelobes (~0.998), so a global-max search over the prior box lands on the
  wrong lobe (tens of metres). The phase is perfect; the *inverse problem* is
  multimodal. Resolving it needs diversity that breaks the λ-periodicity —
  multiple frequencies, **more geometry over time (a moving emitter)**, or a
  sub-λ prior. SNR doesn't help (sidelobes are structural). This sets up Exp 5.

## Exp 5 — carrier-phase tracking of a moving (satellite) emitter

The answer to Exp 4's ambiguity: make the emitter move. The target rides the
**ISS TLE** (real two-line set, at the GPS-valid epoch) — a LEO at ~408 km doing
7.36 km/s, moving **221 km across the sky in 30 s**, visible at 12–19° from all
six nodes. Over that arc the range to each node changes by ~175 km ≈ 80 000 VHF
wavelengths. (lipika's truth-export now also records the exact **source
position**, since range-only reconstruction of a distant target from a tight
array is radially ill-conditioned.)

**Milestone A — the moving-emitter carrier phase is faithful.** Tracking the CW
phase continuously and de-rotating by the truth range leaves a residual that is
*purely the receiver clock* −2π·f·Δt(t); the scatter about it is **3–5 mm RMS
(2–5 millicycles)** at every node. So lipika's Doppler/emission-time path
reproduces a fast mover's carrier phase to the millimetre — the moving-emitter
analogue of Exp 4's static 2 mm check.

**This resolves the ambiguity question.** The per-node clock offset comes back
*wrapped mod the 7.28 ns carrier period* (e.g. configured 130 ns → −1.09 ns).
That is the whole point: continuous phase tracking collapses the ambiguity from
one-per-epoch down to a **single constant bias per node**. You never fix
integers; the 221 km of geometry change over the pass is the diversity that
makes the trajectory observable, with that one per-node bias as a nuisance.

**Next (Milestone B, scoped):** turn the tracked phase into the satellite's
orbit. A tight ground array has poor *radial* observability instantaneously, so
this needs a data-driven Doppler tracking loop (an FLL, no truth) plus a
**Keplerian-dynamics-constrained** fit — the orbit model ties range to the
precisely-measured angular track and Doppler, which a free position/velocity arc
would leave at ~km. That is genuinely its own orbit-determination experiment.
