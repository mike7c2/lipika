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

**Milestone B — carrier-phase orbit determination (no truth in the estimator).**
A data-driven Doppler tracker (FFT-peak FLL, no truth) recovers the continuous
phase to **~6 mm**, then a two-body ECEF state `[r, v]` is fit to the multi-node
phase. Multi-seed (from ~2 km / 20 m/s offsets):

| clock handling | recovered position | behaviour |
|---|---:|---|
| absorbed as per-node **bias + rate** | 0.8–1.0 km (scattered) | ill-conditioned |
| **disciplined** (clock removed, bias only) | **246 ± 1 m** (radial 22 m, transverse 245 m), vel 3.2 m/s | **unimodal** |

Two lessons:

- **Motion makes the OD unimodal — the ambiguity is gone.** With the clock
  disciplined, all seeds converge to the *same* orbit (vs Exp 4's λ-lattice of
  equal sidelobes). Continuous phase + 221 km of geometry change turn the
  multimodal static problem into a single basin. *No integers are fixed.*
- **Clock discipline is essential.** Absorbing the per-node clock as a *rate*
  also eats the differential Doppler (FDOA = the velocity/angular signal), and
  the fit degenerates to ~km. Separating clock drift from FDOA — the job of the
  GNSS time reference — is what makes the orbit observable. (Here the clock is
  removed exactly; realistic GNSS discipline must estimate the drift well enough
  to avoid this degeneracy.)

The residual ~246 m is **model-limited**: the fit RMS (9 mcyc) sits above the
~6 mm phase noise because the two-body arc doesn't perfectly match lipika's SGP4
over the window — radial is tight (22 m), the transverse model error dominates.
So, as expected, the *observable* is mm but the *orbit estimate* is hundreds of
metres; tightening it is a force-model / longer-arc / wider-array question, not
an ambiguity one. Refinements: full force model (or shorter arc) for the model
mismatch; carrier-phase GNSS time transfer for tighter clock separation.

## Exp 6 — owning the clock: closing the realistic loop

Exps 3–5 disciplined the node clocks with a *truth* shortcut (the truth-export
delay, or the configured clock removed exactly). Exp 6 replaces that with a real
estimator that recovers each node's clock from the **GPS L1 IQ alone** — built
in the separate `gnss-rx` crate (Rust + pyo3, reusing lipika's tested IS-GPS-200
primitives), driven from `experiments/gps_clock.py`.

**The receiver (truth-free):** acquisition (parallel code/Doppler search) →
block tracking → a **known-position clock solve** (predict each PRN's code phase
from the broadcast ephemeris + surveyed position, the common code-phase residual
is the clock). **Cross-validated against GNSS-SDR** on the same IQ: our recovered
drift agrees with GNSS-SDR's `user_clk_offset` to ~0.05 ns/s, and the configured
offset to 6.5 ns — two independent estimators, neither touching lipika truth.

**The catch, and the fix — the carrier.** Feeding the *code* clock into the
Exp 5 OD only reaches ~4 km: the conditioning-sensitive carrier-phase OD needs
the clock **drift** known to far better than code gives over a 30 s arc. So
gnss-rx also recovers the clock from the GPS **carrier**: de-chirp each PRN by
the predicted geometry (code aligned with the code clock), square out the NAV,
fit the slope → **drift to ~1×10⁻³ ns/s** (matching every node's configured
drift to 3 digits), ~100× tighter than code.

**The loop closes.** With the truth-free code-offset + carrier-drift clock:

| clock fed to the OD | orbit error |
|---|---:|
| none | 23 km |
| code-only (drift ~0.3 ns/s) | 4.3 km |
| **code offset + carrier drift (gnss-rx)** | **243 m** |
| true clock (Exp 5 ideal) | 246 m |

The realistic, fully-owned, independently-validated clock reaches the **246 m
model-limited ideal** — the complete chain (own GPS receiver → carrier-phase
clock transfer → carrier-phase orbit determination) closes with no truth in any
estimator. The binding requirement is now quantified: distributed carrier-phase
coherence needs **carrier-grade** clock transfer (drift ≪ 0.3 ns/s), which
code-based GNSS over a short arc cannot give but the carrier can.
