# Lipika requirements

This is the machine-checked requirements specification for Lipika, in
[mantra](https://github.com/mhatzl/mantra) Markdown form: each requirement is a
heading `` `id`: Title `` followed by its text. Code claims to satisfy a
requirement with `#[mantra_macros::req_satisfied("id")]` and a test claims to
verify one with `#[mantra_macros::req_test("id")]` (or `assert_req!`). mantra
collects the traces together with the test results and line coverage of the
test run, so a requirement is only **Verified** when a test tracing to it
passed *and* executed the code tracing to it. See
[TRACEABILITY.md](TRACEABILITY.md) for how to run it.

IDs are stable and mirror the public names used elsewhere in the repository:
`fid-N` is the README's FID-N, `fn-N` is FN-N, `nf-N` is NF-N, and `eq-eN` is
equation (EN) of [SIGNAL_MODEL.md](../SIGNAL_MODEL.md). Sub-requirements extend
the parent ID with a hyphenated segment (`fn-4-bpsk`) and name the parent
explicitly with `- **Parents:**`; mantra rolls children up in its hierarchy. `nf-1` and
`nf-4` are flagged `- **Manual:** true`: they are verified by a review record
(`reviews/`) rather than by an automated test.

## Fidelity principles

### `fid-1`: Geometric, time-domain propagation

All range and Doppler effects MUST arise from the emission-time (light-time) solve applied to an algebraic transmitter signal `s(t)`. The chain MUST NOT contain an explicit "apply a Doppler shift" step: carrier Doppler, chip-rate Doppler, sample-rate skew and code-phase walk all emerge from evaluating `s` at the solved emission time. Cross-checked externally by the Doppler GNSS-SDR measures on the generated IQ.

### `fid-2`: No accumulating phase error

The carrier phase of every output sample MUST be computed in f64 from the absolute geometry (range) and clock offset at that sample. Phase MUST NOT be accumulated from the previous sample, so the phase error is bounded and does not grow with run length. Per-sample timing arithmetic MUST be done in epoch-relative time so that sub-nanosecond offsets survive f64 rounding at realistic Unix epochs.

### `fid-3`: No biasing trajectory error

Trajectory interpolation MUST stay within a stated error budget (sub-micrometre over a LEO or GPS pass on the default 10 ms grid) and MUST NOT introduce a secular bias relative to the underlying propagator.

### `fid-4`: Physical, signal-independent noise

A band's noise MUST be thermal, `P_noise = k·T·B` with `B` the sample rate and `T` the configured noise temperature, fixed at the receiver and independent of the arriving signals, so that SNR and C/N0 emerge from transmit power against free-space path loss. A scene MAY pin the floor explicitly with `noise_power_w` for testing.

### `fid-5`: Physical, decomposable, coherent clock

Each receiver's clock offset MUST be `δ(t) = c + m·t + ε_wFM(t) + ε_fFM(t)`: a constant offset, a static fractional-frequency drift, a white-FM Wiener walk calibrated to an Allan deviation at τ = 1 s, and a flicker-FM floor flat in τ, each independently configurable and zero when unset. One clock realisation MUST drive all of a receiver's bands so their carrier phases are mutually coherent (a single oscillator).

## Functional requirements

### `fn-1`: Linear accumulation

A band's IQ MUST be the linear sum of each contributing transmitter's contribution plus exactly one noise realisation for the band. Noise MUST NOT be re-realised per transmitter.

### `fn-2`: Multi-band receivers off a shared clock

A receiver MUST support one or more RF bands, and every band of a receiver MUST be sampled against the same clock realisation, so the clock-induced phase on each band is the same `δ(t)` scaled by that band's frequency.

### `fn-3`: Spectral band-overlap gating

A transmitter MUST contribute to a band only if its occupied spectrum `[carrier ± bandwidth/2]` intersects the band's acceptance window `[lo ± bandwidth/2]`; the test MUST be inclusive at the boundary, and a zero-bandwidth (CW) tone MUST still be admitted by a band that straddles it. Partial overlap is admitted unfiltered (a documented limitation, see the README roadmap).

### `fn-4`: Algebraic signal models

Every signal model MUST be a stateless function `s(t) → ℂ` evaluable at any emission time (never a pre-sampled buffer) and MUST expose its nominal RF carrier frequency; it MUST NOT apply the carrier itself.

### `fn-4-cw`: CW model

- **Parents:** ["fn-4"]

The `cw` kind MUST emit a single tone `exp(j·2π·f_offset·t)` at a configurable baseband offset.

### `fn-4-bpsk`: BPSK model

- **Parents:** ["fn-4"]

The `bpsk` kind MUST emit `bits[⌊(t − epoch)·bit_rate⌋]` where the bit stream is a configured sync word followed by a seeded pseudo-random payload, with the bit index taken relative to the scene epoch.

### `fn-4-ca-code`: GPS L1 C/A model

- **Parents:** ["fn-4"]

The `ca_code` kind MUST emit the 1023-chip Gold code of the configured PRN per IS-GPS-200 §3.3.2.1 at 1.023 Mcps, with a navigation-data overlay selectable as `off` (pure code) or `dummy` (an `0x8B` TLM preamble at every 6 s subframe start followed by pseudo-random bits, each 50 bps data bit spanning 20 code periods).

### `fn-4-lnav`: Full LNAV

- **Parents:** ["fn-4"]

With `nav = "full"` the `ca_code` kind MUST emit a valid IS-GPS-200 LNAV stream (TLM/HOW with correct TOW and parity, subframes 1 to 3 carrying a clock model and an ephemeris) whose ephemeris is derived from the trajectory the simulator flies, so that a receiver can compute a position fix. The mapping MUST be exact for a circular-Keplerian source and the osculating element set at `toe` for any other trajectory.

### `fn-5`: Trajectory models

Every trajectory MUST answer an ECEF position (metres) for any UTC time; non-stationary trajectories are sampled once onto a uniform grid and interpolated (see fid-3).

### `fn-5-stationary`: Stationary trajectories

- **Parents:** ["fn-5"]

`stationary_geodetic` (WGS84 latitude, longitude, height) and `stationary_ecef` MUST model an Earth-fixed point, with geodetic input converted to ECEF.

### `fn-5-tle`: TLE trajectory

- **Parents:** ["fn-5"]

The `tle` kind MUST propagate a two-line element set with SGP4 and rotate the TEME result into ECEF.

### `fn-5-keplerian`: Circular Keplerian trajectory

- **Parents:** ["fn-5"]

The `circular_keplerian` kind MUST fly a two-body circular orbit from semi-major axis, inclination, RAAN and true anomaly at epoch.

### `fn-5-rinex-nav`: Broadcast-ephemeris trajectory

- **Parents:** ["fn-5"]

The `rinex_nav` kind MUST parse GPS broadcast ephemeris from a RINEX 3 navigation file and fly the IS-GPS-200 satellite-position algorithm on it, so a `nav = "full"` transmitter carries the very ephemeris that flies its orbit.

### `fn-6`: Self-describing output

Each receiver band MUST be written as raw little-endian `complex64` (interleaved f32 I, f32 Q, no header) to `<receiver>__<band>.iq.c64`, accompanied by a JSON sidecar `<receiver>__<band>.iq.json` that carries epoch, sample rate, LO frequency and sample count at the top level plus a `descriptor` block naming the contributing transmitters (and PRNs) and embedding the full scene configuration. The sidecar MUST be identical whether the IQ was produced in memory or streamed.

### `fn-7`: Determinism

For a given scene and seed the output MUST be bit-for-bit reproducible, and identical regardless of chunk size and thread count. Parallel transmitter contributions MUST be reduced in a fixed order and all randomness MUST derive from the seed.

## Non-functional requirements

### `nf-1`: Throughput

- **Manual:** true

The per-transmitter contribution kernel MUST be covered by a Criterion benchmark (`iq_generation`) so throughput regressions are measurable.

### `nf-2`: Bounded memory

Peak memory MUST be independent of scene duration: IQ is generated and written in fixed-size chunks and the clock is a forward-only generator carrying O(1) state, so no full IQ buffer or clock track is ever materialised.

### `nf-3`: Configurable parallelism

The number of worker threads for the per-chunk transmitter fan-out MUST be configurable (`--threads`, 0 meaning all cores) and MUST NOT change the output.

### `nf-4`: Third-party interoperability

- **Manual:** true

Output MUST load in NumPy, GNU Radio and GNSS-SDR with no bespoke reader, which the headerless interleaved-f32 format of fn-6 provides. End-to-end evidence is the GNSS-SDR validation in the README (acquisition, tracking and position fix); it is a manual check, not an automated test.

## Signal-model equations

These tie the numbered equations of [SIGNAL_MODEL.md](../SIGNAL_MODEL.md) to
the code that implements them (the §9 equation ↔ code map, now checked by
tracey rather than by hand).

### `eq-e1`: Emission-time solve (E1, §2.1)

For each receiver sample the code MUST solve `τ + ‖P_s(τ) − P_o(t_utc)‖ / c = t_utc = t − δ_o(t)` for the emission time `τ` by fixed-point iteration, and evaluate the source at `τ`.

### `eq-e2`: Sagnac correction (E2, §2.2)

Inside the E1 solve the source position MUST be rotated into the reception-time ECEF frame, `P̃_s = R_3(ω_⊕·Δ)·P_s(τ)`, so the modelled delay is the true inertial light-travel time.

### `eq-e3`: Amplitude (E3, §2.3)

Each contribution's amplitude MUST be `A(t) = (λ / 4π·r(t))·√P_tx · √G_ant(θ_s(t))`: free-space path loss in volts times the receive-antenna gain at the source's arrival elevation.

### `eq-e4`: Received phase (E4, §2.4)

Each contribution's phase MUST be `φ = arg s(τ) + 2π·f_tx·τ − 2π·f_lo·t`: the transmit carrier rides on the solved emission time, the LO on the receiver's clock time, and nothing else is added.

### `eq-e5`: Superposition (E5, §3)

The incident field at a band MUST be `r(t) = Σ_i A_i(t)·s_i(τ_i(t))·e^{jφ_i(t)}` over the gated sources, summed before any post-sum operator.

### `eq-e6`: Antenna gain (E6, §4.1)

The receive antenna gain MUST be piecewise-linear in dBi from the horizon gain at 0° to the peak gain at 90° elevation, and the horizon gain minus the front/back ratio below the horizon, with an optional cutoff elevation below which the source is occluded.

### `eq-e7`: Oscillator (E7, §4.2)

The receiver's time offset MUST be `δ_o(t) = c₀ + m·t + ε_wFM(t) + ε_fFM(t)`, with white-FM synthesised as a Wiener walk calibrated to the τ = 1 s Allan deviation and flicker-FM as a flat Allan floor.

### `eq-e8`: IP3 nonlinearity (E8, §4.4)

The post-sum nonlinearity MUST be memoryless third-order, `y ↦ y·(1 − (4/3)·|y|² / P_IP3)`, with `P_IP3` the input third-order intercept power in watts.

### `eq-e9`: ADC collapse (E9, §4.6)

The ADC MUST clip to `±V_fs` and quantise symmetrically mid-tread, `x[n] = Q_b(clip(y(t_n)))` with `Δ = V_fs / (2^(b−1) − 1)`, deterministically per sample.

### `eq-e10`: Full forward model (E10, §5)

A band's samples MUST be `x[n] = Q_b(clip[ N(g·Σ_i A_i·s_i(τ_i)·e^{jφ_i}) + n[n] ])` at `t_n = n / f_s`: superposition, then thermal noise, then the configured post-sum operator chain in order (gain, nonlinearity, ADC), each an identity when unconfigured.
