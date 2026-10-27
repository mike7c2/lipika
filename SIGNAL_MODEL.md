# The signal model — the algebra of the chain

This is the mathematical specification the code implements. The guiding principle:
**everything is a continuous function of time; the only collapse to discrete
samples is the ADC.** The transmitter is a function `s(t)`; propagation, the
oscillator, the antenna and the front-end are *operators* on it; sampling +
quantisation is the single discrete, lossy step. Because every operator is exact
and stateless-in-form, an *estimator* can run the same algebra **backwards** (§6).

The numbered equations **(E1)–(E10)** each cite the code that implements them; §9
collects the equation ↔ code map. Symbols are in the glossary (§8). What the model
deliberately leaves out is in **§7 (notable exclusions)** — read it before trusting
a result that depends on an unmodelled effect.

---

## 1. The transmitter — `s(t)`

A transmitter is a **stateless complex-baseband function** of emission time
(`src/signal/mod.rs::SignalModel`), riding a nominal RF carrier `f_tx`:

$$ s_{\text{RF}}(t) \;=\; \Re\!\Big[\, s(t)\, e^{\,j 2\pi f_{tx} t}\,\Big], \qquad s(t)\in\mathbb{C}. $$

`s(t)` is the *information*; `f_tx` is the carrier. Concrete models define `s(t)`:

| signal | `s(t)` |
|---|---|
| CW tone (offset `f_b`) | $s(t)=e^{\,j2\pi f_b t}$ |
| BPSK (bits `b_k`, rate `R`) | $s(t)=b_{\lfloor R t\rfloor}\in\{\pm1\}$ |
| GPS C/A (PRN, NAV `d`) | $s(t)=c_{\text{PRN}}(t)\,d(t)\in\{\pm1\}$ |

The crucial property: **`s` is never pre-sampled.** Doppler, delay and clock all
act by changing *where* `s` is evaluated, not by resampling a buffer — which is
what keeps the model exact.

---

## 2. Propagation — one source to one receiver

### 2.1 Emission time (the light-time solve)

The receiver timestamps each sample with its **own clock**, which reads `t = n/f_s`
at sample `n`. That clock is imperfect: its error `δ_o(t)` (§4.2) means the **true**
(coordinate) time at which the sample is actually captured is

$$ t_{\text{utc}} \;=\; t - \delta_o(t). $$

**All geometry is evaluated at true time.** A sample whose true reception instant is
`t_utc` was emitted earlier, at the **emission time** `τ`, because the signal
travelled the slant range at `c`. With the source at `P_s` and the receiver at
`P_o`, `τ` solves the implicit equation:

$$ \boxed{\;\tau \;+\; \frac{\big\lVert \mathbf P_s(\tau) - \mathbf P_o(t_{\text{utc}})\big\rVert}{c} \;=\; t_{\text{utc}} \;=\; t - \delta_o(t)\;} \tag{E1} $$

> **code** — the pipeline forms the true-time grid `t_utc_rel = t_local_rel − delta_t`
> (`pipeline.rs`), then `channel::solve_emission_times_for` runs the per-sample loop.

i.e. emission time + light-travel = true reception time. **This is where the clock
enters the geometry**: `τ` is anchored to `t_utc = t − δ_o`, so every quantity derived
from it — the envelope `s(τ)`, the range `r`, the Doppler — carries the clock. The
sample is *indexed* by clock time `t`; the physics happens at `t_utc`. Everything
dynamic — Doppler, code-rate Doppler, sample-rate skew — falls out of `τ` varying
across samples, with no separate "apply Doppler" or "apply clock" step. This is the
exact (*wideband*) treatment of Doppler as a time-scaling of the waveform rather than
the narrowband `e^{j2πf_d t}` frequency-shift approximation [KW65; VanTrees-III;
Bello63]; in GNSS terms the carrier and code Doppler are then automatically the same
`ṙ/c` scaling [KH17]. The light-time solve itself is the standard ephemeris
"reception-case" construction [SPICE96; IERS2010].

### 2.2 Earth rotation (Sagnac)

`P_s(τ)` and `P_o(t_utc)` are ECEF positions at *different* instants, but the geometric
path must be measured in one frame. The source is rotated into the reception-frame
by Earth's rotation over the transit, `Δ = t_utc − τ = r/c` (linearised, exact to sub-mm):

$$ \tilde{\mathbf P}_s = R_3(\omega_\oplus\,\Delta)\,\mathbf P_s(\tau)\;\approx\;\big(P_{s,x}+\omega_\oplus\Delta\,P_{s,y},\;\; P_{s,y}-\omega_\oplus\Delta\,P_{s,x},\;\; P_{s,z}\big). \tag{E2} $$

> **code** — the `wtau`/`sx_rot`/`sy_rot` block inside `channel::solve_emission_times_for`.

The range used everywhere below is this Sagnac-corrected slant range,
$r(t)=\lVert\tilde{\mathbf P}_s-\mathbf P_o\rVert$. Omitting it costs ~tens of ns —
significant for a timing reference [IERS2010]. (The OD side, `common.emission_ranges`, solves
the *same* equation — reception- vs emission-time geometry was the light-time bug
that cost ~50 m of along-track until fixed.)

### 2.3 Amplitude — free-space path loss + antenna

$$ A(t) \;=\; \underbrace{\frac{\lambda}{4\pi\,r(t)}\sqrt{P_{tx}}}_{\text{FSPL (volts)}}\;\cdot\;\underbrace{\sqrt{G_{\text{ant}}\!\big(\theta_s(t)\big)}}_{\text{rx antenna}},\qquad \lambda=\frac{c}{f_{tx}}. \tag{E3} $$

> **code** — `amp_factor`/`amplitude` (FSPL) × `geom_gain` (antenna √G), in
> `channel::accumulate_transmitter_contribution`'s mix loop.

`P_tx` is EIRP (W); the `λ/4πr` factor is the free-space field amplitude; `√` because
amplitude, not power. `G_ant` is the receiver antenna's power gain at the source's
arrival elevation `θ_s(t)` (§4.1). lipika amplitudes are in **√W**, so `|·|²` is
power directly.

### 2.4 The received phase

Two phase processes meet at the mixer: the **arriving carrier**, whose phase was set
at the source at emission time `τ`, and the receiver **LO**. The LO is generated from
the same oscillator that defines the clock, so at sample `n` its phase is exactly
`2πf_lo·t` — clock time, *by definition of the clock*. The baseband phase is the
arriving carrier phase minus the LO phase:

$$ \boxed{\;\varphi \;=\; \underbrace{\arg s(\tau)}_{\text{message}} \;+\; \underbrace{2\pi f_{tx}\,\tau}_{\text{arriving carrier}} \;-\; \underbrace{2\pi f_{lo}\,t}_{\text{LO (clock time)}}\;} \tag{E4} $$

Only two frequencies appear: `f_tx` rides on `τ` (true-time geometry), `f_lo` rides on
`t` (clock time). Expanding with `τ = t_utc − r/c = t − δ_o − r/c` separates the
familiar terms and matches the code exactly:

$$ \varphi \;=\; \arg s(\tau)\;-\;2\pi f_{tx}\Big(\tfrac{r}{c}+\delta_o\Big)\;+\;2\pi(f_{tx}-f_{lo})\,t. $$

> **code** — the `phase = neg_tau_src*(delay+delta_t) + tau_if*t_local` line in
> `channel::accumulate_transmitter_contribution` (`neg_tau_src=−2πf_tx`,
> `tau_if=2π(f_tx−f_lo)`, `delta_t`=`δ_o`, `t_local`=`t`). The envelope `s(τ)` is
> evaluated at `τ` from the true-time solve (E1), so it carries `δ_o` too.

Reading the terms:

- **message** `arg s(τ)` — sampled at the emission time `τ`. Because `τ` is anchored
  to `t_utc = t − δ_o` (E1), the clock shifts the **code phase** here.
- **carrier delay** `−2πf_tx·r/c` — carrier phase over the light path; its
  time-derivative *is* the carrier Doppler, emerging automatically as `τ` walks.
- **clock** `−2πf_tx·δ_o` — the receiver oscillator, via the LO time base. It carries
  the **same sign and the same `δ_o`** as the code-phase shift in `arg s(τ)`, so the
  clock acts as a *common time bias* on code and carrier alike: the apparent delay is
  `t − τ = r/c + δ_o`, i.e. **pseudorange `= r + c·δ_o`** (the standard GNSS
  observation model [ME06; KH17]). This is why clock and LO are one object — one
  `δ_o` drives sample timing and carrier together.
- **IF / de-ramp** `2π(f_tx−f_lo)t` — for an off-tune LO the residual sits at the
  intermediate frequency `f_tx−f_lo`; vanishes on-tune (`f_lo=f_tx`).

(The `−2πf_tx·δ_o` clock term paired with the clock-time de-ramp is algebraically
identical to the `−2πf_lo·δ_o` "LO phase-noise" form with the de-ramp expressed in
true time `t_utc`; they differ only in bookkeeping. The code keeps the de-ramp on
clock time, hence `f_tx`.)

For an on-tune receiver (`f_lo=f_tx`) this collapses to the form in the code:
$$ \varphi=\arg s(\tau)\;-\;2\pi f_{tx}\Big(\tfrac{r(t)}{c}+\delta_o(t)\Big),\qquad \tau = t-\delta_o-\tfrac{r}{c}. $$

So one source contributes the continuous complex sample $A\,|s(\tau)|\,e^{\,j\varphi}$.

---

## 3. Superposition — the incident field

A receiver band sees the **sum over all in-band sources** (the `bands_overlap`
gate decides "in-band"):

$$ \boxed{\; r(t)\;=\;\sum_{i\in\text{sources}} A_i(t)\,s_i\big(\tau_i(t)\big)\,e^{\,j\varphi_i(t)} \;} \tag{E5} $$

> **code** — `pipeline::generate_band` sums each transmitter's
> `accumulate_transmitter_contribution` into the band buffer `acc` (fixed order ⇒
> deterministic). `bands_overlap` gates the source set.

This is still an exact continuous function. Horizon occlusion enters here: a source
below the receiver's local horizon is dropped (`A_i \to 0`).

---

## 4. The receiver front-end — operators on `r(t)`

The front-end is a composition of operators (`src/rx`). Some act **per source**
(woven into §2, because they depend on each source's geometry); the rest act on
the **summed band** `r(t)`. The chain, in order:

### 4.1 Antenna (per source) — `G_ant`
A gain on the arrival direction. The chosen self-phased QFH is azimuthally
symmetric, so gain depends on elevation `θ` only:
$$ G_{\text{ant}}(\theta)\;[\text{dBi}]=\begin{cases} G_{\text{horiz}} + (G_{\text{peak}}-G_{\text{horiz}})\,\theta/90^\circ & \theta\ge0\\[2pt] G_{\text{horiz}} - \text{F/B} & \theta<0\ \text{(back-lobe)} \end{cases} \tag{E6} $$
applied as `√G` to each source's amplitude (§2.3). `θ` is the same arrival angle the
horizon mask uses.

> **code** — `rx::Antenna::amplitude_at`; folded into the `geom_gain` pass of
> `channel::accumulate_transmitter_contribution`.

### 4.2 The oscillator (per source, via `δ_o`) — `θ_{osc}(t)`
The node has **one** oscillator phase process. Its fractional-frequency noise
`y(t)` integrates to the **time error** `δ_o(t)` that appears in `φ` (§2.4):
$$ \delta_o(t)\;=\;\underbrace{c_0}_{\text{offset}} + \underbrace{m\,t}_{\text{drift}} + \underbrace{\varepsilon_{\text{wFM}}(t)}_{\text{white FM (Allan }\sigma_1)} + \underbrace{\varepsilon_{\text{fFM}}(t)}_{\text{flicker FM (Allan floor)}} \tag{E7} $$

> **code** — `clock::NodeClockSpec` (the parameters), `clock::ClockStream` (the
> streaming `δ_o(t)`: deterministic part + white-FM Wiener walk), `clock::Flicker`
> (the flicker-FM floor).

- `c₀` — static bias (per-node, unobservable absolute; a nuisance in the OD).
- `m` — static fractional-frequency offset (ppm-scale drift).
- `ε_wFM` — white-FM **Wiener walk**, calibrated so $\sigma_y(\tau{=}1\text{s})=\sigma_1$; Allan $\propto\tau^{-1/2}$.
- `ε_fFM` — **flicker FM** (1/f), a flat Allan **floor** — the oscillator-grade knob
  (TCXO ≈ 1e‑9…1e‑10, OCXO ≈ 1e‑12…1e‑13). Built as a sum of log-spaced AR(1)
  relaxations, self-calibrated (`clock::Flicker`).

The `c₀ + m·t +` power-law decomposition is the standard frequency-stability model
[Barnes71; IEEE1139]; the white-/flicker-FM synthesis (Allan-correct, causal,
scale-invariant) follows the discrete power-law-noise method of [Kasdin92; Kasdin95].

The **same `δ_o`** sets the sample-timing grid and the LO time base. So
**clock = LO = sample clock**: one oscillator drives the sample timing — hence `τ`,
the code phase and range, via `t_utc = t − δ_o` (§2.1) — **and** the carrier phase,
via the `−2πf_tx·δ_o` term (§2.4). Because it is one `δ_o`, code and carrier see the
*same* clock: a common time bias, pseudorange `= r + c·δ_o`.

### 4.3 RF/IF gain + AGC (post-sum) — `g`
A linear scalar (dB) lifting the weak (FSPL-attenuated) signal into the ADC range:
$y \mapsto g\,y$.

### 4.4 Nonlinearity (post-sum) — IP3
Memoryless third-order LNA/mixer model with input intercept power `P_{\text{IP3}}`:
$$ y \;\mapsto\; y\Big(1-\tfrac{4}{3}\,\frac{|y|^2}{P_{\text{IP3}}}\Big). \tag{E8} $$

> **code** — `rx::PostOp::Nonlinearity`.

Small-signal linear; gain-compresses and produces IM3 at strong input. The binding
spec for terrestrial rejection. (`P_IP3` from the dBm datasheet value, direct since
amplitude² is power.) This memoryless-cubic intercept model and the gain/NF/quantise
front-end chain are the standard behavioural receiver model [Betz16; Borre07].

### 4.5 Thermal noise (post-sum) — `n(t)`
Additive complex AWGN at the receiver floor $P_n=k_B T B$ (`B` = sample rate):
$y \mapsto y + n,\ n\sim\mathcal{CN}(0,P_n)$.

### 4.6 The ADC — **the collapse**
The only discrete, lossy, non-invertible step: sample on the **clock** grid `t_n=n/f_s`
(§2.1 — `t_n` is clock time; the geometry inside `φ`, `s(τ)`, `A` is evaluated at the
true time `t_n − δ_o`, and aperture jitter is a further time-warp of this grid), clip
to the full scale `±V_{fs}`, and quantise to `b` bits:
$$ x[n]\;=\;Q_b\!\big(\mathrm{clip}_{\pm V_{fs}}(y(t_n))\big),\qquad Q_b(v)=\Delta\,\mathrm{round}(v/\Delta),\ \ \Delta=\frac{V_{fs}}{2^{\,b-1}-1}. \tag{E9} $$

> **code** — `rx::PostOp::Adc`.

---

## 5. The complete forward model

Putting it together — what lipika computes for output sample `n` of a band:

$$ \boxed{\; x[n]\;=\;Q_b\!\Bigg(\mathrm{clip}\bigg[\,\mathcal N\!\Big(g\!\sum_i A_i(t_n)\,s_i(\tau_i)\,e^{\,j\varphi_i(t_n)}\Big)\;+\;n[n]\,\bigg]\Bigg)\;} \quad t_n=\tfrac{n}{f_s} \tag{E10} $$

> **code** — `pipeline::generate_band` per band: the `Σ` (E5) + `add_awgn_at_power`
> + `RxChain::apply_post_sum` (gain, NL, ADC), then `sink.push`.

Everything inside `Q_b(clip(·))` is an exact continuous expression evaluated at
`t_n`; `Q_b` is the single collapse. With an unconfigured front-end (`g=1`, no NL,
ideal ADC) this reduces to the plain `Σ` — which is why adding the front-end leaves
old scenes byte-identical.

---

## 6. The inverse — what an estimator does

Because every operator above is exact and the lossy ones (NL, thermal,
quantisation) are well-characterised noise, an estimator can run the **same algebra
backwards** on `x[n]`, inverting the *invertible* operators and treating the rest as
noise:

- **De-rotation** = inverse of §2.4's carrier+IF terms: multiply by
  `e^{-j[2π(f_tx-f_lo)t - 2πf_tx r̂/c]}` using a *predicted* range `r̂` (from
  ephemeris+position, or the orbit fit) → leaves the slow residual `-2πf·δ_o + noise`.
- **Code-phase prediction** = inverse of the §2.1 time-warp: predict `τ` from the
  broadcast ephemeris + surveyed position → the measured-minus-predicted code shift is
  `δ_o` (the forward model put it there via `t_utc = t − δ_o`, so it is recoverable).
- **Clock recovery** = estimating `δ_o(t)` (E7): since the clock is a common time bias
  (pseudorange `= r + c·δ_o`, §2.4) it is common-mode across all visible SVs (geometry
  wiped); code resolves integer-ns, carrier carries the ps fraction.

Because forward and inverse are the *same* `φ`, `r`, `δ`, a simulator and an
estimator built on this algebra cannot drift apart — and a result proved against the
forward geometry (e.g. the light-time term, E1–E2) transfers directly to the inverse.
This analysis-by-synthesis stance — invert the forward model directly rather than
track loosely — is the basis of GNSS direct position estimation and vector tracking
[Closas07]. **Caveat:** this is *internal* consistency, not fidelity. Any effect the forward model
omits (§7) is invisible to this round-trip, so passing it proves the algebra, not the
physics.

---

## 7. Notable exclusions — what the model does *not* capture

The channel is **vacuum truth**: geometry, free-space loss, antenna, one *receiver*
clock, and a characterised front-end. The effects below are real but deliberately
**not** modelled today. They are listed so an estimator knows what its "truth" omits,
and so the self-consistency of §6 is not mistaken for real-world accuracy.

- **Ionospheric propagation.** No group delay / phase advance / TEC. At UHF this is
  the **largest** missing term — tens to >100 ns of group delay for typical TEC, i.e.
  *larger* than the Sagnac correction the geometry does include (it scales as `1/f²`,
  §2.2). The network's TEC / atmospheric-sounding ambitions live exactly here; today
  the ionosphere is an *estimand*, not a simulated channel.
- **Tropospheric delay.** No neutral-atmosphere delay — small versus the ionosphere
  at these frequencies, but nonzero and growing toward the horizon.
- **Transmitter clock.** `s(t)` is phase-perfect: the *source* oscillator is ideal
  (no transmit `δ_o`, no transmit phase noise). A per-source clock is a clean
  extension (a mirror of §4.2), but for downlink / TDOA we currently assume emitters
  are stable. Adding the term doesn't hurt; we just hold it at zero.
- **Polarisation.** No polarisation state or mismatch — the QFH is a scalar gain
  `G_ant(θ)` (§4.1). No circular/linear mismatch loss, no Faraday rotation.
- **Receiver noise figure.** Thermal noise is `k_B T B` at the floor (§4.5) with no
  noise-figure / system-temperature budget, and is added *post-gain* — so absolute
  SNR tracks the configured gain `g` rather than a cascade NF.
- **Relativistic terms.** First-order (classical) Doppler only, via the light-time
  walk of `τ`; second-order time dilation and gravitational frequency shift are
  dropped (~`6e-10` fractional — borderline only for ps-class carrier-phase work).
- **Multipath / scattering.** One direct ray per source; no ground bounce, no local
  scatterers.

None of these conflict with the "evaluate the operator, don't resample" discipline;
they are scoped out, not overlooked. When one is added it slots in as another exact
operator on `s(τ)` or another term in `δ`.

---

## 8. Glossary

| symbol | meaning |
|---|---|
| `t` | node-**clock** time of a sample (the timestamp `n/f_s`; the LO time base) |
| `t_utc` | true (coordinate) reception time `= t − δ_o`; geometry is solved here |
| `τ` | emission time (light-time solve at `t_utc`, §2.1) |
| `s(t)` | transmitter complex baseband (§1) |
| `f_tx`, `f_lo` | source RF carrier; receiver band LO |
| `f_s`, `B` | band sample rate (= noise bandwidth) |
| `P_s, P_o` | source / observer ECEF position |
| `r(t)` | Sagnac-corrected slant range |
| `c`, `λ`, `ω⊕` | speed of light; wavelength `c/f_tx`; Earth rotation rate |
| `P_tx` | transmit EIRP (W); lipika amplitude² is power |
| `G_ant(θ)` | rx antenna power gain at elevation `θ` |
| `δ_o(t)` | node clock/oscillator time error (§4.2) |
| `c₀, m` | clock static offset; fractional-frequency drift |
| `σ_1` | white-FM Allan deviation at τ=1 s |
| `ε_wFM, ε_fFM` | white-FM walk; flicker-FM floor |
| `g`, `P_IP3` | RF/IF gain; input third-order intercept (W) |
| `P_n=k_BTB` | thermal noise power |
| `V_fs`, `b`, `Δ` | ADC full scale; bits; LSB |
| `Q_b`, `n` (index) | `b`-bit quantiser; output sample index |

---

## 9. Equation ↔ code map

Each formula and the function that implements it. Source files are under
`src/`; every function below carries a `(E…)` back-reference in its doc-comment.

The map is machine-checked: each function below carries a `// r[impl eq.eN]`
annotation, `docs/REQUIREMENTS.md` defines `eq.e1`–`eq.e10`, and `tracey query
rule eq.e4` (or `scripts/trace-check.sh`) lists every implementation and test of
an equation — see `docs/TRACEABILITY.md`.

| eq | what | code |
|----|------|------|
| **E1** | emission-time (light-time) solve | `channel::solve_emission_times_for` |
| **E2** | Sagnac (Earth-rotation) correction | same loop (`wtau`/`sx_rot`/`sy_rot`) |
| **E3** | amplitude — FSPL × antenna √G | `channel::accumulate_transmitter_contribution` (`amp_factor`, `geom_gain`) |
| **E4** | received phase `φ(t)` | `channel::accumulate_transmitter_contribution` (the `phase = …` line) |
| **E5** | superposition `r(t)` | `pipeline::generate_band` (sum of contributions) |
| **E6** | antenna gain `G_ant(θ)` | `rx::Antenna::amplitude_at` |
| **E7** | oscillator `δ_o(t)` | `clock::{NodeClockSpec, ClockStream, Flicker}` |
| **E8** | IP3 nonlinearity | `rx::PostOp::Nonlinearity` |
| **E9** | ADC collapse (clip + quantise) | `rx::PostOp::Adc` |
| **E10** | full forward `x[n]` | `pipeline::generate_band` (Σ + noise + `RxChain::apply_post_sum`) |

The reverse links live in the code: each listed function's doc-comment names its
equation, so `grep -rn "SIGNAL_MODEL" src/` walks from code back to the maths.

---

## 10. References & relation to prior art

No part of this model is novel physics. Each ingredient is standard; the contribution
is the *integration* — exact wideband Doppler, a behavioural front end, a physical
clock-noise process, and an explicit forward/inverse symmetry, kept as one continuous
operator algebra whose only discrete step is the ADC. The literature each piece rests
on:

- **Doppler as time-warp / exact (wideband) channel** (§1, §2.1, §2.4)
  - **[Bello63]** P. A. Bello, "Characterization of randomly time-variant linear channels," *IEEE Trans. Communications Systems*, 11(4):360–393, 1963.
  - **[KW65]** E. J. Kelly and R. P. Wishner, "Matched-filter theory for high-velocity, accelerating targets," *IEEE Trans. Military Electronics*, MIL-9(1):56–69, 1965.
  - **[VanTrees-III]** H. L. Van Trees, *Detection, Estimation, and Modulation Theory, Part III: Radar–Sonar Processing*, Wiley, 1971 (repr. 2001) — wideband ambiguity, Doppler as time-scaling.
- **Light-time solve & Earth-rotation (Sagnac)** (§2.1–2.2)
  - **[SPICE96]** C. H. Acton, "Ancillary data services of NASA's Navigation and Ancillary Information Facility (NAIF)," *Planetary and Space Science*, 44(1):65–70, 1996 — the ephemeris "reception-case" light-time construction.
  - **[IERS2010]** G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*, IERS Technical Note No. 36 — Earth-rotation/Sagnac and relativistic propagation terms.
- **GNSS observation model & code/carrier Doppler** (§2.4, §4.2)
  - **[ME06]** P. Misra and P. Enge, *Global Positioning System: Signals, Measurements, and Performance*, 2nd ed., Ganga-Jamuna Press, 2006 — pseudorange `ρ = r + c·δt`.
  - **[KH17]** E. D. Kaplan and C. J. Hegarty (eds.), *Understanding GPS/GNSS: Principles and Applications*, 3rd ed., Artech House, 2017 — carrier and code Doppler as the shared `ṙ/c` scaling.
- **Clock / oscillator noise** (§4.2)
  - **[Barnes71]** J. A. Barnes et al., "Characterization of frequency stability," *IEEE Trans. Instrumentation and Measurement*, IM-20(2):105–120, 1971.
  - **[IEEE1139]** IEEE Std 1139, *Definitions of Physical Quantities for Fundamental Frequency and Time Metrology — Random Instabilities*.
  - **[Kasdin92]** N. J. Kasdin and T. Walter, "Discrete simulation of power law noise," *Proc. IEEE Frequency Control Symposium*, pp. 274–283, 1992.
  - **[Kasdin95]** N. J. Kasdin, "Discrete simulation of colored noise and stochastic processes and 1/f^α power law noise generation," *Proc. IEEE*, 83(5):802–827, 1995.
- **Receiver front-end (behavioural)** (§4.3–4.6)
  - **[Betz16]** J. W. Betz, *Engineering Satellite-Based Navigation and Timing*, IEEE Press/Wiley, 2016 — noise figure, nonlinearity, quantisation effects.
  - **[Borre07]** K. Borre, D. M. Akos, et al., *A Software-Defined GPS and Galileo Receiver*, Birkhäuser, 2007 — IF-sampling front end and quantisation loss.
- **Forward-model inversion (analysis by synthesis)** (§6)
  - **[Closas07]** P. Closas, C. Fernández-Prades, J. A. Fernández-Rubio, "Maximum likelihood estimation of position in GNSS," *IEEE Signal Processing Letters*, 14(5):359–362, 2007 — direct position estimation; the forward model inverted directly.
