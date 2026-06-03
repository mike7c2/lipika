#!/usr/bin/env python3
"""Known-position GPS L1 clock recovery — truth-free, via gnss-rx.

Acquire + track every PRN (gnss_rx, Rust), then for each epoch compare the
measured C/A code phase to the phase *predicted* from the broadcast ephemeris
(RINEX) + the surveyed receiver position (gnss_rx.predicted_code_phase, reusing
lipika's Sagnac-correct emission solve). The common code-phase residual across
satellites is the receiver clock: cp_meas − cp_pred = −Δt·chip_rate. Fit a line
→ offset + drift. Nothing here uses lipika truth (only the broadcast ephemeris
and the known position, both legitimate a-priori); validated against GNSS-SDR.
"""

import os
import numpy as np
import gnss_rx

REPO = "/home/mike/Projects/DGSN_RESSURECT/lipika"
DATA = "/media/mike/36ada8f0-4d34-4343-ba81-487bb65bb372/lipika-test/exp_clock"
IQ = os.path.join(DATA, "glasgow__gnss_l1.iq.c64")
RINEX = os.path.join(REPO, "examples", "gps_synth.rnx")
FS = 2_048_000.0
EPOCH = 1_717_260_840.0
UTC_OFF = 18.0
CHIP = 1_023_000.0
PRNS = [2, 6, 12, 19, 24, 25]
# Configured (for scoring only): Δt(t) = 1 µs + 10 ns/s · t
CFG_OFF, CFG_DRIFT = 1.0e-6, 1.0e-8


def recover_clock(epoch_step_s=0.5, max_seconds=90.0):
    rx = gnss_rx.geodetic_to_ecef(55.86, -4.25, 10.0)  # surveyed node position
    acq = gnss_rx.acquire(IQ, FS, PRNS)
    acqs = [(p, d) for p, d, c, m in acq if m > 10]
    tr = gnss_rx.track(IQ, FS, acqs, epoch_step_s=epoch_step_s, n_accum=4,
                       max_seconds=max_seconds)
    dts = []
    t_ref = None
    for prn, t, cp, dop, cn0 in tr:
        t = np.array(t)
        cp = np.array(cp)
        cp_pred, _rng = gnss_rx.predicted_code_phase(RINEX, prn, tuple(rx),
                                                     list(t), EPOCH, UTC_OFF)
        d = cp - np.array(cp_pred)
        d = (d + 511.5) % 1023.0 - 511.5          # wrap to ±0.5 code period
        dts.append(-d / CHIP)                      # cp_meas − cp_pred = −Δt·chip_rate
        t_ref = t
    dt = np.median(np.vstack(dts), axis=0)         # common clock across SVs
    return t_ref, dt


def main():
    t, dt = recover_clock()
    A = np.vstack([np.ones_like(t), t]).T
    off, drift = np.linalg.lstsq(A, dt, rcond=None)[0]
    scatter = np.std(dt - A @ [off, drift])
    print(f"=== known-position clock recovery (truth-free) ===")
    print(f"  epochs: {len(t)} over {t[-1]:.0f}s")
    print(f"  recovered:  offset {off*1e6:+.4f} us   drift {drift*1e9:+.3f} ns/s")
    print(f"  configured: offset {CFG_OFF*1e6:+.4f} us   drift {CFG_DRIFT*1e9:+.3f} ns/s")
    print(f"  per-epoch scatter: {scatter*1e9:.2f} ns")


if __name__ == "__main__":
    main()
