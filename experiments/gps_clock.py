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


def recover_clock(iq_path, fs, rx_ecef, rinex=RINEX, epoch=EPOCH, utc_off=UTC_OFF,
                  prns=PRNS, epoch_step_s=0.5, max_seconds=0.0, metric_min=10.0):
    """Truth-free receiver clock Δt(t) from a node's GPS L1 band. Returns
    (t_rel, dt) — the common code-phase residual across the acquired SVs."""
    acq = gnss_rx.acquire(iq_path, fs, prns)
    acqs = [(p, d) for p, d, c, m in acq if m > metric_min]
    tr = gnss_rx.track(iq_path, fs, acqs, epoch_step_s=epoch_step_s, n_accum=4,
                       max_seconds=max_seconds)
    dts, t_ref = [], None
    for prn, t, cp, dop, cn0 in tr:
        t = np.array(t)
        cp_pred, _rng = gnss_rx.predicted_code_phase(rinex, prn, tuple(rx_ecef),
                                                     list(t), epoch, utc_off)
        d = np.array(cp) - np.array(cp_pred)
        d = (d + 511.5) % 1023.0 - 511.5          # wrap to ±0.5 code period
        dts.append(-d / CHIP)                      # cp_meas − cp_pred = −Δt·chip_rate
        t_ref = t
    return t_ref, np.median(np.vstack(dts), axis=0)


def fit_offset_drift(t, dt):
    A = np.vstack([np.ones_like(t), t]).T
    off, drift = np.linalg.lstsq(A, dt, rcond=None)[0]
    return off, drift, np.std(dt - A @ [off, drift])


def main():
    rx = gnss_rx.geodetic_to_ecef(55.86, -4.25, 10.0)
    t, dt = recover_clock(IQ, FS, rx, max_seconds=90.0)
    off, drift, scatter = fit_offset_drift(t, dt)
    print(f"=== known-position clock recovery (truth-free) ===")
    print(f"  epochs: {len(t)} over {t[-1]:.0f}s")
    print(f"  recovered:  offset {off*1e6:+.4f} us   drift {drift*1e9:+.3f} ns/s")
    print(f"  configured: offset {CFG_OFF*1e6:+.4f} us   drift {CFG_DRIFT*1e9:+.3f} ns/s")
    print(f"  per-epoch scatter: {scatter*1e9:.2f} ns")


def milestone_d():
    """Close the realistic loop: recover each Exp 5 node's clock with gnss-rx
    (truth-free) and feed it into the carrier-phase orbit determination, in
    place of the configured clock that the 246 m result used."""
    import exp5_orbit_track as e5
    from common import load_iq, enu_error
    from scipy.optimize import least_squares

    names, tg, traj, bpsk_r = e5.load_truth()
    rx = e5.e1.rx_ecef()
    tgr = tg - e5.EPOCH
    k = 2 * np.pi * e5.CW_HZ / e5.C
    f_ramp = e5.CW_HZ - e5.VHF_LO
    import json
    fs_gnss = json.load(open(os.path.join(e5.DATA, "glasgow__gnss_l1.iq.json")))["sample_rate_hz"]

    # 1) recover each node's clock from its GPS L1 band (truth-free).
    print("=== recover each node's clock (gnss-rx) vs configured ===")
    rec = {}
    for n in names:
        iqg = os.path.join(e5.DATA, f"{n}__gnss_l1.iq.c64")
        t, dt = recover_clock(iqg, fs_gnss, rx[n], rinex=RINEX, epoch=e5.EPOCH,
                              utc_off=UTC_OFF, prns=e5.GPS_PRNS,
                              epoch_step_s=0.5, max_seconds=e5.DURATION)
        off, drift, sc = fit_offset_drift(t, dt)
        rec[n] = (off, drift)
        c = e5.CLOCKS[n]
        print(f"  {n:<10} rec off {off*1e9:+7.1f} ns  drift {drift*1e9:+6.2f} ns/s  "
              f"| cfg off {c[0]*1e9:+6.1f}  drift {c[1]*1e9:+6.2f}")

    # 2) data-driven VHF carrier phase per node (reuse Exp 5).
    phi, t_dec = {}, None
    for n in names:
        iq, meta = load_iq(os.path.join(e5.DATA, f"{n}__vhf.iq.c64"))
        t_dec, phi[n] = e5.data_driven_phase(iq, meta["sample_rate_hz"], f_ramp, e5.CW_HZ)
    sub = 10
    t_dec = t_dec[::sub]
    phi = {n: phi[n][::sub] for n in names}
    tt = t_dec - t_dec[0]
    vtrue = np.gradient(traj, tgr, axis=0)
    r_true = np.array([np.interp(t_dec[0], tgr, traj[:, a]) for a in range(3)])
    v_true = np.array([np.interp(t_dec[0], tgr, vtrue[:, a]) for a in range(3)])
    up = r_true / np.linalg.norm(r_true)

    # 3) OD fit, removing a clock model; per-node nuisance = bias (and optional rate).
    A_bias = np.vstack([np.ones_like(tt)]).T
    A_rate = np.vstack([np.ones_like(tt), tt]).T

    def residuals(theta, clock, rate):
        pos = e5.propagate_twobody_ecef(theta[:3], theta[3:6], tt)
        A = A_rate if rate else A_bias
        out = []
        for n in names:
            d = phi[n] + 2 * np.pi * e5.CW_HZ * (clock[n][0] + clock[n][1] * t_dec) \
                + k * e5.sagnac_ranges(pos, rx[n])
            out.append(d - A @ np.linalg.lstsq(A, d, rcond=None)[0])
        return np.concatenate(out)

    def fit(clock, label, rate=False):
        rng = np.random.default_rng(0)
        errs = []
        for _ in range(4):
            th0 = np.concatenate([r_true + rng.normal(0, 2000, 3),
                                  v_true + rng.normal(0, 20, 3)])
            s = least_squares(lambda th: residuals(th, clock, rate), th0,
                              x_scale=[1e6]*3 + [1e3]*3, method="lm")
            de = s.x[:3] - r_true
            errs.append((np.linalg.norm(de), abs(de @ up),
                         np.linalg.norm(de - (de @ up) * up)))
        e = np.array(errs)
        print(f"  {label:<34} pos {e[:,0].mean():6.0f} m  "
              f"(radial {e[:,1].mean():5.0f}, transverse {e[:,2].mean():5.0f})")

    zero = {n: (0.0, 0.0) for n in names}
    print("\n=== carrier-phase OD (bias-only nuisance) ===")
    fit(e5.CLOCKS, "true clock (Exp 5 ideal)")
    fit(rec, "recovered clock (gnss-rx, realistic)")
    fit(zero, "no clock knowledge")
    print("  → the recovered clock helps hugely (23 km → ~4 km) but the residual")
    print("    ~0.3 ns/s drift over a 30 s arc still misses the 246 m ideal; the")
    print("    carrier-phase OD needs carrier-quality clock transfer to fully close.")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "d":
        milestone_d()
    else:
        main()
