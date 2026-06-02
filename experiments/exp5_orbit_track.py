#!/usr/bin/env python3
"""Experiment 5 — carrier-phase orbit tracking of a moving (satellite) emitter.

Exp 4 showed static single-frequency interferometry is ambiguity-limited. The
resolution: make the emitter move. A LEO travels ~7.5 km/s, so over a ~30 s
window the range to each node changes by hundreds of km = millions of VHF
cycles. Tracking that carrier phase *continuously* turns the ambiguous mod-2π
phase into an unambiguous accumulated range-change, and the changing geometry
pins the trajectory — GNSS-disciplined distributed carrier-phase orbit
determination. No integer ever has to be fixed.

Target: the ISS TLE (from examples/iss_uhf_pass.toml) at the GPS-valid epoch,
emitting a co-located BPSK (coarse seed) + CW (carrier phase). The truth
trajectory is reconstructed by range-multilateration of the truth-export ranges
(exact, Sagnac aside), so no Python SGP4 is needed.

This file is built in milestones; run it to (0) generate + check visibility,
(A) validate the tracked phase against truth range-change, (B) fit the orbit.
"""

import os
import json
import numpy as np
import exp1_tdoa as e1
from common import (load_iq, run_lipika, position_from_ranges, elevation_deg,
                    track_carrier_phase, enu_error, C)

REPO = e1.REPO
DATA = "/media/mike/36ada8f0-4d34-4343-ba81-487bb65bb372/lipika-test/exp_orbit"
BUILD = e1.BUILD
EPOCH = e1.EPOCH
DURATION = 30.0
GPS_PRNS = [2, 6, 12, 19, 24, 25]
RNX = "examples/gps_synth.rnx"
L1 = 1_575_420_000.0
VHF_LO = 137_000_000.0
BPSK_HZ = 137_000_000.0
CW_HZ = 137_300_000.0
# ISS two-line set (same as examples/iss_uhf_pass.toml) — a realistic moving LEO.
ISS_L1 = "1 25544U 98067A   24001.50000000  .00010000  00000+0  18000-3 0  9994"
ISS_L2 = "2 25544  51.6400 100.0000 0001000  90.0000 270.0000 15.50000000000009"
CLOCKS = {
    "glasgow": (0.0, 0.0, 1e-11), "edinburgh": (1.30e-7, 1e-9, 1e-11),
    "stirling": (-2.10e-7, -5e-10, 1e-11), "ayr": (3.40e-7, 2e-9, 1e-11),
    "perth": (-1.55e-7, 8e-10, 1e-11), "dumfries": (2.05e-7, -1e-9, 1e-11),
}


def write_scene(path):
    rx = e1.rx_ecef()
    payload_bits = int(e1.BIT_RATE * DURATION)
    trj = f'trajectory = {{ kind = "tle", line1 = "{ISS_L1}", line2 = "{ISS_L2}" }}'
    L = ["[scene]", f"epoch_unix = {EPOCH}", f"duration_s = {DURATION}", ""]
    L += ["[[transmitter]]", 'name = "tgt_bpsk"', f"carrier_hz = {BPSK_HZ}",
          "tx_power_dbw = 0.0", 'kind = "bpsk"', f"bit_rate_bps = {e1.BIT_RATE}",
          f"payload_bits = {payload_bits}", "payload_seed = 20260601",
          f"sync_word = {e1.SYNC}", trj, ""]
    L += ["[[transmitter]]", 'name = "tgt_cw"', f"carrier_hz = {CW_HZ}",
          "tx_power_dbw = 3.0", 'kind = "cw"', "baseband_offset_hz = 0.0", trj, ""]
    for prn in GPS_PRNS:
        L += ["[[transmitter]]", f'name = "gps_prn{prn:02d}"',
              f"carrier_hz = {L1}", "tx_power_dbw = 27.0", 'kind = "ca_code"',
              f"prn = {prn}", f"epoch_unix = {EPOCH}", 'nav = "full"',
              f'trajectory = {{ kind = "rinex_nav", file = "{RNX}", prn = {prn} }}', ""]
    for name, p in rx.items():
        off, fo, ad = CLOCKS[name]
        L += ["[[receiver]]", f'name = "{name}"',
              f'trajectory = {{ kind = "stationary_ecef", x_m = {p[0]:.4f}, '
              f"y_m = {p[1]:.4f}, z_m = {p[2]:.4f} }}",
              f"clock = {{ initial_time_offset_s = {off}, "
              f"frequency_offset_y0 = {fo}, allan_deviation_1s = {ad} }}",
              "[[receiver.band]]", 'band_id = "gnss_l1"', f"lo_carrier_hz = {L1}",
              "sample_rate_hz = 2048000.0", "noise_temperature_k = 290.0",
              "[[receiver.band]]", 'band_id = "vhf"', f"lo_carrier_hz = {VHF_LO}",
              "sample_rate_hz = 1000000.0", "noise_temperature_k = 290.0", ""]
    with open(path, "w") as f:
        f.write("\n".join(L))


def load_truth():
    """Per-node truth-export: t grid, target ranges, and the exact source
    trajectory (recorded by lipika in the truth-export)."""
    rx = e1.rx_ecef()
    names = list(rx.keys())
    tr = {n: json.load(open(os.path.join(DATA, f"{n}__truth.json"))) for n in names}
    tg = np.array(tr[names[0]]["t_utc"])
    g = {t["name"]: t for t in tr[names[0]]["transmitters"]}["tgt_cw"]
    traj = np.array([g["source_x_m"], g["source_y_m"], g["source_z_m"]]).T  # (N,3)

    def rng(node, txname):
        d = {t["name"]: t for t in tr[node]["transmitters"]}[txname]
        return np.array(d["range_m"])
    bpsk_r = {n: rng(n, "tgt_bpsk") for n in names}
    return names, tg, traj, bpsk_r


def main():
    os.makedirs(BUILD, exist_ok=True)
    cfg = os.path.join(BUILD, "orbit.toml")
    write_scene(cfg)
    if not os.path.exists(os.path.join(DATA, "glasgow__truth.json")):
        print("Generating orbit-track scene (to the drive)...")
        run_lipika(REPO, cfg, DATA, threads=0, chunk=2_000_000, truth_hz=50.0)
    else:
        print("Reusing existing generated data.")

    rx = e1.rx_ecef()
    names, tg, traj, bpsk_r = load_truth()

    # --- Step 0: visibility / geometry over the window ---
    print(f"\n=== Visibility (truth trajectory from range-multilateration) ===")
    print(f"window {DURATION:.0f}s, {len(tg)} truth epochs")
    rsat = np.linalg.norm(traj, axis=1)
    print(f"  target geocentric radius: {rsat.mean()/1e3:.1f} km "
          f"(alt ~{(rsat.mean()-6371e3)/1e3:.0f} km)")
    speed = np.linalg.norm(np.diff(traj, axis=0), axis=1) / np.diff(tg)
    print(f"  target speed: {speed.mean()/1e3:.2f} km/s; "
          f"moved {np.linalg.norm(traj[-1]-traj[0])/1e3:.1f} km over window")
    mid = len(tg) // 2
    print(f"  elevation at mid-epoch:")
    for n in names:
        el0 = elevation_deg(rx[n], traj[0])
        elm = elevation_deg(rx[n], traj[mid])
        ele = elevation_deg(rx[n], traj[-1])
        print(f"    {n:<10} start {el0:5.1f}°  mid {elm:5.1f}°  end {ele:5.1f}°")

    # --- Milestone A: track the CW phase, validate against truth range-change ---
    # De-rotating by the truth range, the residual must be exactly the receiver
    # clock −2π·f_cw·Δt(t). Fitting a line to it recovers the clock offset+drift;
    # the scatter about the line is how well the carrier phase tracks the
    # moving emitter's range (the simulator-fidelity metric).
    f_ramp = CW_HZ - VHF_LO
    tgr = tg - EPOCH  # relative grid for interpolation
    lam = C / CW_HZ
    print(f"\n=== Milestone A: CW phase vs truth range (λ={lam:.3f} m) ===")
    print(f"{'node':<10}{'clk offset (carrier)':>22}{'configured':>12}{'resid RMS':>14}")
    for n in names:
        iq, meta = load_iq(os.path.join(DATA, f"{n}__vhf.iq.c64"))
        fs = meta["sample_rate_hz"]
        t_dec, resid = track_carrier_phase(iq, fs, f_ramp, tgr, bpsk_r[n], CW_HZ)
        # resid(t) = −2π f_cw Δt(t) = −2π f_cw (offset + y0·t); fit a line
        A = np.vstack([np.ones_like(t_dec), t_dec]).T
        coef, *_ = np.linalg.lstsq(A, resid, rcond=None)
        offset_est = -coef[0] / (2 * np.pi * CW_HZ)
        scatter = resid - A @ coef
        rms_cyc = np.std(scatter) / (2 * np.pi)
        print(f"{n:<10}{offset_est*1e9:>19.2f} ns{CLOCKS[n][0]*1e9:>10.1f} ns"
              f"{rms_cyc*1e3:>9.2f} mcyc ({rms_cyc*lam*1e3:.2f} mm)")


if __name__ == "__main__":
    main()
