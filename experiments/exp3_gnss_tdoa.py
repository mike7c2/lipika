#!/usr/bin/env python3
"""Experiment 3 — GNSS-disciplined timing rescues TDOA geolocation.

The keystone distributed-receiver result. Exp 1 showed unsynchronised node
clocks wreck TDOA (~33 m → ~79 m horizontal). Here each node additionally
carries a GNSS L1 band that shares its clock — exactly as real hardware does.
We recover each node's clock offset from its GPS signal (code-phase time
transfer against the truth-export geometry), then use it to de-bias the VHF
TDOA. The corrected fix should recover most of the ideal-clock accuracy.

One scene, two bands per receiver:
  * gnss_l1 (1575.42 MHz, 2.048 Msps): 6 GPS SVs, full LNAV, from gps_synth.rnx
  * vhf     (137 MHz, 8 Msps):         the airborne BPSK emitter to geolocate

Because all of a receiver's bands run off one clock, the offset measured on
gnss_l1 IS the offset corrupting the vhf TDOA. Data goes to the mounted drive.
"""

import os
import numpy as np
import exp1_tdoa as e1
from common import (load_iq, run_lipika, fft_xcorr_tdoa, multilaterate_tdoa,
                    enu_error, C)
from gps import time_transfer_offset

REPO = e1.REPO
DATA = "/media/mike/36ada8f0-4d34-4343-ba81-487bb65bb372/lipika-test/exp_gnss_tdoa"
BUILD = e1.BUILD
EPOCH = e1.EPOCH                      # 1717260840, divisible by 6 for full LNAV
DURATION = 2.0
GPS_PRNS = [2, 6, 12, 19, 24, 25]    # the constellation in examples/gps_synth.rnx
RNX = "examples/gps_synth.rnx"
L1 = 1_575_420_000.0

# Same unsynchronised clocks as Exp 1's "clk" run.
CLOCKS = {
    "glasgow":   (0.0,      0.0,    1e-11),
    "edinburgh": (1.30e-7,  1e-9,   1e-11),
    "stirling":  (-2.10e-7, -5e-10, 1e-11),
    "ayr":       (3.40e-7,  2e-9,   1e-11),
    "perth":     (-1.55e-7, 8e-10,  1e-11),
    "dumfries":  (2.05e-7, -1e-9,   1e-11),
}


def write_scene(path):
    rx = e1.rx_ecef()
    em = e1.emitter_ecef()
    payload_bits = int(e1.BIT_RATE * DURATION)
    L = ["[scene]", f"epoch_unix = {EPOCH}", f"duration_s = {DURATION}", ""]
    # VHF emitter (airborne BPSK)
    L += [
        "[[transmitter]]", 'name = "emitter"', f"carrier_hz = {e1.CARRIER}",
        f"tx_power_dbw = {e1.TX_POWER_DBW}", 'kind = "bpsk"',
        f"bit_rate_bps = {e1.BIT_RATE}", f"payload_bits = {payload_bits}",
        "payload_seed = 20260601", f"sync_word = {e1.SYNC}",
        f'trajectory = {{ kind = "stationary_ecef", x_m = {em[0]:.4f}, '
        f"y_m = {em[1]:.4f}, z_m = {em[2]:.4f} }}", "",
    ]
    # GPS constellation (rinex_nav, full LNAV)
    for prn in GPS_PRNS:
        L += [
            "[[transmitter]]", f'name = "gps_prn{prn:02d}"',
            f"carrier_hz = {L1}", "tx_power_dbw = 27.0", 'kind = "ca_code"',
            f"prn = {prn}", f"epoch_unix = {EPOCH}", 'nav = "full"',
            f'trajectory = {{ kind = "rinex_nav", file = "{RNX}", prn = {prn} }}',
            "",
        ]
    # Receivers: gnss_l1 + vhf bands, one shared clock each.
    for name, p in rx.items():
        off, fo, ad = CLOCKS[name]
        L += [
            "[[receiver]]", f'name = "{name}"',
            f'trajectory = {{ kind = "stationary_ecef", x_m = {p[0]:.4f}, '
            f"y_m = {p[1]:.4f}, z_m = {p[2]:.4f} }}",
            f"clock = {{ initial_time_offset_s = {off}, "
            f"frequency_offset_y0 = {fo}, allan_deviation_1s = {ad} }}",
            "[[receiver.band]]", 'band_id = "gnss_l1"', f"lo_carrier_hz = {L1}",
            "sample_rate_hz = 8192000.0", "noise_temperature_k = 290.0",
            "[[receiver.band]]", 'band_id = "vhf"',
            f"lo_carrier_hz = {e1.CARRIER}", "sample_rate_hz = 8000000.0",
            "noise_temperature_k = 290.0", "",
        ]
    with open(path, "w") as f:
        f.write("\n".join(L))


def recover_clock(node, out_dir):
    """Code-phase GNSS time transfer: per-PRN clock offset, then median."""
    iq, meta = load_iq(os.path.join(out_dir, f"{node}__gnss_l1.iq.c64"))
    fs = meta["sample_rate_hz"]
    import json
    truth = json.load(open(os.path.join(out_dir, f"{node}__truth.json")))
    tg = truth["t_utc"]
    by_name = {t["name"]: t for t in truth["transmitters"]}
    ests, quals = [], []
    for prn in GPS_PRNS:
        t = by_name[f"gps_prn{prn:02d}"]
        off, q = time_transfer_offset(
            iq, fs, EPOCH, tg, t["delay_s"], L1, prn,
            sv_epoch=EPOCH, block_ms=1.0, max_seconds=1.0)
        ests.append(off)
        quals.append(q)
    ests = np.array(ests)
    return float(np.median(ests)), ests, np.array(quals)


def main():
    os.makedirs(BUILD, exist_ok=True)
    cfg = os.path.join(BUILD, "gnss_tdoa.toml")
    write_scene(cfg)
    if not os.path.exists(os.path.join(DATA, "glasgow__truth.json")):
        print("Generating GNSS+VHF scene (this is ~1 GB, to the drive)...")
        run_lipika(REPO, cfg, DATA, threads=0, chunk=2_000_000, truth_hz=50.0)
    else:
        print("Reusing existing generated data.")

    rx = e1.rx_ecef()
    names = list(rx.keys())
    ref = names[0]
    em_true = e1.emitter_ecef()

    # 1) GNSS time transfer per node.
    print("\n=== GNSS code-phase time transfer ===")
    print(f"{'node':<12}{'recovered Δt':>15}{'configured Δt':>16}"
          f"{'error (ns)':>12}{'PRN spread(ns)':>16}")
    dt_hat = {}
    for name in names:
        med, ests, quals = recover_clock(name, DATA)
        cfg_off = CLOCKS[name][0] + CLOCKS[name][1] * 0.5  # offset + drift·mid(1s)
        dt_hat[name] = med
        print(f"{name:<12}{med*1e9:>13.1f}ns{cfg_off*1e9:>14.1f}ns"
              f"{(med-cfg_off)*1e9:>12.1f}{ests.std()*1e9:>16.1f}")

    # 2) VHF TDOA: raw, GNSS-corrected, and the ideal-clock reference.
    iqv = {n: load_iq(os.path.join(DATA, f"{n}__vhf.iq.c64"))[0] for n in names}
    fsv = 8_000_000.0
    max_lag = int(150e3 / C * fsv * 1.5)
    raw_rd, corr_rd = np.zeros(len(names)), np.zeros(len(names))
    for k, name in enumerate(names):
        rd, _q = fft_xcorr_tdoa(iqv[name], iqv[ref], fsv, max_lag)
        raw_rd[k] = rd
        # remove the differential clock offset recovered from GNSS
        corr_rd[k] = rd - (dt_hat[name] - dt_hat[ref]) * C

    def solve(rd, label):
        est, rms = multilaterate_tdoa([rx[n] for n in names], rd, ref=0)
        _, _, _, horiz, vert = enu_error(est, em_true)
        print(f"  {label:<26} horizontal {horiz:8.1f} m   "
              f"vertical {vert:+10.1f} m   (resid {rms:5.1f} m)")
        return horiz

    print("\n=== VHF TDOA geolocation ===")
    solve(raw_rd, "raw (unsynchronised)")
    solve(corr_rd, "GNSS-disciplined")


if __name__ == "__main__":
    main()
