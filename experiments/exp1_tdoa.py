#!/usr/bin/env python3
"""Experiment 1 — TDOA multilateration of a fixed emitter.

The headline distributed-receiver capability: several ground receivers observe
one wideband emitter; pairwise cross-correlation of their baseband IQ yields
time-difference-of-arrival, which multilaterates to the emitter position.

Because the scene is stationary, the ground truth is pure ECEF geometry, so the
error we measure is the *estimator's* (correlation resolution + geometry), not
the data's. Two runs:

  * ideal  — perfect clocks: measured TDOA == geometric TDOA. Validates the
             signal geometry and bounds the achievable accuracy.
  * clk    — each node carries an unknown, distinct clock offset. TDOA is biased
             by the differential offset (c·Δt): this is exactly the synchronisation
             problem a distributed network must solve, here quantified.

Data is large, so IQ is written to the mounted drive.
"""

import os
import numpy as np
from common import (geodetic_to_ecef, ecef_to_geodetic, load_iq, run_lipika,
                    fft_xcorr_tdoa, multilaterate_tdoa, enu_error, C)

REPO = "/home/mike/Projects/DGSN_RESSURECT/lipika"
DATA = "/media/mike/36ada8f0-4d34-4343-ba81-487bb65bb372/lipika-test/exp_tdoa"
BUILD = os.path.join(REPO, "experiments", "build")

CARRIER = 137_000_000.0
SAMPLE_RATE = 8_000_000.0
BIT_RATE = 2_000_000.0
DURATION = 1.0
EPOCH = 1_717_260_840.0
TX_POWER_DBW = 0.0  # 1 W EIRP — strong, this experiment tests geometry not sensitivity

# Distributed ground receivers across central Scotland (lat, lon, height_m).
SITES = {
    "glasgow":   (55.86, -4.25, 30.0),
    "edinburgh": (55.95, -3.19, 60.0),
    "stirling":  (56.12, -3.94, 40.0),
    "ayr":       (55.46, -4.63, 20.0),
    "perth":     (56.40, -3.43, 25.0),
    "dumfries":  (55.07, -3.61, 50.0),
}
# Airborne emitter (10 km) over a point inside the receiver hull — elevation
# gives the geometry some vertical leverage.
EMITTER_GEO = (55.67, -3.78, 10_000.0)

SYNC = [1, -1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1]


def emitter_ecef():
    return geodetic_to_ecef(*EMITTER_GEO)


def rx_ecef():
    return {name: geodetic_to_ecef(*geo) for name, geo in SITES.items()}


def write_scene(path, clocks=None):
    """clocks: optional {name: (init_offset_s, freq_offset_y0, allan_1s)}."""
    rx = rx_ecef()
    em = emitter_ecef()
    payload_bits = int(BIT_RATE * DURATION)
    lines = [
        "[scene]",
        f"epoch_unix = {EPOCH}",
        f"duration_s = {DURATION}",
        "",
        "[[transmitter]]",
        'name = "emitter"',
        f"carrier_hz = {CARRIER}",
        f"tx_power_dbw = {TX_POWER_DBW}",
        'kind = "bpsk"',
        f"bit_rate_bps = {BIT_RATE}",
        f"payload_bits = {payload_bits}",
        "payload_seed = 20260601",
        f"sync_word = {SYNC}",
        f"trajectory = {{ kind = \"stationary_ecef\", x_m = {em[0]:.4f}, "
        f"y_m = {em[1]:.4f}, z_m = {em[2]:.4f} }}",
        "",
    ]
    for name, p in rx.items():
        lines += [
            "[[receiver]]",
            f'name = "{name}"',
            f"trajectory = {{ kind = \"stationary_ecef\", x_m = {p[0]:.4f}, "
            f"y_m = {p[1]:.4f}, z_m = {p[2]:.4f} }}",
        ]
        if clocks and name in clocks:
            off, fo, ad = clocks[name]
            lines.append(
                f"clock = {{ initial_time_offset_s = {off}, "
                f"frequency_offset_y0 = {fo}, allan_deviation_1s = {ad} }}")
        lines += [
            "[[receiver.band]]",
            'band_id = "vhf"',
            f"lo_carrier_hz = {CARRIER}",
            f"sample_rate_hz = {SAMPLE_RATE}",
            "noise_temperature_k = 290.0",
            "",
        ]
    with open(path, "w") as f:
        f.write("\n".join(lines))


def analyze(out_dir, label, clocks=None):
    rx = rx_ecef()
    names = list(rx.keys())
    em_true = emitter_ecef()
    ref = names[0]

    iq = {}
    for name in names:
        a, meta = load_iq(os.path.join(out_dir, f"{name}__vhf.iq.c64"))
        iq[name] = a
        fs = meta["sample_rate_hz"]

    max_lag = int(150e3 / C * fs * 1.5)  # ~150 km baseline budget
    range_diff = np.zeros(len(names))
    print(f"\n=== {label} ===")
    print(f"{'pair':<22}{'meas Δr (m)':>14}{'true Δr (m)':>14}"
          f"{'err (m)':>10}{'quality':>10}")
    truth_rd = {}
    for k, name in enumerate(names):
        rd, q = fft_xcorr_tdoa(iq[name], iq[ref], fs, max_lag)
        range_diff[k] = rd
        rt = np.linalg.norm(em_true - rx[name]) - np.linalg.norm(em_true - rx[ref])
        truth_rd[name] = rt
        # expected clock bias on this pair (for the clk run)
        bias = 0.0
        if clocks:
            off_k = clocks.get(name, (0, 0, 0))[0]
            off_r = clocks.get(ref, (0, 0, 0))[0]
            bias = (off_k - off_r) * C
        tag = f"{name}-{ref}"
        print(f"{tag:<22}{rd:>14.2f}{rt:>14.2f}{rd-rt:>10.2f}{q:>10.1f}"
              + (f"   (clk bias {bias:+.1f} m)" if clocks else ""))

    est, rms = multilaterate_tdoa([rx[n] for n in names], range_diff, ref=0)
    err = np.linalg.norm(est - em_true)
    e, n, u, horiz, vert = enu_error(est, em_true)
    lat, lon, h = ecef_to_geodetic(est)
    tlat, tlon, th = ecef_to_geodetic(em_true)
    print(f"  fit residual RMS:  {rms:8.2f} m")
    print(f"  emitter truth:     lat {tlat:.5f}  lon {tlon:.5f}  h {th:8.1f} m")
    print(f"  emitter estimate:  lat {lat:.5f}  lon {lon:.5f}  h {h:8.1f} m")
    print(f"  ENU error:         E {e:+8.1f}  N {n:+8.1f}  U {u:+8.1f} m")
    print(f"  horizontal error:  {horiz:8.1f} m   vertical error: {vert:+8.1f} m")
    print(f"  3D position error: {err:8.1f} m")
    return est, err, rms


def main():
    os.makedirs(BUILD, exist_ok=True)
    # Distinct, unknown per-node clock offsets (tens-hundreds of ns) + realistic
    # OCXO stability. Over 1 s the walk/drift are negligible; the static offset
    # is the dominant, unsynchronised term.
    clocks = {
        "glasgow":   (0.0,      0.0,   1e-11),
        "edinburgh": (1.30e-7,  1e-9,  1e-11),
        "stirling":  (-2.10e-7, -5e-10, 1e-11),
        "ayr":       (3.40e-7,  2e-9,  1e-11),
        "perth":     (-1.55e-7, 8e-10, 1e-11),
        "dumfries":  (2.05e-7, -1e-9,  1e-11),
    }

    ideal_cfg = os.path.join(BUILD, "tdoa_ideal.toml")
    clk_cfg = os.path.join(BUILD, "tdoa_clk.toml")
    write_scene(ideal_cfg, clocks=None)
    write_scene(clk_cfg, clocks=clocks)

    ideal_out = os.path.join(DATA, "ideal")
    clk_out = os.path.join(DATA, "clk")

    print("Generating ideal-clock scene...")
    run_lipika(REPO, ideal_cfg, ideal_out)
    print("Generating clock-offset scene...")
    run_lipika(REPO, clk_cfg, clk_out)

    analyze(ideal_out, "IDEAL CLOCKS (geometry-limited)")
    analyze(clk_out, "UNSYNCHRONISED CLOCKS (differential offset bias)",
            clocks=clocks)


if __name__ == "__main__":
    main()
