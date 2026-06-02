#!/usr/bin/env python3
"""Experiment 4 — VHF carrier-phase interferometry (coherent geolocation).

The payoff of the arc: use the carrier *phase* (not the code envelope) as the
cross-node observable. At VHF (137 MHz, λ≈2.2 m) the phase difference across a
baseline is ~1000× finer than Exp 1's envelope TDOA — the path from tens of
metres to centimetres. This is the ideal-phase version (lipika models perfect
isotropic antennas), so it isolates geometry + timing + ambiguity from hardware
phase errors.

Scene = the Exp 3 scene (6 nodes, GNSS L1 + VHF bands, per-node clocks) plus one
stationary target emitting two co-located signals: a BPSK at 137 MHz (coarse
envelope-TDOA prior) and a CW tone at 140 MHz (the carrier-phase observable, a
+3 MHz tone in the VHF band, clear of the BPSK main lobe).

Pipeline: recover each node's clock from GNSS (Exp 3 method) → measure the CW
carrier phase per node → remove the clock phase → coarse position from the BPSK
→ coherent position-domain search (distributed beamformer). Run with both the
TRUE clocks (validates the method + lipika's phase fidelity → ~cm) and the
GNSS-recovered clocks (the real-world test — does ~few-ns code sync beat the
λ=2.2 m ambiguity sidelobes?).
"""

import os
import json
import numpy as np
import exp1_tdoa as e1
from common import (load_iq, run_lipika, fft_xcorr_tdoa, multilaterate_tdoa,
                    enu_error, tone_phase, coherent_geolocate,
                    geodetic_to_ecef, ecef_to_geodetic, C)
from gps import time_transfer_offset

REPO = e1.REPO
DATA = "/media/mike/36ada8f0-4d34-4343-ba81-487bb65bb372/lipika-test/exp_interf"
BUILD = e1.BUILD
EPOCH = e1.EPOCH
DURATION = 2.0
GPS_PRNS = [2, 6, 12, 19, 24, 25]
RNX = "examples/gps_synth.rnx"
L1 = 1_575_420_000.0
VHF_LO = 137_000_000.0
BPSK_HZ = 137_000_000.0
CW_HZ = 140_000_000.0          # +3 MHz tone in the VHF band
CW_BB = CW_HZ - VHF_LO          # 3 MHz baseband
CLOCKS = e1.__dict__.get("CLOCKS", None) or {
    "glasgow": (0.0, 0.0, 1e-11), "edinburgh": (1.30e-7, 1e-9, 1e-11),
    "stirling": (-2.10e-7, -5e-10, 1e-11), "ayr": (3.40e-7, 2e-9, 1e-11),
    "perth": (-1.55e-7, 8e-10, 1e-11), "dumfries": (2.05e-7, -1e-9, 1e-11),
}


def write_scene(path):
    rx = e1.rx_ecef()
    em = e1.emitter_ecef()
    payload_bits = int(e1.BIT_RATE * DURATION)
    L = ["[scene]", f"epoch_unix = {EPOCH}", f"duration_s = {DURATION}", ""]
    # Target: co-located BPSK (coarse) + CW (carrier phase).
    trj = (f'trajectory = {{ kind = "stationary_ecef", x_m = {em[0]:.4f}, '
           f"y_m = {em[1]:.4f}, z_m = {em[2]:.4f} }}")
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
              "sample_rate_hz = 8192000.0", "noise_temperature_k = 290.0",
              "[[receiver.band]]", 'band_id = "vhf"', f"lo_carrier_hz = {VHF_LO}",
              "sample_rate_hz = 8000000.0", "noise_temperature_k = 290.0", ""]
    with open(path, "w") as f:
        f.write("\n".join(L))


def set_height(ecef, h_m):
    """Return ECEF at the same lat/lon but height h_m (fix the bad coarse
    vertical; the regional array is vertical-GDOP-limited)."""
    lat, lon, _ = ecef_to_geodetic(ecef)
    return geodetic_to_ecef(lat, lon, h_m)


def bpsk_only(iq, fs, cutoff=2.5e6):
    """Brick-wall low-pass to drop the +3 MHz CW tone before envelope TDOA."""
    n = len(iq)
    f = np.fft.fftfreq(n, 1.0 / fs)
    spec = np.fft.fft(iq)
    spec[np.abs(f) > cutoff] = 0.0
    return np.fft.ifft(spec)


def recover_clock(node):
    iq, meta = load_iq(os.path.join(DATA, f"{node}__gnss_l1.iq.c64"))
    fs = meta["sample_rate_hz"]
    truth = json.load(open(os.path.join(DATA, f"{node}__truth.json")))
    tg = truth["t_utc"]
    by = {t["name"]: t for t in truth["transmitters"]}
    # Integrate over the whole record so the recovered offset references the
    # window midpoint — matching the CW phase reference.
    ests = [time_transfer_offset(iq, fs, EPOCH, tg, by[f"gps_prn{prn:02d}"]["delay_s"],
                                 L1, prn, sv_epoch=EPOCH, block_ms=1.0, max_seconds=DURATION)[0]
            for prn in GPS_PRNS]
    return float(np.median(ests))


def main():
    os.makedirs(BUILD, exist_ok=True)
    cfg = os.path.join(BUILD, "interf.toml")
    write_scene(cfg)
    if not os.path.exists(os.path.join(DATA, "glasgow__truth.json")):
        print("Generating interferometry scene (to the drive)...")
        run_lipika(REPO, cfg, DATA, threads=0, chunk=2_000_000, truth_hz=50.0)
    else:
        print("Reusing existing generated data.")

    rx = e1.rx_ecef()
    names = list(rx.keys())
    em_true = e1.emitter_ecef()
    rx_list = [rx[n] for n in names]

    # --- per-node CW carrier phase + GNSS clock ---
    vhf = {n: load_iq(os.path.join(DATA, f"{n}__vhf.iq.c64")) for n in names}
    fsv = vhf[names[0]][1]["sample_rate_hz"]
    phi = {}
    for n in names:
        ph, fe, amp = tone_phase(vhf[n][0], fsv, CW_BB)
        phi[n] = ph
    t_mid = DURATION / 2.0
    dt_gnss = {n: recover_clock(n) for n in names}
    dt_true = {n: CLOCKS[n][0] + CLOCKS[n][1] * t_mid for n in names}  # Δt at midpoint

    # --- raw BPSK envelope range-diffs (CW tone filtered out) ---
    ref = names[0]
    b_ref = bpsk_only(vhf[ref][0], fsv)
    rd = np.zeros(len(names))
    for k, n in enumerate(names):
        rd[k], _ = fft_xcorr_tdoa(bpsk_only(vhf[n][0], fsv), b_ref,
                                  fsv, int(150e3 / C * fsv * 1.5))
    _, _, th_h = ecef_to_geodetic(em_true)
    f_cw = CW_HZ

    def solve(dt, label):
        # coarse prior: clock-correct the envelope range-diffs (Exp 3 style),
        # multilaterate, then fix the GDOP-blind vertical at truth height.
        rd_c = np.array([rd[k] - (dt[n] - dt[ref]) * C for k, n in enumerate(names)])
        coarse, _ = multilaterate_tdoa(rx_list, rd_c, ref=0)
        x0 = set_height(coarse, th_h)
        _, _, _, ch, _ = enu_error(x0, em_true)
        psi = np.array([phi[n] + 2 * np.pi * f_cw * dt[n] for n in names])
        # coarse-to-fine coherent search around the prior, at truth height
        x1, _, *_ = coherent_geolocate(rx_list, psi, f_cw, x0, half_m=45.0, step_m=0.25)
        xh, jh, _, _, _, _ = coherent_geolocate(rx_list, psi, f_cw, x1, half_m=0.6, step_m=0.004)
        e, n_, u, horiz, vert = enu_error(xh, em_true)
        print(f"  {label:<18} coarse {ch:5.1f} m → carrier-phase "
              f"horizontal {horiz*100:7.1f} cm  (E {e*100:+.1f}, N {n_*100:+.1f} cm)  "
              f"coh {jh:.3f}")
        return horiz

    # --- phase plumbing cross-check (truth-clock ψ vs truth-export geometry) ---
    tr0 = json.load(open(os.path.join(DATA, f"{ref}__truth.json")))
    cwname = "tgt_cw"
    by = {t["name"]: t for t in tr0["transmitters"]}
    print("\nphase cross-check (clock-corrected ψ vs −2π·f·delay_truth, ref node):")
    for n in [ref, names[1]]:
        trn = json.load(open(os.path.join(DATA, f"{n}__truth.json")))
        byn = {t["name"]: t for t in trn["transmitters"]}
        dly = byn[cwname]["delay_s"]
        d_mid = dly[len(dly) // 2]  # delay at the record midpoint (constant here)
        psi_n = phi[n] + 2 * np.pi * f_cw * dt_true[n]
        pred = (-2 * np.pi * f_cw * d_mid)
        resid = np.angle(np.exp(1j * (psi_n - pred)))  # wrapped
        print(f"  {n:<10} residual {np.degrees(resid):+7.1f} deg "
              f"({resid/(2*np.pi)/f_cw*1e9:+.2f} ns)")

    print("\n=== Carrier-phase coherent geolocation (VHF, λ=2.19 m) ===")
    solve(dt_true, "ideal clocks")
    solve(dt_gnss, "GNSS-disciplined")
    print("\n(Exp 1 envelope-only TDOA was ~33 m horizontal for comparison.)")


if __name__ == "__main__":
    main()
