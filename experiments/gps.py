"""GPS L1 C/A helpers for the GNSS time-transfer experiment.

Self-contained: a Python port of lipika's C/A Gold-code generator, plus a
code-phase time-transfer that recovers a receiver's clock offset from its GNSS
band. The trick that keeps this simple and unbiased is the truth-export: lipika
hands us the exact geometric delay (Sagnac-correct, real orbit) baked into the
signal, so we predict each satellite's arrival and Doppler from truth and the
*only* unmodelled term left is the receiver clock — which is what we measure.
"""

import numpy as np

CHIP_RATE = 1_023_000.0
CODE_LEN = 1023

# IS-GPS-200 Table 3-Ia G2 phase-tap pairs (1-indexed), PRN 1..32.
_G2_TAPS = [(2, 6), (3, 7), (4, 8), (5, 9), (1, 9), (2, 10), (1, 8), (2, 9),
            (3, 10), (2, 3), (3, 4), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10),
            (1, 4), (2, 5), (3, 6), (4, 7), (5, 8), (6, 9), (1, 3), (4, 6),
            (5, 7), (6, 8), (7, 9), (8, 10), (1, 6), (2, 7), (3, 8), (4, 9)]


def _g1():
    reg = [1] * 10
    out = []
    for _ in range(CODE_LEN):
        out.append(reg[9])
        fb = reg[2] ^ reg[9]
        reg = [fb] + reg[:9]
    return out


def _g2(a, b):
    reg = [1] * 10
    out = []
    for _ in range(CODE_LEN):
        out.append(reg[a - 1] ^ reg[b - 1])
        fb = reg[1] ^ reg[2] ^ reg[5] ^ reg[7] ^ reg[8] ^ reg[9]
        reg = [fb] + reg[:9]
    return out


def ca_code(prn):
    """1023-chip C/A sequence for `prn` as +1/-1 (float32)."""
    g1 = _g1()
    a, b = _G2_TAPS[prn - 1]
    g2 = _g2(a, b)
    return np.array([1 - 2 * (x ^ y) for x, y in zip(g1, g2)], dtype=np.float32)


def time_transfer_offset(iq, fs, epoch, t_grid, delay_grid, carrier_hz, prn,
                         sv_epoch, block_ms=1.0, max_seconds=None):
    """Recover the receiver clock offset (s) from one PRN's C/A in the GNSS band.

    For each ~1 ms block: wipe the predicted carrier Doppler (from the truth
    delay), build the C/A replica at the predicted code phase (Δt = 0), and
    correlate. The replica assumes a perfect clock, so the residual code shift
    of the correlation peak is the receiver clock offset Δt. Blocks are summed
    non-coherently (robust to NAV-bit flips and the residual frequency offset).
    Returns (offset_s, peak_quality)."""
    code = ca_code(prn)
    block_n = int(round(fs * block_ms * 1e-3))  # samples per block
    n_total = len(iq)
    if max_seconds is not None:
        n_total = min(n_total, int(fs * max_seconds))
    n_blocks = n_total // block_n
    lags = 8  # ± samples to search (Δt is well under a chip)
    acc = np.zeros(2 * lags + 1)
    tau_g = np.asarray(delay_grid)
    tg = np.asarray(t_grid)
    two_pi_fc = 2.0 * np.pi * carrier_hz
    for b in range(n_blocks):
        s = b * block_n
        idx = s + np.arange(block_n)
        t_n = epoch + idx / fs
        tau_n = np.interp(t_n, tg, tau_g)
        # carrier-Doppler wipe (bring code to a stationary phase within the block)
        rx_w = iq[s:s + block_n] * np.exp(1j * two_pi_fc * tau_n)
        # replica at predicted code phase
        chip = (t_n - tau_n - sv_epoch) * CHIP_RATE
        replica = code[np.mod(np.floor(chip).astype(np.int64), CODE_LEN)]
        # circular correlation, small-lag window
        c = np.fft.ifft(np.fft.fft(rx_w) * np.conj(np.fft.fft(replica)))
        mag = np.abs(c)
        window = np.concatenate([mag[:lags + 1], mag[-lags:]])  # lags 0..L, -L..-1
        acc += window
    # acc index 0..lags = lag 0..+L ; lags+1..2L = lag -L..-1
    order = np.concatenate([np.arange(0, lags + 1), np.arange(-lags, 0)])
    peak = np.argmax(acc)
    pk_lag = order[peak]
    # parabolic interpolation around the integer peak in lag space
    reidx = {l: i for i, l in enumerate(order)}
    a0 = acc[reidx[pk_lag]]
    am = acc[reidx.get(pk_lag - 1, peak)]
    ap = acc[reidx.get(pk_lag + 1, peak)]
    denom = (am - 2 * a0 + ap)
    delta = 0.5 * (am - ap) / denom if denom != 0 else 0.0
    lag = pk_lag + delta
    quality = a0 / (np.median(acc) + 1e-30)
    return lag / fs, quality
