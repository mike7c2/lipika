"""Shared helpers for lipika-data experiments.

Everything here reads the raw `complex64 + .json` sidecar that lipika writes,
and works in the WGS84 ECEF frame. Ground truth for stationary scenes is pure
geometry, so experiments use `stationary_ecef` entries whose coordinates are
computed *here* — lipika consumes them verbatim, so the truth positions and the
positions in the signal are identical by construction (no conversion ambiguity).
"""

import json
import os
import subprocess
import numpy as np

C = 299_792_458.0  # m/s, matches lipika's SPEED_OF_LIGHT_M_S

# WGS84
_A = 6_378_137.0
_F = 1.0 / 298.257_223_563
_E2 = _F * (2.0 - _F)


def geodetic_to_ecef(lat_deg, lon_deg, h_m):
    """WGS84 geodetic -> ECEF (metres). Used only to place nodes at realistic
    coordinates; consistency with lipika is guaranteed because we pass the
    resulting ECEF straight into the scene as `stationary_ecef`."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    n = _A / np.sqrt(1.0 - _E2 * np.sin(lat) ** 2)
    x = (n + h_m) * np.cos(lat) * np.cos(lon)
    y = (n + h_m) * np.cos(lat) * np.sin(lon)
    z = (n * (1.0 - _E2) + h_m) * np.sin(lat)
    return np.array([x, y, z])


def ecef_to_geodetic(p):
    """ECEF -> (lat_deg, lon_deg, h_m), Bowring's method (for reporting)."""
    x, y, z = p
    lon = np.arctan2(y, x)
    r = np.hypot(x, y)
    b = _A * (1.0 - _F)
    ep2 = (_A**2 - b**2) / b**2
    th = np.arctan2(_A * z, b * r)
    lat = np.arctan2(z + ep2 * b * np.sin(th) ** 3,
                     r - _E2 * _A * np.cos(th) ** 3)
    n = _A / np.sqrt(1.0 - _E2 * np.sin(lat) ** 2)
    h = r / np.cos(lat) - n
    return np.degrees(lat), np.degrees(lon), h


def load_iq(path):
    """Load a lipika band: returns (iq complex64 array, sidecar dict)."""
    iq = np.fromfile(path, dtype=np.complex64)
    with open(path.replace(".iq.c64", ".iq.json")) as f:
        meta = json.load(f)
    return iq, meta


def run_lipika(repo, config_path, out_dir, threads=0, chunk=2_000_000, truth_hz=0.0):
    """Invoke the release lipika binary on a scene. Streams to `out_dir`.
    truth_hz > 0 also writes per-receiver geometric truth sidecars."""
    binary = os.path.join(repo, "target", "release", "lipika")
    cmd = [binary, "--config", config_path, "--output", out_dir,
           "--threads", str(threads), "--chunk-samples", str(chunk)]
    if truth_hz > 0:
        cmd += ["--truth-hz", str(truth_hz)]
    subprocess.run(cmd, check=True, cwd=repo)


def fft_xcorr_tdoa(a, b, fs, max_lag_samples):
    """Cross-correlate complex baseband `a` against reference `b`.

    Returns (range_diff_m, peak_quality). `a[n] ≈ b[n - D]` (a delayed by D
    samples relative to b) peaks at lag +D, so range_diff = D·c/fs estimates
    r_a − r_b. Sub-sample peak via quadratic interpolation on |corr|.
    Circular FFT correlation is exact for |lag| ≪ N, which holds here
    (max range difference ⇒ a few thousand samples ≪ millions)."""
    n = min(len(a), len(b))
    a = a[:n]
    b = b[:n]
    fa = np.fft.fft(a)
    fb = np.fft.fft(b)
    corr = np.fft.ifft(fa * np.conj(fb))
    mag = np.abs(corr)
    # Candidate lags: 0..max in corr[0..K], negative in corr[N-K..N-1].
    k = int(max_lag_samples)
    idx = np.r_[0:k + 1, n - k:n]
    sub = mag[idx]
    pk = np.argmax(sub)
    peak_idx = idx[pk]
    # quadratic interpolation around the integer peak (wrap-safe)
    im1 = mag[(peak_idx - 1) % n]
    i0 = mag[peak_idx]
    ip1 = mag[(peak_idx + 1) % n]
    denom = (im1 - 2 * i0 + ip1)
    delta = 0.5 * (im1 - ip1) / denom if denom != 0 else 0.0
    lag = peak_idx if peak_idx <= n // 2 else peak_idx - n
    lag = lag + delta
    # peak quality: peak height over median sidelobe
    quality = i0 / (np.median(mag) + 1e-30)
    return lag / fs * C, quality


def _tdoa_residual(x, rx, others, ref, range_diff):
    dref = x - rx[ref]
    rref = np.linalg.norm(dref)
    res, rows = [], []
    for k in others:
        dk = x - rx[k]
        rk = np.linalg.norm(dk)
        res.append(range_diff[k] - (rk - rref))
        rows.append(dk / rk - dref / rref)
    return np.array(res), np.array(rows)


def multilaterate_tdoa(rx_ecef, range_diff, ref=0, x0=None, iters=200):
    """Solve emitter ECEF from TDOA range-differences relative to rx `ref`.

    range_diff[k] = |x − rx_k| − |x − rx_ref|. Levenberg-Marquardt least
    squares — damping is essential here because a regional ground array has
    near-singular vertical geometry, which plain Gauss-Newton sends to
    infinity. Returns (x_ecef, rms_residual_m)."""
    rx = np.asarray(rx_ecef, float)
    others = [k for k in range(len(rx)) if k != ref]
    x = rx.mean(axis=0) if x0 is None else np.array(x0, float)
    res, jac = _tdoa_residual(x, rx, others, ref, range_diff)
    cost = float(res @ res)
    lam = 1e-3
    for _ in range(iters):
        a = jac.T @ jac
        g = jac.T @ res
        step = np.linalg.solve(a + lam * np.diag(np.diag(a) + 1e-9), g)
        xn = x + step
        rn, jn = _tdoa_residual(xn, rx, others, ref, range_diff)
        cn = float(rn @ rn)
        if cn < cost:
            x, res, jac, cost = xn, rn, jn, cn
            lam = max(lam * 0.5, 1e-12)
            if np.linalg.norm(step) < 1e-4:
                break
        else:
            lam *= 4.0
            if lam > 1e12:
                break
    rms = float(np.sqrt(np.mean(np.square(res))))
    return x, rms


def enu_basis(p_ecef):
    """Local East/North/Up unit vectors (ECEF columns) at an ECEF point."""
    lat, lon, _ = ecef_to_geodetic(p_ecef)
    lat, lon = np.radians(lat), np.radians(lon)
    e = np.array([-np.sin(lon), np.cos(lon), 0.0])
    n = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    u = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    return e, n, u


def position_from_ranges(rx_ecef, ranges, x0, iters=80):
    """Least-squares source position from absolute ranges rᵢ = |x − rxᵢ|
    (Levenberg-Marquardt). Used to reconstruct the truth trajectory from the
    truth-export ranges and to seed the orbit fit. Returns (x_ecef, rms_m)."""
    rx = np.asarray(rx_ecef, float)
    r = np.asarray(ranges, float)
    x = np.array(x0, float)
    lam = 1e-2
    def res_jac(x):
        d = x[None, :] - rx
        rr = np.linalg.norm(d, axis=1)
        return r - rr, -d / rr[:, None]
    res, jac = res_jac(x)
    cost = res @ res
    for _ in range(iters):
        a = jac.T @ jac
        step = np.linalg.solve(a + lam * np.diag(np.diag(a) + 1e-9), jac.T @ res)
        xn = x + step
        rn, jn = res_jac(xn)
        if rn @ rn < cost:
            x, res, jac, cost = xn, rn, jn, rn @ rn
            lam = max(lam * 0.5, 1e-12)
            if np.linalg.norm(step) < 1e-4:
                break
        else:
            lam *= 4.0
    return x, float(np.sqrt(np.mean(res**2)))


def elevation_deg(rx_ecef, src_ecef):
    """Elevation angle (deg) of `src` above the local horizon at `rx`."""
    rx = np.asarray(rx_ecef, float)
    up = rx / np.linalg.norm(rx)
    los = np.asarray(src_ecef, float) - rx
    los /= np.linalg.norm(los)
    return float(np.degrees(np.arcsin(np.clip(los @ up, -1, 1))))


def enu_error(est_ecef, truth_ecef):
    """Decompose (est − truth) into local East/North/Up at truth.
    Returns (east_m, north_m, up_m, horizontal_m, vertical_m)."""
    e, n, u = enu_basis(truth_ecef)
    d = np.asarray(est_ecef) - np.asarray(truth_ecef)
    east, north, up = d @ e, d @ n, d @ u
    return east, north, up, float(np.hypot(east, north)), float(up)


def tone_phase(iq, fs, f_baseband, search_hz=3000.0):
    """Phase (rad) of a complex tone near `f_baseband`, referenced to the
    record's **midpoint** in time.

    Refines the tone frequency by an FFT peak + parabolic interpolation, then
    de-rotates and averages. Referencing to the midpoint (centred time axis)
    makes the recovered phase insensitive *to first order* to a residual
    frequency error — a small f-error otherwise drifts the phase linearly across
    the window and biases a t=0 reference by ~π·Δf·T. The midpoint is therefore
    the natural common reference; the caller must evaluate clock Δt and geometric
    delay at the same midpoint time. Returns (phase_rad, f_est_hz, amplitude)."""
    n = len(iq)
    t = (np.arange(n) - n / 2.0) / fs  # centred on the record midpoint
    spec = np.fft.fft(iq)
    freqs = np.fft.fftfreq(n, 1.0 / fs)
    mask = np.abs(freqs - f_baseband) < search_hz
    cand = np.where(mask)[0]
    pk = cand[np.argmax(np.abs(spec[cand]))]
    a0, am, ap = np.abs(spec[pk]), np.abs(spec[pk - 1]), np.abs(spec[(pk + 1) % n])
    denom = am - 2 * a0 + ap
    delta = 0.5 * (am - ap) / denom if denom != 0 else 0.0
    f_est = freqs[pk] + delta * (fs / n)
    z = np.mean(iq * np.exp(-1j * 2 * np.pi * f_est * t))
    return float(np.angle(z)), float(f_est), float(np.abs(z))


def track_carrier_phase(iq, fs, f_ramp, t_grid, range_pred, f_cw, decim_hz=50.0):
    """Continuously track a moving emitter's CW carrier phase.

    De-rotates each ~1/decim_hz block by the predicted baseband model
    (the `f_ramp = f_cw − LO` offset plus the predicted carrier `−k·range_pred`),
    block-averages, and unwraps. With a good `range_pred` the residual is slowly
    varying (just the receiver clock), so the block phase is clean and the unwrap
    is unambiguous — turning the moving emitter's phase into a continuous,
    integer-free observable. Returns (t_dec, residual_phase_rad) where the
    residual is `−k·(range_true − range_pred) − 2π·f_cw·Δt` (k = 2π·f_cw/c)."""
    n = len(iq)
    k = 2.0 * np.pi * f_cw / C
    block = int(round(fs / decim_hz))
    nb = n // block
    t_dec = np.zeros(nb)
    ph = np.zeros(nb)
    for b in range(nb):
        s = b * block
        t = (s + np.arange(block)) / fs
        rp = np.interp(t, t_grid, range_pred)
        derot = iq[s:s + block] * np.exp(1j * (k * rp - 2.0 * np.pi * f_ramp * t))
        z = np.mean(derot)
        ph[b] = np.angle(z)
        t_dec[b] = t.mean()
    return t_dec, np.unwrap(ph)


def estimate_doppler_track(iq, fs, f_ramp, win=8192, hop=4096):
    """Data-driven Doppler(t) of a CW tone (no truth). De-rotates the known
    `f_ramp` offset, then takes the FFT peak (parabolic-interpolated) in each
    window → a coarse instantaneous-frequency track. Returns (t_centres,
    doppler_hz). A smooth polynomial fit of this is precise enough to de-rotate
    the carrier for continuous phase tracking."""
    n = len(iq)
    t = np.arange(n) / fs
    base = iq * np.exp(-1j * 2 * np.pi * f_ramp * t)
    starts = np.arange(0, n - win, hop)
    tc = np.zeros(len(starts))
    fd = np.zeros(len(starts))
    freqs = np.fft.fftfreq(win, 1.0 / fs)
    for j, s in enumerate(starts):
        seg = base[s:s + win]
        sp = np.abs(np.fft.fft(seg))
        pk = np.argmax(sp)
        a0, am, ap = sp[pk], sp[pk - 1], sp[(pk + 1) % win]
        d = 0.5 * (am - ap) / (am - 2 * a0 + ap) if (am - 2 * a0 + ap) != 0 else 0.0
        fd[j] = freqs[pk] + d * (fs / win)
        tc[j] = (s + win / 2) / fs
    return tc, fd


def propagate_twobody_ecef(r0, v0, t, mu=3.986004418e14, omega=7.292_115_146_7e-5):
    """Propagate an ECEF state [r0, v0] (at t[0]) under two-body gravity plus the
    rotating-frame (Coriolis + centrifugal) terms, by RK4 on the uniform grid
    `t`. Returns positions (len(t), 3). The dynamics constraint is what makes the
    satellite's range observable to a tight (radially-blind) ground array."""
    w = np.array([0.0, 0.0, omega])
    def acc(r, v):
        rn = np.linalg.norm(r)
        return -mu * r / rn**3 - 2 * np.cross(w, v) - np.cross(w, np.cross(w, r))
    out = np.zeros((len(t), 3))
    r, v = np.array(r0, float), np.array(v0, float)
    out[0] = r
    for i in range(1, len(t)):
        h = t[i] - t[i - 1]
        k1v, k1a = v, acc(r, v)
        k2v, k2a = v + 0.5 * h * k1a, acc(r + 0.5 * h * k1v, v + 0.5 * h * k1a)
        k3v, k3a = v + 0.5 * h * k2a, acc(r + 0.5 * h * k2v, v + 0.5 * h * k2a)
        k4v, k4a = v + h * k3a, acc(r + h * k3v, v + h * k3a)
        r = r + h / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        v = v + h / 6.0 * (k1a + 2 * k2a + 2 * k3a + k4a)
        out[i] = r
    return out


def sagnac_ranges(pts, rx, omega=7.292_115_146_7e-5):
    """Sagnac-corrected ranges from source positions `pts` (M,3) to receiver
    `rx`, matching lipika's emission solve (linearised Earth-rotation term)."""
    d0 = pts - rx[None, :]
    r0 = np.linalg.norm(d0, axis=1)
    wtau = omega * r0 / C
    xr = pts[:, 0] + wtau * pts[:, 1]
    yr = pts[:, 1] - wtau * pts[:, 0]
    return np.sqrt((xr - rx[0])**2 + (yr - rx[1])**2 + (pts[:, 2] - rx[2])**2)


def coherent_geolocate(rx_ecef, psi, f_cw, x0, half_m=8.0, step_m=0.02):
    """Carrier-phase coherent geolocation (distributed beamformer).

    `psi[i]` is node i's clock-corrected carrier phase ≈ −2π·f_cw·|x−rxᵢ|/c
    (+ a global constant). Grid-search the horizontal plane around `x0` (height
    fixed at x0) maximising the coherent sum `|Σ exp(j(psiᵢ + 2π·f_cw·rᵢ/c))|`;
    the free global phase makes the peak occur at the true geometry. Returns
    (x_ecef, coherence_0to1, east_off_m, north_off_m, J_map, grid)."""
    rxs = np.asarray(rx_ecef, float)
    x0 = np.asarray(x0, float)
    e, n, _ = enu_basis(x0)
    k = 2.0 * np.pi * f_cw / C
    omega = 7.292_115_146_7e-5  # Earth rotation rate, matches lipika's Sagnac term
    grid = np.arange(-half_m, half_m + step_m, step_m)
    de, dn = np.meshgrid(grid, grid, indexing="ij")
    # Candidate positions: (G, G, 3)
    pts = x0[None, None, :] + de[..., None] * e[None, None, :] + dn[..., None] * n[None, None, :]
    acc = np.zeros(de.shape, dtype=complex)
    for i in range(len(rxs)):
        # Sagnac-corrected range, matching lipika's emission solve: rotate the
        # source into the reception-frame by Ω·(range/c) (linearised). Skipping
        # this leaves ~tenths-of-a-metre / tens-of-degrees of VHF phase error.
        d0 = pts - rxs[i][None, None, :]
        r0 = np.linalg.norm(d0, axis=2)
        wtau = omega * r0 / C
        xr = pts[..., 0] + wtau * pts[..., 1]
        yr = pts[..., 1] - wtau * pts[..., 0]
        r = np.sqrt((xr - rxs[i, 0])**2 + (yr - rxs[i, 1])**2 + (pts[..., 2] - rxs[i, 2])**2)
        acc += np.exp(1j * (psi[i] + k * r))
    j_map = np.abs(acc) / len(rxs)
    ij = np.unravel_index(np.argmax(j_map), j_map.shape)
    x_best = pts[ij]
    return x_best, float(j_map[ij]), float(de[ij]), float(dn[ij]), j_map, grid
