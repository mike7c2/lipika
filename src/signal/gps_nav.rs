//! GPS L1 LNAV message generation (IS-GPS-200).
//!
//! Builds a valid 50 bps navigation bit stream — TLM/HOW words with correct
//! TOW and parity, and subframes 1–3 carrying a clock model and ephemeris
//! derived from the satellite's own trajectory — so a receiver can compute a
//! position fix. The ephemeris is generated to be *consistent* with the orbit
//! the simulator flies: for a circular-Keplerian source the mapping is exact;
//! for an arbitrary trajectory (e.g. SGP4 from a TLE) it is the osculating
//! element set at `toe`.

use crate::geometry::kepler::CircularKeplerian;
use crate::geometry::Trajectory;

/// Unix time of the GPS epoch (1980-01-06 00:00:00 UTC).
pub const GPS_EPOCH_UNIX: f64 = 315_964_800.0;
pub const SECONDS_PER_WEEK: f64 = 604_800.0;
/// WGS-84 gravitational parameter used by IS-GPS-200 (matches `kepler.rs`).
pub const GPS_MU: f64 = 3.986_005e14;
/// WGS-84 Earth rotation rate used by IS-GPS-200 (matches `kepler.rs`).
pub const GPS_OMEGA_E_DOT: f64 = 7.292_115_146_7e-5;

const PI: f64 = std::f64::consts::PI;

/// GPS week number + seconds-of-week.
#[derive(Debug, Clone, Copy)]
pub struct GpsTime {
    pub week: u32,
    pub sow: f64,
}

impl GpsTime {
    /// Convert a UTC unix time to GPS time. `utc_offset_s` is GPS−UTC (leap
    /// seconds; 18 for the current era).
    pub fn from_unix(unix: f64, utc_offset_s: f64) -> Self {
        let gps = unix - GPS_EPOCH_UNIX + utc_offset_s;
        let week = (gps / SECONDS_PER_WEEK).floor();
        let sow = gps - week * SECONDS_PER_WEEK;
        Self { week: week as u32, sow }
    }
}

/// Round a seconds-of-week value to the broadcast `toe`/`toc` 16 s grid.
fn round_to_16s(sow: f64) -> f64 { (sow / 16.0).round() * 16.0 }

/// Wrap an angle (rad) to [−π, π).
fn wrap_pi(x: f64) -> f64 {
    let twopi = 2.0 * PI;
    let mut r = (x + PI).rem_euclid(twopi) - PI;
    if r >= PI { r -= twopi; }
    r
}

/// Broadcast ephemeris + clock parameters for one satellite, in SI/radians
/// (encoded to the message's semicircle/scaled-integer fields at emit time).
#[derive(Debug, Clone, Copy)]
pub struct Ephemeris {
    pub prn: u8,
    pub week: u32,
    pub toe_s: f64,
    pub toc_s: f64,
    pub sqrt_a: f64,
    pub e: f64,
    pub m0_rad: f64,
    pub delta_n_rad_s: f64,
    pub i0_rad: f64,
    pub idot_rad_s: f64,
    pub omega0_rad: f64,
    pub omega_dot_rad_s: f64,
    pub argp_rad: f64,
    pub cuc: f64,
    pub cus: f64,
    pub crc: f64,
    pub crs: f64,
    pub cic: f64,
    pub cis: f64,
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    pub tgd: f64,
    pub iode: u8,
    pub iodc: u16,
    pub health: u8,
    pub ura_index: u8,
}

impl Ephemeris {
    /// Exact ephemeris for a circular-Keplerian orbit. `toe` is placed at the
    /// scene epoch (`scene_epoch_unix`); the closed-form mapping reproduces the
    /// trajectory's ECEF positions under the IS-GPS-200 SV-position algorithm
    /// (given the shared μ and Ω̇ₑ constants).
    // r[impl fn.4.lnav] exact ephemeris for a circular-Keplerian source
    pub fn from_circular_keplerian(
        orbit: &CircularKeplerian,
        prn: u8,
        scene_epoch_unix: f64,
        gps: GpsTime,
    ) -> Self {
        let a = orbit.semi_major_axis_m;
        let n = (GPS_MU / a.powi(3)).sqrt();
        // `toe` is transmitted with a 16 s LSB, so it must sit on the 16 s grid;
        // M0 has to be referenced to *that* quantised toe, or the receiver
        // (which uses the decoded toe) places the SV n·Δtoe along its orbit off.
        let toe = round_to_16s(gps.sow);
        let dt0 = scene_epoch_unix - orbit.epoch_unix;
        let dt_toe = dt0 + (toe - gps.sow); // orbit epoch → toe
        // e = 0 ⇒ mean anomaly = true anomaly.
        let m0 = wrap_pi(orbit.true_anomaly_deg_at_epoch.to_radians() + n * dt_toe);
        // Ωk = Ω0 − Ω̇ₑ·t_sow is toe-independent; Ω0 anchors the ECEF node.
        let omega0 = wrap_pi(
            orbit.raan_deg.to_radians() + GPS_OMEGA_E_DOT * gps.sow - GPS_OMEGA_E_DOT * dt0,
        );
        Self {
            prn,
            week: gps.week,
            toe_s: toe,
            toc_s: toe,
            sqrt_a: a.sqrt(),
            e: 0.0,
            m0_rad: m0,
            delta_n_rad_s: 0.0,
            i0_rad: orbit.inclination_deg.to_radians(),
            idot_rad_s: 0.0,
            omega0_rad: omega0,
            omega_dot_rad_s: 0.0,
            argp_rad: 0.0,
            cuc: 0.0, cus: 0.0, crc: 0.0, crs: 0.0, cic: 0.0, cis: 0.0,
            af0: 0.0, af1: 0.0, af2: 0.0, tgd: 0.0,
            iode: 0, iodc: 0, health: 0, ura_index: 0,
        }
    }

    /// Osculating ephemeris from an ECEF state vector at `toe`. Used for
    /// arbitrary trajectories (e.g. SGP4/TLE). Reproduces the position exactly
    /// at `toe`; away from `toe` it follows a two-body propagation, so the
    /// residual vs. the true orbit is the (informative) error metric.
    // r[impl fn.4.lnav] osculating ephemeris at toe for any trajectory
    pub fn from_trajectory(
        traj: &dyn Trajectory,
        prn: u8,
        toe_unix: f64,
        gps_at_toe: GpsTime,
    ) -> Self {
        // Snap toe to the 16 s transmit grid and evaluate the state there, so
        // the osculating elements are consistent with the decoded toe.
        let toe_sow = round_to_16s(gps_at_toe.sow);
        let toe_unix = toe_unix + (toe_sow - gps_at_toe.sow);
        let r = traj.position_ecef(toe_unix);
        let v = traj.velocity_ecef(toe_unix);
        // Rotate ECEF → GPS quasi-inertial frame at toe (undo Ω̇ₑ·toe), adding
        // Earth rotation to the velocity. Then classical osculating elements
        // give the broadcast set directly (Ω0 = inertial node).
        let phi = GPS_OMEGA_E_DOT * toe_sow;
        let (s, c) = phi.sin_cos();
        let r_i = [r[0] * c - r[1] * s, r[0] * s + r[1] * c, r[2]];
        // v_inertial = R_z(phi)·(v_ecef + ω × r_ecef)
        let v_rot = [
            v[0] - GPS_OMEGA_E_DOT * r[1],
            v[1] + GPS_OMEGA_E_DOT * r[0],
            v[2],
        ];
        let v_i = [v_rot[0] * c - v_rot[1] * s, v_rot[0] * s + v_rot[1] * c, v_rot[2]];

        let el = classical_elements(r_i, v_i, GPS_MU);
        Self {
            prn,
            week: gps_at_toe.week,
            toe_s: toe_sow,
            toc_s: toe_sow,
            sqrt_a: el.a.sqrt(),
            e: el.e,
            m0_rad: wrap_pi(el.m),
            delta_n_rad_s: 0.0,
            i0_rad: el.i,
            idot_rad_s: 0.0,
            omega0_rad: wrap_pi(el.raan),
            omega_dot_rad_s: 0.0,
            argp_rad: wrap_pi(el.argp),
            cuc: 0.0, cus: 0.0, crc: 0.0, crs: 0.0, cic: 0.0, cis: 0.0,
            af0: 0.0, af1: 0.0, af2: 0.0, tgd: 0.0,
            iode: 0, iodc: 0, health: 0, ura_index: 0,
        }
    }
}

struct Classical {
    a: f64,
    e: f64,
    i: f64,
    raan: f64,
    argp: f64,
    m: f64,
}

/// Classical osculating elements from an inertial position/velocity state.
fn classical_elements(r: [f64; 3], v: [f64; 3], mu: f64) -> Classical {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    };
    let norm = |a: [f64; 3]| dot(a, a).sqrt();

    let rmag = norm(r);
    let vmag2 = dot(v, v);
    let h = cross(r, v);
    let hmag = norm(h);
    // Node vector n = k × h.
    let nvec = [-h[1], h[0], 0.0];
    let nmag = (nvec[0] * nvec[0] + nvec[1] * nvec[1]).sqrt();
    // Eccentricity vector.
    let vxh = cross(v, h);
    let evec = [
        vxh[0] / mu - r[0] / rmag,
        vxh[1] / mu - r[1] / rmag,
        vxh[2] / mu - r[2] / rmag,
    ];
    let e = norm(evec);
    let energy = vmag2 / 2.0 - mu / rmag;
    let a = -mu / (2.0 * energy);
    let i = (h[2] / hmag).acos();

    let mut raan = if nmag > 1e-12 { (nvec[0] / nmag).acos() } else { 0.0 };
    if nvec[1] < 0.0 { raan = 2.0 * PI - raan; }

    // Argument of perigee.
    let mut argp = if nmag > 1e-12 && e > 1e-12 {
        (dot(nvec, evec) / (nmag * e)).clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    if evec[2] < 0.0 { argp = 2.0 * PI - argp; }

    // True anomaly.
    let mut nu = if e > 1e-12 {
        (dot(evec, r) / (e * rmag)).clamp(-1.0, 1.0).acos()
    } else {
        // Circular: use argument of latitude from node.
        let mut u = if nmag > 1e-12 {
            (dot(nvec, r) / (nmag * rmag)).clamp(-1.0, 1.0).acos()
        } else {
            (r[0] / rmag).acos()
        };
        if r[2] < 0.0 { u = 2.0 * PI - u; }
        u
    };
    if e > 1e-12 && dot(r, v) < 0.0 { nu = 2.0 * PI - nu; }

    // Eccentric → mean anomaly.
    let big_e = 2.0 * ((1.0 - e).sqrt() * (nu / 2.0).sin())
        .atan2((1.0 + e).sqrt() * (nu / 2.0).cos());
    let m = big_e - e * big_e.sin();

    Classical { a, e, i, raan, argp, m }
}

// ───────────────────────── LNAV frame assembly ─────────────────────────

const P2_5: f64 = 1.0 / 32.0;
fn p2(n: i32) -> f64 { 2f64.powi(n) }
fn sc(rad: f64) -> f64 { rad / PI } // radians → semicircles

/// Accumulates the 24 source data bits of one word, MSB-first.
struct WordBits {
    bits: Vec<u8>,
}
impl WordBits {
    fn new() -> Self { Self { bits: Vec::with_capacity(24) } }
    fn push_u(&mut self, val: u64, nbits: usize) {
        for i in (0..nbits).rev() { self.bits.push(((val >> i) & 1) as u8); }
    }
    fn push_i(&mut self, val: i64, nbits: usize) {
        let mask = if nbits >= 64 { u64::MAX } else { (1u64 << nbits) - 1 };
        self.push_u((val as u64) & mask, nbits);
    }
    fn pad_to(&mut self, n: usize) { while self.bits.len() < n { self.bits.push(0); } }
    fn finish(self) -> [u8; 24] {
        assert_eq!(self.bits.len(), 24, "word must be exactly 24 data bits");
        self.bits.try_into().unwrap()
    }
}

fn scaled_i(phys: f64, lsb: f64) -> i64 { (phys / lsb).round() as i64 }
fn scaled_u(phys: f64, lsb: f64) -> u64 { (phys / lsb).round() as u64 }

/// Parity bit masks over the source data bits (IS-GPS-200 Table 20-XIV), in
/// the same convention as gps-sdr-sim: the 24 data bits sit in bit positions
/// 29..6 of a 32-bit word, `D29*`/`D30*` in bits 31/30, parity in bits 5..0.
const PARITY_BMASK: [u32; 6] = [
    0x3B1F_3480, 0x1D8F_9A40, 0x2EC7_CD00, 0x1763_E680, 0x2BB1_F340, 0x0B7A_89C0,
];

/// 30-bit GPS word (complemented data + parity) from a 32-bit source word.
/// Faithful port of gps-sdr-sim's `computeChecksum`: parity is computed over
/// the **source** bits, the data is complemented by `D30*` separately, and for
/// words 2 and 10 (`nib`) data bits 23–24 are solved so parity bits 29,30 = 0.
fn compute_checksum(source: u32, nib: bool) -> u32 {
    let mut d = source & 0x3FFF_FFC0;
    let d29 = (source >> 31) & 1;
    let d30 = (source >> 30) & 1;
    if nib {
        if (d30 + (PARITY_BMASK[4] & d).count_ones()) % 2 == 1 { d ^= 1 << 6; }
        if (d29 + (PARITY_BMASK[5] & d).count_ones()) % 2 == 1 { d ^= 1 << 7; }
    }
    let mut big = d;
    if d30 == 1 { big ^= 0x3FFF_FFC0; }
    big |= ((d29 + (PARITY_BMASK[0] & d).count_ones()) % 2) << 5;
    big |= ((d30 + (PARITY_BMASK[1] & d).count_ones()) % 2) << 4;
    big |= ((d29 + (PARITY_BMASK[2] & d).count_ones()) % 2) << 3;
    big |= ((d30 + (PARITY_BMASK[3] & d).count_ones()) % 2) << 2;
    big |= ((d30 + (PARITY_BMASK[4] & d).count_ones()) % 2) << 1;
    big |= (d29 + (PARITY_BMASK[5] & d).count_ones()) % 2;
    big & 0x3FFF_FFFF
}

/// Encode one 30-bit word from 24 source data bits and the previous word's
/// `D29*`/`D30*`. `solve_tbits` is set for words 2 and 10.
fn encode_word(d_src: &[u8; 24], d29p: u8, d30p: u8, solve_tbits: bool) -> [u8; 30] {
    // Pack into the 32-bit source word: D29* at bit31, D30* at bit30,
    // d1..d24 at bits 29..6.
    let mut source: u32 = (d29p as u32) << 31 | (d30p as u32) << 30;
    for (i, &b) in d_src.iter().enumerate() {
        source |= (b as u32 & 1) << (29 - i);
    }
    let w = compute_checksum(source, solve_tbits);
    let mut out = [0u8; 30];
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = ((w >> (29 - i)) & 1) as u8;
    }
    out
}

/// TLM (word 1): preamble 0x8B + (zero) message + integrity/reserved.
fn tlm_word(d29p: u8, d30p: u8) -> ([u8; 30], (u8, u8)) {
    let mut w = WordBits::new();
    w.push_u(0x8B, 8);  // preamble
    w.push_u(0, 14);    // TLM message
    w.push_u(0, 2);     // integrity status + reserved
    let bits = w.finish();
    let enc = encode_word(&bits, d29p, d30p, false);
    (enc, (enc[28], enc[29]))
}

/// HOW (word 2): TOW count of the next subframe, flags, subframe id, solved
/// parity bits.
fn how_word(tow_count: u32, subframe_id: u8, d29p: u8, d30p: u8) -> ([u8; 30], (u8, u8)) {
    let mut w = WordBits::new();
    w.push_u(tow_count as u64 & 0x1_FFFF, 17); // 17-bit TOW count
    w.push_u(0, 1);  // alert
    w.push_u(0, 1);  // anti-spoof
    w.push_u(subframe_id as u64, 3);
    w.push_u(0, 2);  // t-bits (solved)
    let bits = w.finish();
    let enc = encode_word(&bits, d29p, d30p, true);
    (enc, (enc[28], enc[29]))
}

/// Assemble one 300-bit subframe (10 words). `data_words` supplies the 24
/// source data bits for words 3–10 (word 10 gets its t-bits solved).
fn build_subframe(
    subframe_id: u8,
    tow_count: u32,
    data_words: [[u8; 24]; 8],
    d29p_in: u8,
    d30p_in: u8,
) -> (Vec<u8>, (u8, u8)) {
    let mut out = Vec::with_capacity(300);
    let (w1, p1) = tlm_word(d29p_in, d30p_in);
    out.extend_from_slice(&w1);
    let (w2, mut p) = how_word(tow_count, subframe_id, p1.0, p1.1);
    out.extend_from_slice(&w2);
    for (k, dw) in data_words.iter().enumerate() {
        let solve = k == 7; // word 10
        let enc = encode_word(dw, p.0, p.1, solve);
        out.extend_from_slice(&enc);
        p = (enc[28], enc[29]);
    }
    (out, p)
}

/// Split a 32-bit value into (8 MSB, 24 LSB).
fn split_32(v: i64) -> (u64, u64) {
    let u = (v as u64) & 0xFFFF_FFFF;
    (u >> 24, u & 0xFF_FFFF)
}

fn subframe1_words(eph: &Ephemeris) -> [[u8; 24]; 8] {
    let tgd = scaled_i(eph.tgd, p2(-31));
    let toc = scaled_u(eph.toc_s, 16.0);
    let af2 = scaled_i(eph.af2, p2(-55));
    let af1 = scaled_i(eph.af1, p2(-43));
    let af0 = scaled_i(eph.af0, p2(-31));

    // Word 3
    let mut w3 = WordBits::new();
    w3.push_u(eph.week as u64 & 0x3FF, 10);
    w3.push_u(0, 2); // codes on L2
    w3.push_u(eph.ura_index as u64, 4);
    w3.push_u(eph.health as u64, 6);
    w3.push_u((eph.iodc as u64 >> 8) & 0x3, 2); // IODC MSB
    let w3 = w3.finish();
    // Words 4,5,6 reserved
    let zeros = { let mut w = WordBits::new(); w.pad_to(24); w.finish() };
    // Word 7: reserved(16) + TGD(8)
    let mut w7 = WordBits::new();
    w7.push_u(0, 16);
    w7.push_i(tgd, 8);
    let w7 = w7.finish();
    // Word 8: IODC LSB(8) + toc(16)
    let mut w8 = WordBits::new();
    w8.push_u(eph.iodc as u64 & 0xFF, 8);
    w8.push_u(toc, 16);
    let w8 = w8.finish();
    // Word 9: af2(8) + af1(16)
    let mut w9 = WordBits::new();
    w9.push_i(af2, 8);
    w9.push_i(af1, 16);
    let w9 = w9.finish();
    // Word 10: af0(22) + 2 t-bits
    let mut w10 = WordBits::new();
    w10.push_i(af0, 22);
    w10.push_u(0, 2);
    let w10 = w10.finish();

    // Words 3–10 of subframe 1: word 3, reserved words 4/5/6, then 7–10.
    [w3, zeros, zeros, zeros, w7, w8, w9, w10]
}

fn subframe2_words(eph: &Ephemeris) -> [[u8; 24]; 8] {
    let crs = scaled_i(eph.crs, P2_5);
    let dn = scaled_i(sc(eph.delta_n_rad_s), p2(-43));
    let m0 = scaled_i(sc(eph.m0_rad), p2(-31));
    let cuc = scaled_i(eph.cuc, p2(-29));
    let ecc = scaled_u(eph.e, p2(-33)) as i64;
    let cus = scaled_i(eph.cus, p2(-29));
    let sqa = scaled_u(eph.sqrt_a, p2(-19)) as i64;
    let toe = scaled_u(eph.toe_s, 16.0);

    let (m0_hi, m0_lo) = split_32(m0);
    let (e_hi, e_lo) = split_32(ecc);
    let (a_hi, a_lo) = split_32(sqa);

    // W3: IODE(8) + Crs(16)
    let mut w3 = WordBits::new();
    w3.push_u(eph.iode as u64, 8);
    w3.push_i(crs, 16);
    let w3 = w3.finish();
    // W4: dn(16) + M0 MSB(8)
    let mut w4 = WordBits::new();
    w4.push_i(dn, 16);
    w4.push_u(m0_hi, 8);
    let w4 = w4.finish();
    // W5: M0 LSB(24)
    let mut w5 = WordBits::new();
    w5.push_u(m0_lo, 24);
    let w5 = w5.finish();
    // W6: Cuc(16) + e MSB(8)
    let mut w6 = WordBits::new();
    w6.push_i(cuc, 16);
    w6.push_u(e_hi, 8);
    let w6 = w6.finish();
    // W7: e LSB(24)
    let mut w7 = WordBits::new();
    w7.push_u(e_lo, 24);
    let w7 = w7.finish();
    // W8: Cus(16) + sqrtA MSB(8)
    let mut w8 = WordBits::new();
    w8.push_i(cus, 16);
    w8.push_u(a_hi, 8);
    let w8 = w8.finish();
    // W9: sqrtA LSB(24)
    let mut w9 = WordBits::new();
    w9.push_u(a_lo, 24);
    let w9 = w9.finish();
    // W10: toe(16) + fit(1) + AODO(5) + 2 t-bits
    let mut w10 = WordBits::new();
    w10.push_u(toe, 16);
    w10.push_u(0, 1);
    w10.push_u(0, 5);
    w10.push_u(0, 2);
    let w10 = w10.finish();

    [w3, w4, w5, w6, w7, w8, w9, w10]
}

fn subframe3_words(eph: &Ephemeris) -> [[u8; 24]; 8] {
    let cic = scaled_i(eph.cic, p2(-29));
    let omg0 = scaled_i(sc(eph.omega0_rad), p2(-31));
    let cis = scaled_i(eph.cis, p2(-29));
    let i0 = scaled_i(sc(eph.i0_rad), p2(-31));
    let crc = scaled_i(eph.crc, P2_5);
    let argp = scaled_i(sc(eph.argp_rad), p2(-31));
    let omgdot = scaled_i(sc(eph.omega_dot_rad_s), p2(-43));
    let idot = scaled_i(sc(eph.idot_rad_s), p2(-43));

    let (o_hi, o_lo) = split_32(omg0);
    let (i_hi, i_lo) = split_32(i0);
    let (w_hi, w_lo) = split_32(argp);

    // W3: Cic(16) + Omega0 MSB(8)
    let mut w3 = WordBits::new();
    w3.push_i(cic, 16);
    w3.push_u(o_hi, 8);
    let w3 = w3.finish();
    // W4: Omega0 LSB(24)
    let mut w4 = WordBits::new();
    w4.push_u(o_lo, 24);
    let w4 = w4.finish();
    // W5: Cis(16) + i0 MSB(8)
    let mut w5 = WordBits::new();
    w5.push_i(cis, 16);
    w5.push_u(i_hi, 8);
    let w5 = w5.finish();
    // W6: i0 LSB(24)
    let mut w6 = WordBits::new();
    w6.push_u(i_lo, 24);
    let w6 = w6.finish();
    // W7: Crc(16) + argp MSB(8)
    let mut w7 = WordBits::new();
    w7.push_i(crc, 16);
    w7.push_u(w_hi, 8);
    let w7 = w7.finish();
    // W8: argp LSB(24)
    let mut w8 = WordBits::new();
    w8.push_u(w_lo, 24);
    let w8 = w8.finish();
    // W9: OmegaDot(24)
    let mut w9 = WordBits::new();
    w9.push_i(omgdot, 24);
    let w9 = w9.finish();
    // W10: IODE(8) + IDOT(14) + 2 t-bits
    let mut w10 = WordBits::new();
    w10.push_u(eph.iode as u64, 8);
    w10.push_i(idot, 14);
    w10.push_u(0, 2);
    let w10 = w10.finish();

    [w3, w4, w5, w6, w7, w8, w9, w10]
}

fn zero_data_words() -> [[u8; 24]; 8] {
    let z = { let mut w = WordBits::new(); w.pad_to(24); w.finish() };
    [z; 8]
}

/// Build the ±1 LNAV bit stream for `duration_s` of a satellite, with the
/// first subframe starting at the scene epoch. The scene epoch should fall on
/// a 6 s GPS boundary (true for any `epoch_unix` divisible by 6 with an
/// integer leap-second offset) so subframe edges align with the TOW Z-count.
// r[impl fn.4.lnav] TLM/HOW + subframes 1–3 with parity
pub fn lnav_bits(eph: &Ephemeris, gps_epoch: GpsTime, duration_s: f64) -> Vec<i8> {
    let base_z = (gps_epoch.sow / 6.0).round() as u32;
    let n_subframes = (duration_s / 6.0).ceil() as u32 + 1;
    let mut bits: Vec<i8> = Vec::with_capacity((n_subframes as usize) * 300);
    let mut d29p = 0u8;
    let mut d30p = 0u8;
    for g in 0..n_subframes {
        let subframe_id = ((g % 5) + 1) as u8;
        // HOW carries the TOW count of the *next* subframe's leading edge.
        let tow_count = base_z + g + 1;
        let words = match subframe_id {
            1 => subframe1_words(eph),
            2 => subframe2_words(eph),
            3 => subframe3_words(eph),
            _ => zero_data_words(),
        };
        let (sf, p) = build_subframe(subframe_id, tow_count, words, d29p, d30p);
        d29p = p.0;
        d30p = p.1;
        // 0/1 → ±1 (data bit 1 → −1, matching the ±1 C/A chip convention).
        for b in sf {
            bits.push(1 - 2 * b as i8);
        }
    }
    bits
}

// ───────────────────── SV position oracle (IS-GPS-200 Table 20-IV) ─────────

/// Satellite ECEF position (m) computed from broadcast ephemeris at GPS
/// seconds-of-week `t_sow`, per IS-GPS-200 §20.3.3.4.3. Used to validate that
/// a derived ephemeris reproduces the simulated trajectory.
pub fn sv_position_ecef(eph: &Ephemeris, t_sow: f64) -> [f64; 3] {
    let a = eph.sqrt_a * eph.sqrt_a;
    let n0 = (GPS_MU / a.powi(3)).sqrt();
    let n = n0 + eph.delta_n_rad_s;
    let mut tk = t_sow - eph.toe_s;
    if tk > 302_400.0 { tk -= SECONDS_PER_WEEK; }
    if tk < -302_400.0 { tk += SECONDS_PER_WEEK; }

    let mk = eph.m0_rad + n * tk;
    // Solve Kepler's equation for eccentric anomaly.
    let mut ek = mk;
    for _ in 0..10 {
        ek = mk + eph.e * ek.sin();
    }
    let nu = ((1.0 - eph.e * eph.e).sqrt() * ek.sin()).atan2(ek.cos() - eph.e);
    let phi = nu + eph.argp_rad;
    let (s2, c2) = (2.0 * phi).sin_cos();
    let du = eph.cus * s2 + eph.cuc * c2;
    let dr = eph.crs * s2 + eph.crc * c2;
    let di = eph.cis * s2 + eph.cic * c2;
    let u = phi + du;
    let r = a * (1.0 - eph.e * ek.cos()) + dr;
    let i = eph.i0_rad + di + eph.idot_rad_s * tk;
    let xp = r * u.cos();
    let yp = r * u.sin();
    let omega = eph.omega0_rad + (eph.omega_dot_rad_s - GPS_OMEGA_E_DOT) * tk
        - GPS_OMEGA_E_DOT * eph.toe_s;
    let (so, co) = omega.sin_cos();
    let (si, ci) = i.sin_cos();
    [
        xp * co - yp * ci * so,
        xp * so + yp * ci * co,
        yp * si,
    ]
}

// ───────────────────── RINEX broadcast ephemeris ──────────────────────────

/// A satellite trajectory backed directly by a GPS broadcast `Ephemeris`
/// (e.g. parsed from a RINEX nav file). `position_ecef` is the IS-GPS-200
/// SV-position algorithm, so when this drives a `nav = "full"` transmitter the
/// broadcast message carries the *same* ephemeris and the receiver reconstructs
/// the orbit exactly — no osculating fit, no residual.
// r[impl fn.5.rinex-nav] fly the broadcast ephemeris (IS-GPS-200 SV position)
pub struct GpsBroadcastTrajectory {
    eph: Ephemeris,
    /// UTC unix time of `toe`, so a query time maps to GPS seconds-of-week
    /// without needing the leap offset again: `t_sow = toe_s + (t − toe_unix)`.
    toe_unix: f64,
}

impl GpsBroadcastTrajectory {
    pub fn new(eph: Ephemeris, utc_offset_s: f64) -> Self {
        let toe_unix =
            GPS_EPOCH_UNIX + eph.week as f64 * SECONDS_PER_WEEK + eph.toe_s - utc_offset_s;
        Self { eph, toe_unix }
    }
    /// The broadcast ephemeris this trajectory evaluates — fed straight into
    /// the LNAV builder so signal geometry and broadcast message agree exactly.
    pub fn ephemeris(&self) -> &Ephemeris { &self.eph }
}

impl Trajectory for GpsBroadcastTrajectory {
    fn position_ecef(&self, t_unix: f64) -> [f64; 3] {
        sv_position_ecef(&self.eph, self.eph.toe_s + (t_unix - self.toe_unix))
    }
}

/// Read a 19-character RINEX float field at byte `start` ('D'/'d' exponent
/// markers normalised to 'E'); missing/short fields read as 0.
fn rinex_field(line: &str, start: usize) -> f64 {
    let end = (start + 19).min(line.len());
    if start >= line.len() { return 0.0; }
    line[start..end].trim().replace(['D', 'd'], "E").parse().unwrap_or(0.0)
}

/// Parse GPS nav records from RINEX 3 text and return the `Ephemeris` for `prn`
/// whose `toe` is closest to `near_unix`. Angles in RINEX nav are already in
/// radians/SI, mapping straight onto [`Ephemeris`].
// r[impl fn.5.rinex-nav] RINEX 3 GPS nav parser
pub fn parse_rinex3_gps_nav(
    text: &str,
    prn: u8,
    near_unix: f64,
    utc_offset_s: f64,
) -> Option<Ephemeris> {
    let lines: Vec<&str> = text.lines().collect();
    // Skip the header.
    let mut i = lines.iter().position(|l| l.contains("END OF HEADER")).map_or(0, |p| p + 1);
    let mut best: Option<(f64, Ephemeris)> = None;
    while i + 7 < lines.len() {
        let l0 = lines[i];
        // GPS record starts with "Gnn"; other systems (R/E/C/J) are skipped.
        if l0.len() < 3 || l0.as_bytes()[0] != b'G' {
            i += 1;
            continue;
        }
        let rec_prn: u8 = l0[1..3].trim().parse().unwrap_or(0);
        if rec_prn != prn {
            i += 8; // GPS LNAV records are 8 lines
            continue;
        }
        let (l1, l2, l3, l4, l5, l6) =
            (lines[i + 1], lines[i + 2], lines[i + 3], lines[i + 4], lines[i + 5], lines[i + 6]);
        let toe = rinex_field(l3, 4);
        let week = rinex_field(l5, 42);
        let eph = Ephemeris {
            prn,
            week: week as u32,
            toe_s: toe,
            toc_s: toe, // af terms are tiny; toc≈toe is fine
            af0: rinex_field(l0, 23), af1: rinex_field(l0, 42), af2: rinex_field(l0, 61),
            iode: rinex_field(l1, 4) as u8,
            crs: rinex_field(l1, 23),
            delta_n_rad_s: rinex_field(l1, 42),
            m0_rad: rinex_field(l1, 61),
            cuc: rinex_field(l2, 4),
            e: rinex_field(l2, 23),
            cus: rinex_field(l2, 42),
            sqrt_a: rinex_field(l2, 61),
            cic: rinex_field(l3, 23),
            omega0_rad: rinex_field(l3, 42),
            cis: rinex_field(l3, 61),
            i0_rad: rinex_field(l4, 4),
            crc: rinex_field(l4, 23),
            argp_rad: rinex_field(l4, 42),
            omega_dot_rad_s: rinex_field(l4, 61),
            idot_rad_s: rinex_field(l5, 4),
            health: rinex_field(l6, 23) as u8,
            tgd: rinex_field(l6, 42),
            iodc: rinex_field(l6, 61) as u16,
            ura_index: 0,
        };
        let toe_unix = GPS_EPOCH_UNIX + week * SECONDS_PER_WEEK + toe - utc_offset_s;
        let d = (toe_unix - near_unix).abs();
        if best.as_ref().is_none_or(|(bd, _)| d < *bd) {
            best = Some((d, eph));
        }
        i += 8;
    }
    best.map(|(_, e)| e)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parity_ok(word: &[u8; 30], d29p: u8, d30p: u8) -> bool {
        // Recover source data bits (D_n = d_n ⊕ D30*) and recompute the word.
        let mut src = [0u8; 24];
        for n in 0..24 { src[n] = word[n] ^ d30p; }
        let recomputed = encode_word(&src, d29p, d30p, false);
        recomputed == *word
    }

    // r[verify fn.4.lnav]
    #[test]
    fn parity_self_consistent_across_subframe() {
        let eph = sample_ephemeris();
        let (sf, _) = build_subframe(1, 1234, subframe1_words(&eph), 0, 0);
        // Walk the words, verifying each word's parity against the running
        // previous-parity state.
        let mut d29p = 0u8;
        let mut d30p = 0u8;
        for wi in 0..10 {
            let word: [u8; 30] = sf[wi * 30..wi * 30 + 30].try_into().unwrap();
            assert!(parity_ok(&word, d29p, d30p), "parity fail in word {wi}");
            d29p = word[28];
            d30p = word[29];
        }
    }

    // r[verify fn.4.lnav]
    #[test]
    fn words_2_and_10_force_parity_zero() {
        let eph = sample_ephemeris();
        let (sf, _) = build_subframe(2, 999, subframe2_words(&eph), 0, 0);
        let w2: [u8; 30] = sf[30..60].try_into().unwrap();
        assert_eq!((w2[28], w2[29]), (0, 0), "HOW must force parity 29,30 = 0");
        let w10: [u8; 30] = sf[270..300].try_into().unwrap();
        assert_eq!((w10[28], w10[29]), (0, 0), "word 10 must force parity 29,30 = 0");
    }

    // r[verify fn.4.lnav]
    #[test]
    fn tlm_preamble_present() {
        let (w1, _) = tlm_word(0, 0);
        let preamble: Vec<u8> = w1[..8].to_vec();
        assert_eq!(preamble, vec![1, 0, 0, 0, 1, 0, 1, 1], "0x8B preamble");
    }

    // r[verify fn.4.lnav]
    #[test]
    fn gps_time_roundtrip() {
        // 2024-06-01 16:54:00 UTC, a 6 s-aligned epoch.
        let t = GpsTime::from_unix(1_717_260_840.0, 18.0);
        assert!(t.sow >= 0.0 && t.sow < SECONDS_PER_WEEK);
        assert_eq!(t.sow % 6.0, 0.0, "epoch should land on a 6 s boundary");
    }

    fn sample_ephemeris() -> Ephemeris {
        let orbit = CircularKeplerian {
            semi_major_axis_m: 26_560_000.0,
            inclination_deg: 55.0,
            raan_deg: 40.0,
            true_anomaly_deg_at_epoch: 120.0,
            epoch_unix: 1_717_260_840.0,
        };
        let gps = GpsTime::from_unix(1_717_260_840.0, 18.0);
        Ephemeris::from_circular_keplerian(&orbit, 1, 1_717_260_840.0, gps)
    }

    // ---- independent LNAV decoder (the receiver side) ----

    /// Recover the per-word source data bits (un-complemented) from a ±1 stream
    /// and verify every word's parity, as a receiver would.
    fn decode_words(bits: &[i8]) -> Vec<[u8; 24]> {
        let tx: Vec<u8> = bits.iter().map(|&pm| if pm < 0 { 1 } else { 0 }).collect();
        let nwords = tx.len() / 30;
        let mut d29p = 0u8;
        let mut d30p = 0u8;
        let mut out = Vec::new();
        for w in 0..nwords {
            let word = &tx[w * 30..w * 30 + 30];
            let mut src = [0u8; 24];
            for n in 0..24 { src[n] = word[n] ^ d30p; }
            // Parity check (receiver side): recompute and compare.
            let recomputed = encode_word(&src, d29p, d30p, false);
            assert_eq!(&recomputed[..], word, "parity check failed at word {w}");
            out.push(src);
            d29p = word[28];
            d30p = word[29];
        }
        out
    }

    fn read_u(src: &[u8; 24], start: usize, len: usize) -> u64 {
        (0..len).fold(0u64, |a, i| (a << 1) | src[start + i] as u64)
    }
    fn read_i(src: &[u8; 24], start: usize, len: usize) -> i64 {
        let v = read_u(src, start, len);
        let sign = 1u64 << (len - 1);
        if v & sign != 0 { (v as i64) - (1i64 << len) } else { v as i64 }
    }

    // r[verify fn.4.lnav]
    #[test]
    fn lnav_bits_decode_back_to_the_ephemeris() {
        // Pack → emit → decode must round-trip the key ephemeris fields, the
        // subframe IDs and the TOW progression, with valid parity throughout.
        // This is what GNSS-SDR does; if it fails, the receiver rejects the
        // ephemeris (no fix).
        let eph = sample_ephemeris();
        let gps = GpsTime::from_unix(1_717_260_840.0, 18.0);
        let bits = lnav_bits(&eph, gps, 30.0); // ≥1 full frame
        let words = decode_words(&bits);

        let base_z = (gps.sow / 6.0).round() as u64;
        // Walk subframes; check id + TOW, decode SF1/2/3 fields.
        for (sf_idx, chunk) in words.chunks(10).enumerate() {
            if chunk.len() < 10 { break; }
            let how = &chunk[1];
            let tow = read_u(how, 0, 17);
            let sfid = read_u(how, 19, 3);
            assert_eq!(sfid as usize, (sf_idx % 5) + 1, "subframe id mismatch at {sf_idx}");
            assert_eq!(tow, base_z + sf_idx as u64 + 1, "TOW mismatch at {sf_idx}");

            match sfid {
                1 => {
                    let week = read_u(&chunk[2], 0, 10);
                    assert_eq!(week, (eph.week & 0x3FF) as u64, "WN");
                }
                2 => {
                    let m0_hi = read_u(&chunk[3], 16, 8);
                    let m0_lo = read_u(&chunk[4], 0, 24);
                    let m0 = sext32((m0_hi << 24) | m0_lo) as f64 * p2(-31) * PI;
                    let a_hi = read_u(&chunk[7], 16, 8);
                    let a_lo = read_u(&chunk[8], 0, 24);
                    let sqrt_a = ((a_hi << 24) | a_lo) as f64 * p2(-19);
                    let toe = read_u(&chunk[9], 0, 16) as f64 * 16.0;
                    let crs = read_i(&chunk[2], 8, 16) as f64 * P2_5;
                    let _ = crs;
                    assert!((m0 - eph.m0_rad).abs() < 1e-8, "M0 {m0} vs {}", eph.m0_rad);
                    assert!((sqrt_a - eph.sqrt_a).abs() < 1e-3, "sqrtA {sqrt_a} vs {}", eph.sqrt_a);
                    // toe must encode losslessly: the receiver uses the decoded
                    // toe with M0, so any quantisation here mis-places the SV.
                    assert_eq!(toe, eph.toe_s, "decoded toe != ephemeris toe (must be 16 s-aligned)");
                }
                3 => {
                    let o_hi = read_u(&chunk[2], 16, 8);
                    let o_lo = read_u(&chunk[3], 0, 24);
                    let omega0 = sext32((o_hi << 24) | o_lo) as f64 * p2(-31) * PI;
                    let i_hi = read_u(&chunk[4], 16, 8);
                    let i_lo = read_u(&chunk[5], 0, 24);
                    let i0 = sext32((i_hi << 24) | i_lo) as f64 * p2(-31) * PI;
                    assert!((omega0 - eph.omega0_rad).abs() < 1e-8, "Omega0 {omega0} vs {}", eph.omega0_rad);
                    assert!((i0 - eph.i0_rad).abs() < 1e-8, "i0 {i0} vs {}", eph.i0_rad);
                }
                _ => {}
            }
        }
    }

    fn sext32(v: u64) -> i64 {
        let v = v & 0xFFFF_FFFF;
        if v & 0x8000_0000 != 0 { v as i64 - (1i64 << 32) } else { v as i64 }
    }

    // A RINEX 3 GPS nav record (PRN 01), columns laid out per the spec.
    const RINEX3: &str = concat!(
        "     3.04           N: GNSS NAV DATA    G: GPS              RINEX VERSION / TYPE\n",
        "                                                            END OF HEADER\n",
        "G01 2024 06 01 00 00 00 4.629502538592E-04-1.022364418051E-11 0.000000000000E+00\n",
        "     6.100000000000E+01-1.435937500000E+02 4.916182694377E-09-2.838821577537E+00\n",
        "    -7.413327693939E-06 9.470977378078E-03 6.282702088356E-06 5.153664415359E+03\n",
        "     5.184000000000E+05-1.490116119385E-08 1.665292991674E+00-1.117587089539E-07\n",
        "     9.609727413912E-01 2.157812500000E+02 6.967118599040E-01-8.108285197566E-09\n",
        "    -2.500170875643E-10 1.000000000000E+00 2.317000000000E+03 0.000000000000E+00\n",
        "     2.000000000000E+00 0.000000000000E+00-1.117587089539E-08 6.100000000000E+01\n",
        "     5.110800000000E+05 4.000000000000E+00\n",
    );

    // r[verify fn.5.rinex-nav]
    #[test]
    fn parses_rinex3_gps_nav_record() {
        let toe_unix = GPS_EPOCH_UNIX + 2317.0 * SECONDS_PER_WEEK + 518400.0 - 18.0;
        let eph = parse_rinex3_gps_nav(RINEX3, 1, toe_unix, 18.0).expect("PRN 1 record");
        assert_eq!(eph.week, 2317);
        assert_eq!(eph.toe_s, 518400.0);
        assert!((eph.sqrt_a - 5153.664415359).abs() < 1e-6, "sqrtA {}", eph.sqrt_a);
        assert!((eph.e - 9.470977378078e-3).abs() < 1e-12, "e {}", eph.e);
        assert!((eph.i0_rad - 0.9609727413912).abs() < 1e-12, "i0 {}", eph.i0_rad);
        assert!((eph.delta_n_rad_s - 4.916182694377e-9).abs() < 1e-18);
        // SV position from the parsed ephemeris sits at GPS altitude.
        let p = sv_position_ecef(&eph, eph.toe_s);
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!((26.2e6..26.9e6).contains(&r), "radius {r} not GPS-like");
        assert!(parse_rinex3_gps_nav(RINEX3, 7, toe_unix, 18.0).is_none(), "PRN 7 absent");
    }

    // r[verify fn.5.rinex-nav]
    #[test]
    fn broadcast_trajectory_matches_its_own_ephemeris() {
        let toe_unix = GPS_EPOCH_UNIX + 2317.0 * SECONDS_PER_WEEK + 518400.0 - 18.0;
        let eph = parse_rinex3_gps_nav(RINEX3, 1, toe_unix, 18.0).unwrap();
        let traj = GpsBroadcastTrajectory::new(eph, 18.0);
        // The trajectory evaluates the SV-position algorithm on its ephemeris;
        // querying at toe_unix reproduces sv_position at toe (exact NAV path).
        let p = traj.position_ecef(toe_unix);
        let sv = sv_position_ecef(traj.ephemeris(), traj.ephemeris().toe_s);
        for k in 0..3 { assert!((p[k] - sv[k]).abs() < 1e-6, "axis {k}"); }
    }

    // r[verify fn.4.lnav]
    #[test]
    fn decoded_ephemeris_places_satellite_within_metres() {
        // The faithful end-to-end guard: decode the ephemeris from the bits as
        // a receiver does, then run the SV-position algorithm with the DECODED
        // toe/M0 and confirm it reproduces the simulated orbit. A toe/M0
        // quantisation mismatch (toe has a 16 s LSB) shows up as ~km here.
        let orbit = CircularKeplerian {
            semi_major_axis_m: 26_560_000.0, inclination_deg: 55.0, raan_deg: 40.0,
            true_anomaly_deg_at_epoch: 120.0, epoch_unix: 1_717_260_840.0,
        };
        let gps = GpsTime::from_unix(1_717_260_840.0, 18.0);
        let mut dec = Ephemeris::from_circular_keplerian(&orbit, 1, 1_717_260_840.0, gps);
        let bits = lnav_bits(&dec, gps, 30.0);
        let words = decode_words(&bits);
        for chunk in words.chunks(10) {
            if chunk.len() < 10 { continue; }
            match read_u(&chunk[1], 19, 3) {
                2 => {
                    dec.m0_rad = sext32((read_u(&chunk[3],16,8) << 24) | read_u(&chunk[4],0,24)) as f64 * p2(-31) * PI;
                    dec.sqrt_a = (((read_u(&chunk[7],16,8) << 24) | read_u(&chunk[8],0,24)) as f64) * p2(-19);
                    dec.toe_s = read_u(&chunk[9],0,16) as f64 * 16.0;
                }
                3 => {
                    dec.omega0_rad = sext32((read_u(&chunk[2],16,8) << 24) | read_u(&chunk[3],0,24)) as f64 * p2(-31) * PI;
                    dec.i0_rad = sext32((read_u(&chunk[4],16,8) << 24) | read_u(&chunk[5],0,24)) as f64 * p2(-31) * PI;
                }
                _ => {}
            }
        }
        for dt in [0.0, 30.0, 300.0] {
            let truth = orbit.position_ecef(1_717_260_840.0 + dt);
            let sv = sv_position_ecef(&dec, gps.sow + dt);
            let err = ((truth[0]-sv[0]).powi(2)+(truth[1]-sv[1]).powi(2)+(truth[2]-sv[2]).powi(2)).sqrt();
            assert!(err < 5.0, "decoded-ephemeris SV error {err:.1} m at dt={dt} (toe/M0 mismatch?)");
        }
    }

    // r[verify fn.4.lnav]
    #[test]
    fn ephemeris_from_tle_is_consistent_over_window() {
        // A TLE-derived osculating ephemeris reproduces the SGP4 ECEF exactly
        // at toe, and stays close over a fix window (two-body vs SGP4 only).
        use crate::geometry::tle::TleTrajectory;
        let l1 = "1 26407U 00040A   26150.58242625  .00000046  00000+0  00000+0 0  9999";
        let l2 = "2 26407  54.8538 215.8617 0121278 302.3883 230.2243  2.00557636189594";
        let traj = TleTrajectory::from_tle(l1, l2).unwrap();
        // Pick a toe at a 6 s-aligned unix time near the TLE epoch.
        let toe_unix = (traj.epoch_unix() / 6.0).round() * 6.0;
        let gps = GpsTime::from_unix(toe_unix, 18.0);
        let eph = Ephemeris::from_trajectory(&traj, 22, toe_unix, gps);

        // At toe: sub-metre.
        let p0 = traj.position_ecef(toe_unix);
        let s0 = sv_position_ecef(&eph, gps.sow);
        for k in 0..3 {
            assert!((p0[k] - s0[k]).abs() < 1.0, "toe axis {k}: {} vs {}", p0[k], s0[k]);
        }
        // Over 60 s: within a few km (osculating drift), nowhere near the
        // ~90 km a frame mismatch would produce.
        let pt = traj.position_ecef(toe_unix + 60.0);
        let st = sv_position_ecef(&eph, gps.sow + 60.0);
        let err = ((pt[0]-st[0]).powi(2)+(pt[1]-st[1]).powi(2)+(pt[2]-st[2]).powi(2)).sqrt();
        assert!(err < 5000.0, "60 s drift {err} m too large (frame mismatch?)");
    }

    // r[verify fn.4.lnav]
    #[test]
    fn ephemeris_reproduces_circular_trajectory() {
        use crate::geometry::Trajectory;
        let orbit = CircularKeplerian {
            semi_major_axis_m: 26_560_000.0,
            inclination_deg: 55.0,
            raan_deg: 40.0,
            true_anomaly_deg_at_epoch: 120.0,
            epoch_unix: 1_717_260_840.0,
        };
        let gps = GpsTime::from_unix(1_717_260_840.0, 18.0);
        let eph = Ephemeris::from_circular_keplerian(&orbit, 1, 1_717_260_840.0, gps);

        // Across a few minutes, the broadcast ephemeris (via the IS-GPS-200 SV
        // algorithm) must match the trajectory's ECEF to sub-mm.
        for dt in [0.0, 1.0, 30.0, 120.0, 600.0] {
            let t_unix = 1_717_260_840.0 + dt;
            let truth = orbit.position_ecef(t_unix);
            let sv = sv_position_ecef(&eph, gps.sow + dt);
            for k in 0..3 {
                assert!((truth[k] - sv[k]).abs() < 1e-3,
                    "axis {k} at dt={dt}: truth {} vs sv {}", truth[k], sv[k]);
            }
        }
    }
}
