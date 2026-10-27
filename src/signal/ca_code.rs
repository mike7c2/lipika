//! GPS L1 C/A code generator. Port of the Python `signals/gnss.py` module.
//!
//! Reference: IS-GPS-200 §3.3.2.1. Two LFSRs (G1, G2) produce a 1023-chip
//! Gold sequence per PRN at a 1.023 Mcps chip rate.
//!
//! Optionally a 50 bps navigation data stream is BPSK-overlaid on the code
//! (each data bit spans 20 ms = 20 C/A code periods; the emitted chip is
//! `gold_chip × nav_bit`, both ±1). The data is *not* a valid ephemeris —
//! it carries the mandatory 8-bit preamble (`0x8B`) at the start of every
//! 300-bit (6 s) subframe with the remaining 292 bits filled
//! deterministically-pseudo-randomly. That is enough to drive a real
//! receiver's bit synchroniser, Costas loop and telemetry-decoder state
//! machine into a sustained tracking lock (it will fail parity / produce no
//! PVT fix, which is fine for signal-chain validation). Disable it for the
//! pure-code convention a downstream code-phase clock-recovery pipeline may
//! assume.

use ndarray::Array1;
use num_complex::Complex32;

use super::SignalModel;

pub const CA_CHIP_RATE_HZ: f64 = 1_023_000.0;
pub const CA_CODE_LENGTH: usize = 1023;

/// GPS navigation data rate (bits per second).
pub const NAV_BIT_RATE_HZ: f64 = 50.0;
/// Bits per subframe (one subframe = 6 s at 50 bps).
pub const NAV_SUBFRAME_BITS: i64 = 300;
/// The 8-bit telemetry-word preamble `0x8B` = `1000 1011`, as data bits.
const NAV_PREAMBLE_BITS: [u8; 8] = [1, 0, 0, 0, 1, 0, 1, 1];

/// G2 phase tap pairs by PRN (1..32) per IS-GPS-200 Table 3-Ia. Values are
/// 1-indexed in the spec; converted to 0-indexed at access time.
const G2_TAPS: &[(usize, usize)] = &[
    (2, 6), (3, 7), (4, 8), (5, 9), (1, 9), (2, 10), (1, 8), (2, 9), (3, 10),
    (2, 3), (3, 4), (5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (1, 4), (2, 5),
    (3, 6), (4, 7), (5, 8), (6, 9), (1, 3), (4, 6), (5, 7), (6, 8), (7, 9),
    (8, 10), (1, 6), (2, 7), (3, 8), (4, 9),
];

fn generate_g1() -> Vec<i8> {
    let mut reg = [1i8; 10];
    let mut out = Vec::with_capacity(CA_CODE_LENGTH);
    for _ in 0..CA_CODE_LENGTH {
        out.push(reg[9]);
        let fb = reg[2] ^ reg[9];
        for i in (1..10).rev() { reg[i] = reg[i - 1]; }
        reg[0] = fb;
    }
    out
}

fn generate_g2(tap_a: usize, tap_b: usize) -> Vec<i8> {
    let mut reg = [1i8; 10];
    let mut out = Vec::with_capacity(CA_CODE_LENGTH);
    for _ in 0..CA_CODE_LENGTH {
        out.push(reg[tap_a - 1] ^ reg[tap_b - 1]);
        let fb = reg[1] ^ reg[2] ^ reg[5] ^ reg[7] ^ reg[8] ^ reg[9];
        for i in (1..10).rev() { reg[i] = reg[i - 1]; }
        reg[0] = fb;
    }
    out
}

/// Returns the 1023-chip C/A sequence for `prn` as +1/-1.
// r[impl fn.4.ca-code] IS-GPS-200 §3.3.2.1 Gold code
pub fn ca_code(prn: usize) -> Vec<i8> {
    assert!((1..=32).contains(&prn), "PRN must be 1..32");
    let g1 = generate_g1();
    let (a, b) = G2_TAPS[prn - 1];
    let g2 = generate_g2(a, b);
    g1.iter().zip(g2.iter()).map(|(&x, &y)| {
        // XOR in {0,1} terms, then map to +1/-1.
        let bit = x ^ y;
        1 - 2 * bit
    }).collect()
}

/// Deterministic pseudo-random navigation data overlay.
#[derive(Clone, Copy)]
struct NavData {
    /// Per-source seed (folds in the PRN so different SVs carry different
    /// payloads while still sharing the fixed preamble).
    seed: u64,
}

impl NavData {
    /// Data bit (as ±1) at global bit index `bit_idx`. The first 8 bits of
    /// every 300-bit subframe are the fixed `0x8B` preamble; the rest are
    /// deterministic pseudo-random flips (so the Costas loop sees regular
    /// phase inversions and never stalls on a constant run).
    fn bit_pm(&self, bit_idx: i64) -> i8 {
        let pos = bit_idx.rem_euclid(NAV_SUBFRAME_BITS);
        let bit = if pos < NAV_PREAMBLE_BITS.len() as i64 {
            NAV_PREAMBLE_BITS[pos as usize]
        } else {
            // splitmix64 on (seed, bit_idx) → 1 random bit.
            let mut x = self.seed
                ^ (bit_idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            x ^= x >> 30; x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x ^= x >> 27; x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^= x >> 31;
            (x & 1) as u8
        };
        1 - 2 * bit as i8 // {0,1} → {+1,-1}
    }
}

/// Navigation-data overlay mode on the C/A code.
enum Nav {
    /// Pure code, no data.
    Off,
    /// Deterministic preamble + pseudo-random payload (good for lock, fails
    /// parity — no position fix).
    Dummy(NavData),
    /// A precomputed ±1 LNAV bit stream (valid frames + ephemeris → position
    /// fix). Indexed by the running 50 bps bit number.
    Full(Vec<i8>),
}

// r[impl fn.4.ca-code] C/A chips × nav overlay (off / dummy / full)
pub struct CaCodeSource {
    pub carrier_hz: f64,
    pub prn: usize,
    /// Time tag for chip 0, in **seconds relative to the scene epoch**. The
    /// baseband is evaluated at epoch-relative emission times, so keeping this
    /// relative preserves sub-ns code-phase precision (an absolute Unix tag
    /// would quantise to the ~238 ns f64 ULP).
    pub epoch_rel_s: f64,
    code: Vec<i8>,
    nav: Nav,
}

impl CaCodeSource {
    /// Pure C/A code, no navigation data (the code-phase-only convention).
    /// Kept for benches/tests and any pure-code downstream pipeline.
    pub fn new(carrier_hz: f64, prn: usize, epoch_rel_s: f64) -> Self {
        let code = ca_code(prn);
        Self { carrier_hz, prn, epoch_rel_s, code, nav: Nav::Off }
    }

    /// C/A code with a 50 bps navigation data overlay (preamble + pseudo-
    /// random payload). `nav_seed` is folded with the PRN so distinct SVs
    /// carry distinct data while sharing the fixed `0x8B` preamble.
    pub fn with_nav_data(carrier_hz: f64, prn: usize, epoch_rel_s: f64, nav_seed: u64) -> Self {
        let code = ca_code(prn);
        let seed = nav_seed
            .wrapping_add((prn as u64).wrapping_mul(0xD1B5_4A32_D192_ED03));
        Self { carrier_hz, prn, epoch_rel_s, code, nav: Nav::Dummy(NavData { seed }) }
    }

    /// C/A code carrying a fully-formed LNAV bit stream (see
    /// [`super::gps_nav`]). `lnav_bits` is the ±1 sequence at 50 bps starting
    /// at `epoch_rel_s`.
    // r[impl fn.4.lnav]
    pub fn with_full_nav(carrier_hz: f64, prn: usize, epoch_rel_s: f64, lnav_bits: Vec<i8>) -> Self {
        let code = ca_code(prn);
        Self { carrier_hz, prn, epoch_rel_s, code, nav: Nav::Full(lnav_bits) }
    }
}

impl SignalModel for CaCodeSource {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        let mut out = Array1::<Complex32>::zeros(t_emit.len());
        for (i, &t) in t_emit.iter().enumerate() {
            let dt = t - self.epoch_rel_s;
            let chip_phase = dt * CA_CHIP_RATE_HZ;
            // Modulo in the rust-friendly way (handle negatives correctly).
            let idx = chip_phase.floor() as i64;
            let idx_mod = idx.rem_euclid(CA_CODE_LENGTH as i64) as usize;
            let mut chip = self.code[idx_mod];
            let nav_pm = match &self.nav {
                Nav::Off => 1,
                Nav::Dummy(d) => d.bit_pm((dt * NAV_BIT_RATE_HZ).floor() as i64),
                Nav::Full(bits) => {
                    let bit_idx = (dt * NAV_BIT_RATE_HZ).floor() as i64;
                    if bit_idx < 0 {
                        bits.first().copied().unwrap_or(1)
                    } else {
                        let j = bit_idx as usize;
                        bits.get(j).or_else(|| bits.last()).copied().unwrap_or(1)
                    }
                }
            };
            chip *= nav_pm;
            out[i] = Complex32::new(chip as f32, 0.0);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // r[verify fn.4.ca-code]
    #[test]
    fn first_ten_chips_of_prn_1_match_spec() {
        // IS-GPS-200 Table 3-Ia, octal 1440: first 10 chips of PRN 1 are
        // 1100100000 (binary) which maps to -1 -1 +1 +1 -1 +1 +1 +1 +1 +1.
        let code = ca_code(1);
        let expected: [i8; 10] = [-1, -1, 1, 1, -1, 1, 1, 1, 1, 1];
        assert_eq!(&code[..10], &expected, "PRN 1 first 10 chips diverge from IS-GPS-200");
    }

    // r[verify fn.4.ca-code]
    #[test]
    fn pure_code_has_no_nav_overlay() {
        let src = CaCodeSource::new(1_575_420_000.0, 1, 0.0);
        assert!(matches!(src.nav, Nav::Off));
        // First chip is just the gold chip (no sign flip).
        let t = Array1::from_vec(vec![1e-9]);
        let out = src.evaluate_baseband_batch(&t);
        assert_eq!(out[0].re, ca_code(1)[0] as f32);
    }

    // r[verify fn.4.lnav]
    #[test]
    fn full_nav_overlays_the_supplied_bit_stream() {
        // First NAV bit is the 0x8B preamble MSB (1 → −1); the first chip
        // should be gold_chip × (−1).
        let bits = vec![-1i8, 1, 1, 1, -1, 1, -1, -1];
        let src = CaCodeSource::with_full_nav(1_575_420_000.0, 1, 0.0, bits);
        let t = Array1::from_vec(vec![1e-9]);
        let out = src.evaluate_baseband_batch(&t);
        assert_eq!(out[0].re, (ca_code(1)[0] * -1) as f32);
    }

    // r[verify fn.4.ca-code]
    #[test]
    fn nav_preamble_is_0x8b_at_subframe_starts() {
        let nav = NavData { seed: 12345 };
        let expected: [i8; 8] = [-1, 1, 1, 1, -1, 1, -1, -1]; // 0x8B → ±1
        // Subframe 0 and subframe 1 (bit 300) both start with the preamble.
        for sf in 0..3 {
            let base = sf * NAV_SUBFRAME_BITS;
            for (k, &e) in expected.iter().enumerate() {
                assert_eq!(nav.bit_pm(base + k as i64), e,
                    "preamble bit {k} of subframe {sf} wrong");
            }
        }
    }

    // r[verify fn.4.ca-code]
    #[test]
    fn nav_overlay_multiplies_gold_chip() {
        let src = CaCodeSource::with_nav_data(1_575_420_000.0, 1, 0.0, 42);
        let nav_pm = match &src.nav {
            Nav::Dummy(d) => d.bit_pm(0),
            _ => panic!("expected dummy nav"),
        };
        // A time inside the first NAV bit (bit 0, chip 0).
        let t = Array1::from_vec(vec![1e-9]);
        let out = src.evaluate_baseband_batch(&t);
        let expected = (ca_code(1)[0] * nav_pm) as f32;
        assert_eq!(out[0].re, expected);
    }

    // r[verify fn.4.ca-code]
    #[test]
    fn nav_data_flips_within_a_subframe() {
        // The 292 payload bits must not be a constant run, or the Costas
        // loop would have nothing to track.
        let nav = NavData { seed: 999 };
        let payload: Vec<i8> = (8..NAV_SUBFRAME_BITS).map(|i| nav.bit_pm(i)).collect();
        assert!(payload.iter().any(|&b| b == 1));
        assert!(payload.iter().any(|&b| b == -1));
    }
}
