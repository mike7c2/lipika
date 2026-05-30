//! GPS L1 C/A code generator. Port of the Python `signals/gnss.py` module.
//!
//! Reference: IS-GPS-200 §3.3.2.1. Two LFSRs (G1, G2) produce a 1023-chip
//! Gold sequence per PRN at a 1.023 Mcps chip rate. Navigation data is
//! not modelled (per DGSN convention — we only care about code phase).

use ndarray::Array1;
use num_complex::Complex32;

use super::SignalSource;

pub const CA_CHIP_RATE_HZ: f64 = 1_023_000.0;
pub const CA_CODE_LENGTH: usize = 1023;

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

pub struct CaCodeSource {
    pub carrier_hz: f64,
    pub prn: usize,
    /// GPS time tag for chip 0 (UTC seconds). For the synthetic GPS
    /// scene we just use the pass-start epoch.
    pub epoch_unix: f64,
    code: Vec<i8>,
}

impl CaCodeSource {
    pub fn new(carrier_hz: f64, prn: usize, epoch_unix: f64) -> Self {
        let code = ca_code(prn);
        Self { carrier_hz, prn, epoch_unix, code }
    }
}

impl SignalSource for CaCodeSource {
    fn carrier_hz(&self) -> f64 { self.carrier_hz }

    fn evaluate_baseband_batch(&self, t_emit: &Array1<f64>) -> Array1<Complex32> {
        let mut out = Array1::<Complex32>::zeros(t_emit.len());
        for (i, &t) in t_emit.iter().enumerate() {
            let chip_phase = (t - self.epoch_unix) * CA_CHIP_RATE_HZ;
            // Modulo in the rust-friendly way (handle negatives correctly).
            let idx = chip_phase.floor() as i64;
            let idx_mod = idx.rem_euclid(CA_CODE_LENGTH as i64) as usize;
            out[i] = Complex32::new(self.code[idx_mod] as f32, 0.0);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_ten_chips_of_prn_1_match_spec() {
        // IS-GPS-200 Table 3-Ia, octal 1440: first 10 chips of PRN 1 are
        // 1100100000 (binary) which maps to -1 -1 +1 +1 -1 +1 +1 +1 +1 +1.
        let code = ca_code(1);
        let expected: [i8; 10] = [-1, -1, 1, 1, -1, 1, 1, 1, 1, 1];
        assert_eq!(&code[..10], &expected, "PRN 1 first 10 chips diverge from IS-GPS-200");
    }
}
