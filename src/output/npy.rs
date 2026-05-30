//! Write a complex64 ndarray to a NumPy `.npy` file.
//!
//! We hand-roll the .npy v1.0 header (well-documented, stable) rather than
//! fight npyz's typing system over the complex dtype. Output is byte-
//! compatible with what `numpy.save` produces for a complex64 1-D array,
//! readable as `np.load(path)` on the Python side.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use ndarray::Array1;
use num_complex::Complex32;

pub fn write_complex64_npy(path: &Path, iq: &Array1<Complex32>) -> Result<(), anyhow::Error> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    let n = iq.len();

    // Header dict text; trailing newline required.
    let header_body = format!("{{'descr': '<c8', 'fortran_order': False, 'shape': ({n},), }}");
    let prefix_len = 10; // 6 magic + 2 version + 2 header-length u16
    // Pad with spaces so prefix + header is a multiple of 64 (convention).
    let total_aligned = (prefix_len + header_body.len() + 1).next_multiple_of(64);
    let padded_header_len = total_aligned - prefix_len;
    let spaces = padded_header_len - header_body.len() - 1;
    let mut header = header_body;
    header.extend(std::iter::repeat(' ').take(spaces));
    header.push('\n');
    debug_assert_eq!(header.len(), padded_header_len);
    debug_assert_eq!((prefix_len + header.len()) % 64, 0);

    writer.write_all(b"\x93NUMPY")?;
    writer.write_all(&[1u8, 0u8])?;
    writer.write_all(&(header.len() as u16).to_le_bytes())?;
    writer.write_all(header.as_bytes())?;

    let mut buf = [0u8; 8];
    for v in iq.iter() {
        buf[..4].copy_from_slice(&v.re.to_le_bytes());
        buf[4..].copy_from_slice(&v.im.to_le_bytes());
        writer.write_all(&buf)?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn npy_header_is_valid_complex64() {
        let path = std::env::temp_dir().join("dgsn_iq_gen_test.npy");
        let iq = Array1::from(vec![
            Complex32::new(1.0, 2.0),
            Complex32::new(3.0, -4.0),
            Complex32::new(0.5, 0.25),
        ]);
        write_complex64_npy(&path, &iq).unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[..6], b"\x93NUMPY");
        assert_eq!(bytes[6], 1);
        assert_eq!(bytes[7], 0);
        let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        let header = std::str::from_utf8(&bytes[10..10 + header_len]).unwrap();
        assert!(header.contains("'descr': '<c8'"));
        assert!(header.contains("'shape': (3,)"));
        assert_eq!(bytes.len(), 10 + header_len + 24);
        let _ = std::fs::remove_file(&path);
    }
}
