//! 31-bit plaintext encoding: bytes ↔ polynomial coefficients.

use crate::params::{N, BITS_PER_COEFF, SLOT_BYTES};

/// Encode a byte slice (exactly SLOT_BYTES = 15,872 bytes) into N = 4096
/// polynomial coefficients, each holding 31 bits.
pub fn encode(data: &[u8]) -> [u32; N] {
    assert_eq!(data.len(), SLOT_BYTES, "encode: input must be exactly {} bytes", SLOT_BYTES);
    let mut coeffs = [0u32; N];
    let mut bit_pos = 0usize;
    for i in 0..N {
        let mut val = 0u32;
        for b in 0..BITS_PER_COEFF {
            let byte_idx = (bit_pos + b) / 8;
            let bit_idx = (bit_pos + b) % 8;
            if byte_idx < data.len() {
                val |= (((data[byte_idx] >> bit_idx) & 1) as u32) << b;
            }
        }
        coeffs[i] = val;
        bit_pos += BITS_PER_COEFF;
    }
    coeffs
}

/// Decode N = 4096 polynomial coefficients (31 bits each) back into
/// SLOT_BYTES = 15,872 bytes.
pub fn decode(coeffs: &[u32; N]) -> Vec<u8> {
    let mut data = vec![0u8; SLOT_BYTES];
    let mut bit_pos = 0usize;
    for i in 0..N {
        let val = coeffs[i];
        for b in 0..BITS_PER_COEFF {
            let byte_idx = (bit_pos + b) / 8;
            let bit_idx = (bit_pos + b) % 8;
            if byte_idx < data.len() {
                data[byte_idx] |= (((val >> b) & 1) as u8) << bit_idx;
            }
        }
        bit_pos += BITS_PER_COEFF;
    }
    data
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let mut data = vec![0u8; SLOT_BYTES];
        // Fill with a pattern
        for i in 0..data.len() {
            data[i] = (i % 256) as u8;
        }
        let coeffs = encode(&data);
        let decoded = decode(&coeffs);
        assert_eq!(data, decoded, "encode/decode roundtrip failed");
    }

    #[test]
    fn test_encode_zeros() {
        let data = vec![0u8; SLOT_BYTES];
        let coeffs = encode(&data);
        for c in &coeffs {
            assert_eq!(*c, 0);
        }
    }

    #[test]
    fn test_coefficients_in_range() {
        let data = vec![0xFFu8; SLOT_BYTES];
        let coeffs = encode(&data);
        for c in &coeffs {
            assert!(*c < (1u32 << BITS_PER_COEFF), "coefficient out of range: {}", c);
        }
    }
}
