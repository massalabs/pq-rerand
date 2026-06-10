//! 31-bit plaintext encoding: bytes ↔ polynomial coefficients.

use crate::params::{N, BITS_PER_COEFF, SLOT_BYTES};

/// Encode a byte slice (exactly SLOT_BYTES = 15,872 bytes) into N = 4096
/// polynomial coefficients, each holding 31 bits.
///
/// The input is treated as a little-endian bitstream; coefficient `i` holds
/// bits `[31·i, 31·i + 31)`. Implemented with unaligned 64-bit window reads
/// (~31× fewer memory operations than a bit-by-bit loop). The memory-access
/// pattern depends only on the coefficient index, never on the (secret)
/// plaintext content.
pub fn encode(data: &[u8]) -> [u32; N] {
    assert_eq!(data.len(), SLOT_BYTES, "encode: input must be exactly {} bytes", SLOT_BYTES);
    const MASK: u64 = (1u64 << BITS_PER_COEFF) - 1;
    // Pad so the final 8-byte window read stays in bounds.
    let mut padded = [0u8; SLOT_BYTES + 8];
    padded[..SLOT_BYTES].copy_from_slice(data);
    let mut coeffs = [0u32; N];
    for (i, c) in coeffs.iter_mut().enumerate() {
        let bit = i * BITS_PER_COEFF;
        let byte = bit / 8;
        let shift = bit % 8; // ≤ 7, so shift + 31 ≤ 38 < 64
        let w = u64::from_le_bytes(padded[byte..byte + 8].try_into().unwrap());
        *c = ((w >> shift) & MASK) as u32;
    }
    coeffs
}

/// Decode N = 4096 polynomial coefficients (31 bits each) back into
/// SLOT_BYTES = 15,872 bytes. Inverse of [`encode`].
pub fn decode(coeffs: &[u32; N]) -> Vec<u8> {
    let mut padded = vec![0u8; SLOT_BYTES + 8];
    for (i, &c) in coeffs.iter().enumerate() {
        let bit = i * BITS_PER_COEFF;
        let byte = bit / 8;
        let shift = bit % 8;
        let mut w = u64::from_le_bytes(padded[byte..byte + 8].try_into().unwrap());
        w |= (c as u64) << shift;
        padded[byte..byte + 8].copy_from_slice(&w.to_le_bytes());
    }
    padded.truncate(SLOT_BYTES);
    padded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let mut data = vec![0u8; SLOT_BYTES];
        // Fill with a pattern
        for (i, b) in data.iter_mut().enumerate() {
            *b = (i % 256) as u8;
        }
        let coeffs = encode(&data);
        let decoded = decode(&coeffs);
        assert_eq!(data, decoded, "encode/decode roundtrip failed");
    }

    /// Bit-by-bit reference encoder (the original implementation), used to pin
    /// down the bitstream layout of the windowed encoder.
    fn encode_reference(data: &[u8]) -> [u32; N] {
        let mut coeffs = [0u32; N];
        let mut bit_pos = 0usize;
        for c in coeffs.iter_mut() {
            let mut val = 0u32;
            for b in 0..BITS_PER_COEFF {
                let byte_idx = (bit_pos + b) / 8;
                let bit_idx = (bit_pos + b) % 8;
                val |= (((data[byte_idx] >> bit_idx) & 1) as u32) << b;
            }
            *c = val;
            bit_pos += BITS_PER_COEFF;
        }
        coeffs
    }

    #[test]
    fn test_encode_matches_reference() {
        let mut data = vec![0u8; SLOT_BYTES];
        for (i, b) in data.iter_mut().enumerate() {
            *b = ((i * 131 + 17) % 256) as u8;
        }
        assert_eq!(encode(&data), encode_reference(&data));
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
