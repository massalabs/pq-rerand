//! 31-bit plaintext encoding: bytes ↔ polynomial coefficients.

use crate::params::{N, BITS_PER_COEFF, SLOT_BYTES};
use crate::poly::{Secret, secret};
use zeroize::Zeroizing;

const MASK: u64 = (1u64 << BITS_PER_COEFF) - 1;
// `decode` writes the slot out in 32-bit words.
const _: () = assert!((N * BITS_PER_COEFF).is_multiple_of(32));

/// A plaintext slot: the `N` coefficients of a message polynomial.
///
/// [`encode`] produces coefficients below 2^31. Decryption produces coefficients
/// below `t`: one of 2^31 or more shows that the ciphertext did not decrypt to an
/// encoded payload, which [`decode`] reports. The coefficients live on the heap
/// and are zeroized on drop.
pub struct Plaintext(pub(crate) Secret);

/// Position of coefficient `i` in the byte stream: the offset of an 8-byte window
/// and the shift of the coefficient inside it. The window of the last coefficient
/// is aligned on the end of the slot, so no window reaches past it.
#[inline(always)]
fn window(i: usize) -> (usize, usize) {
    let bit = i * BITS_PER_COEFF;
    let byte = (bit / 8).min(SLOT_BYTES - 8);
    (byte, bit - 8 * byte) // shift <= 33, so shift + 31 <= 64
}

/// Encode a slot of SLOT_BYTES = 15,872 bytes into N = 4096
/// polynomial coefficients, each holding 31 bits.
///
/// The input is treated as a little-endian bitstream; coefficient `i` holds
/// bits `[31·i, 31·i + 31)`. Implemented with unaligned 64-bit window reads
/// (~31× fewer memory operations than a bit-by-bit loop). The memory-access
/// pattern depends only on the coefficient index, never on the (secret)
/// plaintext content.
pub fn encode(data: &[u8; SLOT_BYTES]) -> Plaintext {
    let mut coeffs = secret();
    for (i, c) in coeffs.iter_mut().enumerate() {
        let (byte, shift) = window(i);
        let w = u64::from_le_bytes(data[byte..byte + 8].try_into().unwrap());
        *c = ((w >> shift) & MASK) as u32;
    }
    Plaintext(coeffs)
}

/// Decode a plaintext slot back into SLOT_BYTES = 15,872 bytes. Inverse of
/// [`encode`]. The bytes are zeroized on drop.
///
/// Returns `None` if a coefficient is 2^31 or more, which no encoded slot contains:
/// the ciphertext did not decrypt to a payload (wrong key, modified ciphertext,
/// exhausted noise budget). The whole slot is scanned before the verdict, in
/// constant time. Treat `None` exactly like a failed authentication of the
/// payload, and do not let anyone who can submit ciphertexts tell the two apart
/// (not even by timing): the difference would reveal comparisons on secret values.
pub fn decode(plaintext: &Plaintext) -> Option<Zeroizing<Vec<u8>>> {
    let mut bytes = Zeroizing::new(vec![0u8; SLOT_BYTES]);
    let mut seen = 0u32;
    // Bits wait in `acc` until 32 of them can be written out at `at`.
    let (mut acc, mut pending, mut at) = (0u64, 0, 0);
    for &c in plaintext.0.iter() {
        seen |= c;
        acc |= (c as u64 & MASK) << pending;
        pending += BITS_PER_COEFF;
        if pending >= 32 {
            bytes[at..at + 4].copy_from_slice(&(acc as u32).to_le_bytes());
            acc >>= 32;
            pending -= 32;
            at += 4;
        }
    }
    (seen >> BITS_PER_COEFF == 0).then_some(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pattern(mul: usize, add: usize) -> [u8; SLOT_BYTES] {
        std::array::from_fn(|i| ((i * mul + add) % 256) as u8)
    }

    #[test]
    fn test_encode_decode_roundtrip() {
        for data in [pattern(1, 0), pattern(131, 17), [0u8; SLOT_BYTES], [0xFF; SLOT_BYTES]] {
            let decoded = decode(&encode(&data)).expect("an encoded slot must decode");
            assert_eq!(data[..], decoded[..], "encode/decode roundtrip failed");
        }
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
        for data in [pattern(131, 17), [0xFF; SLOT_BYTES]] {
            assert_eq!(*encode(&data).0, encode_reference(&data));
        }
    }

    #[test]
    fn test_encode_zeros() {
        let coeffs = encode(&[0u8; SLOT_BYTES]);
        for c in coeffs.0.iter() {
            assert_eq!(*c, 0);
        }
    }

    #[test]
    fn test_coefficients_in_range() {
        let coeffs = encode(&[0xFFu8; SLOT_BYTES]);
        for c in coeffs.0.iter() {
            assert_eq!(*c, (1u32 << BITS_PER_COEFF) - 1, "coefficient out of range: {}", c);
        }
    }

    /// A coefficient of 2^31 or more is reported, wherever it is, and is not
    /// carried into its neighbor.
    #[test]
    fn test_decode_rejects_out_of_range() {
        let data = pattern(7, 3);
        for i in [0, 1, N / 2, N - 1] {
            for extra in [1u32 << 31, u32::MAX] {
                let mut m = encode(&data);
                m.0[i] |= extra;
                assert!(decode(&m).is_none(), "coefficient {i} accepted");
            }
        }
        // Exactly 2^31, alone in a zero slot: the smallest value that is rejected.
        for i in [0, N - 1] {
            let mut m = Plaintext(secret());
            m.0[i] = 1 << 31;
            assert!(decode(&m).is_none(), "coefficient {i} = 2^31 accepted");
        }
        // The largest canonical value is accepted everywhere.
        let mut largest = Plaintext(secret());
        largest.0.fill((1 << 31) - 1);
        assert!(decode(&largest).is_some());
    }
}
