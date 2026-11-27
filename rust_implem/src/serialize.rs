//! Serialization: CiphertextSlot ↔ bytes.

use crate::params::{N, T, Q2};
use crate::encrypt::CiphertextSlot;

impl CiphertextSlot {
    /// Serialized byte size.
    pub const BYTES: usize = 4 * N * 4; // 65,536

    /// Serialize to bytes: `c0_t[4096] | c0_q2[4096] | c1_t[4096] | c1_q2[4096]`,
    /// each coefficient a little-endian u32 (65,536 bytes in all).
    ///
    /// A slot carries no version byte: it has the format of the key it was made
    /// with, and a serialized key starts with
    /// [`FORMAT_VERSION`](crate::params::FORMAT_VERSION). Nothing in a slot names
    /// that key or its format: a container that keeps slots apart from their key
    /// must record the version itself.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(Self::BYTES);
        for arr in [&self.c0_t, &self.c0_q2, &self.c1_t, &self.c1_q2] {
            for &val in arr.iter() {
                buf.extend_from_slice(&val.to_le_bytes());
            }
        }
        buf
    }

    /// Deserialize from bytes.
    ///
    /// Returns `None` if `data` is not exactly `BYTES` long or if any
    /// coefficient is out of range (not fully reduced modulo `t` / `q2`).
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() != Self::BYTES {
            return None;
        }
        let mut ct = CiphertextSlot::zeroed();
        let arrays: [(&mut [u32; N], u32); 4] = [
            (&mut ct.c0_t, T),
            (&mut ct.c0_q2, Q2),
            (&mut ct.c1_t, T),
            (&mut ct.c1_q2, Q2),
        ];
        let mut offset = 0;
        for (arr, modulus) in arrays {
            for val in arr.iter_mut() {
                let bytes: [u8; 4] = data[offset..offset + 4].try_into().ok()?;
                *val = u32::from_le_bytes(bytes);
                if *val >= modulus {
                    return None;
                }
                offset += 4;
            }
        }
        Some(ct)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serialize_roundtrip() {
        let mut ct = CiphertextSlot::zeroed();
        ct.c0_t[0] = 42;
        ct.c0_q2[100] = 123456;
        ct.c1_t[N - 1] = 999999;
        ct.c1_q2[2048] = 0xDEADBEEF;

        let bytes = ct.to_bytes();
        assert_eq!(bytes.len(), CiphertextSlot::BYTES);

        let ct2 = CiphertextSlot::from_bytes(&bytes).unwrap();
        assert_eq!(ct.c0_t, ct2.c0_t);
        assert_eq!(ct.c0_q2, ct2.c0_q2);
        assert_eq!(ct.c1_t, ct2.c1_t);
        assert_eq!(ct.c1_q2, ct2.c1_q2);
    }

    #[test]
    fn test_deserialize_rejects_bad_length() {
        assert!(CiphertextSlot::from_bytes(&[0u8; 10]).is_none());
        assert!(CiphertextSlot::from_bytes(&vec![0u8; CiphertextSlot::BYTES + 1]).is_none());
    }

    #[test]
    fn test_deserialize_rejects_out_of_range() {
        // A c0_t coefficient ≥ T must be rejected.
        let mut bytes = vec![0u8; CiphertextSlot::BYTES];
        bytes[..4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(CiphertextSlot::from_bytes(&bytes).is_none());

        // In every limb, at both ends: the modulus is rejected and the value just
        // below it is accepted and written back unchanged.
        for (limb, modulus) in [T, Q2, T, Q2].into_iter().enumerate() {
            for index in [0, N - 1] {
                let at = (limb * N + index) * 4;
                let mut bytes = vec![0u8; CiphertextSlot::BYTES];
                bytes[at..at + 4].copy_from_slice(&modulus.to_le_bytes());
                assert!(CiphertextSlot::from_bytes(&bytes).is_none(), "limb {limb}: modulus accepted");
                bytes[at..at + 4].copy_from_slice(&(modulus - 1).to_le_bytes());
                let ct = CiphertextSlot::from_bytes(&bytes).expect("reduced coefficient rejected");
                assert_eq!(ct.to_bytes(), bytes);
            }
        }
    }
}
