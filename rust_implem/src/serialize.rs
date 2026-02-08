//! Serialization: CiphertextSlot ↔ bytes.
//!
//! Format per slot: c0_t[4096] | c0_q2[4096] | c1_t[4096] | c1_q2[4096]
//! Each coefficient is a little-endian u32 (4 bytes).
//! Total per slot: 4 × 4096 × 4 = 65,536 bytes = 64 KiB.

use crate::params::N;
use crate::encrypt::CiphertextSlot;

/// Bytes per serialized slot.
pub const SLOT_CT_BYTES: usize = 4 * N * 4; // 65,536

/// Serialize a ciphertext slot to bytes (little-endian u32 arrays).
pub fn serialize_slot(ct: &CiphertextSlot) -> Vec<u8> {
    let mut buf = Vec::with_capacity(SLOT_CT_BYTES);
    for arr in [&ct.c0_t, &ct.c0_q2, &ct.c1_t, &ct.c1_q2] {
        for &val in arr.iter() {
            buf.extend_from_slice(&val.to_le_bytes());
        }
    }
    buf
}

/// Deserialize a ciphertext slot from bytes.
pub fn deserialize_slot(data: &[u8]) -> CiphertextSlot {
    assert_eq!(data.len(), SLOT_CT_BYTES);
    let mut ct = CiphertextSlot {
        c0_t: [0u32; N],
        c0_q2: [0u32; N],
        c1_t: [0u32; N],
        c1_q2: [0u32; N],
    };
    let arrays: [&mut [u32; N]; 4] = [
        &mut ct.c0_t,
        &mut ct.c0_q2,
        &mut ct.c1_t,
        &mut ct.c1_q2,
    ];
    let mut offset = 0;
    for arr in arrays {
        for i in 0..N {
            let bytes: [u8; 4] = data[offset..offset + 4].try_into().unwrap();
            arr[i] = u32::from_le_bytes(bytes);
            offset += 4;
        }
    }
    ct
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serialize_roundtrip() {
        let mut ct = CiphertextSlot {
            c0_t: [0u32; N],
            c0_q2: [0u32; N],
            c1_t: [0u32; N],
            c1_q2: [0u32; N],
        };
        ct.c0_t[0] = 42;
        ct.c0_q2[100] = 123456;
        ct.c1_t[N - 1] = 999999;
        ct.c1_q2[2048] = 0xDEADBEEF;

        let bytes = serialize_slot(&ct);
        assert_eq!(bytes.len(), SLOT_CT_BYTES);

        let ct2 = deserialize_slot(&bytes);
        assert_eq!(ct.c0_t, ct2.c0_t);
        assert_eq!(ct.c0_q2, ct2.c0_q2);
        assert_eq!(ct.c1_t, ct2.c1_t);
        assert_eq!(ct.c1_q2, ct2.c1_q2);
    }
}
