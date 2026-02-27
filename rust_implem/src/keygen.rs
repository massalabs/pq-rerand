//! Key generation.

use crate::params::{N, T, Q2, SIGMA};
use crate::poly::{Poly, NttContext};
use crate::sampling::{sample_gaussian, sample_uniform};
use rand::Rng;
use zeroize::{Zeroize, Zeroizing, ZeroizeOnDrop};

/// Secret key: the secret polynomial s reduced into both CRT limbs.
///
/// Zeroized on drop to prevent secret material from lingering in memory.
/// Does not implement `Debug` to prevent accidental logging of secrets.
#[derive(Clone, Zeroize, ZeroizeOnDrop)]
pub struct SecretKey {
    pub s_t: [u32; N],
    pub s_q2: [u32; N],
}

impl SecretKey {
    /// Serialized byte size.
    pub const BYTES: usize = 2 * N * 4;

    /// Serialize to bytes (little-endian u32 arrays).
    ///
    /// The returned buffer is zeroized on drop.
    #[must_use]
    pub fn to_bytes(&self) -> Zeroizing<Vec<u8>> {
        let mut buf = Vec::with_capacity(Self::BYTES);
        for arr in [&self.s_t, &self.s_q2] {
            for &val in arr.iter() {
                buf.extend_from_slice(&val.to_le_bytes());
            }
        }
        Zeroizing::new(buf)
    }

    /// Deserialize from bytes.
    ///
    /// Returns `None` if `data` is not exactly `BYTES` long.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() != Self::BYTES {
            return None;
        }
        let mut sk = SecretKey {
            s_t: [0u32; N],
            s_q2: [0u32; N],
        };
        let arrays: [&mut [u32; N]; 2] = [&mut sk.s_t, &mut sk.s_q2];
        let mut offset = 0;
        for arr in arrays {
            for val in arr.iter_mut() {
                *val = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?);
                offset += 4;
            }
        }
        Some(sk)
    }
}

/// Public key: (a, b = a·s + e) in CRT form, stored in NTT domain for efficiency.
#[derive(Clone)]
pub struct PublicKey {
    /// a mod t (NTT domain)
    pub a_ntt_t: [u32; N],
    /// a mod q₂ (NTT domain)
    pub a_ntt_q2: [u32; N],
    /// b mod t (NTT domain)
    pub b_ntt_t: [u32; N],
    /// b mod q₂ (NTT domain)
    pub b_ntt_q2: [u32; N],
}

impl PublicKey {
    /// Serialized byte size.
    pub const BYTES: usize = 4 * N * 4;

    /// Serialize to bytes (little-endian u32 arrays).
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(Self::BYTES);
        for arr in [&self.a_ntt_t, &self.a_ntt_q2, &self.b_ntt_t, &self.b_ntt_q2] {
            for &val in arr.iter() {
                buf.extend_from_slice(&val.to_le_bytes());
            }
        }
        buf
    }

    /// Deserialize from bytes.
    ///
    /// Returns `None` if `data` is not exactly `BYTES` long.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() != Self::BYTES {
            return None;
        }
        let mut pk = PublicKey {
            a_ntt_t: [0u32; N],
            a_ntt_q2: [0u32; N],
            b_ntt_t: [0u32; N],
            b_ntt_q2: [0u32; N],
        };
        let arrays: [&mut [u32; N]; 4] = [
            &mut pk.a_ntt_t,
            &mut pk.a_ntt_q2,
            &mut pk.b_ntt_t,
            &mut pk.b_ntt_q2,
        ];
        let mut offset = 0;
        for arr in arrays {
            for val in arr.iter_mut() {
                *val = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?);
                offset += 4;
            }
        }
        Some(pk)
    }
}

/// Generate a keypair.
pub fn keygen<R: Rng>(rng: &mut R, ctx: &NttContext) -> (SecretKey, PublicKey) {
    // a is uniform per modulus
    let a_t = sample_uniform(rng, T);
    let a_q2 = sample_uniform(rng, Q2);

    // s, e are small signed integers
    let s_int = sample_gaussian(rng, SIGMA);
    let e_int = sample_gaussian(rng, SIGMA);

    let s = Poly::from_signed(&s_int);
    let e = Poly::from_signed(&e_int);

    // Forward NTT of a
    let mut a_ntt_t = a_t;
    let mut a_ntt_q2 = a_q2;
    ctx.tables_t.forward(&mut a_ntt_t);
    ctx.tables_q2.forward(&mut a_ntt_q2);

    // Forward NTT of s (used only for potential future optimization)
    let (_s_ntt_t, _s_ntt_q2) = ctx.forward(&s);

    // b = a*s + e in each limb (NTT domain for a*s, then INTT, then add e, then NTT)
    // Simpler: compute in coefficient domain via ring_mul
    let a_poly = Poly { limb_t: a_t, limb_q2: a_q2 };
    let as_poly = ctx.ring_mul(&a_poly, &s);
    let b_poly = as_poly.add(&e);

    // Store b in NTT domain
    let (b_ntt_t, b_ntt_q2) = ctx.forward(&b_poly);

    let sk = SecretKey {
        s_t: s.limb_t,
        s_q2: s.limb_q2,
    };

    let pk = PublicKey {
        a_ntt_t,
        a_ntt_q2,
        b_ntt_t,
        b_ntt_q2,
    };

    (sk, pk)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn test_keypair() -> (SecretKey, PublicKey) {
        let ctx = NttContext::new();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        keygen(&mut rng, &ctx)
    }

    #[test]
    fn pk_serialization_roundtrip() {
        let (_sk, pk) = test_keypair();
        let bytes = pk.to_bytes();
        assert_eq!(bytes.len(), PublicKey::BYTES);
        let pk2 = PublicKey::from_bytes(&bytes).unwrap();
        assert_eq!(pk.a_ntt_t, pk2.a_ntt_t);
        assert_eq!(pk.a_ntt_q2, pk2.a_ntt_q2);
        assert_eq!(pk.b_ntt_t, pk2.b_ntt_t);
        assert_eq!(pk.b_ntt_q2, pk2.b_ntt_q2);
    }

    #[test]
    fn sk_serialization_roundtrip() {
        let (sk, _pk) = test_keypair();
        let bytes = sk.to_bytes();
        assert_eq!(bytes.len(), SecretKey::BYTES);
        let sk2 = SecretKey::from_bytes(&bytes).unwrap();
        assert_eq!(sk.s_t, sk2.s_t);
        assert_eq!(sk.s_q2, sk2.s_q2);
    }

    #[test]
    fn pk_from_invalid_bytes() {
        assert!(PublicKey::from_bytes(&[0u8; 10]).is_none());
    }

    #[test]
    fn sk_from_invalid_bytes() {
        assert!(SecretKey::from_bytes(&[0u8; 10]).is_none());
    }
}
