//! Key generation.

use crate::params::{N, T, Q2};
use crate::poly::{Poly, NttContext, NTT_CONTEXT};
use crate::sampling::{sample_gaussian_cdt, sample_uniform};
use rand::{CryptoRng, Rng};
use zeroize::{Zeroize, Zeroizing, ZeroizeOnDrop};

/// Secret key: the secret polynomial s reduced into both CRT limbs.
///
/// `s_t`/`s_q2` hold the coefficient-domain limbs (serialized form); `s_ntt_t`/
/// `s_ntt_q2` cache the forward NTT of `s` so that decryption — which only ever
/// needs `s` in the NTT domain — avoids re-transforming the key on every call.
///
/// Zeroized on drop to prevent secret material from lingering in memory.
/// Does not implement `Debug` to prevent accidental logging of secrets.
#[derive(Clone, Zeroize, ZeroizeOnDrop)]
pub struct SecretKey {
    /// `s mod t`, coefficient domain.
    pub s_t: [u32; N],
    /// `s mod q₂`, coefficient domain.
    pub s_q2: [u32; N],
    /// `s mod t`, NTT domain (cached for decryption).
    pub s_ntt_t: [u32; N],
    /// `s mod q₂`, NTT domain (cached for decryption).
    pub s_ntt_q2: [u32; N],
}

impl SecretKey {
    /// Build a secret key from its coefficient-domain limbs, computing and caching
    /// the NTT-domain representation used by decryption.
    fn from_limbs(s_t: [u32; N], s_q2: [u32; N]) -> Self {
        let mut s_ntt_t = s_t;
        let mut s_ntt_q2 = s_q2;
        NTT_CONTEXT.tables_t.forward(&mut s_ntt_t);
        NTT_CONTEXT.tables_q2.forward(&mut s_ntt_q2);
        SecretKey { s_t, s_q2, s_ntt_t, s_ntt_q2 }
    }
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
    /// Returns `None` if `data` is not exactly `BYTES` long or if any limb
    /// coefficient is out of range (not fully reduced modulo `t` / `q₂`).
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() != Self::BYTES {
            return None;
        }
        let mut s_t = [0u32; N];
        let mut s_q2 = [0u32; N];
        let arrays: [(&mut [u32; N], u32); 2] = [(&mut s_t, T), (&mut s_q2, Q2)];
        let mut offset = 0;
        for (arr, modulus) in arrays {
            for val in arr.iter_mut() {
                *val = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?);
                if *val >= modulus {
                    return None;
                }
                offset += 4;
            }
        }
        Some(SecretKey::from_limbs(s_t, s_q2))
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
    /// Returns `None` if `data` is not exactly `BYTES` long or if any limb
    /// coefficient is out of range (not fully reduced modulo `t` / `q₂`).
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
        let arrays: [(&mut [u32; N], u32); 4] = [
            (&mut pk.a_ntt_t, T),
            (&mut pk.a_ntt_q2, Q2),
            (&mut pk.b_ntt_t, T),
            (&mut pk.b_ntt_q2, Q2),
        ];
        let mut offset = 0;
        for (arr, modulus) in arrays {
            for val in arr.iter_mut() {
                *val = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?);
                if *val >= modulus {
                    return None;
                }
                offset += 4;
            }
        }
        Some(pk)
    }
}

/// Generate a keypair.
///
/// The RNG must be cryptographically secure (enforced via the [`CryptoRng`]
/// marker bound). Transient secret material (the signed coefficients of `s` and
/// `e` and their CRT limbs) is zeroized before returning.
pub fn keygen<R: Rng + CryptoRng>(rng: &mut R, ctx: &NttContext) -> (SecretKey, PublicKey) {
    // a is uniform per modulus
    let a_t = sample_uniform(rng, T);
    let a_q2 = sample_uniform(rng, Q2);

    // s, e are small signed integers (constant-time base-width sampler)
    let mut s_int = sample_gaussian_cdt(rng);
    let mut e_int = sample_gaussian_cdt(rng);

    let mut s = Poly::from_signed(&s_int);
    let mut e = Poly::from_signed(&e_int);
    s_int.zeroize();
    e_int.zeroize();

    // Forward NTT of a
    let mut a_ntt_t = a_t;
    let mut a_ntt_q2 = a_q2;
    ctx.tables_t.forward(&mut a_ntt_t);
    ctx.tables_q2.forward(&mut a_ntt_q2);

    // b = a*s + e in each limb (NTT domain for a*s, then INTT, then add e, then NTT)
    // Simpler: compute in coefficient domain via ring_mul
    let a_poly = Poly { limb_t: a_t, limb_q2: a_q2 };
    let mut as_poly = ctx.ring_mul(&a_poly, &s);
    let b_poly = as_poly.add(&e);
    as_poly.limb_t.zeroize();
    as_poly.limb_q2.zeroize();
    e.limb_t.zeroize();
    e.limb_q2.zeroize();

    // Store b in NTT domain
    let (b_ntt_t, b_ntt_q2) = ctx.forward(&b_poly);

    // Cache s in the NTT domain so decryption never re-transforms the key.
    let sk = SecretKey::from_limbs(s.limb_t, s.limb_q2);
    s.limb_t.zeroize();
    s.limb_q2.zeroize();

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
