//! Key generation.

use crate::params::{FORMAT_VERSION, N, T, Q2};
use crate::ntt::{pointwise_add, pointwise_mul, pointwise_sub};
use crate::poly::{NttContext, NTT_CONTEXT, Secret, is_small, secret, set_signed, zeroed};
use crate::sampling::{sample_gaussian_cdt, sample_uniform};
use rand_core::{CryptoRng, Rng};
use zeroize::{Zeroize, Zeroizing, ZeroizeOnDrop};

/// Secret key: the secret polynomial s reduced into both CRT limbs.
///
/// `s_t`/`s_q2` hold the coefficient-domain limbs (serialized form); `s_ntt_t`/
/// `s_ntt_q2` cache the forward NTT of `s` so that decryption — which only ever
/// needs `s` in the NTT domain — avoids re-transforming the key on every call.
///
/// The limbs live on the heap, so moving a key copies no key material, and they
/// are zeroized on drop to prevent secret material from lingering in memory.
/// Does not implement `Debug` to prevent accidental logging of secrets.
pub struct SecretKey {
    /// `s mod t`, coefficient domain.
    pub(crate) s_t: Secret,
    /// `s mod q₂`, coefficient domain.
    pub(crate) s_q2: Secret,
    /// `s mod t`, NTT domain (cached for decryption).
    pub(crate) s_ntt_t: Secret,
    /// `s mod q₂`, NTT domain (cached for decryption).
    pub(crate) s_ntt_q2: Secret,
}

impl Zeroize for SecretKey {
    fn zeroize(&mut self) {
        for limb in [&mut self.s_t, &mut self.s_q2, &mut self.s_ntt_t, &mut self.s_ntt_q2] {
            limb.zeroize();
        }
    }
}

impl ZeroizeOnDrop for SecretKey {}

impl Clone for SecretKey {
    /// Copy the key into new zeroizing storage. Every clone is one more copy of
    /// the key in memory: share a reference where one is enough.
    fn clone(&self) -> Self {
        let mut sk = Self::zeroed();
        sk.s_t.copy_from_slice(&*self.s_t);
        sk.s_q2.copy_from_slice(&*self.s_q2);
        sk.s_ntt_t.copy_from_slice(&*self.s_ntt_t);
        sk.s_ntt_q2.copy_from_slice(&*self.s_ntt_q2);
        sk
    }
}

impl SecretKey {
    /// An all-zero key, to be filled in place.
    fn zeroed() -> Self {
        SecretKey { s_t: secret(), s_q2: secret(), s_ntt_t: secret(), s_ntt_q2: secret() }
    }

    /// Compute and cache the NTT-domain representation used by decryption from
    /// the coefficient-domain limbs.
    fn cache_ntt(&mut self, ctx: &NttContext) {
        self.s_ntt_t.copy_from_slice(&*self.s_t);
        self.s_ntt_q2.copy_from_slice(&*self.s_q2);
        ctx.tables_t.forward(&mut self.s_ntt_t);
        ctx.tables_q2.forward(&mut self.s_ntt_q2);
    }
}

impl SecretKey {
    /// Serialized byte size.
    pub const BYTES: usize = 1 + 2 * N * 4;

    /// Serialize to bytes ([`FORMAT_VERSION`], then little-endian u32 arrays).
    ///
    /// The returned buffer is zeroized on drop.
    #[must_use]
    pub fn to_bytes(&self) -> Zeroizing<Vec<u8>> {
        let mut buf = Vec::with_capacity(Self::BYTES);
        buf.push(FORMAT_VERSION);
        for arr in [&self.s_t, &self.s_q2] {
            for &val in arr.iter() {
                buf.extend_from_slice(&val.to_le_bytes());
            }
        }
        Zeroizing::new(buf)
    }

    /// Deserialize from bytes.
    ///
    /// Returns `None` if `data` is not exactly `BYTES` long, does not start with
    /// [`FORMAT_VERSION`] or does not hold a secret as key generation produces
    /// them: the two limbs must be the reductions modulo `t` and `q2` of one
    /// non-zero polynomial with coefficients of magnitude at most 29 and a squared
    /// norm of at most `2 N sigma^2`. Both limbs are read in full before the
    /// verdict on them, in constant time, and a rejected key is zeroized.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() != Self::BYTES || data[0] != FORMAT_VERSION {
            return None;
        }
        let mut sk = Self::zeroed();
        let (bytes_t, bytes_q2) = data[1..].split_at(N * 4);
        for (limb, bytes) in [(&mut sk.s_t, bytes_t), (&mut sk.s_q2, bytes_q2)] {
            for (val, word) in limb.iter_mut().zip(bytes.as_chunks::<4>().0) {
                *val = u32::from_le_bytes(*word);
            }
        }
        if !is_small(&sk.s_t, &sk.s_q2) {
            return None;
        }
        sk.cache_ntt(&NTT_CONTEXT);
        Some(sk)
    }
}

/// Public key: (a, b = a·s + e) in CRT form, stored in NTT domain for efficiency.
///
/// Public keys are not secret but they must be **authentic**: whoever can replace
/// a public key can read what is then encrypted under it and can link the
/// ciphertexts re-randomized under it. No check on the key alone can tell an
/// honest key from one whose secret someone else knows. Check a public key
/// against its secret key with [`PublicKey::verify_against`] where the secret key
/// is at hand, and authenticate every other public key, including those only
/// used to re-randomize, by other means (for example a hash kept in authenticated
/// storage).
///
/// The serialized form is the byte [`FORMAT_VERSION`], then the four arrays
/// below, each coefficient a little-endian `u32`. The NTT-domain order is part of
/// that format: index `k` holds the evaluation at `psi^(2*bitrev(k)+1)`, where
/// `bitrev` reverses 12 bits and `psi` is the 8192nd root of unity 101,573,657
/// modulo `t`, or 1,953,722,822 modulo `q2`.
///
/// The limbs live on the heap, like those of the other keys and of ciphertexts:
/// returning or moving a key moves pointers, not 64 KiB through the stack.
#[derive(Clone)]
pub struct PublicKey {
    /// a mod t (NTT domain)
    pub(crate) a_ntt_t: Box<[u32; N]>,
    /// a mod q₂ (NTT domain)
    pub(crate) a_ntt_q2: Box<[u32; N]>,
    /// b mod t (NTT domain)
    pub(crate) b_ntt_t: Box<[u32; N]>,
    /// b mod q₂ (NTT domain)
    pub(crate) b_ntt_q2: Box<[u32; N]>,
}

impl PublicKey {
    /// An all-zero key, to be filled in place.
    fn zeroed() -> Self {
        PublicKey { a_ntt_t: zeroed(), a_ntt_q2: zeroed(), b_ntt_t: zeroed(), b_ntt_q2: zeroed() }
    }

    /// Serialized byte size.
    pub const BYTES: usize = 1 + 4 * N * 4;

    /// Serialize to bytes ([`FORMAT_VERSION`], then little-endian u32 arrays).
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(Self::BYTES);
        buf.push(FORMAT_VERSION);
        for arr in [&self.a_ntt_t, &self.a_ntt_q2, &self.b_ntt_t, &self.b_ntt_q2] {
            for &val in arr.iter() {
                buf.extend_from_slice(&val.to_le_bytes());
            }
        }
        buf
    }

    /// Deserialize from bytes.
    ///
    /// Returns `None` if `data` is not exactly `BYTES` long, if it does not start
    /// with [`FORMAT_VERSION`], if any limb coefficient is out of range (not fully
    /// reduced modulo `t` / `q2`), or if the key is degenerate: `a` must be
    /// invertible (no zero coefficient in the NTT domain) and `b` must not be
    /// zero, in both limbs. This rejects in particular the all-zero key and keys
    /// with `a = 0` or `b = 0`, under which anyone can decrypt. Key generation
    /// never outputs such a key (before version 0.3 of the crate it did, for about
    /// one key in 350,000: such a key pair has to be generated again).
    ///
    /// This is a sanity check, not an authentication: see [`PublicKey`].
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() != Self::BYTES || data[0] != FORMAT_VERSION {
            return None;
        }
        let mut pk = PublicKey::zeroed();
        let arrays: [(&mut [u32; N], u32); 4] = [
            (&mut pk.a_ntt_t, T),
            (&mut pk.a_ntt_q2, Q2),
            (&mut pk.b_ntt_t, T),
            (&mut pk.b_ntt_q2, Q2),
        ];
        let mut offset = 1;
        for (arr, modulus) in arrays {
            for val in arr.iter_mut() {
                *val = u32::from_le_bytes(data[offset..offset + 4].try_into().ok()?);
                if *val >= modulus {
                    return None;
                }
                offset += 4;
            }
        }
        let b_is_nonzero = |limb: &[u32; N]| limb.iter().any(|&x| x != 0);
        if !pk.a_is_invertible() || !b_is_nonzero(&pk.b_ntt_t) || !b_is_nonzero(&pk.b_ntt_q2) {
            return None;
        }
        Some(pk)
    }

    /// Whether `a` is invertible: none of its NTT-domain coefficients is zero.
    fn a_is_invertible(&self) -> bool {
        [&self.a_ntt_t, &self.a_ntt_q2].iter().all(|limb| limb.iter().all(|&x| x != 0))
    }

    /// Whether this public key is a key of the secret key `sk`: `b - a*s` must be
    /// an error polynomial as key generation draws them (non-zero, coefficients of
    /// magnitude at most 29, squared norm at most `2 N sigma^2`), the same in both
    /// limbs.
    ///
    /// A public key that passes behaves as the one generated with `sk`: it was not
    /// replaced by the key of another secret, nor by a degenerate key (`a` is
    /// invertible in every `PublicKey`). Run it on a stored public key before
    /// encrypting under it for its own secret key. The scan is constant-time: the
    /// error it recovers is as secret as the key.
    #[must_use]
    pub fn verify_against(&self, sk: &SecretKey) -> bool {
        let ctx = &NTT_CONTEXT;
        let (mut e_t, mut e_q2) = (secret(), secret());
        // a*s - b = -e, back in the coefficient domain.
        e_t.copy_from_slice(&*self.a_ntt_t);
        e_q2.copy_from_slice(&*self.a_ntt_q2);
        pointwise_mul(&mut e_t, &sk.s_ntt_t, T);
        pointwise_mul(&mut e_q2, &sk.s_ntt_q2, Q2);
        pointwise_sub(&mut e_t, &self.b_ntt_t, T);
        pointwise_sub(&mut e_q2, &self.b_ntt_q2, Q2);
        ctx.tables_t.inverse(&mut e_t);
        ctx.tables_q2.inverse(&mut e_q2);
        is_small(&e_t, &e_q2)
    }
}

/// Generate a keypair.
///
/// The RNG must be cryptographically secure (enforced via the [`CryptoRng`]
/// marker bound). The secret `s` and the error `e` only ever exist in zeroizing
/// heap storage; the product `a*s` is formed in the public key itself, where
/// adding `e` overwrites it.
///
/// The generator is read in a fixed order (`a mod t`, `a mod q2`, both drawn again
/// while `a` has a zero NTT coefficient, then `s`, `e`; see
/// [`crate::sampling`]), so the key pair is a deterministic function of the
/// generator stream: a seeded generator with a specified algorithm, such as
/// ChaCha20, derives the same key pair on every platform.
pub fn keygen<R: Rng + CryptoRng>(rng: &mut R, ctx: &NttContext) -> (SecretKey, PublicKey) {
    let mut sk = SecretKey::zeroed();
    let mut pk = PublicKey::zeroed();
    // a is uniform per modulus, among the invertible elements: it is drawn again
    // in the rare case (probability 2^-18.4) where one of its NTT coefficients
    // is zero. Like the rejection sampling itself, this only depends on `a`,
    // which is public.
    loop {
        sample_uniform(rng, T, &mut pk.a_ntt_t);
        sample_uniform(rng, Q2, &mut pk.a_ntt_q2);
        ctx.tables_t.forward(&mut pk.a_ntt_t);
        ctx.tables_q2.forward(&mut pk.a_ntt_q2);
        if pk.a_is_invertible() {
            break;
        }
    }

    // s, e are small signed integers (constant-time base-width sampler)
    let mut small = secret::<i32>();
    let (mut e_t, mut e_q2) = (secret(), secret());
    sample_gaussian_cdt(rng, &mut small);
    set_signed(&mut sk.s_t, &mut sk.s_q2, &small);
    sk.cache_ntt(ctx);
    sample_gaussian_cdt(rng, &mut small);
    set_signed(&mut e_t, &mut e_q2, &small);
    ctx.tables_t.forward(&mut e_t);
    ctx.tables_q2.forward(&mut e_q2);

    // b = a*s + e in each limb, evaluated in the NTT domain (the transform is
    // linear), where the public key is stored.
    pk.b_ntt_t.copy_from_slice(&*pk.a_ntt_t);
    pk.b_ntt_q2.copy_from_slice(&*pk.a_ntt_q2);
    pointwise_mul(&mut pk.b_ntt_t, &sk.s_ntt_t, T);
    pointwise_mul(&mut pk.b_ntt_q2, &sk.s_ntt_q2, Q2);
    pointwise_add(&mut pk.b_ntt_t, &e_t, T);
    pointwise_add(&mut pk.b_ntt_q2, &e_q2, Q2);

    (sk, pk)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sampling::MAX_BASE_SAMPLE;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rand_core::{Infallible, TryCryptoRng, TryRng};

    fn test_keypair() -> (SecretKey, PublicKey) {
        let mut rng = StdRng::seed_from_u64(42);
        keygen(&mut rng, &NTT_CONTEXT)
    }

    /// Replace one coefficient of a serialized limb (the limbs follow the
    /// version byte).
    fn patch(bytes: &mut [u8], limb: usize, index: usize, value: u32) {
        let at = 1 + (limb * N + index) * 4;
        bytes[at..at + 4].copy_from_slice(&value.to_le_bytes());
    }

    fn word(bytes: &[u8], limb: usize, index: usize) -> u32 {
        let at = 1 + (limb * N + index) * 4;
        u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap())
    }

    /// A serialized key starts with the version byte, and no other value of that
    /// byte is accepted.
    #[test]
    fn serialized_keys_carry_the_format_version() {
        let (sk, pk) = test_keypair();
        let (mut sk_bytes, mut pk_bytes) = (sk.to_bytes().to_vec(), pk.to_bytes());
        assert_eq!((FORMAT_VERSION, sk_bytes[0], pk_bytes[0]), (1, 1, 1));
        for version in (0..=u8::MAX).filter(|&v| v != FORMAT_VERSION) {
            (sk_bytes[0], pk_bytes[0]) = (version, version);
            assert!(SecretKey::from_bytes(&sk_bytes).is_none(), "secret key of version {version} accepted");
            assert!(PublicKey::from_bytes(&pk_bytes).is_none(), "public key of version {version} accepted");
        }
        // Without the byte, the limbs alone are not a key.
        assert!(SecretKey::from_bytes(&sk_bytes[1..]).is_none());
        assert!(PublicKey::from_bytes(&pk_bytes[1..]).is_none());
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
        assert_eq!(*sk.s_t, *sk2.s_t);
        assert_eq!(*sk.s_q2, *sk2.s_q2);
        // The cache of the parsed key is the one key generation computed.
        assert_eq!(*sk.s_ntt_t, *sk2.s_ntt_t);
        assert_eq!(*sk.s_ntt_q2, *sk2.s_ntt_q2);
        // A clone is a deep copy, and zeroizing clears every limb.
        let mut sk3 = sk.clone();
        assert_eq!(sk3.to_bytes(), bytes);
        assert_eq!(*sk3.s_ntt_t, *sk.s_ntt_t);
        assert_eq!(*sk3.s_ntt_q2, *sk.s_ntt_q2);
        sk3.zeroize();
        for limb in [&sk3.s_t, &sk3.s_q2, &sk3.s_ntt_t, &sk3.s_ntt_q2] {
            assert!(limb.iter().all(|&x| x == 0));
        }
    }

    #[test]
    fn pk_from_invalid_bytes() {
        assert!(PublicKey::from_bytes(&[]).is_none());
        assert!(PublicKey::from_bytes(&[0u8; 10]).is_none());
        let (_sk, pk) = test_keypair();
        let good = pk.to_bytes();
        assert!(PublicKey::from_bytes(&good[..good.len() - 1]).is_none());
        assert!(PublicKey::from_bytes(&[&good[..], &[0]].concat()).is_none());
        // A coefficient at or above its modulus, in each limb.
        for (limb, modulus) in [(0, T), (1, Q2), (2, T), (3, Q2)] {
            for value in [modulus, u32::MAX] {
                let mut bytes = good.clone();
                patch(&mut bytes, limb, N - 1, value);
                assert!(PublicKey::from_bytes(&bytes).is_none(), "limb {limb}: {value} accepted");
            }
            let mut bytes = good.clone();
            patch(&mut bytes, limb, N - 1, modulus - 1);
            assert!(PublicKey::from_bytes(&bytes).is_some(), "limb {limb}: maximum rejected");
        }
    }

    /// Degenerate keys, under which anyone can decrypt, are rejected: the all-zero
    /// key, a = 0 or b = 0 in either limb, and a zero coefficient anywhere in a.
    #[test]
    fn pk_rejects_degenerate_keys() {
        // Limbs filled with one byte value, after the version byte.
        let filled = |byte: u8| [&[FORMAT_VERSION][..], &vec![byte; PublicKey::BYTES - 1]].concat();
        assert!(PublicKey::from_bytes(&filled(0)).is_none());
        assert!(PublicKey::from_bytes(&filled(1)).is_some());
        let (_sk, pk) = test_keypair();
        let good = pk.to_bytes();
        for limb in 0..4 {
            let mut bytes = good.clone();
            bytes[1..][limb * N * 4..(limb + 1) * N * 4].fill(0);
            assert!(PublicKey::from_bytes(&bytes).is_none(), "zero limb {limb} accepted");
            for index in [0, N - 1] {
                let mut bytes = good.clone();
                patch(&mut bytes, limb, index, 0);
                // A single zero makes a non-invertible; in b it means nothing.
                assert_eq!(PublicKey::from_bytes(&bytes).is_none(), limb < 2, "zero in limb {limb}");
            }
        }
    }

    #[test]
    fn sk_from_invalid_bytes() {
        assert!(SecretKey::from_bytes(&[]).is_none());
        assert!(SecretKey::from_bytes(&[0u8; 10]).is_none());
        let (sk, _pk) = test_keypair();
        let good = sk.to_bytes();
        assert!(SecretKey::from_bytes(&good[..good.len() - 1]).is_none());
        assert!(SecretKey::from_bytes(&[&good[..], &[0u8]].concat()).is_none());
    }

    /// The limbs must be the reductions of one polynomial with small coefficients.
    #[test]
    fn sk_rejects_inconsistent_or_large_limbs() {
        let (sk, _pk) = test_keypair();
        let good = sk.to_bytes().to_vec();
        let max = MAX_BASE_SAMPLE;

        // One limb perturbed, at the first and at the last coefficient.
        for (limb, modulus) in [(0, T), (1, Q2)] {
            for index in [0, N - 1] {
                for delta in [1, 12_345] {
                    let mut bytes = good.clone();
                    let v = (word(&good, limb, index) as u64 + delta) % modulus as u64;
                    patch(&mut bytes, limb, index, v as u32);
                    assert!(SecretKey::from_bytes(&bytes).is_none(), "limb {limb} + {delta} accepted");
                }
            }
        }
        // Both limbs agree on a coefficient: accepted up to the largest magnitude
        // that key generation produces, rejected just beyond it.
        for (value, valid) in [(max, true), (max + 1, false), (1 << 20, false)] {
            for (vt, vq2) in [(value, value), (T - value, Q2 - value)] {
                let mut bytes = good.clone();
                patch(&mut bytes, 0, 5, vt);
                patch(&mut bytes, 1, 5, vq2);
                assert_eq!(SecretKey::from_bytes(&bytes).is_some(), valid, "coefficient {value}");
            }
        }
        // Out-of-range words and uniformly random limbs.
        for (limb, value) in [(0, T), (0, u32::MAX), (1, Q2), (1, u32::MAX)] {
            let mut bytes = good.clone();
            patch(&mut bytes, limb, 9, value);
            assert!(SecretKey::from_bytes(&bytes).is_none());
        }
        let mut rng = StdRng::seed_from_u64(5);
        let mut random = vec![FORMAT_VERSION; SecretKey::BYTES];
        rand::Rng::fill_bytes(&mut rng, &mut random[1..]);
        for i in 0..2 * N {
            let reduced = word(&random, 0, i) % if i < N { T } else { Q2 };
            patch(&mut random, 0, i, reduced);
        }
        assert!(SecretKey::from_bytes(&random).is_none(), "random limbs accepted");
    }

    /// A public key verifies against its own secret key only.
    #[test]
    fn verify_against_accepts_only_the_matching_key() {
        let mut rng = StdRng::seed_from_u64(13);
        let (sk, pk) = keygen(&mut rng, &NTT_CONTEXT);
        let (other_sk, other_pk) = keygen(&mut rng, &NTT_CONTEXT);
        assert!(pk.verify_against(&sk));
        assert!(other_pk.verify_against(&other_sk));
        assert!(!other_pk.verify_against(&sk), "foreign key accepted");
        assert!(!pk.verify_against(&other_sk), "foreign secret accepted");
        // The parsed keys verify too.
        let pk2 = PublicKey::from_bytes(&pk.to_bytes()).unwrap();
        let sk2 = SecretKey::from_bytes(&sk.to_bytes()).unwrap();
        assert!(pk2.verify_against(&sk2));

        // One coefficient of b changed: the error is no longer small.
        let mut tampered = pk.clone();
        tampered.b_ntt_t[0] = crate::ntt::addmod(tampered.b_ntt_t[0], 1, T);
        assert!(!tampered.verify_against(&sk));
        let mut tampered = pk.clone();
        tampered.b_ntt_q2[N - 1] = crate::ntt::addmod(tampered.b_ntt_q2[N - 1], 1, Q2);
        assert!(!tampered.verify_against(&sk));

        // An error that is small but differs between the limbs.
        let mut one = [0u32; N];
        one[0] = 1;
        NTT_CONTEXT.tables_q2.forward(&mut one);
        let mut tampered = pk.clone();
        pointwise_add(&mut tampered.b_ntt_q2, &one, Q2);
        assert!(!tampered.verify_against(&sk));

        // An error just past the bound, the same in both limbs: b + (max + 1 - e_0).
        // The key itself has errors within the bound, so shifting b by a constant
        // moves coefficient 0 of the error only.
        let max = MAX_BASE_SAMPLE as i32;
        let accepted: Vec<i32> = (-2 * max - 1..=2 * max + 1)
            .filter(|&shift| {
                let mut c = [0i32; N];
                c[0] = shift;
                let (mut c_t, mut c_q2) = ([0u32; N], [0u32; N]);
                set_signed(&mut c_t, &mut c_q2, &c);
                NTT_CONTEXT.tables_t.forward(&mut c_t);
                NTT_CONTEXT.tables_q2.forward(&mut c_q2);
                let mut shifted = pk.clone();
                pointwise_add(&mut shifted.b_ntt_t, &c_t, T);
                pointwise_add(&mut shifted.b_ntt_q2, &c_q2, Q2);
                shifted.verify_against(&sk)
            })
            .collect();
        // Exactly the shifts that keep e_0 + shift within [-max, max]: a run of
        // 2 * max + 1 consecutive values that contains 0.
        assert_eq!(accepted.len() as i32, 2 * max + 1);
        assert!(accepted.contains(&0));
        assert_eq!(accepted[accepted.len() - 1] - accepted[0], 2 * max);
    }

    /// Generator whose first `N` 32-bit draws are `prefix[0]` and whose next `N`
    /// are `prefix[1]`, then follows `inner`.
    struct Prefixed {
        prefix: [u32; 2],
        drawn: usize,
        inner: StdRng,
    }

    impl TryRng for Prefixed {
        type Error = Infallible;
        fn try_next_u32(&mut self) -> Result<u32, Infallible> {
            if self.drawn < 2 * N {
                self.drawn += 1;
                return Ok(self.prefix[(self.drawn - 1) / N]);
            }
            self.inner.try_next_u32()
        }
        fn try_next_u64(&mut self) -> Result<u64, Infallible> {
            self.inner.try_next_u64()
        }
        fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), Infallible> {
            self.inner.try_fill_bytes(dst)
        }
    }

    impl TryCryptoRng for Prefixed {}

    /// A candidate a that is not invertible (here a = 0 in both limbs, modulo t
    /// only and modulo q2 only) is discarded and drawn again from the generator.
    #[test]
    fn keygen_discards_non_invertible_a() {
        for prefix in [[0, 0], [0, 1], [1, 0]] {
            let mut rng = Prefixed { prefix, drawn: 0, inner: StdRng::seed_from_u64(77) };
            let (sk, pk) = keygen(&mut rng, &NTT_CONTEXT);
            assert!(pk.a_is_invertible());
            assert!(pk.verify_against(&sk));
            assert!(PublicKey::from_bytes(&pk.to_bytes()).is_some());

            // Same key as a run on the generator alone.
            let (sk2, pk2) = keygen(&mut StdRng::seed_from_u64(77), &NTT_CONTEXT);
            assert_eq!(pk.to_bytes(), pk2.to_bytes());
            assert_eq!(sk.to_bytes(), sk2.to_bytes());
        }
    }

    /// The secret and the error behind a generated key are small, and the key
    /// is what the definition b = a*s + e gives.
    #[test]
    fn keygen_output_is_well_formed() {
        for seed in 0..4 {
            let mut rng = StdRng::seed_from_u64(seed);
            let (sk, pk) = keygen(&mut rng, &NTT_CONTEXT);
            assert!(is_small(&sk.s_t, &sk.s_q2));
            assert!(sk.s_t.iter().any(|&x| x != 0), "s = 0");
            assert!(PublicKey::from_bytes(&pk.to_bytes()).is_some());
            assert!(pk.verify_against(&sk));
        }
    }
}
