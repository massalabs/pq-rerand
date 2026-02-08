//! Key generation.

use crate::params::{N, T, Q2, SIGMA};
use crate::poly::{Poly, NttContext};
use crate::sampling::{sample_gaussian, sample_uniform};
use rand::Rng;

/// Secret key: the secret polynomial s reduced into both CRT limbs.
#[derive(Clone)]
pub struct SecretKey {
    pub s_t: [u32; N],
    pub s_q2: [u32; N],
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
