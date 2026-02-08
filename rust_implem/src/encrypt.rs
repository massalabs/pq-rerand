//! Encryption: base BFV encrypt + noise-flooding smudge.

use crate::params::{N, T, Q2, SIGMA, DELTA_T};
use crate::ntt::{pointwise_add, pointwise_mul, mulmod, addmod};
use crate::poly::{Poly, NttContext};
use crate::keygen::PublicKey;
use crate::sampling::sample_gaussian;
use rand::Rng;

/// A single ciphertext slot: (c0_t, c0_q2, c1_t, c1_q2), all in coefficient domain.
#[derive(Clone, Debug)]
pub struct CiphertextSlot {
    pub c0_t: [u32; N],
    pub c0_q2: [u32; N],
    pub c1_t: [u32; N],
    pub c1_q2: [u32; N],
}

/// Encrypt a single message polynomial M (coefficients in [0, 2^31 - 1]).
/// Performs base encryption + noise-flooding smudge.
pub fn encrypt_slot<R: Rng>(
    rng: &mut R,
    ctx: &NttContext,
    pk: &PublicKey,
    message: &[u32; N],
    sigma_flood: f64,
) -> CiphertextSlot {
    // --- Phase 1: Base BFV encrypt ---
    let r_int = sample_gaussian(rng, SIGMA);
    let e1_int = sample_gaussian(rng, SIGMA);
    let e2_int = sample_gaussian(rng, SIGMA);

    let r = Poly::from_signed(&r_int);
    let e1 = Poly::from_signed(&e1_int);
    let e2 = Poly::from_signed(&e2_int);

    // NTT of r
    let (r_ntt_t, r_ntt_q2) = ctx.forward(&r);

    // c1 = a*r + e1 (in each limb)
    let ar_ntt_t = pointwise_mul(&pk.a_ntt_t, &r_ntt_t, T);
    let ar_ntt_q2 = pointwise_mul(&pk.a_ntt_q2, &r_ntt_q2, Q2);
    let mut c1_t = ar_ntt_t;
    let mut c1_q2 = ar_ntt_q2;
    ctx.tables_t.inverse(&mut c1_t);
    ctx.tables_q2.inverse(&mut c1_q2);
    c1_t = pointwise_add(&c1_t, &e1.limb_t, T);
    c1_q2 = pointwise_add(&c1_q2, &e1.limb_q2, Q2);

    // c0 = b*r + e2 + Δ_t·M (t-limb) or b*r + e2 + 0 (q2-limb)
    let br_ntt_t = pointwise_mul(&pk.b_ntt_t, &r_ntt_t, T);
    let br_ntt_q2 = pointwise_mul(&pk.b_ntt_q2, &r_ntt_q2, Q2);
    let mut c0_t = br_ntt_t;
    let mut c0_q2 = br_ntt_q2;
    ctx.tables_t.inverse(&mut c0_t);
    ctx.tables_q2.inverse(&mut c0_q2);
    c0_t = pointwise_add(&c0_t, &e2.limb_t, T);
    c0_q2 = pointwise_add(&c0_q2, &e2.limb_q2, Q2);

    // Add Δ_t · M to c0_t (message embedding — message lives in t-limb only)
    for i in 0..N {
        let delta_m = mulmod(DELTA_T, message[i], T);
        c0_t[i] = addmod(c0_t[i], delta_m, T);
    }

    // --- Phase 2: Noise-flooding smudge ---
    let rf_int = sample_gaussian(rng, sigma_flood);
    let e1f_int = sample_gaussian(rng, sigma_flood);
    let e2f_int = sample_gaussian(rng, sigma_flood);

    let rf = Poly::from_signed(&rf_int);
    let e1f = Poly::from_signed(&e1f_int);
    let e2f = Poly::from_signed(&e2f_int);

    let (rf_ntt_t, rf_ntt_q2) = ctx.forward(&rf);

    // c0 += b*rf + e2f
    let brf_ntt_t = pointwise_mul(&pk.b_ntt_t, &rf_ntt_t, T);
    let brf_ntt_q2 = pointwise_mul(&pk.b_ntt_q2, &rf_ntt_q2, Q2);
    let mut brf_t = brf_ntt_t;
    let mut brf_q2 = brf_ntt_q2;
    ctx.tables_t.inverse(&mut brf_t);
    ctx.tables_q2.inverse(&mut brf_q2);
    c0_t = pointwise_add(&c0_t, &pointwise_add(&brf_t, &e2f.limb_t, T), T);
    c0_q2 = pointwise_add(&c0_q2, &pointwise_add(&brf_q2, &e2f.limb_q2, Q2), Q2);

    // c1 += a*rf + e1f
    let arf_ntt_t = pointwise_mul(&pk.a_ntt_t, &rf_ntt_t, T);
    let arf_ntt_q2 = pointwise_mul(&pk.a_ntt_q2, &rf_ntt_q2, Q2);
    let mut arf_t = arf_ntt_t;
    let mut arf_q2 = arf_ntt_q2;
    ctx.tables_t.inverse(&mut arf_t);
    ctx.tables_q2.inverse(&mut arf_q2);
    c1_t = pointwise_add(&c1_t, &pointwise_add(&arf_t, &e1f.limb_t, T), T);
    c1_q2 = pointwise_add(&c1_q2, &pointwise_add(&arf_q2, &e1f.limb_q2, Q2), Q2);

    CiphertextSlot { c0_t, c0_q2, c1_t, c1_q2 }
}
