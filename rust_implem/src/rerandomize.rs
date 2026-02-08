//! Public re-randomization: C' = C + Enc(0) with narrow noise.

use crate::params::{T, Q2, SIGMA};
use crate::ntt::{pointwise_add, pointwise_mul};
use crate::poly::{Poly, NttContext};
use crate::keygen::PublicKey;
use crate::encrypt::CiphertextSlot;
use crate::sampling::sample_gaussian;
use rand::Rng;

/// Re-randomize a single ciphertext slot by adding a fresh Enc(0).
pub fn rerandomize_slot<R: Rng>(
    rng: &mut R,
    ctx: &NttContext,
    pk: &PublicKey,
    ct: &CiphertextSlot,
) -> CiphertextSlot {
    let r_int = sample_gaussian(rng, SIGMA);
    let e1_int = sample_gaussian(rng, SIGMA);
    let e2_int = sample_gaussian(rng, SIGMA);

    let r = Poly::from_signed(&r_int);
    let e1 = Poly::from_signed(&e1_int);
    let e2 = Poly::from_signed(&e2_int);

    let (r_ntt_t, r_ntt_q2) = ctx.forward(&r);

    // Enc(0) = (b*r + e2, a*r + e1) in each limb
    let br_ntt_t = pointwise_mul(&pk.b_ntt_t, &r_ntt_t, T);
    let br_ntt_q2 = pointwise_mul(&pk.b_ntt_q2, &r_ntt_q2, Q2);
    let mut br_t = br_ntt_t;
    let mut br_q2 = br_ntt_q2;
    ctx.tables_t.inverse(&mut br_t);
    ctx.tables_q2.inverse(&mut br_q2);

    let ar_ntt_t = pointwise_mul(&pk.a_ntt_t, &r_ntt_t, T);
    let ar_ntt_q2 = pointwise_mul(&pk.a_ntt_q2, &r_ntt_q2, Q2);
    let mut ar_t = ar_ntt_t;
    let mut ar_q2 = ar_ntt_q2;
    ctx.tables_t.inverse(&mut ar_t);
    ctx.tables_q2.inverse(&mut ar_q2);

    // c0' = c0 + b*r + e2
    let c0_t = pointwise_add(&ct.c0_t, &pointwise_add(&br_t, &e2.limb_t, T), T);
    let c0_q2 = pointwise_add(&ct.c0_q2, &pointwise_add(&br_q2, &e2.limb_q2, Q2), Q2);

    // c1' = c1 + a*r + e1
    let c1_t = pointwise_add(&ct.c1_t, &pointwise_add(&ar_t, &e1.limb_t, T), T);
    let c1_q2 = pointwise_add(&ct.c1_q2, &pointwise_add(&ar_q2, &e1.limb_q2, Q2), Q2);

    CiphertextSlot { c0_t, c0_q2, c1_t, c1_q2 }
}
