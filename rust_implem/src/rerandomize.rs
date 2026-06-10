//! Public re-randomization: C' = C + Enc(0) with narrow noise.

use crate::params::{T, Q2};
use crate::ntt::{pointwise_add, pointwise_mul};
use crate::poly::{Poly, NttContext};
use crate::keygen::PublicKey;
use crate::encrypt::CiphertextSlot;
use crate::sampling::sample_gaussian_cdt;
use rand::{CryptoRng, Rng};
use zeroize::Zeroize;

/// Re-randomize a single ciphertext slot by adding a fresh Enc(0).
///
/// The RNG must be cryptographically secure (enforced via the [`CryptoRng`]
/// marker bound). Transient re-randomization randomness is zeroized before
/// returning.
pub fn rerandomize_slot<R: Rng + CryptoRng>(
    rng: &mut R,
    ctx: &NttContext,
    pk: &PublicKey,
    ct: &CiphertextSlot,
) -> CiphertextSlot {
    let mut r_int = sample_gaussian_cdt(rng);
    let mut e1_int = sample_gaussian_cdt(rng);
    let mut e2_int = sample_gaussian_cdt(rng);

    let mut r = Poly::from_signed(&r_int);
    let mut e1 = Poly::from_signed(&e1_int);
    let mut e2 = Poly::from_signed(&e2_int);
    r_int.zeroize();
    e1_int.zeroize();
    e2_int.zeroize();

    let (mut r_ntt_t, mut r_ntt_q2) = ctx.forward(&r);

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

    // Best-effort scrubbing of the transient re-randomization secrets.
    r.limb_t.zeroize();
    r.limb_q2.zeroize();
    e1.limb_t.zeroize();
    e1.limb_q2.zeroize();
    e2.limb_t.zeroize();
    e2.limb_q2.zeroize();
    r_ntt_t.zeroize();
    r_ntt_q2.zeroize();

    CiphertextSlot { c0_t, c0_q2, c1_t, c1_q2 }
}
