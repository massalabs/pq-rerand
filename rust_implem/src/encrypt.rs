//! Encryption: base BFV encrypt + noise-flooding smudge.

use crate::params::{N, T, Q2, DELTA_T};
use crate::ntt::{pointwise_add, pointwise_mul, barrett_mu, barrett_reduce, addmod};
use crate::poly::{Poly, NttContext};
use crate::keygen::PublicKey;
use crate::sampling::{sample_gaussian, sample_gaussian_cdt};
use rand::{CryptoRng, Rng};
use zeroize::Zeroize;

/// A single ciphertext slot: (c0_t, c0_q2, c1_t, c1_q2), all in coefficient domain.
#[derive(Clone, Debug)]
pub struct CiphertextSlot {
    /// `c₀ mod t`, coefficient domain.
    pub c0_t: [u32; N],
    /// `c₀ mod q₂`, coefficient domain.
    pub c0_q2: [u32; N],
    /// `c₁ mod t`, coefficient domain.
    pub c1_t: [u32; N],
    /// `c₁ mod q₂`, coefficient domain.
    pub c1_q2: [u32; N],
}

/// Encrypt a single message polynomial M (coefficients in [0, 2^31 - 1]).
/// Performs base encryption + noise-flooding smudge.
///
/// The RNG must be cryptographically secure (enforced via the [`CryptoRng`]
/// marker bound). Transient encryption randomness is zeroized before returning.
pub fn encrypt_slot<R: Rng + CryptoRng>(
    rng: &mut R,
    ctx: &NttContext,
    pk: &PublicKey,
    message: &[u32; N],
    sigma_flood: f64,
) -> CiphertextSlot {
    // Sample base-width noise (constant-time CDT) and flood-width noise, then add
    // them coefficient-wise. The base encryption and the flooding Enc(0) are both
    // linear in their randomness/noise, so a separate smudge step
    //     C = (b·r+e₂+Δ·M, a·r+e₁) ,  C += (b·r_f+e_{2,f}, a·r_f+e_{1,f})
    // produces exactly the same distribution as a single encryption with combined
    // noise R = r+r_f, E₁ = e₁+e_{1,f}, E₂ = e₂+e_{2,f}. Folding the two together
    // halves the NTT work (one forward of R, one inverse per ciphertext component)
    // while leaving the output distribution — and hence the security analysis —
    // unchanged. Flooding therefore costs only its extra sampling, not extra NTTs.
    let mut r_b = sample_gaussian_cdt(rng);
    let mut e1_b = sample_gaussian_cdt(rng);
    let mut e2_b = sample_gaussian_cdt(rng);
    let mut r_f = sample_gaussian(rng, sigma_flood);
    let mut e1_f = sample_gaussian(rng, sigma_flood);
    let mut e2_f = sample_gaussian(rng, sigma_flood);

    let mut r_int = [0i64; N];
    let mut e1_int = [0i64; N];
    let mut e2_int = [0i64; N];
    for i in 0..N {
        r_int[i] = r_b[i] + r_f[i];
        e1_int[i] = e1_b[i] + e1_f[i];
        e2_int[i] = e2_b[i] + e2_f[i];
    }
    r_b.zeroize();
    e1_b.zeroize();
    e2_b.zeroize();
    r_f.zeroize();
    e1_f.zeroize();
    e2_f.zeroize();

    let mut r = Poly::from_signed(&r_int);
    let mut e1 = Poly::from_signed(&e1_int);
    let mut e2 = Poly::from_signed(&e2_int);
    r_int.zeroize();
    e1_int.zeroize();
    e2_int.zeroize();

    // NTT of the combined randomness R = r + r_f.
    let (mut r_ntt_t, mut r_ntt_q2) = ctx.forward(&r);

    // c1 = a·R + E1 (in each limb)
    let mut c1_t = pointwise_mul(&pk.a_ntt_t, &r_ntt_t, T);
    let mut c1_q2 = pointwise_mul(&pk.a_ntt_q2, &r_ntt_q2, Q2);
    ctx.tables_t.inverse(&mut c1_t);
    ctx.tables_q2.inverse(&mut c1_q2);
    c1_t = pointwise_add(&c1_t, &e1.limb_t, T);
    c1_q2 = pointwise_add(&c1_q2, &e1.limb_q2, Q2);

    // c0 = b·R + E2 + Δ_t·M (t-limb) or b·R + E2 + 0 (q2-limb)
    let mut c0_t = pointwise_mul(&pk.b_ntt_t, &r_ntt_t, T);
    let mut c0_q2 = pointwise_mul(&pk.b_ntt_q2, &r_ntt_q2, Q2);
    ctx.tables_t.inverse(&mut c0_t);
    ctx.tables_q2.inverse(&mut c0_q2);
    c0_t = pointwise_add(&c0_t, &e2.limb_t, T);
    c0_q2 = pointwise_add(&c0_q2, &e2.limb_q2, Q2);

    // Add Δ_t · M to c0_t (message embedding — message lives in t-limb only).
    // Division-free reduction via Barrett.
    let mu_t = barrett_mu(T);
    for i in 0..N {
        let delta_m = barrett_reduce(DELTA_T as u64 * message[i] as u64, T, mu_t);
        c0_t[i] = addmod(c0_t[i], delta_m, T);
    }

    // Best-effort scrubbing of the remaining transient secrets (encryption
    // randomness in coefficient and NTT domain).
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
