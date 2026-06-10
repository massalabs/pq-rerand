//! Decryption using the noise-limb trick.
//!
//! Because Δ·M ≡ 0 (mod q₂), the q₂-limb of v = c0 - c1·s is pure noise.
//! Center-lift it to recover the signed noise, then subtract from the t-limb
//! and multiply by INV_DELTA_T to recover M.

use crate::params::{N, T, Q2, INV_DELTA_T};
use crate::ntt::{barrett_mu, barrett_reduce, csub, submod};
use crate::poly::NttContext;
use crate::keygen::SecretKey;
use crate::encrypt::CiphertextSlot;
use zeroize::Zeroize;

/// Decrypt a single ciphertext slot, returning the message polynomial
/// coefficients in [0, t).
pub fn decrypt_slot(
    ctx: &NttContext,
    sk: &SecretKey,
    ct: &CiphertextSlot,
) -> [u32; N] {
    // Compute c1 * s in each limb via NTT. c1 is in coefficient domain, so it is
    // forward-transformed here; s is already cached in the NTT domain (see
    // `SecretKey`), saving two forward transforms per decryption.
    let mut c1_ntt_t = ct.c1_t;
    let mut c1_ntt_q2 = ct.c1_q2;
    ctx.tables_t.forward(&mut c1_ntt_t);
    ctx.tables_q2.forward(&mut c1_ntt_q2);

    // c1*s in NTT domain
    let mu_t = barrett_mu(T);
    let mu_q2 = barrett_mu(Q2);
    let mut cs_ntt_t = [0u32; N];
    let mut cs_ntt_q2 = [0u32; N];
    for i in 0..N {
        cs_ntt_t[i] = barrett_reduce(c1_ntt_t[i] as u64 * sk.s_ntt_t[i] as u64, T, mu_t);
        cs_ntt_q2[i] = barrett_reduce(c1_ntt_q2[i] as u64 * sk.s_ntt_q2[i] as u64, Q2, mu_q2);
    }

    // INTT to get c1*s in coefficient domain
    ctx.tables_t.inverse(&mut cs_ntt_t);
    ctx.tables_q2.inverse(&mut cs_ntt_q2);

    // v = c0 - c1*s in each limb
    let mut v_t = [0u32; N];
    let mut v_q2 = [0u32; N];
    for i in 0..N {
        v_t[i] = submod(ct.c0_t[i], cs_ntt_t[i], T);
        v_q2[i] = submod(ct.c0_q2[i], cs_ntt_q2[i], Q2);
    }

    // Noise-limb trick: v_q2 is pure noise. The center-lift and final reduction
    // run on secret-dependent values, so they are done branchlessly (constant-time):
    // no data-dependent control flow and no data-indexed memory access.
    let half_q2 = (Q2 / 2) as u64;
    let q2_i = Q2 as i64;
    let t_u = T as u64;
    let t_i = T as i64;
    let mut message = [0u32; N];
    for i in 0..N {
        let v = v_q2[i] as u64;
        // gt_mask = all-ones iff v > q₂/2 (i.e. the center-lift maps v ↦ v − q₂).
        let gt_mask = (half_q2.wrapping_sub(v) >> 63).wrapping_neg();
        let noise = v as i64 - (q2_i & gt_mask as i64); // ∈ [−q₂/2, q₂/2]

        // M = (v_t − noise) · Δ_t⁻¹ mod t. The value (v_t − noise) ∈ (−T, 3T);
        // shift by +T into (0, 3T) and reduce with two branchless subtractions.
        let u = ((v_t[i] as i64 - noise) + t_i) as u64;
        let u = csub(csub(u, t_u), t_u);
        message[i] = barrett_reduce(u * INV_DELTA_T as u64, T, mu_t);
    }

    // Best-effort scrubbing of secret-derived intermediates (c1·s and the
    // decryption noise).
    cs_ntt_t.zeroize();
    cs_ntt_q2.zeroize();
    v_t.zeroize();
    v_q2.zeroize();

    message
}
