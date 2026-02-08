//! Decryption using the noise-limb trick.
//!
//! Because Δ·M ≡ 0 (mod q₂), the q₂-limb of v = c0 - c1·s is pure noise.
//! Center-lift it to recover the signed noise, then subtract from the t-limb
//! and multiply by INV_DELTA_T to recover M.

use crate::params::{N, T, Q2, INV_DELTA_T};
use crate::ntt::{mulmod, submod};
use crate::poly::NttContext;
use crate::keygen::SecretKey;
use crate::encrypt::CiphertextSlot;

/// Decrypt a single ciphertext slot, returning the message polynomial
/// coefficients in [0, t).
pub fn decrypt_slot(
    ctx: &NttContext,
    sk: &SecretKey,
    ct: &CiphertextSlot,
) -> [u32; N] {
    // Compute c1 * s in each limb via NTT
    // c1 is in coefficient domain; we need to NTT it, pointwise-mul with NTT(s), INTT.
    let mut c1_ntt_t = ct.c1_t;
    let mut c1_ntt_q2 = ct.c1_q2;
    ctx.tables_t.forward(&mut c1_ntt_t);
    ctx.tables_q2.forward(&mut c1_ntt_q2);

    let mut s_ntt_t = sk.s_t;
    let mut s_ntt_q2 = sk.s_q2;
    ctx.tables_t.forward(&mut s_ntt_t);
    ctx.tables_q2.forward(&mut s_ntt_q2);

    // c1*s in NTT domain
    let mut cs_ntt_t = [0u32; N];
    let mut cs_ntt_q2 = [0u32; N];
    for i in 0..N {
        cs_ntt_t[i] = mulmod(c1_ntt_t[i], s_ntt_t[i], T);
        cs_ntt_q2[i] = mulmod(c1_ntt_q2[i], s_ntt_q2[i], Q2);
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

    // Noise-limb trick: v_q2 is pure noise.
    // Center-lift to signed, then recover M from v_t.
    let half_q2 = Q2 / 2;
    let mut message = [0u32; N];
    for i in 0..N {
        // Center-lift noise from q₂-limb
        let noise_i: i64 = if v_q2[i] <= half_q2 {
            v_q2[i] as i64
        } else {
            v_q2[i] as i64 - Q2 as i64
        };

        // M = (v_t - noise) * INV_DELTA_T mod t
        let tmp = (v_t[i] as i64 - noise_i).rem_euclid(T as i64);
        message[i] = mulmod(tmp as u32, INV_DELTA_T, T);
    }

    message
}
