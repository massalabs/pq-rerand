//! Encryption: base BFV encrypt + noise-flooding smudge.

use crate::params::{N, T, Q2, DELTA_T, SLOT_BYTES};
use crate::ntt::{pointwise_mul, mul_shoup, shoup, addmod};
use crate::poly::{NTT_CONTEXT, add_signed, secret, set_signed, zeroed};
use crate::encoding::encode;
use crate::keygen::PublicKey;
use crate::sampling::{add_flood, sample_gaussian_cdt, wiped};
use rand_core::{CryptoRng, Rng};

/// A single ciphertext slot: (c0_t, c0_q2, c1_t, c1_q2), all in coefficient domain.
///
/// The limbs live on the heap: returning or moving a slot moves four pointers, not
/// 64 KiB through every stack frame on the way.
#[derive(Clone, Debug)]
pub struct CiphertextSlot {
    /// `c₀ mod t`, coefficient domain.
    pub(crate) c0_t: Box<[u32; N]>,
    /// `c₀ mod q₂`, coefficient domain.
    pub(crate) c0_q2: Box<[u32; N]>,
    /// `c₁ mod t`, coefficient domain.
    pub(crate) c1_t: Box<[u32; N]>,
    /// `c₁ mod q₂`, coefficient domain.
    pub(crate) c1_q2: Box<[u32; N]>,
}

impl CiphertextSlot {
    /// The all-zero slot, to be written in place.
    pub(crate) fn zeroed() -> Self {
        CiphertextSlot { c0_t: zeroed(), c0_q2: zeroed(), c1_t: zeroed(), c1_q2: zeroed() }
    }
}

/// Noise flooding applied by encryption (paper section 6).
///
/// The flood packet makes a fresh ciphertext look, to the secret-key holder, like
/// one that was already re-randomized `KAPPA_F` times. It only hides the number of
/// re-randomizations while that number is small next to `KAPPA_F` and the holder
/// sees a single version of the ciphertext: two versions of one ciphertext differ
/// by exactly the noise added in between. Neither confidentiality nor
/// unlinkability against parties without the secret key depends on it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flooding {
    /// No flood packet (`kappa_f = 0`).
    Off,
    /// One flood packet of width [`SIGMA_FLOOD`](crate::params::SIGMA_FLOOD)
    /// (`kappa_f` = [`KAPPA_F`](crate::params::KAPPA_F)): the analyzed configuration.
    Standard,
}

/// Overwrite `ct` with a fresh encryption of zero under `pk`: `(b*R + E2, a*R + E1)`
/// with base-width noise, widened by a flood packet if `flood` says so.
///
/// This is the whole of re-randomization and, up to the message term, of
/// encryption. The randomness and the products only ever exist in zeroizing heap
/// storage and in `ct` itself, which the caller completes in place. The callers
/// run it through [`wiped`]: a compiler may keep a few coefficients of the noise
/// in the stack frame of this function, which is then overwritten as well.
pub(crate) fn encrypt_zero<R: Rng + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    flood: Flooding,
    ct: &mut CiphertextSlot,
) {
    let ctx = &NTT_CONTEXT;
    // Sample base-width noise (constant-time CDT) and flood-width noise, then add
    // them coefficient-wise. The base encryption and the flooding Enc(0) are both
    // linear in their randomness/noise, so a separate smudge step
    //     C = (b·r+e₂+Δ·M, a·r+e₁) ,  C += (b·r_f+e_{2,f}, a·r_f+e_{1,f})
    // produces exactly the same distribution as a single encryption with combined
    // noise R = r+r_f, E₁ = e₁+e_{1,f}, E₂ = e₂+e_{2,f}. Folding the two together
    // halves the NTT work (one forward of R, one inverse per ciphertext component)
    // while leaving the output distribution — and hence the security analysis —
    // unchanged. Flooding therefore costs only its extra sampling, not extra NTTs.
    let (mut r, mut e1, mut e2) = (secret::<i32>(), secret::<i32>(), secret::<i32>());
    sample_gaussian_cdt(rng, &mut r);
    sample_gaussian_cdt(rng, &mut e1);
    sample_gaussian_cdt(rng, &mut e2);
    if flood == Flooding::Standard {
        add_flood(rng, &mut r);
        add_flood(rng, &mut e1);
        add_flood(rng, &mut e2);
    }

    // NTT of the combined randomness R = r + r_f.
    let (mut r_ntt_t, mut r_ntt_q2) = (secret(), secret());
    set_signed(&mut r_ntt_t, &mut r_ntt_q2, &r);
    ctx.tables_t.forward(&mut r_ntt_t);
    ctx.tables_q2.forward(&mut r_ntt_q2);

    // c1 = a·R + E1 (in each limb)
    ct.c1_t.copy_from_slice(&*pk.a_ntt_t);
    ct.c1_q2.copy_from_slice(&*pk.a_ntt_q2);
    pointwise_mul(&mut ct.c1_t, &r_ntt_t, T);
    pointwise_mul(&mut ct.c1_q2, &r_ntt_q2, Q2);
    ctx.tables_t.inverse(&mut ct.c1_t);
    ctx.tables_q2.inverse(&mut ct.c1_q2);
    add_signed(&mut ct.c1_t, &mut ct.c1_q2, &e1);

    // c0 = b*R + E2 (in each limb)
    ct.c0_t.copy_from_slice(&*pk.b_ntt_t);
    ct.c0_q2.copy_from_slice(&*pk.b_ntt_q2);
    pointwise_mul(&mut ct.c0_t, &r_ntt_t, T);
    pointwise_mul(&mut ct.c0_q2, &r_ntt_q2, Q2);
    ctx.tables_t.inverse(&mut ct.c0_t);
    ctx.tables_q2.inverse(&mut ct.c0_q2);
    add_signed(&mut ct.c0_t, &mut ct.c0_q2, &e2);
}

/// Encrypt a single slot: the `SLOT_BYTES` bytes of `payload`, which the
/// application has authenticated (see the crate documentation). Performs base
/// encryption + noise-flooding smudge, as selected by `flood`.
///
/// The RNG must be cryptographically secure (enforced via the [`CryptoRng`]
/// marker bound). The encryption randomness reveals the payload: like the message
/// polynomial that encodes the payload, it only ever exists in heap storage that
/// is zeroized before returning.
///
/// `pk` must be authentic (see [`PublicKey`]).
pub fn encrypt_slot<R: Rng + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    payload: &[u8; SLOT_BYTES],
    flood: Flooding,
) -> CiphertextSlot {
    let mut ct = CiphertextSlot::zeroed();
    wiped(|| {
        encrypt_zero(rng, pk, flood, &mut ct);
        let message = encode(payload);

        // Add Delta_t * M to c0_t (the message lives in the t limb only), in the
        // wiped frame as well: a compiler may keep coefficients of the message on
        // the stack. Shoup multiplication has no final conditional subtraction
        // that a compiler could turn into a branch on the message.
        const DELTA_SHOUP: u64 = shoup(DELTA_T, T);
        for i in 0..N {
            let delta_m = mul_shoup(message.0[i], DELTA_T, DELTA_SHOUP, T);
            ct.c0_t[i] = addmod(ct.c0_t[i], delta_m, T);
        }
    });
    ct
}
