//! Public re-randomization: C' = C + Enc(0) with narrow noise.

use crate::params::{T, Q2};
use crate::ntt::pointwise_add;
use crate::keygen::PublicKey;
use crate::encrypt::{CiphertextSlot, Flooding, encrypt_zero};
use crate::sampling::wiped;
use rand_core::{CryptoRng, Rng};

/// Re-randomize a single ciphertext slot by adding a fresh Enc(0).
///
/// The RNG must be cryptographically secure (enforced via the [`CryptoRng`]
/// marker bound). The added encryption of zero links the output to the input: it
/// only ever exists in zeroizing heap storage and in the returned ciphertext.
///
/// `pk` must be authentic (see [`PublicKey`]): re-randomizing under a key that
/// someone else chose lets them link the output to the input.
///
/// # Re-randomization budget
///
/// Every call adds noise. Decryption is only guaranteed, with the analyzed failure
/// probability (2^-106 per batch of 1024 slots, for all but a negligible fraction
/// of keys), while a ciphertext has been re-randomized at most
/// [`K_MAX_RERAND`](crate::params::K_MAX_RERAND) times since it was encrypted;
/// beyond that, decryption becomes increasingly likely to fail. The count is not
/// part of the ciphertext (it would break unlinkability), so this function can
/// neither know nor enforce it: the application must keep an upper bound on it
/// (for example age times maximum rate) and refresh the ciphertext, by decrypting
/// and encrypting it again, before the limit.
pub fn rerandomize_slot<R: Rng + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    ct: &CiphertextSlot,
) -> CiphertextSlot {
    // Enc(0) = (b*r + e2, a*r + e1) in each limb, written straight to the output.
    let mut out = CiphertextSlot::zeroed();
    wiped(|| encrypt_zero(rng, pk, Flooding::Off, &mut out));

    // c0' = c0 + b*r + e2
    pointwise_add(&mut out.c0_t, &ct.c0_t, T);
    pointwise_add(&mut out.c0_q2, &ct.c0_q2, Q2);

    // c1' = c1 + a*r + e1
    pointwise_add(&mut out.c1_t, &ct.c1_t, T);
    pointwise_add(&mut out.c1_q2, &ct.c1_q2, Q2);

    out
}
