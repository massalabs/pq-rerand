//! Multi-threaded batch operations.

use crate::params::SLOT_BYTES;
use crate::sampling::wiped;
use crate::keygen::{PublicKey, SecretKey};
use crate::encrypt::{encrypt_slot, CiphertextSlot, Flooding};
use crate::rerandomize::rerandomize_slot;
use crate::decrypt::{Payloads, decrypt_into};
use chacha20::ChaCha20Rng;
use rand_core::{CryptoRng, Rng, SeedableRng};
use core::hint::black_box;
use rayon::prelude::*;
use zeroize::Zeroizing;

/// One independent generator per slot, each keyed with 32 bytes drawn from `rng`.
///
/// The keys are drawn before the parallel section, in slot order. Keys and
/// generators live on the heap and are zeroized on drop, and the stack that their
/// construction used is wiped.
fn slot_rngs<R: Rng + CryptoRng>(rng: &mut R, slots: usize) -> Vec<ChaCha20Rng> {
    wiped(|| {
        let mut seeds = Zeroizing::new(vec![[0u8; 32]; slots]);
        seeds.iter_mut().for_each(|seed| rng.fill_bytes(seed));
        seeds.iter().map(|seed| ChaCha20Rng::from_seed(*seed)).collect()
    })
}

/// Encrypt a batch of payload slots in parallel.
///
/// `rng` must be cryptographically secure, like for [`encrypt_slot`]; `flood`
/// selects the flood packet.
///
/// Like the other batch functions, this one fans the per-slot work out across
/// the current rayon thread pool (the global one, or the pool the call is made
/// in through `ThreadPool::install`). Every slot draws its randomness from its
/// own ChaCha20 generator, keyed with 256 fresh bits of `rng` and zeroized after
/// use, here and in [`rerandomize_batch`].
pub fn encrypt_batch<R: Rng + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    payloads: &[[u8; SLOT_BYTES]],
    flood: Flooding,
) -> Vec<CiphertextSlot> {
    let mut rngs = slot_rngs(rng, payloads.len());
    payloads
        .par_iter()
        .zip(rngs.par_iter_mut())
        .map(|(payload, slot_rng)| encrypt_slot(slot_rng, pk, payload, flood))
        .collect()
}

/// Re-randomize a batch of ciphertext slots in parallel.
///
/// `rng` must be cryptographically secure. The re-randomization budget documented
/// on [`rerandomize_slot`] applies to every slot, and tracking it is the caller's
/// responsibility.
pub fn rerandomize_batch<R: Rng + CryptoRng>(
    rng: &mut R,
    pk: &PublicKey,
    cts: &[CiphertextSlot],
) -> Vec<CiphertextSlot> {
    let mut rngs = slot_rngs(rng, cts.len());
    cts.par_iter()
        .zip(rngs.par_iter_mut())
        .map(|(ct, slot_rng)| rerandomize_slot(slot_rng, pk, ct))
        .collect()
}

/// Decrypt a batch of ciphertext slots in parallel and pass their payloads, in
/// order, to `open`, which authenticates them. Returns what `open` returns, or
/// `None` if a slot did not decrypt to a payload.
///
/// This is [`decrypt_slot`](crate::decrypt_slot) for a payload that
/// spans several slots, with the same rules for `open` and the same single
/// failure: every slot is decrypted whatever the others give, a slot that did
/// not decrypt to a payload reaches `open` as zeros, and the payloads only exist
/// in heap storage that is zeroized before returning. The batch has one
/// verdict: slots that are authenticated one by one go through `decrypt_slot`.
pub fn decrypt_batch<T>(
    sk: &SecretKey,
    cts: &[CiphertextSlot],
    open: impl FnOnce(&mut [[u8; SLOT_BYTES]]) -> Option<T>,
) -> Option<T> {
    let mut payloads = Payloads::new(cts.len());
    let slots = payloads.slots().par_iter_mut().zip(cts.par_iter());
    // The slots that did not decrypt to a payload are counted rather than
    // combined as truth values. Between threads a result travels in an `Option`,
    // whose tag shares the byte of a `bool`: testing the tag would then be a jump
    // on a byte that holds the verdict.
    let invalid: usize = slots.map(|(payload, ct)| usize::from(!decrypt_into(sk, ct, payload))).sum();
    // Read once `open` has accepted, as in `decrypt_slot`.
    open(payloads.slots()).filter(|_| black_box(invalid) == 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::keygen;
    use crate::params::{DELTA_T, NUM_SLOTS, T};
    use rand::rngs::StdRng;

    fn payloads(count: usize) -> Vec<[u8; SLOT_BYTES]> {
        (0..count).map(|s| std::array::from_fn(|i| (i * 7 + s * 13) as u8)).collect()
    }

    #[test]
    fn batch_roundtrip_parallel() {
        let mut rng = StdRng::seed_from_u64(1);
        let (sk, pk) = keygen(&mut rng);

        let data = payloads(32);
        for flood in [Flooding::Standard, Flooding::Off] {
            let cts = encrypt_batch(&mut rng, &pk, &data, flood);
            let cts2 = rerandomize_batch(&mut rng, &pk, &cts);
            let recovered = decrypt_batch(&sk, &cts2, |slots| Some(slots.to_vec()));
            assert!(recovered.as_ref() == Some(&data), "parallel batch roundtrip failed");
        }
        assert!(encrypt_batch(&mut rng, &pk, &[], Flooding::Standard).is_empty());
        // An empty batch has an empty payload, which `open` judges like any other.
        assert_eq!(decrypt_batch(&sk, &[], |slots| Some(slots.len())), Some(0));
    }

    /// One failure for the whole batch: `open` rejecting the payloads, or a slot
    /// that did not decrypt to one, which `open` then sees as zeros next to the
    /// payloads of the other slots, and which fails the batch whatever `open`
    /// answers.
    #[test]
    fn batch_decryption_has_one_failure() {
        let mut rng = StdRng::seed_from_u64(6);
        let (sk, pk) = keygen(&mut rng);
        let data = payloads(4);
        let cts = encrypt_batch(&mut rng, &pk, &data, Flooding::Standard);
        // What `open` is given, and the outcome when it accepts.
        let accepted = |cts: &[CiphertextSlot]| {
            let mut given = Vec::new();
            let outcome = decrypt_batch(&sk, cts, |slots| {
                given = slots.to_vec();
                Some(())
            });
            (outcome, given)
        };
        assert!(accepted(&cts) == (Some(()), data.clone()));
        assert_eq!(decrypt_batch(&sk, &cts, |_| None::<()>), None);

        // Every slot in turn: the first coefficient of its message, m, is moved to
        // t - 1, which is not below 2^31, by adding Delta_t * (t - 1 - m) to c0 in
        // the t limb.
        for bad in 0..cts.len() {
            let (mut cts, mut data) = (cts.clone(), data.clone());
            let m = u32::from_le_bytes(data[bad][..4].try_into().unwrap()) & 0x7fff_ffff;
            let shift = (DELTA_T as u64 * (T - 1 - m) as u64 % T as u64) as u32;
            cts[bad].c0_t[0] = crate::ntt::addmod(cts[bad].c0_t[0], shift, T);
            data[bad] = [0; SLOT_BYTES];
            assert!(accepted(&cts) == (None, data), "slot {bad}");
        }
    }

    /// Slots do not share randomness, and neither do two calls made with one
    /// generator: no two c1 components are equal, which a reused generator
    /// state would make them.
    #[test]
    fn batch_slots_and_calls_use_independent_randomness() {
        let mut rng = StdRng::seed_from_u64(2);
        let (_sk, pk) = keygen(&mut rng);
        let data = payloads(4);

        let mut cts = encrypt_batch(&mut rng, &pk, &data, Flooding::Standard);
        cts.extend(encrypt_batch(&mut rng, &pk, &data, Flooding::Standard));
        for i in 0..cts.len() {
            for j in 0..i {
                assert_ne!(cts[i].c1_t, cts[j].c1_t, "slots {i} and {j} share randomness");
                assert_ne!(cts[i].c1_q2, cts[j].c1_q2, "slots {i} and {j} share randomness");
            }
        }

        // Re-randomization: the encryptions of zero added to the slots all differ,
        // within a call and across two calls on the same input.
        let delta = |out: &CiphertextSlot, input: &CiphertextSlot| {
            let mut d = *out.c1_q2;
            crate::ntt::pointwise_sub(&mut d, &input.c1_q2, crate::params::Q2);
            d
        };
        let first = rerandomize_batch(&mut rng, &pk, &cts);
        let second = rerandomize_batch(&mut rng, &pk, &cts);
        let mut deltas: Vec<_> = first.iter().zip(&cts).map(|(o, i)| delta(o, i)).collect();
        deltas.extend(second.iter().zip(&cts).map(|(o, i)| delta(o, i)));
        for i in 0..deltas.len() {
            for j in 0..i {
                assert_ne!(deltas[i], deltas[j], "re-randomizations {i} and {j} share randomness");
            }
        }
    }

    /// The nominal batch runs on small worker stacks: nothing large is moved
    /// through the frames of the recursive work splitting.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "slow, and frames are larger, without optimizations")]
    fn nominal_batch_fits_small_worker_stacks() {
        let mut rng = StdRng::seed_from_u64(5);
        let (sk, pk) = keygen(&mut rng);
        let data = payloads(NUM_SLOTS);

        let pool = rayon::ThreadPoolBuilder::new().num_threads(16).stack_size(256 * 1024).build().unwrap();
        let recovered = pool.install(|| {
            let cts = encrypt_batch(&mut rng, &pk, &data, Flooding::Standard);
            let cts = rerandomize_batch(&mut rng, &pk, &cts);
            decrypt_batch(&sk, &cts, |slots| Some(slots.to_vec()))
        });
        assert!(recovered == Some(data));
    }

    /// Slot `i` uses the generator keyed with the `i`-th 32 bytes of the caller's
    /// generator, whatever the thread that runs it, and `flood` reaches every slot.
    #[test]
    fn batch_matches_sequential_slots() {
        let mut rng = StdRng::seed_from_u64(3);
        let (_sk, pk) = keygen(&mut rng);
        let data = payloads(5);
        // The generators that a batch call gives to its slots when seeded with 4.
        let slot_rngs = || {
            let mut seeder = StdRng::seed_from_u64(4);
            (0..data.len()).map(move |_| {
                let mut seed = [0u8; 32];
                seeder.fill_bytes(&mut seed);
                ChaCha20Rng::from_seed(seed)
            })
        };
        let bytes = |cts: &[CiphertextSlot]| cts.iter().map(CiphertextSlot::to_bytes).collect::<Vec<_>>();

        let mut cts = Vec::new();
        for flood in [Flooding::Off, Flooding::Standard] {
            let slots = data.iter().zip(slot_rngs());
            let expected: Vec<_> = slots.map(|(d, mut rng)| encrypt_slot(&mut rng, &pk, d, flood)).collect();
            cts = encrypt_batch(&mut StdRng::seed_from_u64(4), &pk, &data, flood);
            assert!(bytes(&cts) == bytes(&expected), "encrypt_batch, {flood:?}");
        }
        let slots = cts.iter().zip(slot_rngs());
        let expected: Vec<_> = slots.map(|(ct, mut rng)| rerandomize_slot(&mut rng, &pk, ct)).collect();
        let rerandomized = rerandomize_batch(&mut StdRng::seed_from_u64(4), &pk, &cts);
        assert!(bytes(&rerandomized) == bytes(&expected), "rerandomize_batch");
    }
}
