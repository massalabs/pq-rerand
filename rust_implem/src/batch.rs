//! Multi-threaded batch operations.
//!
//! Batch encryption, re-randomization and decryption are embarrassingly parallel
//! across ciphertext slots. These helpers fan the per-slot work out across the
//! current rayon thread pool (the global one, or the pool the call is made in
//! through `ThreadPool::install`). Every slot draws its randomness from its own
//! ChaCha20 generator, keyed with 256 fresh bits of the caller's generator and
//! zeroized after use.

use crate::poly::NttContext;
use crate::sampling::wiped;
use crate::encoding::Plaintext;
use crate::keygen::{PublicKey, SecretKey};
use crate::encrypt::{encrypt_slot, CiphertextSlot, Flooding};
use crate::rerandomize::rerandomize_slot;
use crate::decrypt::decrypt_slot;
use chacha20::ChaCha20Rng;
use rand_core::{CryptoRng, Rng, SeedableRng};
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

/// Encrypt a batch of plaintext slots in parallel.
///
/// `rng` must be cryptographically secure, like for [`encrypt_slot`]; `flood`
/// selects the flood packet.
pub fn encrypt_batch<R: Rng + CryptoRng>(
    rng: &mut R,
    ctx: &NttContext,
    pk: &PublicKey,
    messages: &[Plaintext],
    flood: Flooding,
) -> Vec<CiphertextSlot> {
    let mut rngs = slot_rngs(rng, messages.len());
    messages
        .par_iter()
        .zip(rngs.par_iter_mut())
        .map(|(m, slot_rng)| encrypt_slot(slot_rng, ctx, pk, m, flood))
        .collect()
}

/// Re-randomize a batch of ciphertext slots in parallel.
///
/// `rng` must be cryptographically secure. The re-randomization budget documented
/// on [`rerandomize_slot`] applies to every slot, and tracking it is the caller's
/// responsibility.
pub fn rerandomize_batch<R: Rng + CryptoRng>(
    rng: &mut R,
    ctx: &NttContext,
    pk: &PublicKey,
    cts: &[CiphertextSlot],
) -> Vec<CiphertextSlot> {
    let mut rngs = slot_rngs(rng, cts.len());
    cts.par_iter()
        .zip(rngs.par_iter_mut())
        .map(|(ct, slot_rng)| rerandomize_slot(slot_rng, ctx, pk, ct))
        .collect()
}

/// Decrypt a batch of ciphertext slots in parallel.
///
/// The caveats of [`decrypt_slot`] apply to every slot.
pub fn decrypt_batch(
    ctx: &NttContext,
    sk: &SecretKey,
    cts: &[CiphertextSlot],
) -> Vec<Plaintext> {
    cts.par_iter().map(|ct| decrypt_slot(ctx, sk, ct)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::{decode, encode};
    use crate::keygen::keygen;
    use crate::params::{NUM_SLOTS, SLOT_BYTES};
    use crate::poly::NTT_CONTEXT;
    use crate::serialize::serialize_slot;
    use rand::rngs::StdRng;

    fn payloads(count: usize) -> Vec<[u8; SLOT_BYTES]> {
        (0..count).map(|s| std::array::from_fn(|i| (i * 7 + s * 13) as u8)).collect()
    }

    #[test]
    fn batch_roundtrip_parallel() {
        let ctx = &NTT_CONTEXT;
        let mut rng = StdRng::seed_from_u64(1);
        let (sk, pk) = keygen(&mut rng, ctx);

        let data = payloads(32);
        let messages: Vec<Plaintext> = data.iter().map(encode).collect();
        for flood in [Flooding::Standard, Flooding::Off] {
            let cts = encrypt_batch(&mut rng, ctx, &pk, &messages, flood);
            let cts2 = rerandomize_batch(&mut rng, ctx, &pk, &cts);
            let recovered = decrypt_batch(ctx, &sk, &cts2);
            assert_eq!(recovered.len(), data.len());
            for (m, d) in recovered.iter().zip(&data) {
                assert_eq!(decode(m).expect("valid slot")[..], d[..], "parallel batch roundtrip failed");
            }
        }
        assert!(encrypt_batch(&mut rng, ctx, &pk, &[], Flooding::Standard).is_empty());
    }

    /// Slots do not share randomness, and neither do two calls made with one
    /// generator: no two c1 components are equal, which a reused generator
    /// state would make them.
    #[test]
    fn batch_slots_and_calls_use_independent_randomness() {
        let ctx = &NTT_CONTEXT;
        let mut rng = StdRng::seed_from_u64(2);
        let (_sk, pk) = keygen(&mut rng, ctx);
        let messages: Vec<Plaintext> = payloads(4).iter().map(encode).collect();

        let mut cts = encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard);
        cts.extend(encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard));
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
        let first = rerandomize_batch(&mut rng, ctx, &pk, &cts);
        let second = rerandomize_batch(&mut rng, ctx, &pk, &cts);
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
        let ctx = &NTT_CONTEXT;
        let mut rng = StdRng::seed_from_u64(5);
        let (sk, pk) = keygen(&mut rng, ctx);
        let data = payloads(NUM_SLOTS);
        let messages: Vec<Plaintext> = data.iter().map(encode).collect();

        let pool = rayon::ThreadPoolBuilder::new().num_threads(16).stack_size(256 * 1024).build().unwrap();
        let recovered = pool.install(|| {
            let cts = encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard);
            let cts = rerandomize_batch(&mut rng, ctx, &pk, &cts);
            decrypt_batch(ctx, &sk, &cts)
        });
        assert_eq!(recovered.len(), NUM_SLOTS);
        for (m, d) in recovered.iter().zip(&data) {
            assert_eq!(decode(m).expect("valid slot")[..], d[..]);
        }
    }

    /// Slot `i` uses the generator keyed with the `i`-th 32 bytes of the caller's
    /// generator, whatever the thread that runs it, and `flood` reaches every slot.
    #[test]
    fn batch_matches_sequential_slots() {
        let ctx = &NTT_CONTEXT;
        let mut rng = StdRng::seed_from_u64(3);
        let (_sk, pk) = keygen(&mut rng, ctx);
        let messages: Vec<Plaintext> = payloads(5).iter().map(encode).collect();
        // The generators that a batch call gives to its slots when seeded with 4.
        let slot_rngs = || {
            let mut seeder = StdRng::seed_from_u64(4);
            (0..messages.len()).map(move |_| {
                let mut seed = [0u8; 32];
                seeder.fill_bytes(&mut seed);
                ChaCha20Rng::from_seed(seed)
            })
        };
        let bytes = |cts: &[CiphertextSlot]| cts.iter().map(serialize_slot).collect::<Vec<_>>();

        let mut cts = Vec::new();
        for flood in [Flooding::Off, Flooding::Standard] {
            let slots = messages.iter().zip(slot_rngs());
            let expected: Vec<_> = slots.map(|(m, mut rng)| encrypt_slot(&mut rng, ctx, &pk, m, flood)).collect();
            cts = encrypt_batch(&mut StdRng::seed_from_u64(4), ctx, &pk, &messages, flood);
            assert!(bytes(&cts) == bytes(&expected), "encrypt_batch, {flood:?}");
        }
        let slots = cts.iter().zip(slot_rngs());
        let expected: Vec<_> = slots.map(|(ct, mut rng)| rerandomize_slot(&mut rng, ctx, &pk, ct)).collect();
        let rerandomized = rerandomize_batch(&mut StdRng::seed_from_u64(4), ctx, &pk, &cts);
        assert!(bytes(&rerandomized) == bytes(&expected), "rerandomize_batch");
    }
}
