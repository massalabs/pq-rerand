//! Multi-threaded batch operations.
//!
//! Batch encryption, re-randomization and decryption are embarrassingly parallel
//! across ciphertext slots. These helpers fan the per-slot work out across a rayon
//! thread pool. Per-slot RNGs (ChaCha-based `StdRng`, a CSPRNG) are derived
//! deterministically from a 256-bit base seed with per-slot domain separation, so
//! parallel runs are reproducible given the seed. Callers must supply a fresh,
//! uniformly random 256-bit seed per batch operation (e.g. from [`rand::rngs::OsRng`]).

use crate::params::N;
use crate::poly::NttContext;
use crate::keygen::{PublicKey, SecretKey};
use crate::encrypt::{encrypt_slot, CiphertextSlot};
use crate::rerandomize::rerandomize_slot;
use crate::decrypt::decrypt_slot;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rayon::prelude::*;
use std::sync::OnceLock;

/// Per-slot work uses large stack frames (several length-N arrays), so ensure the
/// rayon workers have a generous stack. Installs a global pool on first use; if the
/// caller already configured one, this is a no-op.
fn ensure_pool() {
    static INIT: OnceLock<()> = OnceLock::new();
    INIT.get_or_init(|| {
        let _ = rayon::ThreadPoolBuilder::new()
            .stack_size(16 * 1024 * 1024)
            .build_global();
    });
}

/// Derive a per-slot 256-bit seed from the base seed and the slot index.
///
/// Splices a mixed slot counter into bytes 0..8 of the base seed. Distinct slots
/// therefore use distinct ChaCha keys; the mixing constant (SplitMix64 increment)
/// only spreads the index bits and carries no security weight.
#[inline]
fn slot_seed(base_seed: &[u8; 32], i: usize) -> [u8; 32] {
    let mut seed = *base_seed;
    let mixed = (i as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    for (b, m) in seed[..8].iter_mut().zip(mixed.to_le_bytes()) {
        *b ^= m;
    }
    seed
}

/// Encrypt a batch of message polynomials in parallel.
///
/// `base_seed` must be a fresh, uniformly random 256-bit value.
pub fn encrypt_batch(
    ctx: &NttContext,
    pk: &PublicKey,
    messages: &[[u32; N]],
    sigma_flood: f64,
    base_seed: &[u8; 32],
) -> Vec<CiphertextSlot> {
    ensure_pool();
    messages
        .par_iter()
        .enumerate()
        .map(|(i, m)| {
            let mut rng = StdRng::from_seed(slot_seed(base_seed, i));
            encrypt_slot(&mut rng, ctx, pk, m, sigma_flood)
        })
        .collect()
}

/// Re-randomize a batch of ciphertext slots in parallel.
///
/// `base_seed` must be a fresh, uniformly random 256-bit value.
pub fn rerandomize_batch(
    ctx: &NttContext,
    pk: &PublicKey,
    cts: &[CiphertextSlot],
    base_seed: &[u8; 32],
) -> Vec<CiphertextSlot> {
    ensure_pool();
    cts.par_iter()
        .enumerate()
        .map(|(i, ct)| {
            let mut rng = StdRng::from_seed(slot_seed(base_seed, i));
            rerandomize_slot(&mut rng, ctx, pk, ct)
        })
        .collect()
}

/// Decrypt a batch of ciphertext slots in parallel.
pub fn decrypt_batch(
    ctx: &NttContext,
    sk: &SecretKey,
    cts: &[CiphertextSlot],
) -> Vec<[u32; N]> {
    ensure_pool();
    cts.par_iter().map(|ct| decrypt_slot(ctx, sk, ct)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SIGMA_FLOOD;
    use crate::keygen::keygen;

    #[test]
    fn batch_roundtrip_parallel() {
        let ctx = NttContext::new();
        let mut rng = StdRng::seed_from_u64(1);
        let (sk, pk) = keygen(&mut rng, &ctx);

        let messages: Vec<[u32; N]> = (0..32)
            .map(|s| {
                let mut m = [0u32; N];
                for (i, c) in m.iter_mut().enumerate() {
                    *c = ((i + s) as u32) % (1 << 31);
                }
                m
            })
            .collect();

        let cts = encrypt_batch(&ctx, &pk, &messages, SIGMA_FLOOD, &[42u8; 32]);
        let cts2 = rerandomize_batch(&ctx, &pk, &cts, &[43u8; 32]);
        let recovered = decrypt_batch(&ctx, &sk, &cts2);
        assert_eq!(messages, recovered, "parallel batch roundtrip failed");
    }

    #[test]
    fn slot_seeds_are_distinct() {
        let base = [7u8; 32];
        let s0 = slot_seed(&base, 0);
        let s1 = slot_seed(&base, 1);
        assert_ne!(s0, s1);
        // Only the first 8 bytes carry the slot index.
        assert_eq!(s0[8..], base[8..]);
    }
}
