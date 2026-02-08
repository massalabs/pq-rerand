//! Integration tests: encrypt → decrypt = identity, rerand preserves plaintext.

use pq_rerand::params::*;
use pq_rerand::poly::NttContext;
use pq_rerand::keygen::keygen;
use pq_rerand::encrypt::encrypt_slot;
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::rerandomize::rerandomize_slot;
use pq_rerand::encoding::{encode, decode};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn setup() -> (pq_rerand::keygen::SecretKey, pq_rerand::keygen::PublicKey, NttContext, StdRng) {
    let ctx = NttContext::new();
    let mut rng = StdRng::seed_from_u64(42);
    let (sk, pk) = keygen(&mut rng, &ctx);
    (sk, pk, ctx, rng)
}

#[test]
fn test_encrypt_decrypt_zeros() {
    let (sk, pk, ctx, mut rng) = setup();
    let message = [0u32; N];
    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    let recovered = decrypt_slot(&ctx, &sk, &ct);
    assert_eq!(message, recovered, "decrypt of zeros failed");
}

#[test]
fn test_encrypt_decrypt_ones() {
    let (sk, pk, ctx, mut rng) = setup();
    let message = [1u32; N];
    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    let recovered = decrypt_slot(&ctx, &sk, &ct);
    assert_eq!(message, recovered, "decrypt of ones failed");
}

#[test]
fn test_encrypt_decrypt_random() {
    let (sk, pk, ctx, mut rng) = setup();
    // Random message with 31-bit coefficients
    let mut message = [0u32; N];
    for i in 0..N {
        message[i] = rand::Rng::gen_range(&mut rng, 0..(1u32 << BITS_PER_COEFF));
    }
    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    let recovered = decrypt_slot(&ctx, &sk, &ct);
    assert_eq!(message, recovered, "decrypt of random message failed");
}

#[test]
fn test_encrypt_decrypt_max_values() {
    let (sk, pk, ctx, mut rng) = setup();
    // Max 31-bit value
    let message = [(1u32 << BITS_PER_COEFF) - 1; N];
    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    let recovered = decrypt_slot(&ctx, &sk, &ct);
    assert_eq!(message, recovered, "decrypt of max values failed");
}

#[test]
fn test_rerand_preserves_plaintext() {
    let (sk, pk, ctx, mut rng) = setup();
    let mut message = [0u32; N];
    for i in 0..N {
        message[i] = rand::Rng::gen_range(&mut rng, 0..(1u32 << BITS_PER_COEFF));
    }
    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    let ct2 = rerandomize_slot(&mut rng, &ctx, &pk, &ct);
    let recovered = decrypt_slot(&ctx, &sk, &ct2);
    assert_eq!(message, recovered, "rerand broke plaintext");
}

#[test]
fn test_multiple_rerands() {
    let (sk, pk, ctx, mut rng) = setup();
    let mut message = [0u32; N];
    for i in 0..N {
        message[i] = rand::Rng::gen_range(&mut rng, 0..(1u32 << BITS_PER_COEFF));
    }
    let mut ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    // Apply 100 re-randomizations
    for _ in 0..100 {
        ct = rerandomize_slot(&mut rng, &ctx, &pk, &ct);
    }
    let recovered = decrypt_slot(&ctx, &sk, &ct);
    assert_eq!(message, recovered, "100 rerands broke plaintext");
}

#[test]
fn test_encoding_roundtrip_via_crypto() {
    let (sk, pk, ctx, mut rng) = setup();
    // Create a realistic plaintext slot from bytes
    let mut data = vec![0u8; SLOT_BYTES];
    for i in 0..data.len() {
        data[i] = rand::Rng::gen(&mut rng);
    }
    let message = encode(&data);
    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);
    let ct2 = rerandomize_slot(&mut rng, &ctx, &pk, &ct);
    let recovered = decrypt_slot(&ctx, &sk, &ct2);
    let decoded = decode(&recovered);
    assert_eq!(data, decoded, "full encode→encrypt→rerand→decrypt→decode roundtrip failed");
}
