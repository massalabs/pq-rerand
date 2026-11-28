//! Integration tests: encrypt → decrypt = identity, rerand preserves plaintext.

use chacha20::ChaCha20Rng;
use pq_rerand::params::*;
use pq_rerand::{decrypt_slot, encrypt_slot, keygen, rerandomize_slot};
use pq_rerand::{CiphertextSlot, Flooding, PublicKey, SecretKey};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

fn setup() -> (SecretKey, PublicKey, StdRng) {
    let mut rng = StdRng::seed_from_u64(42);
    let (sk, pk) = keygen(&mut rng);
    (sk, pk, rng)
}

/// The payload that `ct` decrypts to, accepted as it is. Only a test may do
/// that: an application authenticates the payload here.
fn payload(sk: &SecretKey, ct: &CiphertextSlot) -> Option<[u8; SLOT_BYTES]> {
    decrypt_slot(sk, ct, |payload| Some(*payload))
}

fn random_slot(rng: &mut StdRng) -> [u8; SLOT_BYTES] {
    let mut data = [0u8; SLOT_BYTES];
    rng.fill_bytes(&mut data);
    data
}

/// Encrypt `data`, decrypt it and compare.
fn roundtrip(data: &[u8; SLOT_BYTES], flood: Flooding) {
    let (sk, pk, mut rng) = setup();
    let ct = encrypt_slot(&mut rng, &pk, data, flood);
    let recovered = payload(&sk, &ct).expect("decryption failed");
    assert_eq!(data[..], recovered[..]);
}

#[test]
fn test_encrypt_decrypt_zeros() {
    roundtrip(&[0u8; SLOT_BYTES], Flooding::Standard);
}

#[test]
fn test_encrypt_decrypt_ones() {
    roundtrip(&[1u8; SLOT_BYTES], Flooding::Standard);
}

#[test]
fn test_encrypt_decrypt_random() {
    let data = random_slot(&mut StdRng::seed_from_u64(7));
    roundtrip(&data, Flooding::Standard);
}

#[test]
fn test_encrypt_decrypt_max_values() {
    // Every coefficient at the maximum 31-bit value
    roundtrip(&[0xFFu8; SLOT_BYTES], Flooding::Standard);
}

#[test]
fn test_rerand_preserves_plaintext() {
    let (sk, pk, mut rng) = setup();
    let data = random_slot(&mut rng);
    let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
    let ct2 = rerandomize_slot(&mut rng, &pk, &ct);
    let recovered = payload(&sk, &ct2).expect("rerand broke plaintext");
    assert_eq!(data[..], recovered[..], "rerand broke plaintext");
}

#[test]
fn test_multiple_rerands() {
    let (sk, pk, mut rng) = setup();
    let data = random_slot(&mut rng);
    let mut ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
    // Apply 100 re-randomizations
    for _ in 0..100 {
        ct = rerandomize_slot(&mut rng, &pk, &ct);
    }
    let recovered = payload(&sk, &ct).expect("100 rerands broke plaintext");
    assert_eq!(data[..], recovered[..], "100 rerands broke plaintext");
}

#[test]
fn test_encrypt_decrypt_no_flood() {
    // Flooding is optional (paper section 6) and is disabled explicitly.
    let (sk, pk, mut rng) = setup();
    let data = random_slot(&mut rng);
    let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Off);
    let ct2 = rerandomize_slot(&mut rng, &pk, &ct);
    let recovered = payload(&sk, &ct2).expect("no-flood roundtrip failed");
    assert_eq!(data[..], recovered[..], "no-flood roundtrip failed");
}

/// The single-slot operations run on a thread that asks for a stack of 16 KiB,
/// which the C library may round up: keys, ciphertexts, randomness and the
/// payload of a decryption live on the heap.
#[test]
#[cfg_attr(debug_assertions, ignore = "frames are larger without optimizations")]
fn test_single_slot_operations_fit_a_small_stack() {
    let data = Box::new(random_slot(&mut StdRng::seed_from_u64(7)));
    let thread = std::thread::Builder::new().stack_size(16 * 1024).spawn(move || {
        let (sk, pk, mut rng) = setup();
        let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
        let ct = rerandomize_slot(&mut rng, &pk, &ct);
        decrypt_slot(&sk, &ct, |payload| (payload == &*data).then_some(())).is_some()
    });
    assert!(thread.unwrap().join().unwrap());
}

#[test]
fn test_serialization_roundtrip_via_crypto() {
    let (sk, pk, mut rng) = setup();
    let data = random_slot(&mut rng);
    let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
    let ct2 = CiphertextSlot::from_bytes(&ct.to_bytes()).expect("valid ciphertext rejected");
    let recovered = payload(&sk, &ct2).expect("valid slot");
    assert_eq!(data[..], recovered[..], "serialize -> deserialize -> decrypt failed");
}

#[test]
fn test_key_serialization_rejects_out_of_range() {
    let (sk, pk, _rng) = setup();

    // The first coefficient follows the version byte.
    let mut pk_bytes = pk.to_bytes();
    pk_bytes[1..5].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(PublicKey::from_bytes(&pk_bytes).is_none(), "out-of-range pk accepted");

    let sk_bytes = sk.to_bytes();
    let mut bad = sk_bytes.to_vec();
    bad[1..5].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(SecretKey::from_bytes(&bad).is_none(), "out-of-range sk accepted");
}

/// Keys that went through their serialized form work like the originals, and a
/// stored public key is checked against its secret key.
#[test]
fn test_roundtrip_with_many_keys() {
    for seed in 0..8 {
        let mut rng = StdRng::seed_from_u64(seed);
        let (sk, pk) = keygen(&mut rng);
        let sk = SecretKey::from_bytes(&sk.to_bytes()).expect("valid secret key rejected");
        let pk = PublicKey::from_bytes(&pk.to_bytes()).expect("valid public key rejected");
        assert!(pk.verify_against(&sk));

        let data = random_slot(&mut rng);
        let flood = if seed % 2 == 0 { Flooding::Standard } else { Flooding::Off };
        let mut ct = encrypt_slot(&mut rng, &pk, &data, flood);
        for _ in 0..3 {
            ct = rerandomize_slot(&mut rng, &pk, &ct);
        }
        let recovered = payload(&sk, &ct).expect("valid slot");
        assert_eq!(data[..], recovered[..], "seed {seed}");
    }
}

/// Decrypting under another key yields unrelated coefficients, a share of which
/// is no valid encoding: the caller gets no payload, or one that cannot authenticate.
#[test]
fn test_wrong_key_does_not_decrypt() {
    let (sk, pk, mut rng) = setup();
    let (other_sk, other_pk) = keygen(&mut rng);
    assert!(!pk.verify_against(&other_sk));
    assert!(!other_pk.verify_against(&sk));
    let data = random_slot(&mut rng);
    let mut rejected = 0;
    for _ in 0..64 {
        let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
        match payload(&other_sk, &ct) {
            None => rejected += 1,
            Some(garbage) => assert_ne!(data[..], garbage[..]),
        }
    }
    // A slot of uniform coefficients has one of 2^31 or more with probability 0.145.
    assert!((1..32).contains(&rejected), "{rejected} of 64 wrong-key decryptions rejected");
}

/// Random and perturbed inputs to the three parsers: an input of the wrong length,
/// of another version or with an unreduced coefficient is rejected, and whatever
/// parses as a ciphertext decrypts without a panic (in debug builds: without an
/// arithmetic overflow).
#[test]
fn test_parsers_on_random_inputs() {
    let (sk, pk, mut rng) = setup();
    let data = random_slot(&mut rng);
    let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
    // In all three formats, the limbs alternate between t and q2. They follow the
    // version byte of a key (a slot has none): `header` below.
    let reduced = |data: &[u8]| {
        let words = data.as_chunks::<4>().0.iter().map(|w| u32::from_le_bytes(*w));
        words.enumerate().all(|(i, w)| w < [T, Q2][i / N % 2])
    };
    type Parser = fn(&[u8]) -> bool;
    let parsers: [(Vec<u8>, Parser); 3] = [
        (pk.to_bytes(), |data| PublicKey::from_bytes(data).is_some()),
        (sk.to_bytes().to_vec(), |data| SecretKey::from_bytes(data).is_some()),
        (ct.to_bytes(), |data| CiphertextSlot::from_bytes(data).is_some()),
    ];
    let mut accepted = 0;
    for round in 0..300 {
        for (good, parses) in &parsers {
            let header = good.len() % 4;
            let mut data = good.clone();
            // Up to two words replaced by a value around a modulus or a random one.
            for _ in 0..rng.next_u32() % 3 {
                let at = header + rng.next_u32() as usize % (data.len() / 4) * 4;
                let word = [0, 1, T - 1, T, Q2 - 1, Q2, u32::MAX, rng.next_u32()][rng.next_u32() as usize % 8];
                data[at..at + 4].copy_from_slice(&word.to_le_bytes());
            }
            // One round in eight replaces the version byte of a key.
            if round % 8 == 5 {
                data[..header].fill(rng.next_u32() as u8);
            }
            if round % 3 == 0 {
                data.resize(good.len() + rng.next_u32() as usize % 9 - 4, 0);
            }
            if round % 50 == 49 {
                rng.fill_bytes(&mut data);
            }
            let valid = parses(&data);
            let well_formed = || data[..header] == good[..header] && reduced(&data[header..]);
            assert!(!valid || (data.len() == good.len() && well_formed()), "round {round}: invalid input accepted");
            accepted += valid as usize;
            if let Some(ct) = CiphertextSlot::from_bytes(&data) {
                let _ = payload(&sk, &ct);
            }
        }
    }
    // The perturbations leave a share of the inputs valid.
    assert!((100..800).contains(&accepted), "{accepted} inputs accepted");
}

/// FNV-1a, to pin serialized values in the vectors below.
fn digest(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| (h ^ b as u64).wrapping_mul(0x0100_0000_01b3))
}

/// Known answers. Key generation, encryption and re-randomization are
/// deterministic functions of the generator stream, on every platform: a key pair
/// derived from a ChaCha20 seed never changes.
#[test]
fn test_known_answers() {
    // (seed byte, digests of the public key and of the secret key after their version
    // byte, of the ciphertext and of the re-randomized ciphertext)
    let vectors: [(u8, [u64; 4]); 3] = [
        (1, [0x535a36f8fb8f7970, 0xef883001f8d26859, 0xad160b239bbf48d9, 0x305b63e4b2eb74b7]),
        (2, [0x6f6693337c982f14, 0x4b97040803cac428, 0xb08b752027c8fa2b, 0xe370ac6aee44651d]),
        (255, [0x531fa3634e572be1, 0x6d2e3075808998d9, 0x59df3719c8cad91d, 0x4b5d75d692a7406a]),
    ];
    let got = vectors.map(|(seed, _)| {
        let mut rng = ChaCha20Rng::from_seed([seed; 32]);
        let (sk, pk) = keygen(&mut rng);
        let data: [u8; SLOT_BYTES] = std::array::from_fn(|i| (i as u8).wrapping_mul(seed).wrapping_add(seed));
        let ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
        let rerandomized = rerandomize_slot(&mut rng, &pk, &ct);
        let digests = [
            digest(&pk.to_bytes()[1..]),
            digest(&sk.to_bytes()[1..]),
            digest(&ct.to_bytes()),
            digest(&rerandomized.to_bytes()),
        ];
        (seed, digests)
    });
    assert_eq!(got, vectors, "{got:#018x?}");
}
