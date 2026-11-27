//! Memory-residue scan, for Linux.
//!
//! A secret that stays in memory after its use can be read later from a memory
//! image, a swap or hibernation file, or a crash report. This program runs the
//! operations on known secrets, drops everything, and then searches all writable
//! memory of the process (heap, statics and the stacks of all threads, read
//! through `/proc/self/mem`) for these secrets and for values that reveal them:
//!   * the keys and the output of the generators involved, and the keys of the
//!     per-slot generators of the batch layer;
//!   * the secret key and the error of its public key, the randomness of an
//!     encryption and of a re-randomization, and the plaintext, in the forms in
//!     which the library holds them (as sampled, reduced, transformed), and the
//!     payload;
//!   * what the library computes from them on the way to its results: the
//!     products of a secret with a part of the public key or of a ciphertext,
//!     the encryption of zero that a message is added to or that a
//!     re-randomization adds, and the noisy plaintext that a decryption starts
//!     from.
//!
//! Nothing that an operation leaves is overwritten before the search: freed
//! memory is never handed out again, and every operation runs on threads of its
//! own, which are kept until the end.
//!
//! Run with:
//!   cargo run --release --example residue
//! on every target the crate is deployed to, with the flags of the deployed
//! build: what a compiler copies to the stack is a property of a build, not of
//! the source. On x86 it also depends on the CPU, as the generator selects its
//! code at run time. The exit status is zero only if nothing was found.

// The transform of the library, compiled into this program as well, to compute
// the transformed values that are searched for.
#[allow(dead_code)]
#[path = "../src/ntt.rs"]
mod ntt;
mod params {
    pub use pq_rerand::params::*;
}

use chacha20::ChaCha20Rng;
use pq_rerand::batch::{decrypt_batch, encrypt_batch, rerandomize_batch};
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::encoding::{decode, encode, Plaintext};
use pq_rerand::encrypt::{encrypt_slot, CiphertextSlot, Flooding};
use pq_rerand::keygen::{keygen, PublicKey, SecretKey};
use pq_rerand::params::{BITS_PER_COEFF, DELTA_T, N, Q2, SLOT_BYTES, T};
use pq_rerand::poly::NTT_CONTEXT;
use pq_rerand::rerandomize::rerandomize_slot;
use pq_rerand::sampling::{add_flood, sample_gaussian_cdt};
use pq_rerand::serialize::{deserialize_slot, serialize_slot};
use rand_core::{Rng, SeedableRng};
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::{HashMap, HashSet};
use std::hint::black_box;
use std::io::{Read, Seek, SeekFrom};
use std::ops::Range;
use zeroize::{Zeroize, Zeroizing};

/// Freed memory is never handed out again: a buffer that was freed without being
/// zeroized stays as it was left, instead of being overwritten by the next
/// allocation.
struct NoReuse;

unsafe impl GlobalAlloc for NoReuse {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn dealloc(&self, _: *mut u8, _: Layout) {}
}

#[global_allocator]
static ALLOCATOR: NoReuse = NoReuse;

/// Generator keys are searched word by word, as a key may be spilled that way:
/// one counts as found when three of its eight words are (single words also
/// match unrelated data). Keys and generator output are searched in windows of 8
/// bytes as well, one at every word. Other secrets are searched in windows of 16
/// or 32 bytes that hold at least three distinct non-zero words: shorter or
/// plainer windows also match unrelated small numbers.
const WORD: usize = 4;
const KEY: usize = 32;
const ARENA: usize = 1 << 22;
/// Memory is read in pieces of this size, which overlap by one page.
const PIECE: usize = 1 << 20;
const PAGE: usize = 4096;
/// Slots per batch call.
const SLOTS: usize = 16;
/// Stack of the threads that run the operations.
const STACK: usize = 512 * 1024;

static TABLES: [ntt::NttTables; 2] = [ntt::NttTables::new(T), ntt::NttTables::new(Q2)];

struct Needle {
    name: String,
    bytes: Range<usize>,
    window: usize,
}

struct Scanner {
    /// Reference copies of the secrets. The scan skips this buffer and its read
    /// buffer; both are allocated once and never move.
    arena: Vec<u8>,
    buffer: Vec<u8>,
    needles: Vec<Needle>,
}

impl Scanner {
    fn new() -> Self {
        Scanner { arena: Vec::with_capacity(ARENA), buffer: vec![0; PIECE], needles: Vec::new() }
    }

    /// Register a secret, to be searched for in windows of `window` bytes. The
    /// bytes arrive one by one, so that no copy of them is made on the way.
    fn add(&mut self, name: &str, window: usize, bytes: impl IntoIterator<Item = u8>) -> usize {
        let start = self.arena.len();
        self.arena.extend(bytes);
        assert!(self.arena.len() <= ARENA, "the arena must not move");
        self.needles.push(Needle { name: name.to_string(), bytes: start..self.arena.len(), window });
        self.needles.len() - 1
    }

    /// The reference copy of a secret.
    fn bytes(&self, needle: usize) -> &[u8] {
        &self.arena[self.needles[needle].bytes.clone()]
    }

    /// Search the writable memory of the process. Returns, per needle, the number
    /// of its windows (or keys) that were found, or nothing where the memory of
    /// the process cannot be read this way.
    fn scan(&mut self) -> Option<Vec<usize>> {
        let maps = std::fs::read_to_string("/proc/self/maps").ok()?;
        let mut memory = std::fs::File::open("/proc/self/mem").ok()?;
        let address = |v: &Vec<u8>| v.as_ptr() as usize..v.as_ptr() as usize + v.capacity();
        let skipped = [address(&self.arena), address(&self.buffer)];

        // Window -> (needle, index of the window in the needle).
        let mut windows: HashMap<&[u8], (usize, usize)> = HashMap::new();
        for (n, needle) in self.needles.iter().enumerate() {
            // Windows of 8 bytes start at every word, the others follow one another,
            // from the first whole word: after the version byte of a serialized key.
            let step = if needle.window == 8 { WORD } else { needle.window };
            let bytes = &self.arena[needle.bytes.clone()];
            for (i, window) in bytes[bytes.len() % WORD..].windows(needle.window).step_by(step).enumerate() {
                let words: HashSet<_> = window.as_chunks::<WORD>().0.iter().filter(|w| **w != [0; WORD]).collect();
                if needle.window < 16 || words.len() >= 3 {
                    windows.insert(window, (n, i));
                }
            }
        }

        let mut found = vec![HashSet::new(); self.needles.len()];
        for line in maps.lines() {
            let mut fields = line.split_whitespace();
            let (range, permissions) = (fields.next().unwrap(), fields.next().unwrap());
            if !permissions.starts_with("rw") {
                continue;
            }
            let (start, end) = range.split_once('-').unwrap();
            let (mut at, end) = (usize::from_str_radix(start, 16).unwrap(), usize::from_str_radix(end, 16).unwrap());
            loop {
                let piece = &mut self.buffer[..PIECE.min(end - at)];
                // Page by page: a page that cannot be read, such as the guard page
                // of a thread stack, counts as zeros.
                for (i, page) in piece.chunks_mut(PAGE).enumerate() {
                    let start = SeekFrom::Start((at + i * PAGE) as u64);
                    if memory.seek(start).and_then(|_| memory.read_exact(page)).is_err() {
                        page.fill(0);
                    }
                }
                for offset in 0..piece.len() {
                    let inside = |r: &Range<usize>| at + offset + KEY > r.start && at + offset < r.end;
                    if skipped.iter().any(inside) {
                        continue;
                    }
                    for width in [WORD, 8, 16, 32] {
                        if let Some(&(n, i)) = piece.get(offset..offset + width).and_then(|w| windows.get(w)) {
                            found[n].insert(i);
                        }
                    }
                }
                if at + piece.len() == end {
                    break;
                }
                at += PIECE - PAGE;
            }
        }

        let per_key = KEY / WORD;
        let count = |(needle, hits): (&Needle, &HashSet<usize>)| match needle.window {
            // Keys found: those with at least three of their words in memory.
            WORD => (0..needle.bytes.len() / KEY)
                .filter(|k| hits.iter().filter(|&&i| i / per_key == *k).count() >= 3)
                .count(),
            _ => hits.len(),
        };
        Some(self.needles.iter().zip(&found).map(count).collect())
    }
}

/// Run `operation` on new threads, which are kept, and far below the first frames
/// of the thread that starts it: what the operation leaves on a stack is out of
/// the reach of everything that runs afterwards.
fn isolated<T: Send>(kept: &mut Vec<ThreadPool>, threads: usize, operation: impl FnOnce() -> T + Send) -> T {
    #[inline(never)]
    fn call<T>(operation: impl FnOnce() -> T) -> T {
        operation()
    }
    #[inline(never)]
    fn deep<T>(operation: impl FnOnce() -> T) -> T {
        let mut padding = [0u8; 64 * 1024];
        black_box(&mut padding);
        let out = call(operation);
        // Used after the call, which is therefore not a tail call.
        black_box(&mut padding);
        out
    }
    let pool = ThreadPoolBuilder::new().num_threads(threads).stack_size(STACK).build().expect("threads");
    let out = pool.install(|| deep(operation));
    kept.push(pool);
    out
}

/// Byte `i` of the test value number `tag`: generator keys and the payload.
fn byte(tag: u8, i: usize) -> u8 {
    let x = (tag as u32 + 1).wrapping_mul(0x9E37_79B1) ^ (i as u32 + 1).wrapping_mul(0x85EB_CA77);
    ((x ^ (x >> 15)).wrapping_mul(0x2C1B_3C6D) >> 24) as u8
}

/// Overwrite the stack below the caller, where this program itself handled
/// secrets. The library has to do the same for the generators it runs.
#[inline(never)]
fn wipe_stack() {
    black_box(&mut [0u8; 64 * 1024]);
}

/// The generator with key number `tag`, on the heap.
fn generator(tag: u8) -> Box<ChaCha20Rng> {
    #[inline(never)]
    fn build(tag: u8) -> Box<ChaCha20Rng> {
        Box::new(ChaCha20Rng::from_seed(std::array::from_fn(|i| byte(tag, i))))
    }
    let rng = build(tag);
    wipe_stack();
    rng
}

/// `N` zero coefficients on the heap, zeroized on drop: this program must not
/// leave copies of the secrets behind either.
fn zeros<T: Copy + Default + Zeroize>() -> Box<Zeroizing<[T; N]>> {
    Box::new(Zeroizing::new([T::default(); N]))
}

/// The `i`-th array of `N` little-endian words in `bytes`, after the version
/// byte if `bytes` is a serialized key.
fn limb(bytes: &[u8], i: usize) -> Box<Zeroizing<[u32; N]>> {
    let mut limb = zeros();
    let words = bytes[bytes.len() % WORD + 4 * N * i..].as_chunks::<WORD>().0;
    limb.iter_mut().zip(words).for_each(|(coefficient, word)| *coefficient = u32::from_le_bytes(*word));
    limb
}

/// The bytes of an array of words.
fn le(words: &[u32; N]) -> impl Iterator<Item = u8> + '_ {
    words.iter().flat_map(|word| word.to_le_bytes())
}

/// Register the keys of the five generators and what the operations draw from
/// them: the output of generators 1 to 3 (key generation, encryption,
/// re-randomization), the randomness `R`, `E1`, `E2` of the encryption and of the
/// re-randomization with the products `a*R` and `b*R` under the public key `pk`,
/// and the per-slot keys of the batch calls made with generators 4 and 5.
#[inline(never)]
fn register_randomness(scanner: &mut Scanner, pk: &PublicKey) {
    let pk_bytes = pk.to_bytes();
    for tag in 1..=5 {
        scanner.add(&format!("key of generator {tag}"), WORD, (0..KEY).map(|i| byte(tag, i)));
        scanner.add(&format!("key of generator {tag}, 8 bytes"), 8, (0..KEY).map(|i| byte(tag, i)));
    }
    for tag in 1..=3 {
        // As much as the largest operation draws: a flooded encryption takes 15 N
        // words of 64 bits.
        let mut rng = generator(tag);
        scanner.add(&format!("output of generator {tag}"), 8, (0..30 * N).flat_map(|_| rng.next_u32().to_le_bytes()));
    }
    for (tag, flood, what) in [(2, true, "encryption"), (3, false, "re-randomization")] {
        let mut rng = generator(tag);
        let mut samples = [zeros::<i32>(), zeros(), zeros()];
        samples.iter_mut().for_each(|s| sample_gaussian_cdt(&mut *rng, s));
        if flood {
            samples.iter_mut().for_each(|s| add_flood(&mut *rng, s));
        }
        let window = if flood { 16 } else { 32 };
        for (s, name) in samples.iter().zip(["R", "E1", "E2"]) {
            scanner.add(&format!("{what}: {name}"), window, s.iter().flat_map(|x| x.to_le_bytes()));
        }
        for (i, tables) in TABLES.iter().enumerate() {
            let p = tables.modulus;
            let mut r = zeros();
            r.iter_mut().zip(samples[0].iter()).for_each(|(r, &x)| *r = (x as i64).rem_euclid(p as i64) as u32);
            scanner.add(&format!("{what}: R mod {p}"), window, le(&r));
            tables.forward(&mut r);
            scanner.add(&format!("{what}: NTT(R) mod {p}"), 16, le(&r));
            // The products a*R and b*R, transformed and in the coefficient domain.
            for (j, name) in [(i, "a"), (2 + i, "b")] {
                let mut product = limb(&pk_bytes, j);
                ntt::pointwise_mul(&mut product, &r, p);
                scanner.add(&format!("{what}: NTT({name}*R) mod {p}"), 16, le(&product));
                tables.inverse(&mut product);
                scanner.add(&format!("{what}: {name}*R mod {p}"), 16, le(&product));
            }
        }
    }
    for tag in [4, 5] {
        // Slot `i` is keyed with the `i`-th 32 bytes of the caller's generator.
        let mut rng = generator(tag);
        let mut keys = Zeroizing::new(vec![0u8; SLOTS * KEY]);
        keys.as_chunks_mut::<KEY>().0.iter_mut().for_each(|key| rng.fill_bytes(key));
        scanner.add(&format!("slot keys of generator {tag}"), WORD, keys.iter().copied());
        scanner.add(&format!("slot keys of generator {tag}, 8 bytes"), 8, keys.iter().copied());
    }
}

/// Register the secret `s` and the error `e = b - a*s` of a key pair, in the forms
/// that the library computes. Returns the needles of `s` and of its transform
/// modulo `q2`.
#[inline(never)]
fn register_key(scanner: &mut Scanner, sk: &SecretKey, pk: &PublicKey) -> [usize; 2] {
    let (sk_bytes, pk_bytes) = (sk.to_bytes(), pk.to_bytes());
    let mut needles = [scanner.add("key: s", 32, sk_bytes.iter().copied()); 2];
    for (i, tables) in TABLES.iter().enumerate() {
        let p = tables.modulus;
        let (mut s, mut e) = (limb(&sk_bytes, i), limb(&pk_bytes, i));
        tables.forward(&mut s);
        needles[1] = scanner.add(&format!("key: NTT(s) mod {p}"), 16, le(&s));
        // The product a*s, then a*s - b = -e, transformed and in the coefficient
        // domain.
        ntt::pointwise_mul(&mut e, &s, p);
        scanner.add(&format!("key: NTT(a*s) mod {p}"), 16, le(&e));
        ntt::pointwise_sub(&mut e, &limb(&pk_bytes, 2 + i), p);
        scanner.add(&format!("key: NTT(-e) mod {p}"), 16, le(&e));
        scanner.add(&format!("key: NTT(e) mod {p}"), 16, e.iter().flat_map(|&x| ((p - x) % p).to_le_bytes()));
        tables.inverse(&mut e);
        scanner.add(&format!("key: -e mod {p}"), 32, le(&e));
        scanner.add(&format!("key: e mod {p}"), 32, e.iter().flat_map(|&x| ((p - x) % p).to_le_bytes()));
        if p == T {
            // Both as sampled: the small signed values behind the limbs.
            let signed = |x: u32| if x > p / 2 { (p - x) as i32 } else { -(x as i32) };
            scanner.add("key: e as sampled", 32, e.iter().flat_map(|&x| signed(x).to_le_bytes()));
            let s = limb(&sk_bytes, i);
            scanner.add("key: s as sampled", 32, s.iter().flat_map(|&x| (-signed(x)).to_le_bytes()));
        }
    }
    needles
}

/// Register what the decryption of `ct` computes: the product `c1*s`, transformed
/// and in the coefficient domain, and the noisy plaintext `c0 - c1*s`.
#[inline(never)]
fn register_product(scanner: &mut Scanner, key: usize, ct: &CiphertextSlot) {
    let ct_bytes = serialize_slot(ct);
    for (i, tables) in TABLES.iter().enumerate() {
        let p = tables.modulus;
        let (mut s, mut c1, mut v) = (limb(scanner.bytes(key), i), limb(&ct_bytes, 2 + i), limb(&ct_bytes, i));
        tables.forward(&mut s);
        tables.forward(&mut c1);
        ntt::pointwise_mul(&mut c1, &s, p);
        scanner.add(&format!("decryption: NTT(c1*s) mod {p}"), 16, le(&c1));
        tables.inverse(&mut c1);
        scanner.add(&format!("decryption: c1*s mod {p}"), 16, le(&c1));
        ntt::pointwise_sub(&mut v, &c1, p);
        scanner.add(&format!("decryption: c0 - c1*s mod {p}"), 16, le(&v));
    }
}

/// Register the encryption of zero that the re-randomization of `fresh` into `ct`
/// added: it links the two.
#[inline(never)]
fn register_added(scanner: &mut Scanner, fresh: &CiphertextSlot, ct: &CiphertextSlot) {
    let (fresh_bytes, ct_bytes) = (serialize_slot(fresh), serialize_slot(ct));
    for i in 0..4 {
        let mut added = limb(&ct_bytes, i);
        ntt::pointwise_sub(&mut added, &limb(&fresh_bytes, i), TABLES[i % 2].modulus);
        scanner.add(&format!("re-randomization: Enc(0), limb {i}"), 16, le(&added));
    }
}

/// A stored key in use: parsed from its bytes, moved around, then dropped.
struct Session {
    sk: SecretKey,
}

fn main() {
    let ctx = &NTT_CONTEXT;
    let mut scanner = Scanner::new();
    let mut kept = Vec::new();
    let payload = scanner.add("payload", 32, (0..SLOT_BYTES).map(|i| byte(0, i)));
    let bit = |i: usize| (byte(0, i / 8) >> (i % 8)) as u32 & 1;
    let coefficient = |c: usize| (0..BITS_PER_COEFF).map(|b| bit(c * BITS_PER_COEFF + b) << b).sum::<u32>();
    scanner.add("plaintext: M", 16, (0..N).flat_map(|c| coefficient(c).to_le_bytes()));
    let embedded = |c: usize| (DELTA_T as u64 * coefficient(c) as u64 % T as u64) as u32;
    scanner.add("plaintext: Delta_t * M mod t", 16, (0..N).flat_map(|c| embedded(c).to_le_bytes()));

    let (sk, pk) = isolated(&mut kept, 1, || keygen(&mut *generator(1), ctx));
    drop(isolated(&mut kept, 1, || sk.to_bytes()));
    register_randomness(&mut scanner, &pk);
    let [key, transformed] = register_key(&mut scanner, &sk, &pk);
    wipe_stack();

    // Control: the scan sees the key, in both of its forms, while it is in use.
    let Some(seen) = scanner.scan() else {
        println!("The memory of the process cannot be read through /proc/self/mem: nothing was checked.");
        std::process::exit(2);
    };
    println!("control: {} windows of s and {} of NTT(s) found while the key is in use\n", seen[key], seen[transformed]);
    assert!(seen[key] > 0 && seen[transformed] > 0, "the scan does not see a key in use");

    // One slot: encrypt, re-randomize, check the public key, drop the key.
    let slot: &[u8; SLOT_BYTES] = scanner.bytes(payload).try_into().unwrap();
    let message = isolated(&mut kept, 1, || encode(slot));
    let fresh = isolated(&mut kept, 1, || encrypt_slot(&mut *generator(2), ctx, &pk, &message, Flooding::Standard));
    let ct = isolated(&mut kept, 1, || rerandomize_slot(&mut *generator(3), ctx, &pk, &fresh));
    // The flood sampler on its own, outside the wiped frame of an encryption.
    isolated(&mut kept, 1, || add_flood(&mut *generator(3), &mut zeros()));
    assert!(isolated(&mut kept, 1, || pk.verify_against(&sk)));
    drop(isolated(&mut kept, 1, || sk.clone()));
    drop((sk, message));

    // The life of a stored key: open a session from the stored bytes, read, move
    // the session, read again, lock.
    let stored = scanner.bytes(key);
    let session = Session { sk: isolated(&mut kept, 1, || SecretKey::from_bytes(stored)).expect("valid key") };
    let read = |kept: &mut Vec<ThreadPool>, session: &Session| {
        let message = isolated(kept, 1, || decrypt_slot(ctx, &session.sk, &ct));
        isolated(kept, 1, || decode(&message)).expect("valid slot")[..] == slot[..]
    };
    assert!(read(&mut kept, &session));
    let session = Box::new(session);
    assert!(read(&mut kept, &session));

    // Batch calls, each on four threads.
    let messages: Vec<Plaintext> = isolated(&mut kept, 1, || (0..SLOTS).map(|_| encode(slot)).collect());
    let cts = isolated(&mut kept, 4, || encrypt_batch(&mut *generator(4), ctx, &pk, &messages, Flooding::Standard));
    let cts = isolated(&mut kept, 4, || rerandomize_batch(&mut *generator(5), ctx, &pk, &cts));
    let decrypted = isolated(&mut kept, 4, || decrypt_batch(ctx, &session.sk, &cts));
    assert!(isolated(&mut kept, 1, || decrypted.iter().all(|m| decode(m).expect("valid slot")[..] == slot[..])));

    // A slot that decode rejects: coefficient 0 raised to 2^31.
    let mut tampered = serialize_slot(&ct);
    let c0 = u32::from_le_bytes(tampered[..4].try_into().unwrap()) as u64;
    let shift = DELTA_T as u64 * ((1 << 31) - coefficient(0) as u64) % T as u64;
    tampered[..4].copy_from_slice(&(((c0 + shift) % T as u64) as u32).to_le_bytes());
    let tampered = deserialize_slot(&tampered).expect("reduced");
    let rejected = isolated(&mut kept, 1, || decrypt_slot(ctx, &session.sk, &tampered));
    assert!(isolated(&mut kept, 1, || decode(&rejected)).is_none());
    drop((session, messages, decrypted, rejected));

    // The encryption of zero that the message was added to: its first limb.
    let c0 = limb(&serialize_slot(&fresh), 0);
    let zero = (0..N).flat_map(|c| ntt::submod(c0[c], embedded(c), T).to_le_bytes());
    scanner.add("encryption: c0 - Delta_t * M mod t", 16, zero);
    register_product(&mut scanner, key, &ct);
    register_added(&mut scanner, &fresh, &ct);
    wipe_stack();
    let found = scanner.scan().expect("readable a moment ago");
    for (needle, count) in scanner.needles.iter().zip(&found) {
        let unit = if needle.window == WORD { "keys" } else { "windows" };
        println!("{:<42} {count:>5} {unit} found{}", needle.name, if *count > 0 { "  <-- RESIDUE" } else { "" });
    }
    let clean = found.iter().all(|&count| count == 0);
    println!("\n{}", if clean { "No residue found." } else { "Secrets were left in memory." });
    std::process::exit(!clean as i32);
}
