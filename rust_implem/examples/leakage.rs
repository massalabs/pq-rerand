//! Constant-time / timing-leakage harness (TVLA / dudect-style).
//!
//! For each target we measure timer ticks (`rdtsc`) over two input classes
//! (fixed vs. random) and apply Welch's t-test, as in the dudect methodology
//! [Reparaz–Balasch–Tunstall, DATE'17] and Goodwill et al.'s TVLA fixed-vs-random
//! test. A small |t| (conventionally `|t| < 4.5`) means no input-dependent timing
//! was detected; a large |t| flags leakage. We report the raw statistic and
//! versions with the slow-outlier tail cropped (interrupts / migrations inflate it).
//!
//! The class of every measurement is drawn at random, so that drift (frequency
//! scaling, temperature) cannot line up with a class, and the two classes only
//! differ in the content of one input:
//!   * the generator stream, for the samplers, key generation (with the public
//!     polynomial `a` held fixed: its rejection sampling takes a time that depends
//!     on `a`, which is public), encryption and re-randomization: the fixed class
//!     always replays one stream, the random class a fresh one;
//!   * the secret key, for decryption and for secret-key parsing;
//!   * the plaintext, for encryption, and the ciphertext, for re-randomization
//!     and decryption.
//!
//! Whatever the class, the same work is done before a measurement and the
//! operation runs on the same memory: an input that lives in memory is first
//! assembled in one working buffer from a fixed value and a random one, which are
//! both read in full and of which a mask keeps one. Only the content of that
//! buffer differs between the classes: the memory that was accessed, and with it
//! the state of the caches, does not depend on the class. (Inputs read from one
//! memory pool per class already differ in timing by that fact alone, at a level
//! that a few hundred thousand measurements detect.)
//!
//! Two controls calibrate the harness: one gives both classes the same inputs
//! and must not be flagged (false-positive check), the other runs a loop with one
//! data-dependent branch per coefficient and must be flagged (sensitivity check).
//!
//! (We test the full operations rather than micro-benchmarking individual modular
//! primitives: a ~300-cycle primitive sits below reliable `rdtsc` resolution, where
//! the t-test reflects measurement quantization rather than the operation.)
//!
//! Run with: `cargo run --release --example leakage [scale]`, where `scale`
//! multiplies the number of measurements (default 1).
//! Results are statistical and machine-dependent; run multiple times under low
//! load, on every target the crate is deployed to.

use pq_rerand::params::{N, T, SLOT_BYTES};
use pq_rerand::poly::NTT_CONTEXT;
use pq_rerand::keygen::{keygen, SecretKey};
use pq_rerand::encoding::encode;
use pq_rerand::encrypt::{encrypt_slot, Flooding};
use pq_rerand::rerandomize::rerandomize_slot;
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::sampling::{add_flood, sample_gaussian_cdt};
use pq_rerand::serialize::{deserialize_slot, serialize_slot};
use rand::{Rng, SeedableRng, TryCryptoRng, TryRng};
use rand::rngs::StdRng;
use std::convert::Infallible;
use std::hint::black_box;

/// Read a high-resolution time-stamp counter.
///
/// * x86_64: serialized `rdtsc`.
/// * aarch64 (mobile / Apple silicon / ARM servers): the EL0-readable virtual
///   counter `cntvct_el0`, whose resolution depends on the platform (from one
///   nanosecond to tens of them; increase the scale if the t-test looks noisy).
/// * other targets: nanoseconds from a monotonic clock.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
fn rdtsc() -> u64 {
    // Serialize, then read the timestamp counter.
    unsafe {
        core::arch::x86_64::_mm_lfence();
        let t = core::arch::x86_64::_rdtsc();
        core::arch::x86_64::_mm_lfence();
        t
    }
}

#[cfg(target_arch = "aarch64")]
#[inline(always)]
fn rdtsc() -> u64 {
    let t: u64;
    unsafe {
        core::arch::asm!("isb", "mrs {t}, cntvct_el0", t = out(reg) t, options(nostack, preserves_flags));
    }
    t
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
#[inline(always)]
fn rdtsc() -> u64 {
    use std::time::Instant;
    use std::sync::OnceLock;
    static START: OnceLock<Instant> = OnceLock::new();
    START.get_or_init(Instant::now).elapsed().as_nanos() as u64
}

/// Ticks taken by one call of `f`.
#[inline(always)]
fn ticks<T>(f: impl FnOnce() -> T) -> f64 {
    let start = rdtsc();
    black_box(f());
    rdtsc().wrapping_sub(start) as f64
}

fn mean_var(xs: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let mean = xs.iter().sum::<f64>() / n;
    let var = xs.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / (n - 1.0);
    (mean, var)
}

/// Welch's two-sample t-statistic.
fn welch_t(a: &[f64], b: &[f64]) -> f64 {
    let (ma, va) = mean_var(a);
    let (mb, vb) = mean_var(b);
    (ma - mb) / (va / a.len() as f64 + vb / b.len() as f64).sqrt()
}

/// dudect-style cropped t: drop the slow tail using a *common* threshold (the
/// `keep`-quantile of the pooled measurements), so the same cutoff is applied to
/// both classes (unbiased), then run Welch's test on what remains.
fn cropped_t(c0: &[f64], c1: &[f64], keep: f64) -> f64 {
    let mut all: Vec<f64> = c0.iter().chain(c1.iter()).copied().collect();
    all.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let thr = all[(((all.len() as f64) * keep) as usize).min(all.len() - 1)];
    let a: Vec<f64> = c0.iter().copied().filter(|&x| x <= thr).collect();
    let b: Vec<f64> = c1.iter().copied().filter(|&x| x <= thr).collect();
    welch_t(&a, &b)
}

/// Take `n` measurements of `measure(class)`, drawing the class of each one at
/// random, and report the t-statistics between the two classes.
fn run(name: &str, n: usize, mut measure: impl FnMut(usize) -> f64) {
    let mut select = StdRng::seed_from_u64(0xC1A55);
    let (mut c0, mut c1) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for _ in 0..n {
        let class = (select.next_u32() & 1) as usize;
        let t = measure(class);
        if class == 0 { c0.push(t) } else { c1.push(t) }
    }
    let t_raw = welch_t(&c0, &c1).abs();
    let t_c90 = cropped_t(&c0, &c1, 0.90).abs();
    let t_c95 = cropped_t(&c0, &c1, 0.95).abs();
    let t = t_raw.max(t_c90).max(t_c95);
    let verdict = if t < 4.5 { "no leakage detected" } else { "LEAKAGE DETECTED" };
    let (m0, _) = mean_var(&c0);
    let (m1, _) = mean_var(&c1);
    println!(
        "{name:<34} |t|_raw={t_raw:7.2}  |t|_crop90={t_c90:7.2}  |t|_crop95={t_c95:7.2}  \
         mean0={m0:.0} mean1={m1:.0}  -> {verdict}",
    );
}

/// Generator that answers 32-bit draws from one stream and 64-bit draws from
/// another. Key generation draws the public polynomial `a` in 32-bit words and
/// its secrets in 64-bit words, so this holds `a` fixed while the secrets vary.
struct Split {
    words32: StdRng,
    words64: StdRng,
}

impl TryRng for Split {
    type Error = Infallible;
    fn try_next_u32(&mut self) -> Result<u32, Infallible> {
        Ok(self.words32.next_u32())
    }
    fn try_next_u64(&mut self) -> Result<u64, Infallible> {
        Ok(self.words64.next_u64())
    }
    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), Infallible> {
        self.words64.fill_bytes(dst);
        Ok(())
    }
}

// Only so that the type is accepted by the functions under test.
impl TryCryptoRng for Split {}

/// POSITIVE CONTROL: one data-dependent branch per coefficient, as a compiler
/// could emit for a conditional subtraction -> variable time.
#[inline(never)]
fn branchy_sum(buf: &[u32; N]) -> u32 {
    let mut acc = 0u32;
    for &x in buf.iter() {
        if x >= T {
            acc = black_box(acc.wrapping_add(x - T));
        } else {
            acc = acc.wrapping_add(x);
        }
    }
    acc
}

const POOL: usize = 64;

/// Fill `work` with `fixed` (class 0) or with `random` (class 1). Both are read in
/// full whatever the class, and a mask keeps one of them.
fn select(work: &mut [u8], fixed: &[u8], random: &[u8], class: usize) {
    let mask = black_box((class as u8).wrapping_neg());
    for ((w, &f), &r) in work.iter_mut().zip(fixed).zip(random) {
        *w = f ^ ((f ^ r) & mask);
    }
}

fn random_slot(rng: &mut StdRng) -> [u8; SLOT_BYTES] {
    let mut data = [0u8; SLOT_BYTES];
    rng.fill_bytes(&mut data);
    data
}

fn main() {
    let scale: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(1);
    let ctx = &NTT_CONTEXT;
    // Fixed test inputs: reproducible measurements, not a way to seed real keys.
    let mut rng = StdRng::seed_from_u64(1);
    let (sk, pk) = keygen(&mut rng, ctx);
    let fixed_data = random_slot(&mut rng);
    let fixed_msg = encode(&fixed_data);
    let fixed_ct = encrypt_slot(&mut rng, ctx, &pk, &fixed_msg, Flooding::Standard);

    // Generator streams: the fixed class always replays one, the random class
    // takes a fresh one (which is drawn in both cases).
    let fixed_seed = [0x5au8; 32];
    let mut fresh = StdRng::seed_from_u64(0xF2E5);
    let mut seed_for = |class: usize| -> [u8; 32] {
        let mut seed = [0u8; 32];
        fresh.fill_bytes(&mut seed);
        if class == 0 { fixed_seed } else { seed }
    };

    // Inputs that live in memory, in serialized form: the fixed value and a pool
    // of random ones.
    let fixed_key = sk.to_bytes();
    let random_keys: Vec<Vec<u8>> = (0..POOL).map(|_| keygen(&mut rng, ctx).0.to_bytes().to_vec()).collect();
    let random_payloads: Vec<[u8; SLOT_BYTES]> = (0..POOL).map(|_| random_slot(&mut rng)).collect();
    let fixed_ct_bytes = serialize_slot(&fixed_ct);
    let encrypted = |d| serialize_slot(&encrypt_slot(&mut rng, ctx, &pk, &encode(d), Flooding::Standard));
    let random_cts: Vec<Vec<u8>> = random_payloads.iter().map(encrypted).collect();
    // For the branch control only: one pool of inputs per class.
    let below: [u32; N] = std::array::from_fn(|_| rng.next_u32() % T);
    let anywhere: Vec<[u32; N]> = (0..POOL).map(|_| std::array::from_fn(|_| rng.next_u32())).collect();
    let words = [vec![below; POOL], anywhere];

    // The working buffers that the inputs are assembled in.
    let mut samples = [0i32; N];
    let mut work_words = [0u32; N];
    let mut work_key = fixed_key.to_vec();
    let mut work_payload = fixed_data;
    let mut work_ct = fixed_ct_bytes.clone();

    println!("TVLA fixed-vs-random timing-leakage test (Welch t; |t| < 4.5 => no leakage detected)\n");
    let mut i = 0;
    let mut next = || {
        i = (i + 1) % POOL;
        i
    };

    println!("Controls (the first must pass, the second must be flagged):");
    run("control: same inputs", 20_000 * scale, |_| {
        let mut g = StdRng::from_seed(fixed_seed);
        ticks(|| sample_gaussian_cdt(&mut g, &mut samples))
    });
    run("control: one branch per coeff", 20_000 * scale, |class| {
        work_words = words[class][next()];
        ticks(|| branchy_sum(black_box(&work_words)))
    });

    println!("\nSecret randomness (fixed vs. random generator stream):");
    run("base sampler", 20_000 * scale, |class| {
        let mut g = StdRng::from_seed(seed_for(class));
        ticks(|| sample_gaussian_cdt(&mut g, &mut samples))
    });
    run("flood sampler", 20_000 * scale, |class| {
        let mut g = StdRng::from_seed(seed_for(class));
        ticks(|| add_flood(&mut g, &mut samples))
    });
    run("keygen (fixed a)", 5_000 * scale, |class| {
        let mut g = Split { words32: StdRng::from_seed(fixed_seed), words64: StdRng::from_seed(seed_for(class)) };
        ticks(|| keygen(&mut g, ctx))
    });
    run("encrypt_slot", 5_000 * scale, |class| {
        let mut g = StdRng::from_seed(seed_for(class));
        ticks(|| encrypt_slot(&mut g, ctx, &pk, &fixed_msg, Flooding::Standard))
    });
    run("rerandomize_slot", 5_000 * scale, |class| {
        let mut g = StdRng::from_seed(seed_for(class));
        ticks(|| rerandomize_slot(&mut g, ctx, &pk, &fixed_ct))
    });

    println!("\nSecret key (fixed vs. random key):");
    run("decrypt_slot", 20_000 * scale, |class| {
        select(&mut work_key, &fixed_key, &random_keys[next()], class);
        let key = SecretKey::from_bytes(&work_key).expect("valid key");
        ticks(|| decrypt_slot(ctx, &key, &fixed_ct))
    });
    run("SecretKey::from_bytes", 20_000 * scale, |class| {
        select(&mut work_key, &fixed_key, &random_keys[next()], class);
        ticks(|| SecretKey::from_bytes(&work_key))
    });

    println!("\nPlaintext and ciphertext (fixed vs. random):");
    run("encrypt_slot (plaintext)", 5_000 * scale, |class| {
        let mut g = StdRng::from_seed(fixed_seed);
        select(&mut work_payload, &fixed_data, &random_payloads[next()], class);
        let msg = encode(&work_payload);
        ticks(|| encrypt_slot(&mut g, ctx, &pk, &msg, Flooding::Standard))
    });
    run("rerandomize_slot (ciphertext)", 5_000 * scale, |class| {
        let mut g = StdRng::from_seed(fixed_seed);
        select(&mut work_ct, &fixed_ct_bytes, &random_cts[next()], class);
        let ct = deserialize_slot(&work_ct).expect("valid ciphertext");
        ticks(|| rerandomize_slot(&mut g, ctx, &pk, &ct))
    });
    run("decrypt_slot (ciphertext)", 20_000 * scale, |class| {
        select(&mut work_ct, &fixed_ct_bytes, &random_cts[next()], class);
        let ct = deserialize_slot(&work_ct).expect("valid ciphertext");
        ticks(|| decrypt_slot(ctx, &sk, &ct))
    });

    println!(
        "\nNotes: the second control is expected to show LEAKAGE; the first control and\n\
         the constant-time targets are expected to stay below the 4.5 threshold.\n\
         Re-run under low load."
    );
}
