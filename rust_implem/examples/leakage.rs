//! Constant-time / timing-leakage harness (TVLA / dudect-style).
//!
//! For each target we measure CPU-cycle counts (`rdtsc`) over two input classes
//! (fixed vs. random) and apply Welch's t-test, as in the dudect methodology
//! [Reparaz–Balasch–Tunstall, DATE'17] and Goodwill et al.'s TVLA fixed-vs-random
//! test. A small |t| (conventionally `|t| < 4.5`) means no input-dependent timing
//! was detected; a large |t| flags leakage. We report the raw statistic and a
//! version with the slow-outlier tail cropped (interrupts / migrations inflate it).
//!
//! Targets (fixed-vs-random over the *data* input). To avoid a cache-locality
//! confound, both classes stream through equally-sized pools (identical memory
//! footprint and access pattern); only the data content differs. Following
//! dudect convention, the "fixed" class is an arbitrary fixed *random* constant:
//!   * `encrypt_slot` — fixed vs. random *plaintext* (tests message-dependence
//!     of the encryption path).
//!   * `rerandomize` — fixed- vs. random-plaintext *input ciphertexts* (public
//!     path; tests operand-dependence).
//!   * `decrypt_slot` — fixed- vs. random-plaintext ciphertexts; the
//!     secret-dependent path exercising the whole branchless arithmetic layer
//!     (Barrett reduction, csub/submod, center-lift).
//!   * `leaky_reduce` — POSITIVE CONTROL: a branchy reduction that *should* be
//!     flagged, proving the harness can detect leakage.
//!
//! Two additional diagnostics document a *hardware* (not software) effect:
//!   * `encrypt (all-zero)` — all-zero vs. random plaintext. On some CPUs
//!     (observed on Zen 4) this can be flagged at a ~0.02% timing difference
//!     even though the code is branchless, because operating on all-zero data
//!     costs marginally less energy/time at the silicon level (data-dependent
//!     bit-toggle effects). The signal is at the edge of detectability at this
//!     sample count and does not appear on every run.
//!   * `hw zero-load (baseline)` — a pure load+XOR loop with *no* crypto
//!     arithmetic over the same two pools. It is robustly flagged on such CPUs,
//!     proving the all-zero signal originates in the hardware data path, not in
//!     this crate's control flow or memory addressing.
//!
//! In the intended deployment this degenerate-input channel is absent by
//! construction: payloads are AEAD ciphertexts, which are indistinguishable
//! from uniformly random bytes, so an all-zero (or otherwise low-entropy)
//! plaintext never occurs.
//!
//! (We test the full operations rather than micro-benchmarking individual modular
//! primitives: a ~300-cycle primitive sits below reliable `rdtsc` resolution, where
//! the t-test reflects measurement quantization rather than the operation.)
//!
//! Run with: `cargo run --release --example leakage`
//! Results are statistical and machine-dependent; run multiple times under low load.

use pq_rerand::params::{N, T, SIGMA_FLOOD, BITS_PER_COEFF};
use pq_rerand::poly::NttContext;
use pq_rerand::keygen::{keygen, SecretKey, PublicKey};
use pq_rerand::encrypt::{encrypt_slot, CiphertextSlot};
use pq_rerand::rerandomize::rerandomize_slot;
use pq_rerand::decrypt::decrypt_slot;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use std::hint::black_box;

/// Read a high-resolution time-stamp counter.
///
/// * x86_64: serialized `rdtsc`.
/// * aarch64 (mobile / Apple silicon / ARM servers): the EL0-readable virtual
///   counter `cntvct_el0` (coarser than `rdtsc`; increase sample counts if the
///   t-test looks noisy).
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

fn report(name: &str, c0: &[f64], c1: &[f64]) {
    let t_raw = welch_t(c0, c1).abs();
    let t_c90 = cropped_t(c0, c1, 0.90).abs();
    let t_c95 = cropped_t(c0, c1, 0.95).abs();
    let t = t_raw.max(t_c90).max(t_c95);
    let verdict = if t < 4.5 { "no leakage detected" } else { "LEAKAGE DETECTED" };
    let (m0, _) = mean_var(c0);
    let (m1, _) = mean_var(c1);
    println!(
        "{name:<22} |t|_raw={t_raw:7.2}  |t|_crop90={t_c90:7.2}  |t|_crop95={t_c95:7.2}  \
         mean0={m0:.0}c mean1={m1:.0}c  -> {verdict}",
    );
}

fn rand_message(rng: &mut StdRng) -> [u32; N] {
    let mut m = [0u32; N];
    for c in m.iter_mut() {
        *c = rng.gen::<u32>() & ((1u32 << BITS_PER_COEFF) - 1);
    }
    m
}

/// POSITIVE CONTROL: data-dependent number of loop iterations -> variable time.
#[inline(never)]
fn leaky_reduce(mut x: u64) -> u64 {
    let p = T as u64;
    while x >= p {
        x -= p;
    }
    x
}

const W: usize = 256;
const POOL: usize = 512;

/// Positive control: matched-pool driver with the branchy reducer whose
/// iteration count depends on the operand value.
fn leak_leaky(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(0x1EA4);
    let p = T as u64;
    // class 0: values in [0, p) -> 0 loop iterations.
    let pool0: Vec<[u64; W]> = (0..POOL).map(|_| [p / 2; W]).collect();
    // class 1: values in [0, 64p) -> 0..64 iterations (data-dependent).
    let pool1: Vec<[u64; W]> = (0..POOL)
        .map(|_| { let mut b = [0u64; W]; b.iter_mut().for_each(|x| *x = rng.gen::<u64>() % (64 * p)); b })
        .collect();
    let run = |buf: &[u64; W]| -> f64 {
        let s = rdtsc();
        let mut acc = 0u64;
        for &x in buf.iter() { acc ^= leaky_reduce(black_box(x)); }
        black_box(acc);
        (rdtsc() - s) as f64
    };
    let mut c0 = Vec::with_capacity(n);
    let mut c1 = Vec::with_capacity(n);
    for i in 0..n {
        c0.push(run(&pool0[i % POOL]));
        c1.push(run(&pool1[i % POOL]));
    }
    (c0, c1)
}

/// Measure `decrypt_slot` for fixed-plaintext vs. random-plaintext ciphertexts.
///
/// To avoid a cache-locality confound, *both* classes stream through equally-sized
/// ciphertext pools (same memory footprint and access pattern); the only difference
/// is the secret-dependent content (a fixed plaintext, each freshly re-encrypted with
/// independent noise, vs. random plaintexts).
fn leak_decrypt(
    ctx: &NttContext,
    sk: &SecretKey,
    pk: &PublicKey,
    n: usize,
) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(0xDEC0);
    let fixed_msg = rand_message(&mut rng); // arbitrary fixed constant (dudect)
    let pool0: Vec<CiphertextSlot> = (0..POOL)
        .map(|_| encrypt_slot(&mut rng, ctx, pk, &fixed_msg, SIGMA_FLOOD))
        .collect();
    let pool1: Vec<CiphertextSlot> = (0..POOL)
        .map(|_| { let m = rand_message(&mut rng); encrypt_slot(&mut rng, ctx, pk, &m, SIGMA_FLOOD) })
        .collect();
    let run = |ct: &CiphertextSlot| -> f64 {
        let s = rdtsc();
        black_box(decrypt_slot(ctx, sk, black_box(ct)));
        (rdtsc() - s) as f64
    };
    let mut c0 = Vec::with_capacity(n);
    let mut c1 = Vec::with_capacity(n);
    for i in 0..n {
        c0.push(run(&pool0[i % POOL]));
        c1.push(run(&pool1[i % POOL]));
    }
    (c0, c1)
}

/// Measure `encrypt_slot` for fixed vs. random plaintexts. Both classes stream
/// through equally-sized message pools (matched footprint); only the content
/// differs. Encryption noise is drawn from one shared RNG advanced identically
/// across both classes, so its variance is unbiased. `fixed_msg` selects the
/// fixed class: a fixed random constant (dudect convention) or the degenerate
/// all-zero plaintext (hardware-effect diagnostic; see module docs).
fn leak_encrypt(
    ctx: &NttContext,
    pk: &PublicKey,
    n: usize,
    fixed_msg: &[u32; N],
) -> (Vec<f64>, Vec<f64>) {
    let mut mrng = StdRng::seed_from_u64(0xE2C0);
    let pool0: Vec<[u32; N]> = (0..POOL).map(|_| *fixed_msg).collect();
    let pool1: Vec<[u32; N]> = (0..POOL).map(|_| rand_message(&mut mrng)).collect();
    let mut erng = StdRng::seed_from_u64(0xE2C1);
    let run = |msg: &[u32; N], r: &mut StdRng| -> f64 {
        let s = rdtsc();
        black_box(encrypt_slot(r, ctx, pk, black_box(msg), SIGMA_FLOOD));
        (rdtsc() - s) as f64
    };
    let (mut c0, mut c1) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for i in 0..n {
        c0.push(run(&pool0[i % POOL], &mut erng));
        c1.push(run(&pool1[i % POOL], &mut erng));
    }
    (c0, c1)
}

/// Measure `rerandomize_slot` for fixed- vs. random-plaintext input ciphertexts.
/// Both classes stream through equally-sized ciphertext pools (matched footprint);
/// the re-randomization noise comes from one shared RNG advanced identically.
fn leak_rerand(ctx: &NttContext, pk: &PublicKey, n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut erng = StdRng::seed_from_u64(0x4E4D);
    let mut mrng = StdRng::seed_from_u64(0x4E4F);
    let fixed_msg = rand_message(&mut mrng); // arbitrary fixed constant (dudect)
    let pool0: Vec<CiphertextSlot> = (0..POOL)
        .map(|_| encrypt_slot(&mut erng, ctx, pk, &fixed_msg, SIGMA_FLOOD))
        .collect();
    let pool1: Vec<CiphertextSlot> = (0..POOL)
        .map(|_| { let m = rand_message(&mut erng); encrypt_slot(&mut erng, ctx, pk, &m, SIGMA_FLOOD) })
        .collect();
    let mut rrng = StdRng::seed_from_u64(0x4E4E);
    let run = |ct: &CiphertextSlot, r: &mut StdRng| -> f64 {
        let s = rdtsc();
        black_box(rerandomize_slot(r, ctx, pk, black_box(ct)));
        (rdtsc() - s) as f64
    };
    let (mut c0, mut c1) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for i in 0..n {
        c0.push(run(&pool0[i % POOL], &mut rrng));
        c1.push(run(&pool1[i % POOL], &mut rrng));
    }
    (c0, c1)
}

/// Hardware baseline: a pure streaming load+XOR over all-zero vs. random pools.
/// Contains no crypto arithmetic at all; any flagged difference is a property of
/// the CPU's data path (data-dependent bit-toggle energy/time), not of this crate.
fn leak_hw_zero_load(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(0x0B5E);
    let pool0: Vec<[u32; N]> = (0..POOL).map(|_| [0u32; N]).collect();
    let pool1: Vec<[u32; N]> = (0..POOL).map(|_| rand_message(&mut rng)).collect();
    let run = |buf: &[u32; N]| -> f64 {
        let s = rdtsc();
        let mut acc = 0u32;
        for _ in 0..8 {
            for &x in buf.iter() {
                acc ^= black_box(x);
            }
        }
        let e = rdtsc();
        black_box(acc);
        (e - s) as f64
    };
    let (mut c0, mut c1) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for i in 0..n {
        c0.push(run(&pool0[i % POOL]));
        c1.push(run(&pool1[i % POOL]));
    }
    (c0, c1)
}

fn main() {
    let ctx = NttContext::new();
    let mut rng = StdRng::seed_from_u64(1);
    let (sk, pk) = keygen(&mut rng, &ctx);

    println!("TVLA fixed-vs-random timing-leakage test (Welch t; |t| < 4.5 => no leakage detected)\n");

    // Warm up caches / frequency.
    let _ = leak_leaky(5_000);

    let fixed_random = rand_message(&mut rng);
    let (a, b) = leak_encrypt(&ctx, &pk, 30_000, &fixed_random);
    report("encrypt_slot", &a, &b);

    let (a, b) = leak_rerand(&ctx, &pk, 30_000);
    report("rerandomize", &a, &b);

    let (a, b) = leak_decrypt(&ctx, &sk, &pk, 80_000);
    report("decrypt_slot", &a, &b);

    let (a, b) = leak_leaky(400_000);
    report("leaky_reduce (control)", &a, &b);

    println!("\nHardware data-path diagnostics (see module docs):");
    let (a, b) = leak_encrypt(&ctx, &pk, 30_000, &[0u32; N]);
    report("encrypt (all-zero msg)", &a, &b);

    let (a, b) = leak_hw_zero_load(200_000);
    report("hw zero-load baseline", &a, &b);

    println!(
        "\nNotes: the positive control is expected to show LEAKAGE; the constant-time\n\
         targets are expected to stay below the 4.5 threshold. On CPUs with\n\
         data-dependent bit-toggle timing (observed on Zen 4), the zero-load\n\
         baseline is robustly flagged and the all-zero encrypt diagnostic may be\n\
         flagged at a ~0.02% difference: a degenerate all-zero input can be\n\
         distinguishable at the hardware level even for branchless code (the\n\
         zero-load baseline contains no crypto arithmetic). AEAD-wrapped payloads\n\
         (the intended deployment) are uniformly random, so this input never occurs.\n\
         Re-run under low load."
    );
}
