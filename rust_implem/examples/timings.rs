//! Standalone timing harness for the comparison table.
//!
//! Run with:
//!   RUSTFLAGS="-C target-cpu=native" cargo run --release --example timings
//! (target-cpu=native lets the compiler auto-vectorize the branchless
//! Barrett/NTT kernels with the host SIMD ISA; omit it for a portable build.)
//!
//! Reports median wall-clock timings for single-slot and 1024-slot batch
//! operations, both single-threaded and multi-threaded (rayon).

use pq_rerand::params::NUM_SLOTS;
use pq_rerand::poly::NTT_CONTEXT;
use pq_rerand::keygen::keygen;
use pq_rerand::encoding::{encode, Plaintext};
use pq_rerand::encrypt::{encrypt_slot, Flooding};
use pq_rerand::rerandomize::rerandomize_slot;
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::batch::{encrypt_batch, rerandomize_batch, decrypt_batch};
use std::time::{Duration, Instant};

fn median(mut v: Vec<Duration>) -> Duration {
    v.sort();
    v[v.len() / 2]
}

fn time_n(iters: usize, mut f: impl FnMut()) -> Duration {
    let mut samples = Vec::with_capacity(iters);
    for _ in 0..iters {
        let t = Instant::now();
        f();
        samples.push(t.elapsed());
    }
    median(samples)
}

fn ms(d: Duration) -> f64 { d.as_secs_f64() * 1e3 }

fn main() {
    let threads = rayon::current_num_threads();
    let ctx = &NTT_CONTEXT;
    // Seeded from the operating system, as every generator given to this crate
    // must be.
    let mut rng = rand::rng();
    let (sk, pk) = keygen(&mut rng, ctx);

    let message = encode(&std::array::from_fn(|i| i as u8));
    let messages: Vec<Plaintext> = (0..NUM_SLOTS).map(|s| encode(&std::array::from_fn(|i| (i + s) as u8))).collect();

    // Warm up (caches, rayon pool spin-up).
    let ct0 = encrypt_slot(&mut rng, ctx, &pk, &message, Flooding::Standard);
    let _ = decrypt_slot(ctx, &sk, &ct0);
    let _ = encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard);

    // ---- Single-slot ----
    let enc = time_n(200, || {
        let _ = encrypt_slot(&mut rng, ctx, &pk, &message, Flooding::Standard);
    });
    let ct = encrypt_slot(&mut rng, ctx, &pk, &message, Flooding::Standard);
    let rer = time_n(400, || {
        let _ = rerandomize_slot(&mut rng, ctx, &pk, &ct);
    });
    let dec = time_n(400, || {
        let _ = decrypt_slot(ctx, &sk, &ct);
    });

    // ---- Batch (B = 1024), single-threaded sequential ----
    let enc_seq = time_n(7, || {
        let _: Vec<_> = messages.iter()
            .map(|m| encrypt_slot(&mut rng, ctx, &pk, m, Flooding::Standard)).collect();
    });
    let cts = encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard);
    let rer_seq = time_n(7, || {
        let _: Vec<_> = cts.iter()
            .map(|c| rerandomize_slot(&mut rng, ctx, &pk, c)).collect();
    });
    let dec_seq = time_n(7, || {
        let _: Vec<_> = cts.iter().map(|c| decrypt_slot(ctx, &sk, c)).collect();
    });

    // ---- Batch (B = 1024), multi-threaded (rayon) ----
    let enc_par = time_n(7, || {
        let _ = encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard);
    });
    let rer_par = time_n(7, || {
        let _ = rerandomize_batch(&mut rng, ctx, &pk, &cts);
    });
    let dec_par = time_n(7, || {
        let _ = decrypt_batch(ctx, &sk, &cts);
    });

    println!("threads (rayon)        : {threads}");
    println!();
    println!("=== Single ciphertext (median) ===");
    println!("  encrypt (base+flood) : {:.3} ms", ms(enc));
    println!("  re-randomize         : {:.3} ms", ms(rer));
    println!("  decrypt              : {:.3} ms", ms(dec));
    println!();
    println!("=== Batch B=1024, single-threaded (median) ===");
    println!("  encrypt              : {:.3} s  ({:.3} ms/ct)", enc_seq.as_secs_f64(), ms(enc_seq) / NUM_SLOTS as f64);
    println!("  re-randomize         : {:.3} s  ({:.3} ms/ct)", rer_seq.as_secs_f64(), ms(rer_seq) / NUM_SLOTS as f64);
    println!("  decrypt              : {:.3} s  ({:.3} ms/ct)", dec_seq.as_secs_f64(), ms(dec_seq) / NUM_SLOTS as f64);
    println!();
    println!("=== Batch B=1024, multi-threaded rayon ({threads} threads, median) ===");
    println!("  encrypt              : {:.3} s  ({:.3} ms/ct)", enc_par.as_secs_f64(), ms(enc_par) / NUM_SLOTS as f64);
    println!("  re-randomize         : {:.3} s  ({:.3} ms/ct)", rer_par.as_secs_f64(), ms(rer_par) / NUM_SLOTS as f64);
    println!("  decrypt              : {:.3} s  ({:.3} ms/ct)", dec_par.as_secs_f64(), ms(dec_par) / NUM_SLOTS as f64);
}
