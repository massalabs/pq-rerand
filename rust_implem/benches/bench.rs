use criterion::{criterion_group, criterion_main, Criterion, black_box};
use pq_rerand::params::*;
use pq_rerand::poly::NttContext;
use pq_rerand::keygen::keygen;
use pq_rerand::encrypt::encrypt_slot;
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::rerandomize::rerandomize_slot;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::time::Duration;

fn bench_single_slot(c: &mut Criterion) {
    let ctx = NttContext::new();
    let mut rng = StdRng::seed_from_u64(12345);
    let (sk, pk) = keygen(&mut rng, &ctx);

    let mut message = [0u32; N];
    for (i, c) in message.iter_mut().enumerate() {
        *c = (i as u32) % (1 << BITS_PER_COEFF);
    }

    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);

    c.bench_function("keygen", |b| {
        b.iter(|| keygen(&mut rng, &ctx))
    });

    c.bench_function("encrypt_slot (base+smudge)", |b| {
        b.iter(|| {
            encrypt_slot(&mut rng, &ctx, &pk, black_box(&message), SIGMA_FLOOD)
        })
    });

    c.bench_function("rerandomize_slot", |b| {
        b.iter(|| {
            rerandomize_slot(&mut rng, &ctx, &pk, black_box(&ct))
        })
    });

    c.bench_function("decrypt_slot", |b| {
        b.iter(|| {
            decrypt_slot(&ctx, &sk, black_box(&ct))
        })
    });
}

fn bench_full_message(c: &mut Criterion) {
    let ctx = NttContext::new();
    let mut rng = StdRng::seed_from_u64(99999);
    let (sk, pk) = keygen(&mut rng, &ctx);

    let mut message = [0u32; N];
    for (i, c) in message.iter_mut().enumerate() { *c = (i as u32) % (1 << BITS_PER_COEFF); }
    let slots: Vec<_> = (0..NUM_SLOTS)
        .map(|_| encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD))
        .collect();

    // Heavy 15.5 MiB batch ops (single-threaded); fewer samples to keep wall-clock
    // reasonable. Multi-threaded (rayon) batch timings are produced by the
    // `timings` example, which configures a large-stack pool.
    let mut group = c.benchmark_group("full_15.5MiB_1024_slots");
    group.sample_size(10).measurement_time(Duration::from_secs(8));

    group.bench_function("rerandomize (single-threaded)", |b| {
        b.iter(|| {
            let _: Vec<_> = slots.iter()
                .map(|ct| rerandomize_slot(&mut rng, &ctx, &pk, black_box(ct)))
                .collect();
        })
    });
    group.bench_function("decrypt (single-threaded)", |b| {
        b.iter(|| {
            let _: Vec<_> = slots.iter()
                .map(|ct| decrypt_slot(&ctx, &sk, black_box(ct)))
                .collect();
        })
    });
    group.finish();
}

criterion_group!(benches, bench_single_slot, bench_full_message);
criterion_main!(benches);
