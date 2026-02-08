use criterion::{criterion_group, criterion_main, Criterion, black_box};
use pq_rerand::params::*;
use pq_rerand::poly::NttContext;
use pq_rerand::keygen::keygen;
use pq_rerand::encrypt::encrypt_slot;
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::rerandomize::rerandomize_slot;
use rand::SeedableRng;
use rand::rngs::StdRng;

fn bench_single_slot(c: &mut Criterion) {
    let ctx = NttContext::new();
    let mut rng = StdRng::seed_from_u64(12345);
    let (sk, pk) = keygen(&mut rng, &ctx);

    let mut message = [0u32; N];
    for i in 0..N {
        message[i] = (i as u32) % (1 << BITS_PER_COEFF);
    }

    let ct = encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD);

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

    // Encrypt 1024 slots
    let mut message = [0u32; N];
    for i in 0..N { message[i] = (i as u32) % (1 << BITS_PER_COEFF); }
    let slots: Vec<_> = (0..NUM_SLOTS)
        .map(|_| encrypt_slot(&mut rng, &ctx, &pk, &message, SIGMA_FLOOD))
        .collect();

    c.bench_function("rerandomize_full_15.5MiB (1024 slots, single-threaded)", |b| {
        b.iter(|| {
            let _: Vec<_> = slots.iter()
                .map(|ct| rerandomize_slot(&mut rng, &ctx, &pk, ct))
                .collect();
        })
    });

    c.bench_function("decrypt_full_15.5MiB (1024 slots, single-threaded)", |b| {
        b.iter(|| {
            let _: Vec<_> = slots.iter()
                .map(|ct| decrypt_slot(&ctx, &sk, ct))
                .collect();
        })
    });
}

criterion_group!(benches, bench_single_slot, bench_full_message);
criterion_main!(benches);
