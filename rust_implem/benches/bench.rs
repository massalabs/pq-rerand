use criterion::{criterion_group, criterion_main, Criterion, black_box};
use pq_rerand::params::*;
use pq_rerand::{decrypt_slot, encrypt_slot, keygen, rerandomize_slot, Flooding};
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::time::Duration;

/// Stands for the authenticated decryption of an application: it is handed the
/// payload and accepts it.
fn open(payload: &mut [u8; SLOT_BYTES]) -> Option<u8> {
    Some(black_box(payload)[0])
}

fn bench_single_slot(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(12345);
    let (sk, pk) = keygen(&mut rng);

    let payload: [u8; SLOT_BYTES] = std::array::from_fn(|i| i as u8);

    let ct = encrypt_slot(&mut rng, &pk, &payload, Flooding::Standard);

    c.bench_function("keygen", |b| {
        b.iter(|| keygen(&mut rng))
    });

    c.bench_function("encrypt_slot (base+smudge)", |b| {
        b.iter(|| {
            encrypt_slot(&mut rng, &pk, black_box(&payload), Flooding::Standard)
        })
    });

    c.bench_function("rerandomize_slot", |b| {
        b.iter(|| {
            rerandomize_slot(&mut rng, &pk, black_box(&ct))
        })
    });

    c.bench_function("decrypt_slot", |b| {
        b.iter(|| {
            decrypt_slot(&sk, black_box(&ct), open)
        })
    });
}

fn bench_full_message(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(99999);
    let (sk, pk) = keygen(&mut rng);

    let payload: [u8; SLOT_BYTES] = std::array::from_fn(|i| i as u8);
    let slots: Vec<_> = (0..NUM_SLOTS)
        .map(|_| encrypt_slot(&mut rng, &pk, &payload, Flooding::Standard))
        .collect();

    // Heavy 15.5 MiB batch ops (single-threaded); fewer samples to keep wall-clock
    // reasonable. Multi-threaded (rayon) batch timings are produced by the
    // `timings` example.
    let mut group = c.benchmark_group("full_15.5MiB_1024_slots");
    group.sample_size(10).measurement_time(Duration::from_secs(8));

    group.bench_function("rerandomize (single-threaded)", |b| {
        b.iter(|| {
            let _: Vec<_> = slots.iter()
                .map(|ct| rerandomize_slot(&mut rng, &pk, black_box(ct)))
                .collect();
        })
    });
    group.bench_function("decrypt (single-threaded)", |b| {
        b.iter(|| {
            let _: Vec<_> = slots.iter()
                .map(|ct| decrypt_slot(&sk, black_box(ct), open))
                .collect();
        })
    });
    group.finish();
}

criterion_group!(benches, bench_single_slot, bench_full_message);
criterion_main!(benches);
