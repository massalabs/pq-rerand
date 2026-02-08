# pq-rerand

An efficient, large-block, post-quantum publicly re-randomizable encryption scheme based on Ring-LWE.

**Public re-randomization** allows anyone holding the public key to transform a ciphertext into a fresh-looking encryption of the same plaintext — without learning anything about the plaintext and without increasing the ciphertext size. The original and re-randomized ciphertexts are computationally unlinkable under the Decision Ring-LWE assumption.

## Why this scheme?

Classical ElGamal supports unlimited re-randomization with statistical unlinkability, but is broken by quantum computers. Existing lattice-based (Ring-LWE) encryption supports re-randomization in principle (by adding an encryption of zero), but practical implementations typically rely on multi-prime RNS stacks, approximate decryption with rounding, and carry significant complexity.

This construction is engineered around a **two-limb CRT modulus** `q = t · q₂` where both primes fit in 32-bit words. The key insight is that embedding the plaintext as `Δ·M = q₂·M` makes the message vanish modulo `q₂`, so the `q₂`-limb carries **only** the decryption noise. This enables:

- **Exact noise extraction and message recovery** without CRT recombination and without rounding
- **Public re-randomization** via `ct' = ct + Enc(pk, 0)` with narrow Gaussian noise
- **Computational unlinkability** (IND\$) under the Decision Ring-LWE assumption
- **Large plaintext blocks**: ~15.5 KiB per ciphertext slot, ~15.5 MiB per batch of 1024 slots
- **Simple 32-bit NTT** implementation (both CRT primes < 2³², all intermediate products fit in u64)

## Warning

**This implementation is a research prototype for academic evaluation only.**

It has NOT been audited, is NOT constant-time, does NOT protect against side-channel attacks, and MUST NOT be used in any production system. The Gaussian sampler uses floating-point Box-Muller, which is biased in the tails and leaks timing information.

## Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| n | 4,096 | Ring dimension (power of 2) |
| t | 2,147,565,569 | Plaintext modulus, prime ≈ 2³¹, t ≡ 1 (mod 2n) |
| q₂ | 4,294,828,033 | Noise-limb modulus, prime ≈ 2³², q₂ ≡ 1 (mod 2n) |
| q = t·q₂ | ≈ 2⁶³ | Combined modulus (never materialized on the hot path) |
| σ | 3.2 | Narrow discrete Gaussian width |
| σ_flood | ≈ 93,733 | Optional flooding Gaussian width (σ·√κ_f) |
| Plaintext/slot | 15,872 bytes (≈ 15.5 KiB) | 4096 coefficients × 31 bits |
| Ciphertext/slot | 65,536 bytes (64 KiB) | 4 × 4096 × 4 bytes (u32 arrays) |
| Batch (1024 slots) | ~15.5 MiB plaintext / 64 MiB ciphertext | |

## Security

- **HE Security Standard v1.1**: for n = 4096 and σ ≈ 3.2, the 128-bit classical threshold is log₂q ≤ 109 and the 128-bit quantum threshold is log₂q ≤ 101. Our log₂q = 63 uses only 58% of the classical budget and 62% of the quantum budget.
- **Lattice estimator** (core-SVP sieving model): best attack requires BKZ block size β ≈ 1349, costing ≥ 2³⁹⁴ classical / ≥ 2³⁵⁷ quantum operations.

## Re-randomization budget

With a per-batch failure probability of 2⁻¹⁰⁶ and B = 1024 slots:

| Configuration | k_max | Real-time equivalent |
|---------------|-------|----------------------|
| No flooding (κ_f = 0) | ≈ 1.5 × 10¹⁰ | ~476 years at 1 rerand/s |
| With flooding (κ_f = 8.58 × 10⁸) | ≈ 1.4 × 10¹⁰ additional | ~448 years at 1 rerand/s |

## Benchmark results

Single-threaded, `--release` mode, no SIMD intrinsics. Measured on AMD Ryzen 7 (Zen 4, 3.8 GHz base / 5.0 GHz boost), 32 GiB DDR5, Linux 6.17, rustc 1.91.1.

| Operation | Per ciphertext slot | Full batch (1024 slots) |
|-----------|--------------------:|------------------------:|
| Encrypt (base + flood) | 2.6 ms | 2.7 s |
| Re-randomize | 1.3 ms | 1.4 s |
| Decrypt | 1.2 ms | 1.2 s |

Batch operations are embarrassingly parallel across slots (not yet parallelized in this reference implementation).

## Repository structure

```
pq-rerand/
├── rust_implem/          # Rust reference implementation (~1050 LOC)
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs        # Crate root
│   │   ├── params.rs     # Scheme parameters and constants
│   │   ├── ntt.rs        # Number Theoretic Transform (negacyclic)
│   │   ├── poly.rs       # CRT polynomial types and ring arithmetic
│   │   ├── sampling.rs   # Gaussian and uniform sampling
│   │   ├── encoding.rs   # 31-bit plaintext encoding (bytes ↔ coefficients)
│   │   ├── keygen.rs     # Key generation
│   │   ├── encrypt.rs    # Encryption with optional noise flooding
│   │   ├── rerandomize.rs# Public re-randomization
│   │   ├── decrypt.rs    # Decryption via noise-limb trick
│   │   └── serialize.rs  # Ciphertext serialization
│   ├── benches/
│   │   └── bench.rs      # Criterion benchmarks
│   └── tests/
│       └── correctness.rs# Integration tests (22 tests total)
├── paper/                # LaTeX source of the accompanying paper
├── tools/                # Python scripts (noise simulation, figures, lattice estimator)
└── README.md
```

## Building and testing

```bash
cd rust_implem
cargo build --release
cargo test --release
cargo bench
```

## How it works

1. **KeyGen**: sample secret `s` and error `e` from a narrow Gaussian; compute `b = a·s + e` in both CRT limbs.
2. **Encrypt**: sample randomness `(r, e₁, e₂)`, compute `c₁ = a·r + e₁` and `c₀ = b·r + e₂ + Δ_t·M` (t-limb) / `c₀ = b·r + e₂` (q₂-limb). Optionally add a one-time flooding encryption of zero.
3. **ReRand**: sample fresh `(r', e₁', e₂')` and add `Enc(pk, 0)` to the ciphertext.
4. **Decrypt**: compute `v = c₀ - c₁·s` in each limb. The q₂-limb `v_{q₂}` is pure noise (since `Δ·M ≡ 0 mod q₂`). Center-lift to recover the signed noise `ν`, then recover `M = (v_t - ν) · Δ_t⁻¹ mod t`.
