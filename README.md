# pq-rerand

An efficient, large-block, post-quantum publicly re-randomizable encryption scheme based on Ring-LWE.

This repository is the reference implementation of the scheme published in the *Journal of Cryptographic Engineering*:

> Vodenicarevic, D., Fleiser, A., Seznec, P. et al. **Two-limb CRT Ring-LWE encryption with exact decryption and public re-randomization.** *J Cryptogr Eng* 16, 13 (2026). [https://link.springer.com/article/10.1007/s13389-026-00401-2](https://link.springer.com/article/10.1007/s13389-026-00401-2)

**Public re-randomization** allows anyone holding the public key to transform a ciphertext into a fresh-looking encryption of the same plaintext — without learning anything about the plaintext and without increasing the ciphertext size. The original and re-randomized ciphertexts are computationally unlinkable under the Decision Ring-LWE assumption.

## Why this scheme?

Classical ElGamal supports unlimited re-randomization with statistical unlinkability, but is broken by quantum computers. Existing lattice-based (Ring-LWE) encryption supports re-randomization in principle (by adding an encryption of zero), but practical implementations typically rely on multi-prime RNS stacks, approximate decryption with rounding, and carry significant complexity.

This construction is engineered around a **two-limb CRT modulus** `q = t · q₂` where both primes fit in 32-bit words. The key insight is that embedding the plaintext as `Δ·M = q₂·M` makes the message vanish modulo `q₂`, so the `q₂`-limb carries **only** the decryption noise. This enables:

- **Exact noise extraction and message recovery** without CRT recombination and without rounding
- **Public re-randomization** via `ct' = ct + Enc(pk, 0)` with narrow Gaussian noise
- **Computational unlinkability** (IND\$) under the Decision Ring-LWE assumption
- **Large plaintext blocks**: ~15.5 KiB per ciphertext slot, ~15.5 MiB per batch of 1024 slots
- **Simple 32-bit NTT** implementation (both CRT primes < 2³², all intermediate products fit in u64)

## Security & side-channel posture

**This implementation has not yet been externally audited.** It is written to production-quality standards and is intended to be audit-ready:

- The base-width Gaussian sampler (keygen, base encryption, re-randomization) is a **constant-time CDT** discrete-Gaussian sampler (fixed-length branchless table scan).
- The flood-width sampler is a **branchless constant-time inverse-CDF** sampler (no rejection loop, no data-dependent branch or table indexing).
- The NTT / modular-arithmetic layer — including the secret-dependent decryption path (center-lift and final reduction) — is **branchless and division-free by construction**: control flow and memory-access patterns are independent of operand values.
- All RNG inputs are bounded by `rand::CryptoRng`, enforcing a CSPRNG at the type level. Secret keys are zeroized on drop; transient secrets are zeroized best-effort. The library is `#![forbid(unsafe_code)]`.
- A TVLA/dudect-style fixed-vs-random timing-leakage harness (`cargo run --release --example leakage`) reports no measurable leakage (Welch |t| < 4.5) on the encryption, re-randomization and decryption paths, while flagging a deliberately variable-time control (|t| > 3000). The harness also includes a hardware diagnostic: on CPUs with data-dependent bit-toggle timing (observed on Zen 4), *all-zero* vs random data content is distinguishable at the ~0.02% level even for a pure load loop with no crypto arithmetic; AEAD-wrapped payloads (the intended deployment) are uniformly random, so such degenerate inputs cannot occur.

Branchless source plus one leakage run is not a substitute for microarchitectural verification: production deployments should repeat the timing-leakage assessment on their target platform (desktop and mobile) and an external audit is recommended before production use.

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

Measured on AMD Ryzen 7 260 (8C/16T, Zen 4, 3.8 GHz base / 5.0 GHz boost), 32 GiB DDR5, Linux 7.0, rustc 1.94.1, `--release` with `RUSTFLAGS="-C target-cpu=native"` (auto-vectorized; no hand-written SIMD intrinsics). Criterion medians.

**Fixed-frequency config** (reproducible: `cpupower frequency-set -g performance`, turbo boost disabled, all cores pinned at 3.8 GHz):

| Operation | 1 thread | 16 threads (per slot, B=1024) |
|-----------|---------:|------------------------------:|
| Key generation | 0.92 ms | — |
| Encrypt (base + folded flood) | 0.80 ms | 0.096 ms |
| Re-randomize | 0.51 ms | 0.070 ms |
| Decrypt | 0.21 ms | 0.025 ms |

With turbo boost enabled (OS-default scaling): encrypt 0.63 ms, re-randomize 0.41 ms, decrypt 0.16 ms (single thread).

Versus the original reference implementation at the same fixed frequency (encrypt 3.24 ms, re-randomize 1.62 ms, decrypt 1.47 ms), the optimized arithmetic (Barrett reduction, Shoup-multiply butterflies with precomputed bit-reversed twiddles, folded flooding, cached key NTT) is **3.1–7.1× faster** single-threaded; rayon batch operations add a further ~7–8× across 16 hardware threads. A full 15.5 MiB batch (1024 slots) re-randomizes in 0.072 s and decrypts in 0.026 s on 16 threads.

Reproduce with:

```bash
sudo cpupower frequency-set -g performance
echo 0 | sudo tee /sys/devices/system/cpu/cpufreq/boost
RUSTFLAGS="-C target-cpu=native" cargo bench               # criterion
RUSTFLAGS="-C target-cpu=native" cargo run --release --example timings
RUSTFLAGS="-C target-cpu=native" cargo run --release --example leakage
```

## Repository structure

```
pq-rerand/
├── rust_implem/               # Rust implementation (~1550 LOC)
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs             # Crate root (forbid(unsafe_code), security notes)
│   │   ├── params.rs          # Scheme parameters and constants
│   │   ├── ntt.rs             # Negacyclic NTT (Barrett/Shoup, branchless, const tables)
│   │   ├── poly.rs            # CRT polynomial types and ring arithmetic
│   │   ├── sampling.rs        # Constant-time CDT + inverse-CDF Gaussian samplers
│   │   ├── encoding.rs        # 31-bit plaintext encoding (bytes ↔ coefficients)
│   │   ├── keygen.rs          # Key generation (zeroized secret key, cached NTT)
│   │   ├── encrypt.rs         # Encryption with folded noise flooding
│   │   ├── rerandomize.rs     # Public re-randomization
│   │   ├── decrypt.rs         # Branchless decryption via noise-limb trick
│   │   ├── serialize.rs       # Validated ciphertext (de)serialization
│   │   └── batch.rs           # Multi-threaded (rayon) batch operations
│   ├── benches/
│   │   └── bench.rs           # Criterion benchmarks
│   ├── examples/
│   │   ├── leakage.rs         # TVLA/dudect timing-leakage harness (x86_64 + aarch64)
│   │   └── timings.rs         # Standalone timing harness (single + batch, 1/N threads)
│   └── tests/
│       └── correctness.rs     # Integration tests
├── tools/                     # Python scripts for analysis and figures
│   ├── make_figures.py        # Generate publication figures (requires numpy, matplotlib, scipy)
│   ├── security_estimate.py   # HE Standard v1.1 security cross-check (pure Python)
│   └── lattice_estimate_sage.py  # Lattice-estimator script (requires SageMath)
├── LICENSE                    # MIT
└── README.md
```

## Building and testing

```bash
cd rust_implem
cargo build --release    # portable build
cargo test --release
cargo clippy --all-targets -- -D warnings
RUSTFLAGS="-C target-cpu=native" cargo bench   # benchmarks with host SIMD
```

The crate is portable across desktop and mobile targets (pure safe Rust, no architecture-specific intrinsics in the library; verified to build for `aarch64-unknown-linux-gnu`).

## Running the analysis tools

```bash
# Generate publication figures (outputs to paper/figures/)
pip install numpy matplotlib scipy
python3 tools/make_figures.py

# Quick security check (pure Python, no dependencies)
python3 tools/security_estimate.py

# Full lattice-estimator analysis (requires SageMath + lattice-estimator)
sage tools/lattice_estimate_sage.py
```

## How it works

1. **KeyGen**: sample secret `s` and error `e` from a narrow Gaussian; compute `b = a·s + e` in both CRT limbs.
2. **Encrypt**: sample randomness `(r, e₁, e₂)`, compute `c₁ = a·r + e₁` and `c₀ = b·r + e₂ + Δ_t·M` (t-limb) / `c₀ = b·r + e₂` (q₂-limb). Optionally add a one-time flooding encryption of zero.
3. **ReRand**: sample fresh `(r', e₁', e₂')` and add `Enc(pk, 0)` to the ciphertext.
4. **Decrypt**: compute `v = c₀ - c₁·s` in each limb. The q₂-limb `v_{q₂}` is pure noise (since `Δ·M ≡ 0 mod q₂`). Center-lift to recover the signed noise `ν`, then recover `M = (v_t - ν) · Δ_t⁻¹ mod t`.

## Citing this work

The scheme, its security analysis, and the parameter selection are described in the peer-reviewed article, published in the *Journal of Cryptographic Engineering* (Springer): [https://link.springer.com/article/10.1007/s13389-026-00401-2](https://link.springer.com/article/10.1007/s13389-026-00401-2) (DOI: [10.1007/s13389-026-00401-2](https://doi.org/10.1007/s13389-026-00401-2)).

If you use this work, please cite:

```bibtex
@article{vodenicarevic2026pqrerand,
  author  = {Vodenicarevic, Damir and Fleiser, Andrei and Seznec, Pierre and
             Mayen Naranjo, Karen and Foucher, Lucas and Besan{\c{c}}on, L{\'e}o and
             Alabarbe, Thybault and Morcillo, Jean-Fran{\c{c}}ois and
             Reynes, Benjamin and Urvoy, Lilian},
  title   = {Two-limb {CRT} {Ring-LWE} encryption with exact decryption and
             public re-randomization},
  journal = {Journal of Cryptographic Engineering},
  volume  = {16},
  pages   = {13},
  year    = {2026},
  doi     = {10.1007/s13389-026-00401-2},
  url     = {https://link.springer.com/article/10.1007/s13389-026-00401-2}
}
```

## License

This repository is released under the [MIT License](LICENSE).

