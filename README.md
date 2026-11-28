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

- The base-width Gaussian sampler (keygen, base encryption, re-randomization) is a **constant-time CDT** discrete-Gaussian sampler (fixed-length branchless table scan). Its table is a constant: every nonzero magnitude has its probability in the discrete Gaussian rounded to the 63 bits of the uniform draw, and zero has the rest. The sampled distribution is within statistical distance 2^-60.96 of the discrete Gaussian and is the same on every platform; `tools/lean/cdt_table.lean` rebuilds the table in Lean 4, compares it with the Rust constant and bounds that distance.
- The flood-width sampler adds up eight samples of a second table, as the digits of a number in base 4, in **constant time and without floating point** (no rejection loop, no branch or memory address that depends on the data). No value is more likely than under the discrete Gaussian by a factor above 1 + 6.5e-8, and values far in the tail are less likely; a unit test computes the exact distribution of a sample from the table.
- The NTT / modular-arithmetic layer — including the secret-dependent decryption path (center-lift and final reduction) — is **branchless and division-free by construction**: control flow and memory-access patterns are independent of operand values.
- All RNG inputs are bounded by `rand_core::CryptoRng`, enforcing a CSPRNG at the type level (the bound cannot tell how the generator was seeded: seed it from the operating system). The samplers read the generator in a fixed, documented order, so a key pair derived from a seeded generator is the same on every platform. The library is `#![forbid(unsafe_code)]`.
- Secret keys, plaintexts and transient secrets (encryption and re-randomization randomness, decryption intermediates) are held in heap storage that is zeroized on drop and are computed there in place: these arrays are not moved or returned through the stack. The per-slot generators of the batch layer are zeroized too, and the part of the stack that a generator worked on is wiped after each sampling call. Keys and ciphertexts keep their limbs on the heap, so in release builds no operation needs more than 16 KiB of stack, and a full batch runs on worker stacks of 256 KiB. These measures are best effort: a compiler may still spill individual values, and CPU registers keep parts of the last arrays they processed.
- A TVLA/dudect-style fixed-vs-random timing-leakage harness (`cargo run --release --example leakage`) measures both samplers, key generation, encryption, re-randomization, decryption and secret-key parsing over their secret inputs (generator stream, secret key) as well as over plaintexts and ciphertexts, and the two causes of a failed decryption against each other, drawing the class of each measurement at random and reading the same memory whatever the class. On the machine of the benchmarks below, in portable and native builds (eight runs of 10,000 to 200,000 measurements per test), no sampler and no operation is flagged: the largest |t| is 3.0 against a threshold of 4.5, while a control with one data-dependent branch per coefficient gives |t| above 2700.
- A secret-taint harness (`examples/ctgrind.rs`, run under Valgrind) marks the secrets as undefined and checks that no conditional branch and no memory address depends on them, in the samplers and in every operation; only the final verdict of `SecretKey::from_bytes`, and of a decryption once its payload was accepted, may. A decryption that fails causes no report, whether the slot held no payload or the application rejected it. It runs in CI on x86-64 (baseline and x86-64-v3 builds), 32-bit x86 and aarch64. It cannot run under emulation, nor on a build that uses AVX-512: on 32-bit ARM, CI checks instead that the compiled crate calls no rounding, conversion or division routine of the run-time libraries, that it contains no floating-point arithmetic, comparison or conversion instruction, and that the flood sampler calls no library routine but `memset`. The 32-bit ARM and Android builds were otherwise checked by reading their assembly.
- A memory-residue scan (`examples/residue.rs`) runs the operations on known secrets, drops everything, and searches the heap and the stacks of all threads of the process for the keys and the output of the generators, for the secrets in the forms in which the library holds them (keys, errors, encryption randomness and plaintexts, as sampled, reduced or transformed) and for what it computes from them on the way to its results (products with public values, the encryption of zero that a message is added to or that a re-randomization adds, the noisy plaintext of a decryption). Freed memory is not reused and every operation runs on threads of its own, so nothing is overwritten before the search. It runs in CI on x86-64 (baseline and x86-64-v3 builds), aarch64, 32-bit ARM and 32-bit x86.

These properties hold for optimized builds without overflow checks (the release profile): debug builds branch on secret values and must not process real keys. They are properties of a build, not of the source alone, and on x86 of the CPU as well (the generator selects its code at run time; CI also runs the memory scan emulating a CPU without AVX2): repeat the taint check and the memory scan, like the timing assessment, on every target platform.

Branchless source plus one leakage run is not a substitute for microarchitectural verification: production deployments should repeat the timing-leakage assessment on their target platform (desktop and mobile) and an external audit is recommended before production use.

## Using it safely

The interface is key generation, encryption, re-randomization and decryption of a slot (a payload of 15,872 bytes), the same on a batch of slots in parallel (`encrypt_batch`, `rerandomize_batch`, `decrypt_batch`), and `to_bytes` / `from_bytes` on `SecretKey`, `PublicKey` and `CiphertextSlot`:

```rust
use pq_rerand::{decrypt_slot, encrypt_slot, keygen, rerandomize_slot, Flooding};

let (sk, pk) = keygen(&mut rng);
let ct = encrypt_slot(&mut rng, &pk, &payload, Flooding::Standard);
let ct = rerandomize_slot(&mut rng, &pk, &ct);    // anyone, with the public key alone
let message = decrypt_slot(&sk, &ct, aead_open);  // `None` for every failure
```

The RLWE layer is malleable by design (that is what makes public re-randomization possible) and only CPA-secure. An application must:

- **Authenticate the payload**: encrypt an AEAD ciphertext. Decryption takes the function that authenticates the payload (`aead_open` above) and returns nothing but what that function returns. It must return `None` unless the payload authenticates, in a time that does not depend on the payload, and must authenticate every byte it is given.
- **Hide decryption outcomes** from anyone who can submit or modify ciphertexts, including through errors, retries, replies, later writes and timing. Each observed success or failure of a modified ciphertext is a comparison on a secret value, and some 16,000 of them recover the secret key, also when honest parties re-randomize in between. The usual defense, re-encrypting the result to check it, is not available to a re-randomizable scheme. Decryption has a single failure, `None`, for a slot that did not decrypt to a payload and for a payload that was rejected, after the same work: it is that failure itself that must stay hidden.
- **Treat unauthenticated decryption output as key material**: a crafted ciphertext can decrypt to a function of the secret key. The payload only exists in zeroizing storage for the duration of the call, and the function that authenticates it may work in that buffer, as an AEAD that decrypts in place does: it must not log or return any of the payload unless it authenticates, and must leave no copy of it behind.
- **Count failed decryptions per key** and rotate the key after a small number of them.
- **Authenticate public keys**, including those only used to re-randomize: whoever replaces a public key can read what is encrypted under it and link what is re-randomized under it. `PublicKey::verify_against` checks a stored public key against its secret key.
- **Track the re-randomization budget** (below): the library cannot count re-randomizations.

Any party handling a ciphertext can also add an offset to it that stays below the decryption threshold: the plaintext still decrypts, later re-randomizations do not remove the offset and the secret-key holder can read it. Unlinkability holds against observers without the secret key; a relay colluding with the recipient can mark a ciphertext unless every hop proves that it added a well-formed encryption of zero. The crate documentation states these requirements in full.

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

The count is deliberately absent from the ciphertext, so nothing in the library tracks or enforces this budget: the application keeps an upper bound on it (for example age times maximum rate) and refreshes a ciphertext, by decrypting and encrypting it again, before the limit (`params::K_MAX_RERAND`). Flooding hides the number of re-randomizations from the secret-key holder only while that number is small next to `kappa_f` and the holder sees a single version of the ciphertext.

## Benchmark results

Measured on AMD Ryzen 7 260 (8C/16T, Zen 4, 3.8 GHz base / 5.0 GHz boost), 32 GiB DDR5, Linux 7.0, rustc 1.100.0, `--release` with `RUSTFLAGS="-C target-cpu=native"` (auto-vectorized; no hand-written SIMD intrinsics). Criterion medians.

**Fixed-frequency config** (reproducible: `cpupower frequency-set -g performance`, turbo boost disabled, all cores pinned at 3.8 GHz):

| Operation | 1 thread | 16 threads (per slot, B=1024) |
|-----------|---------:|------------------------------:|
| Key generation | 0.33 ms | - |
| Encrypt (base + folded flood) | 0.61 ms | 0.083 ms |
| Re-randomize | 0.30 ms | 0.044 ms |
| Decrypt | 0.17 ms | 0.023 ms |

Encryption and decryption include the encoding and the decoding of the payload. A portable build (no `target-cpu` flag) takes 0.43 ms, 0.79 ms, 0.45 ms and 0.23 ms on one thread at the same fixed frequency. With turbo boost enabled (OS-default scaling), the native build takes 0.49 ms to encrypt, 0.23 ms to re-randomize and 0.13 ms to decrypt (single thread).

Versus the original reference implementation at the same fixed frequency (encrypt 3.24 ms, re-randomize 1.62 ms, decrypt 1.47 ms), the optimized arithmetic (Barrett reduction, Shoup-multiply butterflies with precomputed bit-reversed twiddles, folded flooding, cached key NTT) is **5.3 to 8.6 times faster** single-threaded; rayon batch operations add a further factor of about 7 across 16 hardware threads. A full 15.5 MiB batch (1024 slots) re-randomizes in 0.045 s and decrypts in 0.023 s on 16 threads.

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
├── .github/workflows/ci.yml   # Tests (x86-64, aarch64, 32-bit ARM and x86), minimum Rust version, target builds, taint, residue and sampler-table checks
├── rust_implem/               # Rust implementation (~800 LOC of library code, plus tests)
│   ├── Cargo.toml
│   ├── src/
│   │   ├── lib.rs             # Crate root (public interface, forbid(unsafe_code), security notes)
│   │   ├── params.rs          # Scheme parameters and constants
│   │   ├── ntt.rs             # Negacyclic NTT (Barrett/Shoup, branchless, const tables)
│   │   ├── poly.rs            # NTT context, zeroizing storage and CRT-limb helpers
│   │   ├── sampling.rs        # Constant-time table samplers (base and flood widths)
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
│   │   ├── ctgrind.rs         # Secret-taint harness for Valgrind (x86, x86_64, aarch64)
│   │   ├── leakage.rs         # TVLA/dudect timing-leakage harness (x86_64 + aarch64)
│   │   ├── residue.rs         # Memory-residue scan of the heap and all stacks (Linux)
│   │   └── timings.rs         # Standalone timing harness (single + batch, 1/N threads)
│   └── tests/
│       └── correctness.rs     # Integration tests
├── tools/                     # Analysis scripts and figures
│   ├── make_figures.py        # Generate publication figures (requires numpy, matplotlib, scipy)
│   ├── security_estimate.py   # HE Standard v1.1 security cross-check (pure Python)
│   ├── lattice_estimate_sage.py  # Lattice-estimator script (requires SageMath)
│   └── lean/                  # Base sampler table, rebuilt and measured in Lean 4
│       ├── cdt_table.lean     # Builds the table, compares it with the Rust constant, bounds its distance from the discrete Gaussian
│       ├── cdt_table_out.txt  # Its output: the table and the distances of this and other constructions
│       └── lean-toolchain     # Version of Lean
├── LICENSE                    # MIT
└── README.md
```

## Building and testing

```bash
cd rust_implem
cargo build --release    # portable build
cargo test               # debug build: arithmetic overflow is checked
cargo test --release
cargo clippy --all-targets -- -D warnings
RUSTFLAGS="-C target-cpu=native" cargo bench   # benchmarks with host SIMD
cargo build --release --example ctgrind && valgrind target/release/examples/ctgrind   # secret-taint check
cargo run --release --example residue          # memory-residue scan (Linux)
```

Rust 1.88 or later is required. The crate is portable across desktop and mobile targets (pure safe Rust, no architecture-specific intrinsics in the library): the test suite, including its known-answer vectors, gives the same results on x86-64, aarch64, 32-bit ARM and 32-bit x86, and the library also builds for Android, iOS and WebAssembly targets.

## Running the analysis tools

```bash
# Generate publication figures (outputs to paper/figures/)
# Idealized model: it samples rounded Gaussians with NumPy, not with this crate.
pip install numpy matplotlib scipy
python3 tools/make_figures.py

# Quick security check (pure Python, no dependencies)
python3 tools/security_estimate.py

# Full lattice-estimator analysis (requires SageMath + lattice-estimator)
sage tools/lattice_estimate_sage.py

# Rebuild the table of the base sampler, compare it with the Rust constant and
# bound its distance from the discrete Gaussian (requires Lean 4 through elan)
(cd tools/lean && lean ./cdt_table.lean)
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

