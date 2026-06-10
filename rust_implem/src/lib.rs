//! Post-quantum re-randomizable encryption (BFV/LPR + CRT noise-limb trick).
//!
//! # Security and side-channel posture
//!
//! This crate is written to be production-quality but has **not yet been
//! externally audited**. Its side-channel posture:
//!
//! * The base-width Gaussian sampler (keygen / base encryption /
//!   re-randomization) is a constant-time CDT sampler
//!   ([`sampling::sample_gaussian_cdt`]).
//! * The flood-width sampler is a branchless constant-time inverse-CDF sampler
//!   ([`sampling::sample_gaussian`]).
//! * The NTT / modular-arithmetic layer — including the secret-dependent
//!   decryption path (center-lift and final reduction) — is branchless and
//!   division-free by construction: control flow and memory-access patterns are
//!   independent of operand values (see [`ntt::csub`], [`ntt::barrett_reduce`]).
//! * All RNG inputs are bounded by [`rand::CryptoRng`], so a CSPRNG is enforced
//!   at the type level.
//! * Secret keys are zeroized on drop, and transient secrets (encryption /
//!   re-randomization randomness, decryption noise) are zeroized best-effort.
//!
//! A TVLA-style fixed-vs-random timing-leakage harness (`examples/leakage.rs`)
//! reports no measurable leakage (Welch |t| < 4.5) on the encryption,
//! re-randomization or decryption paths, while a branchy control reducer in the
//! same harness is flagged at |t| > 3000 (sensitivity check).
//!
//! Note that branchless source plus one leakage run is not a substitute for
//! microarchitectural verification: compilers may in principle reintroduce
//! branches, so deployments should repeat the timing-leakage assessment on the
//! target platform (desktop and mobile), and an external audit is recommended
//! before production use.
//!
//! # Authenticity
//!
//! The RLWE layer is CPA-secure and intentionally malleable (that is what makes
//! public re-randomization possible). Applications must provide authenticity by
//! encrypting an AEAD ciphertext as the payload and must surface a single
//! generic failure signal at decryption (see the accompanying paper).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod params;
pub mod ntt;
pub mod poly;
pub mod sampling;
pub mod encoding;
pub mod keygen;
pub mod encrypt;
pub mod rerandomize;
pub mod decrypt;
pub mod serialize;
pub mod batch;
