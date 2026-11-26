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
//!   ([`sampling::add_flood`]).
//! * The NTT / modular-arithmetic layer — including the secret-dependent
//!   decryption path (center-lift and final reduction) — is branchless and
//!   division-free by construction: control flow and memory-access patterns are
//!   independent of operand values. Its selections go through
//!   [`core::hint::select_unpredictable`], which asks the compiler not to turn
//!   them into conditional branches. That is a request, not a guarantee: the
//!   compiled code is checked as described below.
//! * All RNG inputs are bounded by [`rand_core::CryptoRng`], so a CSPRNG is enforced
//!   at the type level. The bound says nothing about how the generator was
//!   seeded: seed it from the operating system. Its state gives away the
//!   randomness it produced, so prefer a generator that zeroizes its state on
//!   drop (the batch layer does this for its per-slot generators).
//! * Secret keys, plaintexts and transient secrets (encryption /
//!   re-randomization randomness, decryption noise) are held in heap storage that
//!   is zeroized on drop and are computed there in place: these arrays are not
//!   moved or returned through the stack. The part of the stack that a generator
//!   worked on is wiped after each sampling call. Both measures are best effort:
//!   a compiler may still spill individual values to the stack, and parts of
//!   these arrays remain in CPU registers after an operation, which the operating
//!   system saves to the stack of the thread when it delivers a signal.
//!
//! These properties hold for optimized builds without overflow checks, such as
//! the release profile. Debug builds, and builds with `overflow-checks` enabled,
//! branch on secret values: they must not process real keys.
//!
//! A TVLA-style fixed-vs-random timing-leakage harness (`examples/leakage.rs`)
//! measures the samplers and the main operations over their secret inputs, and a
//! secret-taint harness (`examples/ctgrind.rs`, run under Valgrind) checks that no
//! branch or memory address depends on a secret, except for the final verdict of
//! [`keygen::SecretKey::from_bytes`] and of [`encoding::decode`]. A memory scan
//! (`examples/residue.rs`) searches the heap and the stacks of all threads for
//! secrets after their use.
//!
//! Note that branchless source plus one leakage run is not a substitute for
//! microarchitectural verification: compilers may in principle reintroduce
//! branches, so deployments should repeat these assessments on the
//! target platform (desktop and mobile), and an external audit is recommended
//! before production use.
//!
//! The flood sampler computes in floating point: it needs a target with hardware
//! floating point, and relies on the division and square root of the CPU taking
//! a time that does not depend on their operands.
//!
//! # Requirements on applications
//!
//! The RLWE layer is CPA-secure and intentionally malleable (that is what makes
//! public re-randomization possible). It is a transport for a payload that the
//! application authenticates, and it only protects that payload and its own
//! secret key if the application meets the following requirements.
//!
//! * **Authenticate the payload.** Encrypt an AEAD ciphertext as the payload, and
//!   verify it before using or releasing anything that was decrypted.
//! * **Hide decryption outcomes.** Whether a ciphertext decrypted and verified
//!   must not be observable by anyone who can submit or modify ciphertexts: not
//!   through an error, a reply, a retry, a fallback to another replica, a later
//!   write, or timing. Adding a chosen offset to a ciphertext leaves the plaintext
//!   intact up to a threshold that depends on the secret noise, so each observed
//!   outcome is a comparison on a secret value: some 16,000 of them, about four per
//!   coefficient, recover the secret key, also when honest parties re-randomize
//!   the ciphertext in between. A check that re-encrypts the result, the usual
//!   defense, is not available to a re-randomizable scheme. Give callers a single
//!   generic failure, and treat [`encoding::decode`] returning `None` exactly
//!   like a payload that fails to authenticate.
//! * **Treat unauthenticated decryption output as key material.** A crafted
//!   ciphertext can make [`decrypt::decrypt_slot`] return a function of the secret
//!   key, such as its negation. Keep the output in the zeroizing containers it
//!   comes in, never log or return it, and drop it when authentication fails.
//! * **Count failures and rotate keys.** Count failed decryptions per key, and
//!   replace a key (re-encrypting what it protects) after a small number of them.
//! * **Authenticate public keys**, including those only used to re-randomize
//!   (see [`keygen::PublicKey`]).
//! * **Track the re-randomization budget**
//!   (see [`rerandomize::rerandomize_slot`]).
//!
//! One limit remains whatever the application does. Any party that handles a
//! ciphertext can add an offset to it that stays below the decryption threshold:
//! the plaintext still decrypts and authenticates, later re-randomizations do not
//! remove the offset, and the secret-key holder can read it in the noise. A relay
//! and a recipient that collude can use it to mark a ciphertext and recognize it.
//! Unlinkability holds against observers without the secret key. Where relays
//! must not be able to mark ciphertexts, each hop has to prove that it added a
//! well-formed encryption of zero.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod params;
mod ntt;
pub mod poly;
pub mod sampling;
pub mod encoding;
pub mod keygen;
pub mod encrypt;
pub mod rerandomize;
pub mod decrypt;
pub mod serialize;
pub mod batch;
