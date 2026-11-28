//! Post-quantum re-randomizable encryption (BFV/LPR + CRT noise-limb trick).
//!
//! The interface is [`fn@keygen`], [`encrypt_slot`], [`rerandomize_slot`] and
//! [`decrypt_slot`], their versions for a batch of slots ([`encrypt_batch`],
//! [`rerandomize_batch`], [`decrypt_batch`]) and the serialized forms of
//! [`SecretKey`], [`PublicKey`] and [`CiphertextSlot`] (`to_bytes` and
//! `from_bytes`). A slot holds a payload of [`params::SLOT_BYTES`] bytes.
//!
//! ```
//! use pq_rerand::{decrypt_slot, encrypt_slot, keygen, rerandomize_slot, Flooding};
//! # use pq_rerand::params::SLOT_BYTES;
//! # // Stand-ins for the AEAD of the application, which this crate does not provide.
//! # fn aead_seal(message: &[u8]) -> [u8; SLOT_BYTES] {
//! #     std::array::from_fn(|i| message.get(i).copied().unwrap_or(0))
//! # }
//! # fn aead_open(payload: &mut [u8; SLOT_BYTES]) -> Option<Vec<u8>> {
//! #     Some(payload[..5].to_vec())
//! # }
//!
//! let mut rng = rand::rng();
//! let (sk, pk) = keygen(&mut rng);
//!
//! // The payload is a ciphertext of the AEAD of the application.
//! let payload = aead_seal(b"hello");
//! let ct = encrypt_slot(&mut rng, &pk, &payload, Flooding::Standard);
//!
//! // Anyone can re-randomize it, with the public key alone.
//! let ct = rerandomize_slot(&mut rng, &pk, &ct);
//!
//! // Decryption passes the payload to the authenticated decryption of the
//! // application and returns what that returns: `None` for every failure.
//! let message = decrypt_slot(&sk, &ct, aead_open);
//! assert_eq!(message.as_deref(), Some(&b"hello"[..]));
//! ```
//!
//! # Security and side-channel posture
//!
//! This crate is written to be production-quality but has **not yet been
//! externally audited**. Its side-channel posture:
//!
//! * The base-width Gaussian sampler (keygen / base encryption /
//!   re-randomization) is a constant-time CDT sampler.
//! * The flood-width sampler adds up eight samples of a second table, in constant
//!   time like the base sampler and without floating point.
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
//! [`SecretKey::from_bytes`] and of a decryption whose payload was accepted. A
//! memory scan (`examples/residue.rs`) searches the heap and the stacks of all
//! threads for secrets after their use.
//!
//! Note that branchless source plus one leakage run is not a substitute for
//! microarchitectural verification: compilers may in principle reintroduce
//! branches, so deployments should repeat these assessments on the
//! target platform (desktop and mobile), and an external audit is recommended
//! before production use.
//!
//! # Requirements on applications
//!
//! The RLWE layer is CPA-secure and intentionally malleable (that is what makes
//! public re-randomization possible). It is a transport for a payload that the
//! application authenticates, and it only protects that payload and its own
//! secret key if the application meets the following requirements.
//!
//! * **Authenticate the payload.** Encrypt an AEAD ciphertext as the payload.
//!   Decryption ([`decrypt_slot`], [`decrypt_batch`]) takes the function that
//!   authenticates it, `open`, and returns nothing but what `open` returns. `open`
//!   must return `None` unless the payload authenticates, in a time that does not
//!   depend on the payload, and must authenticate every byte that it is given: a
//!   byte that it ignores, such as padding added after the AEAD, can be changed
//!   without failing it.
//! * **Hide decryption outcomes.** Whether a ciphertext decrypted and verified
//!   must not be observable by anyone who can submit or modify ciphertexts: not
//!   through an error, a reply, a retry, a fallback to another replica, a later
//!   write, or timing. Adding a chosen offset to a ciphertext leaves the plaintext
//!   intact up to a threshold that depends on the secret noise, so each observed
//!   outcome is a comparison on a secret value: some 16,000 of them, about four per
//!   coefficient, recover the secret key, also when honest parties re-randomize
//!   the ciphertext in between. A check that re-encrypts the result, the usual
//!   defense, is not available to a re-randomizable scheme. Decryption has a
//!   single failure, `None`, for a slot that did not decrypt to a payload and for
//!   a payload that `open` rejects, after the same work. It is that `None`
//!   itself, however it shows, that must stay hidden.
//! * **Treat unauthenticated decryption output as key material.** A crafted
//!   ciphertext can make the payload that `open` receives a function of the secret
//!   key, such as a constant minus the key. That payload only exists in zeroizing
//!   storage for the duration of the call, and `open` may work in that buffer, as
//!   an AEAD that decrypts in place does: it must not log or return any of the
//!   payload unless it authenticates, and must leave no copy of it behind. An AEAD
//!   call that returns the plaintext in a buffer of its own may first copy the
//!   payload there, and free that copy without clearing it when the payload does
//!   not authenticate.
//! * **Count failures and rotate keys.** Count failed decryptions per key, and
//!   replace a key (re-encrypting what it protects) after a small number of them.
//! * **Authenticate public keys**, including those only used to re-randomize
//!   (see [`PublicKey`]).
//! * **Track the re-randomization budget** (see [`rerandomize_slot`]).
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
mod poly;
mod sampling;
mod encoding;
mod keygen;
mod encrypt;
mod rerandomize;
mod decrypt;
mod serialize;
mod batch;

pub use batch::{decrypt_batch, encrypt_batch, rerandomize_batch};
pub use decrypt::decrypt_slot;
pub use encrypt::{encrypt_slot, CiphertextSlot, Flooding};
pub use keygen::{keygen, PublicKey, SecretKey};
pub use rerandomize::rerandomize_slot;
