//! Secret-taint harness, to be run under Valgrind (the "ctgrind" method).
//!
//! Memcheck tracks, bit by bit, which data is undefined, and reports every
//! conditional jump and every memory address that depends on undefined data. This
//! program marks the secrets as undefined (every 64-bit word of the generator,
//! which the samplers turn into secrets, and the payload) and runs every
//! operation: a report is a control-flow or memory-access decision that depends
//! on a secret, exactly what constant-time code must not contain. Public keys and
//! ciphertexts, which derive from secrets, stay marked although an observer sees
//! them anyway: that only makes the check stricter. The 32-bit words of the
//! generator stay defined: they only produce the public polynomial `a`.
//!
//! Each operation is checked against the number of reports it may cause: none,
//! or one to four jumps for the two that end in a verdict on secret data
//! (`SecretKey::from_bytes` and `decode` say whether their input is well formed).
//! A control that branches on a secret must be reported, which shows that the
//! tool sees the secrets on this target.
//!
//! Run with:
//!   cargo build --release --example ctgrind
//!   valgrind target/release/examples/ctgrind
//! on every target the crate is deployed to (x86, x86_64 and aarch64 are supported
//! here; a build that uses AVX-512 does not run under Valgrind), with the flags of
//! the deployed build. The exit status is zero only if every check passed under
//! Valgrind. This complements the timing harness (`examples/leakage.rs`): it sees
//! every jump and address, but not instructions whose latency depends on their
//! operands.

use pq_rerand::batch::{decrypt_batch, encrypt_batch, rerandomize_batch};
use pq_rerand::decrypt::decrypt_slot;
use pq_rerand::encoding::{decode, encode};
use pq_rerand::encrypt::{encrypt_slot, Flooding};
use pq_rerand::keygen::{keygen, SecretKey};
use pq_rerand::params::{N, SLOT_BYTES};
use pq_rerand::poly::NTT_CONTEXT;
use pq_rerand::rerandomize::rerandomize_slot;
use pq_rerand::sampling::{add_flood, sample_gaussian_cdt};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng, TryCryptoRng, TryRng};
use std::convert::Infallible;
use std::hint::black_box;

/// A Valgrind client request: a no-op that returns 0 outside Valgrind.
#[cfg(target_arch = "x86_64")]
fn client_request(request: usize, arg1: usize, arg2: usize) -> usize {
    let args = [request, arg1, arg2, 0, 0, 0];
    let mut result = 0usize;
    unsafe {
        core::arch::asm!(
            "rol rdi, 3", "rol rdi, 13", "rol rdi, 61", "rol rdi, 51", "xchg rbx, rbx",
            in("rax") args.as_ptr(), inout("rdx") result, options(nostack),
        );
    }
    result
}

#[cfg(target_arch = "x86")]
fn client_request(request: usize, arg1: usize, arg2: usize) -> usize {
    let args = [request, arg1, arg2, 0, 0, 0];
    let mut result = 0usize;
    unsafe {
        core::arch::asm!(
            "rol edi, 3", "rol edi, 13", "rol edi, 29", "rol edi, 19", "xchg ebx, ebx",
            in("eax") args.as_ptr(), inout("edx") result, options(nostack),
        );
    }
    result
}

#[cfg(target_arch = "aarch64")]
fn client_request(request: usize, arg1: usize, arg2: usize) -> usize {
    let args = [request, arg1, arg2, 0, 0, 0];
    let mut result = 0usize;
    unsafe {
        core::arch::asm!(
            "ror x12, x12, #3", "ror x12, x12, #13", "ror x12, x12, #51", "ror x12, x12, #61",
            "orr x10, x10, x10",
            in("x4") args.as_ptr(), inout("x3") result, options(nostack),
        );
    }
    result
}

#[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
fn client_request(_: usize, _: usize, _: usize) -> usize {
    0
}

const RUNNING_ON_VALGRIND: usize = 0x1001;
const COUNT_ERRORS: usize = 0x1201;
const MAKE_MEM_UNDEFINED: usize = 0x4d43_0001;
const MAKE_MEM_DEFINED: usize = 0x4d43_0002;

/// Mark a value as secret: Memcheck treats its bytes as undefined from here on.
fn secret<T: ?Sized>(value: &mut T) {
    client_request(MAKE_MEM_UNDEFINED, value as *mut T as *mut u8 as usize, size_of_val(value));
}

/// Mark a value as public again.
fn public<T: ?Sized>(value: &mut T) {
    client_request(MAKE_MEM_DEFINED, value as *mut T as *mut u8 as usize, size_of_val(value));
}

/// Generator whose 64-bit words and bytes are secret.
struct Tainted(StdRng);

impl TryRng for Tainted {
    type Error = Infallible;
    fn try_next_u32(&mut self) -> Result<u32, Infallible> {
        Ok(self.0.next_u32())
    }
    fn try_next_u64(&mut self) -> Result<u64, Infallible> {
        let mut word = self.0.next_u64();
        secret(&mut word);
        Ok(word)
    }
    fn try_fill_bytes(&mut self, dst: &mut [u8]) -> Result<(), Infallible> {
        self.0.fill_bytes(dst);
        secret(dst);
        Ok(())
    }
}

// Only so that the type is accepted by the functions under test.
impl TryCryptoRng for Tainted {}

/// Run `f` and compare the number of reports it caused with the allowed range.
fn check<T>(failures: &mut u32, name: &str, allowed: std::ops::RangeInclusive<usize>, f: impl FnOnce() -> T) -> T {
    let before = client_request(COUNT_ERRORS, 0, 0);
    let out = f();
    let reports = client_request(COUNT_ERRORS, 0, 0) - before;
    let ok = allowed.contains(&reports);
    *failures += !ok as u32;
    let expected = match (*allowed.start(), *allowed.end()) {
        (0, 0) => "none allowed".to_string(),
        (low, usize::MAX) => format!("at least {low} expected"),
        (low, high) => format!("{low} to {high} allowed"),
    };
    println!("{name:<28} {reports:>6} reports ({expected}) -> {}", if ok { "ok" } else { "FAILED" });
    out
}

/// What a compiler must not do: one branch per coefficient on a secret.
#[inline(never)]
fn branchy_sum(values: &[i32; N]) -> i32 {
    let mut acc = 0i32;
    for &v in values.iter() {
        if v < 0 {
            acc = black_box(acc.wrapping_sub(v));
        } else {
            acc = acc.wrapping_add(v);
        }
    }
    acc
}

/// The two functions that end in a verdict on secret data branch on it: a few
/// jumps at most, against thousands for a scan that stops at the first bad value.
const VERDICT: std::ops::RangeInclusive<usize> = 1..=4;
const NONE: std::ops::RangeInclusive<usize> = 0..=0;

fn main() {
    let on_valgrind = client_request(RUNNING_ON_VALGRIND, 0, 0) != 0;
    let ctx = &NTT_CONTEXT;
    let mut rng = Tainted(StdRng::seed_from_u64(1));
    let mut failures = 0;
    let f = &mut failures;

    let mut samples = [0i32; N];
    check(f, "base sampler", NONE, || sample_gaussian_cdt(&mut rng, &mut samples));
    check(f, "control: secret branch", N / 4..=usize::MAX, || black_box(branchy_sum(&samples)));
    check(f, "flood sampler", NONE, || add_flood(&mut rng, &mut samples));

    let (sk, pk) = check(f, "keygen", NONE, || keygen(&mut rng, ctx));
    let mut matches = check(f, "PublicKey::verify_against", NONE, || pk.verify_against(&sk));
    public(&mut matches);
    assert!(matches);
    let sk_bytes = check(f, "SecretKey::to_bytes", NONE, || sk.to_bytes());
    let parsed = check(f, "SecretKey::from_bytes", VERDICT, || SecretKey::from_bytes(&sk_bytes));
    assert!(parsed.is_some());

    let mut data = [0x5au8; SLOT_BYTES];
    secret(&mut data);
    let message = check(f, "encode", NONE, || encode(&data));
    let ct = check(f, "encrypt_slot", NONE, || encrypt_slot(&mut rng, ctx, &pk, &message, Flooding::Standard));
    // Kept until the end, like the other results: memory that is freed while
    // marked is handed out again with its marks.
    let _unflooded = check(f, "encrypt_slot (no flood)", NONE, || encrypt_slot(&mut rng, ctx, &pk, &message, Flooding::Off));
    let ct2 = check(f, "rerandomize_slot", NONE, || rerandomize_slot(&mut rng, ctx, &pk, &ct));
    let recovered = check(f, "decrypt_slot", NONE, || decrypt_slot(ctx, &sk, &ct2));
    let mut decoded = check(f, "decode", VERDICT, || decode(&recovered)).expect("decryption failed");
    public(&mut decoded[..]);
    assert!(decoded.iter().all(|&b| b == 0x5a));

    let messages = [encode(&data), encode(&data), encode(&data)];
    let cts = check(f, "encrypt_batch", NONE, || encrypt_batch(&mut rng, ctx, &pk, &messages, Flooding::Standard));
    let cts2 = check(f, "rerandomize_batch", NONE, || rerandomize_batch(&mut rng, ctx, &pk, &cts));
    let recovered = check(f, "decrypt_batch", NONE, || decrypt_batch(ctx, &sk, &cts2));
    assert_eq!(recovered.len(), messages.len());

    if !on_valgrind {
        println!("\nNot running under Valgrind: nothing was checked.");
        std::process::exit(2);
    }
    println!("\n{}", if failures == 0 { "All checks passed." } else { "Secret-dependent branches or addresses found." });
    std::process::exit(failures.min(1) as i32);
}
