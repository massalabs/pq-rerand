//! Number Theoretic Transform (NTT) for negacyclic convolution.
//!
//! Optimized implementation:
//!   * **Barrett reduction** ([`barrett_reduce`]) replaces the per-multiply `% p`
//!     division on the hot path (one 128-bit multiply, a shift and one
//!     conditional subtraction).
//!   * **Precomputed twiddle tables** in bit-reversed order: the merged-negacyclic
//!     Cooley–Tukey forward (Longa–Naehrig) and Gentleman–Sande inverse fold the
//!     ψ^i pre/post-scaling into the twiddles and need no separate bit-reversal,
//!     and no `powmod` inside the butterfly loops.
//!   * Branchless, slice-based pointwise kernels that auto-vectorize under
//!     `-C target-cpu=native` (AVX2/AVX-512).
//!
//! The forward transform maps natural order → bit-reversed order; the inverse maps
//! bit-reversed order → natural order. Forward/pointwise/inverse are always paired,
//! and pointwise products are order-agnostic, so the convention does not matter to
//! the arithmetic. It is part of the serialized form of public keys, which are
//! stored in the NTT domain (see `PublicKey`).

use crate::params::N;
use core::hint::select_unpredictable;

// ─────────────────────────────────────────────────────────────────────────────
// Constant-time modular layer.
//
// All reductions below are *branchless* and use no data-indexed memory accesses,
// so their control flow and memory-access pattern are independent of the operand
// values. Their selections go through `select_unpredictable`: mask arithmetic
// alone does not keep a compiler from lowering a selection to a conditional
// branch, whereas the hint asks for a conditional move and still vectorizes.
// Combined with the data-independent NTT addressing (loop bounds depend
// only on N) and the constant-time CDT sampler, this gives an arithmetic layer
// whose timing does not depend on secret data. (This is constant-time *by
// construction*; it is not formally verified against a specific microarchitecture,
// and production use should still confirm with a timing-leakage test such as
// dudect.)
// ─────────────────────────────────────────────────────────────────────────────

/// Branchless conditional subtract: returns `r - p` if `r >= p`, else `r`.
/// No data-dependent branch.
#[inline(always)]
pub fn csub(r: u64, p: u64) -> u64 {
    select_unpredictable(r >= p, r.wrapping_sub(p), r)
}

/// Modular addition: (a + b) mod m, for `a, b < m`. Branchless / constant-time.
#[inline(always)]
pub fn addmod(a: u32, b: u32, m: u32) -> u32 {
    csub(a as u64 + b as u64, m as u64) as u32
}

/// Modular subtraction: (a - b) mod m, for `a, b < m`. Branchless / constant-time.
#[inline(always)]
pub fn submod(a: u32, b: u32, m: u32) -> u32 {
    let d = a.wrapping_sub(b);
    select_unpredictable(a < b, d.wrapping_add(m), d)
}

/// Modular exponentiation: base^exp mod m. Variable-time: only used to build the
/// (public) twiddle tables at compile time.
const fn powmod(base: u32, mut exp: u64, m: u32) -> u32 {
    let mut result: u64 = 1;
    let mut b = base as u64 % m as u64;
    let modulus = m as u64;
    while exp > 0 {
        if exp & 1 == 1 { result = result * b % modulus; }
        b = b * b % modulus;
        exp >>= 1;
    }
    result as u32
}

/// Barrett constant μ = ⌊2⁶⁴ / p⌋ for a modulus `p < 2³²`.
#[inline(always)]
pub const fn barrett_mu(p: u32) -> u64 {
    ((1u128 << 64) / p as u128) as u64
}

/// Barrett reduction of `x < 2⁶⁴` modulo `p < 2³²`, given `mu = barrett_mu(p)`.
///
/// Computes `q = floor(x*mu / 2^64)`, which is `floor(x/p)` or one less (so
/// `q*p <= x`, no overflow). The remainder `r = x - q*p` therefore lies in
/// `[0, 2p)` and one branchless conditional subtraction finishes the reduction.
/// Constant-time.
#[inline(always)]
pub fn barrett_reduce(x: u64, p: u32, mu: u64) -> u32 {
    let q = ((x as u128 * mu as u128) >> 64) as u64;
    let p64 = p as u64;
    let r = x.wrapping_sub(q.wrapping_mul(p64));
    csub(r, p64) as u32
}

/// Shoup constant for a *fixed* multiplier `w < p`: `w' = ⌊w·2⁶⁴ / p⌋`.
#[inline(always)]
pub const fn shoup(w: u32, p: u32) -> u64 {
    (((w as u128) << 64) / p as u128) as u64
}

/// Shoup multiplication: `(a*w) mod p` for a fixed multiplier `w` with precomputed
/// `w_shoup = shoup(w, p)`, where `a, w < p < 2^32` and `p` is prime. One 64-bit
/// high-multiply and two 64-bit multiplies, with no final correction: unlike
/// Barrett, the quotient is exact. Constant-time.
#[inline(always)]
pub const fn mul_shoup(a: u32, w: u32, w_shoup: u64, p: u32) -> u32 {
    // With w*2^64 = w_shoup*p + e (0 <= e < p) and a*w_shoup = q*2^64 + d
    // (0 <= d < 2^64), the remainder r = a*w - q*p equals (d*p + a*e) / 2^64, so
    // 0 <= r < p + a*p/2^64 < p + 1. It cannot be p either: p would then divide
    // a*w, which for a prime p means a = 0 or w = 0, and both give r = 0.
    let q = ((a as u128 * w_shoup as u128) >> 64) as u64;
    (a as u64).wrapping_mul(w as u64).wrapping_sub(q.wrapping_mul(p as u64)) as u32
}

/// Bit-reverse the low `bits` bits of `x`.
const fn bit_reverse(mut x: usize, bits: usize) -> usize {
    let mut r = 0;
    let mut i = 0;
    while i < bits {
        r = (r << 1) | (x & 1);
        x >>= 1;
        i += 1;
    }
    r
}

/// Precomputed NTT tables for a specific prime.
///
/// Each twiddle is stored together with its Shoup constant so the butterflies use
/// [`mul_shoup`] (exact, no conditional subtraction) instead of general Barrett.
pub struct NttTables {
    /// The NTT prime `p` these tables were built for.
    pub modulus: u32,
    n_inv: u32,           // N^{-1} mod p
    n_inv_shoup: u64,     // shoup(n_inv, p)
    // Twiddles in bit-reversed order: psi_rev[k] = psi^{bitrev(k)},
    // psi_inv_rev[k] = psi^{-bitrev(k)}, each with its Shoup constant.
    psi_rev: [u32; N],
    psi_rev_shoup: [u64; N],
    psi_inv_rev: [u32; N],
    psi_inv_rev_shoup: [u64; N],
}

/// Find a primitive 2N-th root of unity modulo p.
const fn find_psi(p: u32) -> u32 {
    let two_n = 2 * N as u64;
    let exp = (p as u64 - 1) / two_n;
    let mut g = 2u32;
    while g < 200 {
        let w = powmod(g, exp, p);
        if w != 0 && w != 1 {
            // Check order is exactly 2N: w^N ≡ -1 (mod p)
            let w_n = powmod(w, N as u64, p);
            if w_n == p - 1 { return w; }
        }
        g += 1;
    }
    panic!("No 2N-th root of unity found");
}

impl NttTables {
    /// Build the twiddle tables for prime `p` (const-evaluable: tables for the
    /// scheme's two primes are computed at compile time).
    pub const fn new(p: u32) -> Self {
        let psi = find_psi(p);
        let psi_inv = powmod(psi, p as u64 - 2, p);
        let n_inv = powmod(N as u32, p as u64 - 2, p);
        let log_n = N.trailing_zeros() as usize;

        let mut psi_rev = [0u32; N];
        let mut psi_rev_shoup = [0u64; N];
        let mut psi_inv_rev = [0u32; N];
        let mut psi_inv_rev_shoup = [0u64; N];
        let mut k = 0;
        while k < N {
            let e = bit_reverse(k, log_n) as u64;
            let pr = powmod(psi, e, p);
            let pir = powmod(psi_inv, e, p);
            psi_rev[k] = pr;
            psi_rev_shoup[k] = shoup(pr, p);
            psi_inv_rev[k] = pir;
            psi_inv_rev_shoup[k] = shoup(pir, p);
            k += 1;
        }

        NttTables {
            modulus: p,
            n_inv,
            n_inv_shoup: shoup(n_inv, p),
            psi_rev,
            psi_rev_shoup,
            psi_inv_rev,
            psi_inv_rev_shoup,
        }
    }

    /// Forward negacyclic NTT (natural order → bit-reversed order).
    /// Merged Cooley–Tukey (Longa–Naehrig): no separate bit-reversal pass.
    pub fn forward(&self, a: &mut [u32; N]) {
        let p = self.modulus;
        let mut t = N;
        let mut m = 1;
        while m < N {
            t >>= 1;
            let mut i = 0;
            while i < m {
                let j1 = 2 * i * t;
                let s = self.psi_rev[m + i];
                let s_shoup = self.psi_rev_shoup[m + i];
                let mut j = j1;
                while j < j1 + t {
                    let u = a[j];
                    let v = mul_shoup(a[j + t], s, s_shoup, p);
                    a[j] = addmod(u, v, p);
                    a[j + t] = submod(u, v, p);
                    j += 1;
                }
                i += 1;
            }
            m <<= 1;
        }
    }

    /// Inverse negacyclic NTT (bit-reversed order → natural order).
    /// Gentleman–Sande, with final N^{-1} scaling.
    pub fn inverse(&self, a: &mut [u32; N]) {
        let p = self.modulus;
        let mut t = 1;
        let mut m = N;
        while m > 1 {
            let h = m >> 1;
            let mut j1 = 0;
            let mut i = 0;
            while i < h {
                let s = self.psi_inv_rev[h + i];
                let s_shoup = self.psi_inv_rev_shoup[h + i];
                let mut j = j1;
                while j < j1 + t {
                    let u = a[j];
                    let w = a[j + t];
                    a[j] = addmod(u, w, p);
                    a[j + t] = mul_shoup(submod(u, w, p), s, s_shoup, p);
                    j += 1;
                }
                j1 += 2 * t;
                i += 1;
            }
            t <<= 1;
            m >>= 1;
        }
        let n_inv = self.n_inv;
        let n_inv_shoup = self.n_inv_shoup;
        for x in a.iter_mut() {
            *x = mul_shoup(*x, n_inv, n_inv_shoup, p);
        }
    }
}

// The pointwise kernels work in place: returning an array by value would leave a
// copy of a possibly secret operand in a dead stack frame.

/// Pointwise multiplication `a *= b` of two NTT-domain polynomials (Barrett
/// reduction).
pub fn pointwise_mul(a: &mut [u32; N], b: &[u32; N], p: u32) {
    let mu = barrett_mu(p);
    for i in 0..N {
        a[i] = barrett_reduce(a[i] as u64 * b[i] as u64, p, mu);
    }
}

/// Pointwise addition `a += b`. Branchless conditional subtract (constant-time,
/// auto-vectorizable).
pub fn pointwise_add(a: &mut [u32; N], b: &[u32; N], p: u32) {
    for i in 0..N {
        a[i] = addmod(a[i], b[i], p);
    }
}

/// Pointwise subtraction `a -= b`. Branchless (constant-time, auto-vectorizable).
pub fn pointwise_sub(a: &mut [u32; N], b: &[u32; N], p: u32) {
    for i in 0..N {
        a[i] = submod(a[i], b[i], p);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{T, Q2};

    #[test]
    fn test_barrett_matches_div() {
        for &p in &[T, Q2] {
            let mu = barrett_mu(p);
            let top = u64::MAX / p as u64 * p as u64; // largest multiple of p below 2^64
            let cases = [
                0u64, 1, (p as u64) - 1, p as u64, p as u64 + 1,
                (p as u64 - 1) * (p as u64 - 1),
                1234567u64 * 7654321u64,
                (p as u64 - 2) * (p as u64 - 3),
                top - 1, top, top + 1, u64::MAX - 1, u64::MAX,
            ];
            for &x in &cases {
                assert_eq!(barrett_reduce(x, p, mu) as u64, x % p as u64, "p={p} x={x}");
            }
        }
    }

    // Deterministic xorshift for randomized equivalence checks (no test deps).
    fn xs(state: &mut u64) -> u64 {
        let mut x = *state;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        *state = x;
        x
    }

    #[test]
    fn test_shoup_matches_reference() {
        let mut st = 0x0bad_c0de_dead_beefu64;
        for &p in &[T, Q2] {
            for _ in 0..200_000 {
                let a = (xs(&mut st) % p as u64) as u32;
                let w = (xs(&mut st) % p as u64) as u32;
                let ws = shoup(w, p);
                let got = mul_shoup(a, w, ws, p) as u64;
                assert_eq!(got, (a as u64 * w as u64) % p as u64, "p={p} a={a} w={w}");
            }
            // edges
            for &(a, w) in &[(0u32, 0u32), (0, p - 1), (p - 1, 0), (p - 1, p - 1), (1, p - 1)] {
                let ws = shoup(w, p);
                assert_eq!(mul_shoup(a, w, ws, p) as u64, (a as u64 * w as u64) % p as u64);
            }
        }
    }

    #[test]
    fn test_branchless_modops_match_reference() {
        let mut st = 0x1234_5678_9abc_def1u64;
        for &p in &[T, Q2] {
            let mu = barrett_mu(p);
            for _ in 0..200_000 {
                let a = (xs(&mut st) % p as u64) as u32;
                let b = (xs(&mut st) % p as u64) as u32;
                // addmod / submod vs reference
                assert_eq!(addmod(a, b, p) as u64, (a as u64 + b as u64) % p as u64);
                let rs = ((a as i64 - b as i64).rem_euclid(p as i64)) as u64;
                assert_eq!(submod(a, b, p) as u64, rs);
                // barrett product vs reference
                let x = a as u64 * b as u64;
                assert_eq!(barrett_reduce(x, p, mu) as u64, x % p as u64);
                // barrett over the full u64 range
                let x = xs(&mut st);
                assert_eq!(barrett_reduce(x, p, mu) as u64, x % p as u64);
            }
            // csub edge coverage, up to the top of the u64 range
            for r in [0u64, 1, p as u64 - 1, p as u64, p as u64 + 1, 2 * p as u64 - 1, 2 * p as u64, 3 * p as u64 - 1, u64::MAX] {
                let expect = if r >= p as u64 { r - p as u64 } else { r };
                assert_eq!(csub(r, p as u64), expect, "csub r={r} p={p}");
            }
        }
    }

    #[test]
    fn test_ntt_roundtrip_t() {
        let tables = NttTables::new(T);
        let mut a = [0u32; N];
        a[0] = 1; a[1] = 2; a[2] = 3;
        let original = a;
        tables.forward(&mut a);
        tables.inverse(&mut a);
        assert_eq!(a, original);
    }

    #[test]
    fn test_ntt_roundtrip_q2() {
        let tables = NttTables::new(Q2);
        let mut a = [0u32; N];
        a[0] = 42; a[1] = 100; a[N - 1] = 999;
        let original = a;
        tables.forward(&mut a);
        tables.inverse(&mut a);
        assert_eq!(a, original);
    }

    #[test]
    fn test_negacyclic_mul() {
        let tables = NttTables::new(T);
        let mut a = [0u32; N];
        a[0] = 1; a[1] = 1; // 1 + x
        let mut b = [0u32; N];
        b[0] = 1; b[1] = 1; // 1 + x
        tables.forward(&mut a);
        tables.forward(&mut b);
        let mut c = a;
        pointwise_mul(&mut c, &b, T);
        tables.inverse(&mut c);
        // (1+x)^2 = 1 + 2x + x^2
        assert_eq!(c[0], 1);
        assert_eq!(c[1], 2);
        assert_eq!(c[2], 1);
        for (i, &ci) in c.iter().enumerate().skip(3) { assert_eq!(ci, 0, "nonzero at {}", i); }
    }

    #[test]
    fn test_negacyclic_wraparound() {
        let tables = NttTables::new(T);
        let mut a = [0u32; N];
        a[N - 1] = 1; // x^{N-1}
        let mut b = [0u32; N];
        b[1] = 1; // x
        tables.forward(&mut a);
        tables.forward(&mut b);
        let mut c = a;
        pointwise_mul(&mut c, &b, T);
        tables.inverse(&mut c);
        // x^{N-1} * x = x^N = -1 mod (x^N+1)
        assert_eq!(c[0], T - 1); // -1 mod T
        for (i, &ci) in c.iter().enumerate().skip(1) { assert_eq!(ci, 0, "nonzero at {}", i); }
    }

    /// Dense products against the schoolbook negacyclic convolution, and the
    /// pointwise sum and difference against their definitions, in both limbs.
    #[test]
    fn test_dense_mul_matches_schoolbook() {
        let mut st = 0x5eed_0000_0000_0001u64;
        for &p in &[T, Q2] {
            let tables = NttTables::new(p);
            let mut dense = || -> [u32; N] { std::array::from_fn(|_| (xs(&mut st) % p as u64) as u32) };
            let top = [p - 1; N];
            for (a, b) in [(dense(), dense()), (top, top), (top, dense())] {
                // x^N = -1: the products that wrap around are subtracted.
                let (mut pos, mut neg) = (vec![0u128; N], vec![0u128; N]);
                for i in 0..N {
                    for j in 0..N {
                        let prod = (a[i] as u64 * b[j] as u64) as u128;
                        if i + j < N { pos[i + j] += prod } else { neg[i + j - N] += prod }
                    }
                }
                let (mut fa, mut fb) = (a, b);
                tables.forward(&mut fa);
                tables.forward(&mut fb);
                pointwise_mul(&mut fa, &fb, p);
                tables.inverse(&mut fa);
                for k in 0..N {
                    let want = (pos[k] % p as u128 + p as u128 - neg[k] % p as u128) % p as u128;
                    assert_eq!(fa[k] as u128, want, "product mod {p}, coefficient {k}");
                }

                let (mut sum, mut diff) = (a, a);
                pointwise_add(&mut sum, &b, p);
                pointwise_sub(&mut diff, &b, p);
                for i in 0..N {
                    assert_eq!(sum[i] as u64, (a[i] as u64 + b[i] as u64) % p as u64);
                    assert_eq!(diff[i] as u64, (a[i] as u64 + p as u64 - b[i] as u64) % p as u64);
                }
            }
        }
    }

    /// The twiddle tables hold powers of a primitive 2N-th root of unity, in
    /// bit-reversed order, together with their inverses and Shoup constants.
    #[test]
    fn test_tables_are_consistent() {
        for &p in &[T, Q2] {
            let t = NttTables::new(p);
            let psi = t.psi_rev[bit_reverse(1, N.trailing_zeros() as usize)];
            // The root is part of the serialized public-key format.
            assert_eq!(psi, if p == T { 101_573_657 } else { 1_953_722_822 });
            assert_eq!(powmod(psi, N as u64, p), p - 1, "psi^N must be -1");
            assert_eq!((t.n_inv as u64 * N as u64) % p as u64, 1);
            assert_eq!(t.n_inv_shoup, shoup(t.n_inv, p));
            for k in 0..N {
                let e = bit_reverse(k, N.trailing_zeros() as usize) as u64;
                assert_eq!(t.psi_rev[k], powmod(psi, e, p));
                assert_eq!((t.psi_rev[k] as u64 * t.psi_inv_rev[k] as u64) % p as u64, 1);
                assert_eq!(t.psi_rev_shoup[k], shoup(t.psi_rev[k], p));
                assert_eq!(t.psi_inv_rev_shoup[k], shoup(t.psi_inv_rev[k], p));
            }
        }
    }

    /// Output order of the forward transform, which public keys are serialized in:
    /// index k holds the evaluation at psi^(2*bitrev(k)+1).
    #[test]
    fn test_forward_order() {
        let mut st = 0x0123_4567_89ab_cdefu64;
        for &p in &[T, Q2] {
            let tables = NttTables::new(p);
            let psi = find_psi(p);
            let a: [u32; N] = std::array::from_fn(|_| (xs(&mut st) % p as u64) as u32);
            let mut fa = a;
            tables.forward(&mut fa);
            for k in [0, 1, 2, 3, 1000, 2047, 2048, 4095] {
                let x = powmod(psi, 2 * bit_reverse(k, 12) as u64 + 1, p) as u64;
                let horner = a.iter().rev().fold(0u64, |acc, &c| (acc * x + c as u64) % p as u64);
                assert_eq!(fa[k] as u64, horner, "p={p} k={k}");
            }
        }
    }
}
