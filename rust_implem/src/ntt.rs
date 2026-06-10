//! Number Theoretic Transform (NTT) for negacyclic convolution.
//!
//! Optimized implementation:
//!   * **Barrett reduction** ([`barrett_reduce`]) replaces the per-multiply `% p`
//!     division on the hot path (one 128-bit multiply + shift + ≤2 subtractions).
//!   * **Precomputed twiddle tables** in bit-reversed order: the merged-negacyclic
//!     Cooley–Tukey forward (Longa–Naehrig) and Gentleman–Sande inverse fold the
//!     ψ^i pre/post-scaling into the twiddles and need no separate bit-reversal,
//!     and no `powmod` inside the butterfly loops.
//!   * Branchless, slice-based pointwise kernels that auto-vectorize under
//!     `-C target-cpu=native` (AVX2/AVX-512).
//!
//! The forward transform maps natural order → bit-reversed order; the inverse maps
//! bit-reversed order → natural order. Forward/pointwise/inverse are always paired,
//! and pointwise products are order-agnostic, so the convention is internal.

use crate::params::N;

// ─────────────────────────────────────────────────────────────────────────────
// Constant-time modular layer.
//
// All reductions below are *branchless* and use no data-indexed memory accesses,
// so their control flow and memory-access pattern are independent of the operand
// values. Combined with the data-independent NTT addressing (loop bounds depend
// only on N) and the constant-time CDT sampler, this gives an arithmetic layer
// whose timing does not depend on secret data. (This is constant-time *by
// construction*; it is not formally verified against a specific microarchitecture,
// and production use should still confirm with a timing-leakage test such as
// dudect.)
// ─────────────────────────────────────────────────────────────────────────────

/// Branchless conditional subtract: returns `r - p` if `r >= p`, else `r`.
/// Correct for `r < 2^63` and `p < 2^32`. No data-dependent branch.
#[inline(always)]
pub const fn csub(r: u64, p: u64) -> u64 {
    let d = r.wrapping_sub(p);
    // If r < p the subtraction borrows and bit 63 of `d` is set; build an all-ones
    // mask in that case and add `p` back. If r >= p, mask is 0.
    let mask = (d >> 63).wrapping_neg();
    d.wrapping_add(p & mask)
}

/// Modular multiplication: (a * b) mod m, via u64 division. Kept for non-hot-path,
/// non-secret callers (e.g. compile-time table construction). The hot path uses
/// the branchless [`barrett_reduce`] instead.
#[inline(always)]
pub const fn mulmod(a: u32, b: u32, m: u32) -> u32 {
    ((a as u64 * b as u64) % m as u64) as u32
}

/// Modular addition: (a + b) mod m. Branchless / constant-time.
#[inline(always)]
pub const fn addmod(a: u32, b: u32, m: u32) -> u32 {
    csub(a as u64 + b as u64, m as u64) as u32
}

/// Modular subtraction: (a - b) mod m. Branchless / constant-time.
#[inline(always)]
pub const fn submod(a: u32, b: u32, m: u32) -> u32 {
    let d = (a as u64).wrapping_sub(b as u64);
    let mask = (d >> 63).wrapping_neg();
    d.wrapping_add(m as u64 & mask) as u32
}

/// Modular exponentiation: base^exp mod m.
pub const fn powmod(base: u32, mut exp: u64, m: u32) -> u32 {
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
/// Computes `q = ⌊x·μ / 2⁶⁴⌋ ≤ ⌊x/p⌋` (so `q·p ≤ x`, no overflow), then
/// `r = x − q·p ∈ [0, 3p)` and finishes with two branchless conditional
/// subtractions. Constant-time.
#[inline(always)]
pub const fn barrett_reduce(x: u64, p: u32, mu: u64) -> u32 {
    let q = ((x as u128 * mu as u128) >> 64) as u64;
    let p64 = p as u64;
    let r = x.wrapping_sub(q.wrapping_mul(p64));
    csub(csub(r, p64), p64) as u32
}

/// Shoup constant for a *fixed* multiplier `w < p`: `w' = ⌊w·2⁶⁴ / p⌋`.
#[inline(always)]
pub const fn shoup(w: u32, p: u32) -> u64 {
    (((w as u128) << 64) / p as u128) as u64
}

/// Shoup multiplication: `(a·w) mod p` for a fixed multiplier `w` with precomputed
/// `w_shoup = shoup(w, p)`, where `a, w < p < 2³¹`. One 64-bit high-multiply, two
/// 64-bit multiplies and a single branchless conditional subtraction — cheaper than
/// general Barrett (which needs two). Constant-time.
#[inline(always)]
pub const fn mul_shoup(a: u32, w: u32, w_shoup: u64, p: u32) -> u32 {
    // q ≈ ⌊a·w / p⌋, so r = a·w − q·p ∈ [0, 2p); one csub finishes the reduction.
    let q = ((a as u128 * w_shoup as u128) >> 64) as u64;
    let r = (a as u64).wrapping_mul(w as u64).wrapping_sub(q.wrapping_mul(p as u64));
    csub(r, p as u64) as u32
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
/// [`mul_shoup`] (single conditional subtraction) instead of general Barrett.
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

/// Pointwise multiplication of two NTT-domain polynomials (Barrett reduction).
pub fn pointwise_mul(a: &[u32; N], b: &[u32; N], p: u32) -> [u32; N] {
    let mu = barrett_mu(p);
    let mut c = [0u32; N];
    for i in 0..N {
        c[i] = barrett_reduce(a[i] as u64 * b[i] as u64, p, mu);
    }
    c
}

/// Pointwise addition. Branchless conditional subtract (constant-time,
/// auto-vectorizable).
pub fn pointwise_add(a: &[u32; N], b: &[u32; N], p: u32) -> [u32; N] {
    let p64 = p as u64;
    let mut c = [0u32; N];
    for i in 0..N {
        c[i] = csub(a[i] as u64 + b[i] as u64, p64) as u32;
    }
    c
}

/// Pointwise subtraction. Branchless (constant-time, auto-vectorizable).
pub fn pointwise_sub(a: &[u32; N], b: &[u32; N], p: u32) -> [u32; N] {
    let mut c = [0u32; N];
    for i in 0..N {
        c[i] = submod(a[i], b[i], p);
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{T, Q2};

    /// Verify NTT tables are computed at compile time.
    const _: NttTables = NttTables::new(T);
    const _: NttTables = NttTables::new(Q2);

    #[test]
    fn test_barrett_matches_div() {
        for &p in &[T, Q2] {
            let mu = barrett_mu(p);
            let cases = [
                0u64, 1, (p as u64) - 1, p as u64, p as u64 + 1,
                (p as u64 - 1) * (p as u64 - 1),
                1234567u64 * 7654321u64,
                (p as u64 - 2) * (p as u64 - 3),
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
            }
            // csub edge coverage across [0, 3p)
            for r in [0u64, 1, p as u64 - 1, p as u64, p as u64 + 1, 2 * p as u64 - 1, 2 * p as u64, 3 * p as u64 - 1] {
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
        let mut c = pointwise_mul(&a, &b, T);
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
        let mut c = pointwise_mul(&a, &b, T);
        tables.inverse(&mut c);
        // x^{N-1} * x = x^N = -1 mod (x^N+1)
        assert_eq!(c[0], T - 1); // -1 mod T
        for (i, &ci) in c.iter().enumerate().skip(1) { assert_eq!(ci, 0, "nonzero at {}", i); }
    }
}
