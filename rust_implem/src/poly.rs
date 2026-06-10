//! Polynomial types and operations in CRT domain.
//!
//! A polynomial in R_q = Z_q[x]/(x^N+1) is stored as two CRT limbs:
//! one mod t and one mod q₂. All arithmetic is limbwise.

use crate::ntt::{NttTables, pointwise_mul, pointwise_add, pointwise_sub};
use crate::params::{N, T, Q2};

/// Branchless reduction of a small signed coefficient into `[0, m)`.
///
/// Precondition: `|c| < m` (always satisfied by the Gaussian widths used here —
/// base σ = 3.2 and flood σ_f ≈ 93 733 are both astronomically smaller than the
/// 31-/32-bit moduli). Adds `m` iff `c` is negative, using an arithmetic-shift
/// mask, so there is no data-dependent branch and no division.
#[inline(always)]
fn reduce_signed_small(c: i64, m: u32) -> u32 {
    let mask = c >> 63; // -1 (all ones) if c < 0, else 0
    (c + (m as i64 & mask)) as u32
}

/// A polynomial in R_q represented in CRT form (coefficient domain).
/// `limb_t[i]` is the i-th coefficient mod t; `limb_q2[i]` mod q₂.
#[derive(Clone, Debug)]
pub struct Poly {
    /// Coefficients modulo `t`.
    pub limb_t: [u32; N],
    /// Coefficients modulo `q₂`.
    pub limb_q2: [u32; N],
}

impl Poly {
    /// The zero polynomial.
    pub fn zero() -> Self {
        Poly {
            limb_t: [0u32; N],
            limb_q2: [0u32; N],
        }
    }

    /// Create a polynomial from small signed integer coefficients.
    /// Reduces each coefficient into both moduli.
    pub fn from_signed(coeffs: &[i64; N]) -> Self {
        let mut limb_t = [0u32; N];
        let mut limb_q2 = [0u32; N];
        for i in 0..N {
            limb_t[i] = reduce_signed_small(coeffs[i], T);
            limb_q2[i] = reduce_signed_small(coeffs[i], Q2);
        }
        Poly { limb_t, limb_q2 }
    }

    /// Pointwise addition in both limbs.
    pub fn add(&self, other: &Poly) -> Poly {
        Poly {
            limb_t: pointwise_add(&self.limb_t, &other.limb_t, T),
            limb_q2: pointwise_add(&self.limb_q2, &other.limb_q2, Q2),
        }
    }

    /// Pointwise subtraction in both limbs.
    pub fn sub(&self, other: &Poly) -> Poly {
        Poly {
            limb_t: pointwise_sub(&self.limb_t, &other.limb_t, T),
            limb_q2: pointwise_sub(&self.limb_q2, &other.limb_q2, Q2),
        }
    }
}

/// Precomputed NTT context holding tables for both primes.
pub struct NttContext {
    /// Tables for the plaintext-limb prime `t`.
    pub tables_t: NttTables,
    /// Tables for the noise-limb prime `q₂`.
    pub tables_q2: NttTables,
}

/// Shared NTT context (twiddle tables for both primes). `NttContext::new` is a
/// `const fn`, so this static is built at compile time (no runtime table setup).
/// Used where a context is needed outside the hot path (e.g. key deserialization).
pub static NTT_CONTEXT: NttContext = NttContext::new();

impl Default for NttContext {
    fn default() -> Self {
        Self::new()
    }
}

impl NttContext {
    /// Build the context (const-evaluable; see [`NTT_CONTEXT`]).
    pub const fn new() -> Self {
        NttContext {
            tables_t: NttTables::new(T),
            tables_q2: NttTables::new(Q2),
        }
    }

    /// Forward NTT of a CRT polynomial → pair of NTT-domain arrays.
    pub fn forward(&self, p: &Poly) -> ([u32; N], [u32; N]) {
        let mut ft = p.limb_t;
        let mut fq = p.limb_q2;
        self.tables_t.forward(&mut ft);
        self.tables_q2.forward(&mut fq);
        (ft, fq)
    }

    /// Inverse NTT of NTT-domain arrays → CRT polynomial.
    pub fn inverse(&self, ntt_t: &mut [u32; N], ntt_q2: &mut [u32; N]) -> Poly {
        self.tables_t.inverse(ntt_t);
        self.tables_q2.inverse(ntt_q2);
        Poly {
            limb_t: *ntt_t,
            limb_q2: *ntt_q2,
        }
    }

    /// Ring multiplication: a * b in R_q (via NTT in both limbs).
    pub fn ring_mul(&self, a: &Poly, b: &Poly) -> Poly {
        let (a_ntt_t, a_ntt_q2) = self.forward(a);
        let (b_ntt_t, b_ntt_q2) = self.forward(b);
        let mut c_ntt_t = pointwise_mul(&a_ntt_t, &b_ntt_t, T);
        let mut c_ntt_q2 = pointwise_mul(&a_ntt_q2, &b_ntt_q2, Q2);
        self.inverse(&mut c_ntt_t, &mut c_ntt_q2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_signed_small_matches_rem_euclid() {
        for &m in &[T, Q2] {
            for c in -1_000_000i64..=1_000_000 {
                assert_eq!(reduce_signed_small(c, m) as i64, c.rem_euclid(m as i64), "c={c} m={m}");
            }
            // larger magnitudes still within (-m, m)
            for &c in &[-(m as i64 - 1), -(m as i64) / 2, m as i64 / 2, m as i64 - 1] {
                assert_eq!(reduce_signed_small(c, m) as i64, c.rem_euclid(m as i64), "c={c} m={m}");
            }
        }
    }

    #[test]
    fn test_poly_add_sub() {
        let mut a_coeffs = [0i64; N];
        let mut b_coeffs = [0i64; N];
        a_coeffs[0] = 100;
        a_coeffs[1] = 200;
        b_coeffs[0] = 50;
        b_coeffs[1] = 300;

        let a = Poly::from_signed(&a_coeffs);
        let b = Poly::from_signed(&b_coeffs);

        let sum = a.add(&b);
        assert_eq!(sum.limb_t[0], 150);
        assert_eq!(sum.limb_t[1], 500);

        let diff = a.sub(&b);
        assert_eq!(diff.limb_t[0], 50);
        // 200 - 300 mod T = T - 100
        assert_eq!(diff.limb_t[1], T - 100);
    }

    #[test]
    fn test_ring_mul() {
        let ctx = NttContext::new();
        let mut a_coeffs = [0i64; N];
        let mut b_coeffs = [0i64; N];
        a_coeffs[0] = 1;
        a_coeffs[1] = 1; // 1 + x
        b_coeffs[0] = 1;
        b_coeffs[1] = 1; // 1 + x

        let a = Poly::from_signed(&a_coeffs);
        let b = Poly::from_signed(&b_coeffs);
        let c = ctx.ring_mul(&a, &b);

        // (1+x)^2 = 1 + 2x + x^2
        assert_eq!(c.limb_t[0], 1);
        assert_eq!(c.limb_t[1], 2);
        assert_eq!(c.limb_t[2], 1);
        assert_eq!(c.limb_q2[0], 1);
        assert_eq!(c.limb_q2[1], 2);
        assert_eq!(c.limb_q2[2], 1);
    }
}
