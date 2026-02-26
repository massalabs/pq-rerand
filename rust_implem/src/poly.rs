//! Polynomial types and operations in CRT domain.
//!
//! A polynomial in R_q = Z_q[x]/(x^N+1) is stored as two CRT limbs:
//! one mod t and one mod q₂. All arithmetic is limbwise.

use crate::ntt::{NttTables, pointwise_mul, pointwise_add, pointwise_sub};
use crate::params::{N, T, Q2};

/// A polynomial in R_q represented in CRT form (coefficient domain).
/// `limb_t[i]` is the i-th coefficient mod t; `limb_q2[i]` mod q₂.
#[derive(Clone, Debug)]
pub struct Poly {
    pub limb_t: [u32; N],
    pub limb_q2: [u32; N],
}

impl Poly {
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
            limb_t[i] = coeffs[i].rem_euclid(T as i64) as u32;
            limb_q2[i] = coeffs[i].rem_euclid(Q2 as i64) as u32;
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

/// A polynomial in NTT evaluation form (one limb).
#[derive(Clone)]
pub struct NttPoly {
    pub coeffs: [u32; N],
    pub modulus: u32,
}

/// Precomputed NTT context holding tables for both primes.
pub struct NttContext {
    pub tables_t: NttTables,
    pub tables_q2: NttTables,
}

impl NttContext {
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

    /// Ring multiply-accumulate: result += a * b (in NTT domain for efficiency).
    /// Takes a and b already in NTT form.
    pub fn ntt_mul_acc(
        acc_t: &mut [u32; N],
        acc_q2: &mut [u32; N],
        a_ntt_t: &[u32; N],
        a_ntt_q2: &[u32; N],
        b_ntt_t: &[u32; N],
        b_ntt_q2: &[u32; N],
    ) {
        for i in 0..N {
            let pt = (a_ntt_t[i] as u64 * b_ntt_t[i] as u64) % T as u64;
            acc_t[i] = ((acc_t[i] as u64 + pt) % T as u64) as u32;
            let pq = (a_ntt_q2[i] as u64 * b_ntt_q2[i] as u64) % Q2 as u64;
            acc_q2[i] = ((acc_q2[i] as u64 + pq) % Q2 as u64) as u32;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
