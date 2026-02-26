//! Number Theoretic Transform (NTT) for negacyclic convolution.
//!
//! Uses the "pre-multiply by ψ^i" approach:
//! 1. Multiply a[i] by ψ^i (converts negacyclic to cyclic)
//! 2. Standard radix-2 DIT NTT using ω = ψ² (N-th root of unity)
//! 3. Pointwise multiply
//! 4. Standard INTT using ω⁻¹
//! 5. Multiply result[i] by ψ^{-i} and scale by N⁻¹

use crate::params::N;

/// Modular multiplication: (a * b) mod m, via u64.
#[inline(always)]
pub const fn mulmod(a: u32, b: u32, m: u32) -> u32 {
    ((a as u64 * b as u64) % m as u64) as u32
}

/// Modular addition: (a + b) mod m.
#[inline(always)]
pub const fn addmod(a: u32, b: u32, m: u32) -> u32 {
    let s = a as u64 + b as u64;
    (if s >= m as u64 { s - m as u64 } else { s }) as u32
}

/// Modular subtraction: (a - b) mod m.
#[inline(always)]
pub const fn submod(a: u32, b: u32, m: u32) -> u32 {
    if a >= b { a - b } else { (a as u64 + m as u64 - b as u64) as u32 }
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

/// Bit-reverse an index in log2(N) bits.
fn bit_reverse(mut x: usize, bits: usize) -> usize {
    let mut r = 0;
    for _ in 0..bits {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    r
}

/// Precomputed NTT tables for a specific prime.
pub struct NttTables {
    pub modulus: u32,
    omega: u32,      // psi^2 = primitive N-th root of unity
    omega_inv: u32,  // omega^{-1} mod p
    n_inv: u32,      // N^{-1} mod p
    // Precomputed: psi_table[i] = psi^i for i = 0..N-1
    psi_table: [u32; N],
    psi_inv_table: [u32; N],
}

impl NttTables {
    pub const fn new(p: u32) -> Self {
        let psi = find_psi(p);
        let psi_inv = powmod(psi, p as u64 - 2, p);
        let omega = mulmod(psi, psi, p); // ω = ψ²
        let omega_inv = powmod(omega, p as u64 - 2, p);
        let n_inv = powmod(N as u32, p as u64 - 2, p);

        // psi_table[i] = psi^i
        let mut psi_table = [0u32; N];
        let mut psi_inv_table = [0u32; N];
        psi_table[0] = 1;
        psi_inv_table[0] = 1;
        let mut i = 1;
        while i < N {
            psi_table[i] = mulmod(psi_table[i - 1], psi, p);
            psi_inv_table[i] = mulmod(psi_inv_table[i - 1], psi_inv, p);
            i += 1;
        }

        NttTables {
            modulus: p, omega, omega_inv, n_inv,
            psi_table, psi_inv_table,
        }
    }

    /// Forward negacyclic NTT.
    /// Transforms coefficient-domain polynomial to evaluation form.
    pub fn forward(&self, a: &mut [u32; N]) {
        let p = self.modulus;
        // Step 1: Pre-multiply by ψ^i (converts negacyclic to cyclic)
        for i in 0..N {
            a[i] = mulmod(a[i], self.psi_table[i], p);
        }
        // Step 2: Standard Cooley-Tukey DIT radix-2 NTT
        self.dit_ntt_cyclic(a);
    }

    /// Standard DIT NTT for cyclic convolution using omega.
    fn dit_ntt_cyclic(&self, a: &mut [u32; N]) {
        let p = self.modulus;
        let log_n = N.trailing_zeros() as usize;

        // Bit-reverse permutation
        for i in 0..N {
            let j = bit_reverse(i, log_n);
            if i < j { a.swap(i, j); }
        }

        // Butterfly stages
        let mut len = 2;
        while len <= N {
            let half = len / 2;
            let step = N / len;
            // twiddle at position j: omega^(j * step) where omega is N-th root
            for start in (0..N).step_by(len) {
                let mut w: u64 = 1;
                let omega_step = powmod(self.omega, step as u64, p);
                for j in 0..half {
                    let u = a[start + j];
                    let v = mulmod(a[start + j + half], w as u32, p);
                    a[start + j] = addmod(u, v, p);
                    a[start + j + half] = submod(u, v, p);
                    w = w * omega_step as u64 % p as u64;
                }
            }
            len *= 2;
        }
    }

    /// Inverse negacyclic NTT.
    /// Transforms evaluation form back to coefficient domain.
    pub fn inverse(&self, a: &mut [u32; N]) {
        let p = self.modulus;
        // Step 1: Standard INTT (DIT using omega_inv)
        self.dit_intt_cyclic(a);
        // Step 2: Scale by N^{-1} and post-multiply by ψ^{-i}
        for i in 0..N {
            a[i] = mulmod(a[i], self.n_inv, p);
            a[i] = mulmod(a[i], self.psi_inv_table[i], p);
        }
    }

    /// Standard DIT INTT for cyclic convolution using omega_inv.
    fn dit_intt_cyclic(&self, a: &mut [u32; N]) {
        let p = self.modulus;
        let log_n = N.trailing_zeros() as usize;

        // Bit-reverse permutation
        for i in 0..N {
            let j = bit_reverse(i, log_n);
            if i < j { a.swap(i, j); }
        }

        // Butterfly stages (same as forward but with omega_inv)
        let mut len = 2;
        while len <= N {
            let half = len / 2;
            let step = N / len;
            for start in (0..N).step_by(len) {
                let omega_inv_step = powmod(self.omega_inv, step as u64, p);
                let mut w: u64 = 1;
                for j in 0..half {
                    let u = a[start + j];
                    let v = mulmod(a[start + j + half], w as u32, p);
                    a[start + j] = addmod(u, v, p);
                    a[start + j + half] = submod(u, v, p);
                    w = w * omega_inv_step as u64 % p as u64;
                }
            }
            len *= 2;
        }
    }
}

/// Pointwise multiplication of two NTT-domain polynomials.
pub fn pointwise_mul(a: &[u32; N], b: &[u32; N], p: u32) -> [u32; N] {
    let mut c = [0u32; N];
    for i in 0..N { c[i] = mulmod(a[i], b[i], p); }
    c
}

/// Pointwise addition.
pub fn pointwise_add(a: &[u32; N], b: &[u32; N], p: u32) -> [u32; N] {
    let mut c = [0u32; N];
    for i in 0..N { c[i] = addmod(a[i], b[i], p); }
    c
}

/// Pointwise subtraction.
pub fn pointwise_sub(a: &[u32; N], b: &[u32; N], p: u32) -> [u32; N] {
    let mut c = [0u32; N];
    for i in 0..N { c[i] = submod(a[i], b[i], p); }
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
        for i in 3..N { assert_eq!(c[i], 0, "nonzero at {}", i); }
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
        for i in 1..N { assert_eq!(c[i], 0, "nonzero at {}", i); }
    }
}
