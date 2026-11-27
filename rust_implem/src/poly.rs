//! Storage and helpers for polynomials in CRT domain.
//!
//! A polynomial in R_q = Z_q[x]/(x^N+1) is stored as two CRT limbs:
//! one mod t and one mod q₂. All arithmetic is limbwise.

use crate::ntt::{NttTables, addmod};
use crate::params::{N, T, Q2, SIGMA};
use crate::sampling::MAX_BASE_SAMPLE;
use core::hint::select_unpredictable;
use core::ops::{Deref, DerefMut};
use zeroize::Zeroize;

/// Heap storage for `N` secret-dependent coefficients, zeroized on drop.
///
/// Secret arrays are only ever written in place through such a buffer: unlike
/// an array local or an array returned by value, it has a single location that
/// moves do not copy, so no copy is left behind in dead stack frames.
pub(crate) struct Secret<T: Zeroize = u32>(Box<[T; N]>);

/// `N` zero coefficients, allocated on the heap without a detour through the
/// stack.
pub(crate) fn zeroed<T: Copy + Default>() -> Box<[T; N]> {
    match vec![T::default(); N].into_boxed_slice().try_into() {
        Ok(coeffs) => coeffs,
        Err(_) => unreachable!(),
    }
}

/// A zeroed [`Secret`] buffer.
pub(crate) fn secret<T: Copy + Default + Zeroize>() -> Secret<T> {
    Secret(zeroed())
}

impl<T: Zeroize> Deref for Secret<T> {
    type Target = [T; N];
    fn deref(&self) -> &[T; N] {
        &self.0
    }
}

impl<T: Zeroize> DerefMut for Secret<T> {
    fn deref_mut(&mut self) -> &mut [T; N] {
        &mut self.0
    }
}

impl<T: Zeroize> Drop for Secret<T> {
    fn drop(&mut self) {
        self.0.zeroize();
    }
}

/// Branchless reduction of a small signed coefficient into `[0, m)`.
///
/// Requires `|c| < m`, which holds for every `i32` since both moduli exceed
/// 2^31. Adds `m` iff `c` is negative, with no data-dependent branch and no
/// division.
#[inline(always)]
fn reduce_signed_small(c: i32, m: u32) -> u32 {
    select_unpredictable(c < 0, (c as u32).wrapping_add(m), c as u32)
}

/// Write the CRT limbs of the small signed polynomial `e` to `(t, q2)`.
pub(crate) fn set_signed(t: &mut [u32; N], q2: &mut [u32; N], e: &[i32; N]) {
    for i in 0..N {
        t[i] = reduce_signed_small(e[i], T);
        q2[i] = reduce_signed_small(e[i], Q2);
    }
}

/// Add the small signed polynomial `e` to the CRT limbs `(t, q2)`.
pub(crate) fn add_signed(t: &mut [u32; N], q2: &mut [u32; N], e: &[i32; N]) {
    for i in 0..N {
        t[i] = addmod(t[i], reduce_signed_small(e[i], T), T);
        q2[i] = addmod(q2[i], reduce_signed_small(e[i], Q2), Q2);
    }
}

/// Largest squared norm of the secret and of the error of a good key, `2 N sigma^2`
/// (paper Lemma 4.4): the re-randomization budget is derived for such keys, and
/// key generation exceeds it with probability below 2^-900.
const MAX_SQUARED_NORM: u32 = (2.0 * N as f64 * SIGMA * SIGMA) as u32;

/// Whether `(t, q2)` are the reduced CRT limbs of one polynomial that is small as
/// the secret and the error of a key are: not zero, all coefficients of magnitude
/// at most [`MAX_BASE_SAMPLE`] and a squared norm of at most [`MAX_SQUARED_NORM`].
/// Scans every coefficient unconditionally: constant-time.
pub(crate) fn is_small(t: &[u32; N], q2: &[u32; N]) -> bool {
    let mut bad = false;
    let mut squared_norm = 0u32;
    for i in 0..N {
        // The only small value that t[i] can stand for: itself, or t[i] - T.
        let c = t[i].wrapping_sub(select_unpredictable(t[i] > MAX_BASE_SAMPLE, T, 0));
        bad |= (c as i32).unsigned_abs() > MAX_BASE_SAMPLE;
        bad |= t[i] != reduce_signed_small(c as i32, T);
        bad |= q2[i] != reduce_signed_small(c as i32, Q2);
        // Exact while the coefficients are within their bound, and the sum then
        // stays far below 2^32; meaningless (and ignored) otherwise.
        squared_norm = squared_norm.wrapping_add(c.wrapping_mul(c));
    }
    !bad & (1..=MAX_SQUARED_NORM).contains(&squared_norm)
}

/// Precomputed NTT context holding tables for both primes.
pub(crate) struct NttContext {
    /// Tables for the plaintext-limb prime `t`.
    pub(crate) tables_t: NttTables,
    /// Tables for the noise-limb prime `q₂`.
    pub(crate) tables_q2: NttTables,
}

/// Shared NTT context (twiddle tables for both primes). `NttContext::new` is a
/// `const fn`, so this static is built at compile time (no runtime table setup).
pub(crate) static NTT_CONTEXT: NttContext = NttContext::new();

impl NttContext {
    /// Build the context, at compile time (see [`NTT_CONTEXT`]).
    const fn new() -> Self {
        NttContext {
            tables_t: NttTables::new(T),
            tables_q2: NttTables::new(Q2),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_signed_small_matches_rem_euclid() {
        for &m in &[T, Q2] {
            for c in -1_000_000i32..=1_000_000 {
                assert_eq!(reduce_signed_small(c, m) as i64, (c as i64).rem_euclid(m as i64), "c={c} m={m}");
            }
            // the whole i32 range is within (-m, m)
            for &c in &[i32::MIN, i32::MIN + 1, -(1 << 30), 1 << 30, i32::MAX] {
                assert_eq!(reduce_signed_small(c, m) as i64, (c as i64).rem_euclid(m as i64), "c={c} m={m}");
            }
        }
    }

    #[test]
    fn test_set_and_add_signed() {
        let mut e = [0i32; N];
        e[0] = 100;
        e[1] = -300;
        e[2] = i32::MIN;
        let (mut t, mut q2) = ([7u32; N], [7u32; N]);
        set_signed(&mut t, &mut q2, &e);
        assert_eq!((t[0], q2[0]), (100, 100));
        assert_eq!((t[1], q2[1]), (T - 300, Q2 - 300));
        assert_eq!((t[2], q2[2]), (T - (1 << 31), Q2 - (1 << 31)));
        assert_eq!((t[3], q2[3]), (0, 0));

        // 100 + 100, -300 - 300, and a wrap past the modulus in both directions.
        add_signed(&mut t, &mut q2, &e);
        assert_eq!((t[0], q2[0]), (200, 200));
        assert_eq!((t[1], q2[1]), (T - 600, Q2 - 600));
        (t[4], q2[4]) = (T - 1, Q2 - 1);
        (t[5], q2[5]) = (1, 1);
        e[4] = 5;
        e[5] = -5;
        add_signed(&mut t, &mut q2, &e);
        assert_eq!((t[4], q2[4]), (4, 4));
        assert_eq!((t[5], q2[5]), (T - 4, Q2 - 4));
    }

    #[test]
    fn test_is_small() {
        let max = MAX_BASE_SAMPLE as i32;
        let mut e = [0i32; N];
        for (i, c) in (-max..=max).enumerate() {
            e[i] = c; // every value in [-max, max], once
        }
        let (mut t, mut q2) = ([0u32; N], [0u32; N]);
        set_signed(&mut t, &mut q2, &e);
        assert!(is_small(&t, &q2));

        // One coefficient just past the bound, on either side.
        for c in [max + 1, -max - 1] {
            e[N - 1] = c;
            set_signed(&mut t, &mut q2, &e);
            assert!(!is_small(&t, &q2), "{c} accepted");
        }
        e[N - 1] = 0;
        set_signed(&mut t, &mut q2, &e);

        // Limbs that disagree, at the first and at the last coefficient.
        for i in [0, N - 1] {
            let (mut bad_t, mut bad_q2) = (t, q2);
            bad_q2[i] = addmod(q2[i], 1, Q2);
            assert!(!is_small(&t, &bad_q2));
            bad_t[i] = addmod(t[i], 1, T);
            assert!(!is_small(&bad_t, &q2));
        }

        // Unreduced or large limbs.
        for (vt, vq2) in [(T, 0), (T + 5, 5), (5, Q2 + 5), (u32::MAX, u32::MAX), (T / 2, T / 2), (1 << 31, 1 << 31)] {
            let (mut bad_t, mut bad_q2) = (t, q2);
            (bad_t[7], bad_q2[7]) = (vt, vq2);
            assert!(!is_small(&bad_t, &bad_q2), "({vt}, {vq2}) accepted");
        }

        // Small coefficients but too large a norm: N * 5^2 exceeds 2 N sigma^2 and
        // N * 4^2 does not.
        assert_eq!(MAX_SQUARED_NORM, 83_886);
        for (c, small) in [(4, true), (-4, true), (5, false), (-5, false)] {
            set_signed(&mut t, &mut q2, &[c; N]);
            assert_eq!(is_small(&t, &q2), small, "all coefficients {c}");
        }
        // The bound itself is accepted and the next value is not:
        // 83_886 = 99 * 29^2 + 25^2 + 1^2 + 1^2.
        let mut e = [0i32; N];
        e[..99].fill(-max);
        (e[99], e[N - 2], e[N - 1]) = (25, 1, -1);
        set_signed(&mut t, &mut q2, &e);
        assert!(is_small(&t, &q2));
        e[N - 3] = 1;
        set_signed(&mut t, &mut q2, &e);
        assert!(!is_small(&t, &q2));

        // The zero polynomial is not a secret; a single unit coefficient is small.
        let mut e = [0i32; N];
        set_signed(&mut t, &mut q2, &e);
        assert!(!is_small(&t, &q2));
        e[N - 1] = -1;
        set_signed(&mut t, &mut q2, &e);
        assert!(is_small(&t, &q2));
    }
}
