//! Gaussian and uniform sampling.
//!
//! Two error widths are used by the scheme:
//!
//! * **Base width** (`σ = SIGMA`, used by keygen, base encryption and
//!   re-randomization): sampled with a **constant-time** cumulative-distribution-table
//!   (CDT) discrete Gaussian sampler ([`sample_gaussian_cdt`]). The per-coefficient
//!   lookup scans the whole table unconditionally with branchless comparisons, so its
//!   running time and memory-access pattern are independent of the secret value.
//! * **Flood width** (`σ_f ≈ 93 733`): a CDT is impractical at this width
//!   (>10^5 entries). Two samplers are provided, selected per encryption through
//!   [`FloodSampler`]:
//!   - [`sample_gaussian`] uses **constant-time inverse-transform sampling**: a
//!     53-bit uniform mantissa is mapped through a branchless rational
//!     approximation of the probit ([`inv_norm_cdf`]) and rounded. There is no
//!     rejection loop, no data-dependent branch and no table indexing, so timing
//!     and memory-access patterns are independent of the sampled value. Its
//!     output is a *rounded continuous* Gaussian up to the approximation error
//!     of the probit (relative error below `1.15e-9`), which is enough for the
//!     scheme's correctness and noise budget but is not an exactly specified
//!     distribution.
//!   - [`sample_gaussian_karney`] samples the **exact discrete Gaussian**
//!     `D_{Z,σ}`, `Pr[i] ∝ exp(-i²/(2σ²))`, for a rational `σ`, with Karney's
//!     algorithm (C. F. F. Karney, *Sampling exactly from the normal
//!     distribution*, ACM TOMS 42(1), 2016). It uses integer arithmetic and
//!     uniform random digits only, no floating point and no tables, and its
//!     output distribution is exactly `D_{Z,σ}` given an ideal random source.
//!     It is a rejection sampler and therefore **variable-time**; that is
//!     acceptable for the flood packet, which is public-domain randomness (an
//!     encryption of zero): learning its value through a timing channel would
//!     reduce the ciphertext to the base encryption, which is itself a secure
//!     Ring-LWE encryption, so the sampler's timing carries no secret. Where an
//!     exactly specified flood distribution matters for a security proof (as
//!     for the statistical "flooding" lemmas of deniable-storage applications),
//!     this is the sampler to use.
//!
//! Both flood samplers are tested against the exact probability mass function
//! (chi-square goodness of fit, Karney) and for their moments.

use crate::params::{N, SIGMA, SIGMA_FLOOD, SIGMA_FLOOD_DEN, SIGMA_FLOOD_NUM};
use rand::{CryptoRng, Rng};
use std::sync::OnceLock;

/// Maximum sampled magnitude for the base-width CDT (≈ 12.5·σ).
/// The discrete-Gaussian tail beyond this point has mass < 2^-63 and is therefore
/// below the precision of the 63-bit uniform draw used for the table lookup.
const CDT_LEN: usize = 40;

/// Cumulative distribution table for the centered discrete Gaussian `D_{Z,σ}` at
/// `σ = SIGMA`, folded onto non-negative magnitudes: `cdt[k] = floor(2^63 · P(|x| ≤ k))`.
///
/// The table is built once with floating point. This build is *data-independent*
/// (it takes no secret input), whereas the per-sample lookup in
/// [`sample_gaussian_cdt`] is constant-time.
fn cdt_table() -> &'static [u64; CDT_LEN] {
    static TABLE: OnceLock<[u64; CDT_LEN]> = OnceLock::new();
    TABLE.get_or_init(|| {
        let sigma = SIGMA;
        // Folded weights: w(0) = ρ(0), w(k) = 2·ρ(k) for k ≥ 1, ρ(k) = exp(-k²/2σ²).
        // Normalise over the full integer support (tail beyond CDT_LEN is negligible).
        let mut norm = 0.0f64;
        let mut w = [0.0f64; CDT_LEN];
        for (k, wk) in w.iter_mut().enumerate() {
            let rho = (-((k * k) as f64) / (2.0 * sigma * sigma)).exp();
            *wk = if k == 0 { rho } else { 2.0 * rho };
            norm += *wk;
        }
        let scale = (1u64 << 63) as f64;
        let mut cdt = [0u64; CDT_LEN];
        let mut acc = 0.0f64;
        for (k, ck) in cdt.iter_mut().enumerate() {
            acc += w[k] / norm;
            // acc ∈ (0, 1]; scaling by 2^63 keeps every threshold ≤ 2^63 so a
            // 63-bit uniform can always be compared against it.
            *ck = (acc * scale) as u64;
        }
        cdt
    })
}

/// Constant-time discrete-Gaussian sampler for the base width `σ = SIGMA`.
///
/// For each coefficient it draws a 63-bit uniform `u`, computes the magnitude as the
/// number of table thresholds `u` meets or exceeds (a fixed-length, branchless scan),
/// and applies an independent uniform sign. The control flow, loop count and table
/// access pattern do not depend on the sampled value.
pub fn sample_gaussian_cdt<R: Rng + CryptoRng>(rng: &mut R) -> [i64; N] {
    let cdt = cdt_table();
    let mut coeffs = [0i64; N];
    for c in coeffs.iter_mut() {
        let u = (rng.gen::<u64>() >> 1) as i64; // 63-bit uniform, always ≥ 0
        let mut mag: i64 = 0;
        for &threshold in cdt.iter() {
            // ge = 1 if u ≥ threshold else 0, branchless.
            // diff = u - threshold ∈ (-2^63, 2^63) fits in i64; sign bit of !diff is
            // 1 exactly when diff ≥ 0.
            let diff = u - threshold as i64;
            let ge = (((!diff) as u64) >> 63) as i64;
            mag += ge;
        }
        // Uniform sign applied unconditionally (negating 0 is harmless, so no bias).
        let sign = 1 - 2 * ((rng.gen::<u8>() & 1) as i64);
        *c = sign * mag;
    }
    coeffs
}

/// Branchless, constant-time natural logarithm for a normal positive `f64`.
///
/// Decomposes `x = m·2^e` by bit manipulation (`e` from the exponent field, `m ∈
/// [1,2)` from the mantissa) and evaluates `ln(m) = 2·atanh((m-1)/(m+1))` with a
/// fixed-degree odd series (`t ≤ 1/3`, truncation error `< 10⁻⁶`). No branches and
/// no table lookups, so it does not inherit libm's data-dependent timing. Valid for
/// finite normal `x > 0`, which is all this module feeds it (`p ∈ [2⁻⁵⁴, 1)`).
#[inline(always)]
fn ct_ln(x: f64) -> f64 {
    const LN2: f64 = core::f64::consts::LN_2;
    let bits = x.to_bits();
    let e = (((bits >> 52) & 0x7ff) as i64 - 1023) as f64;
    let m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let t = (m - 1.0) / (m + 1.0);
    let t2 = t * t;
    // 2t·(1 + t²/3 + t⁴/5 + t⁶/7 + t⁸/9 + t¹⁰/11); t ≤ 1/3 ⇒ error < 10⁻⁷.
    let s = 1.0
        + t2 * (1.0 / 3.0
            + t2 * (1.0 / 5.0 + t2 * (1.0 / 7.0 + t2 * (1.0 / 9.0 + t2 * (1.0 / 11.0)))));
    e * LN2 + 2.0 * t * s
}

/// Branchless approximation of the inverse standard-normal CDF (probit) `Φ⁻¹(p)`
/// for `p ∈ (0,1)`, via Acklam's rational approximation with mask-selected regions
/// (relative error `< 1.15·10⁻⁹`; the constant-time `ln` adds `< 10⁻⁶`).
///
/// All three Acklam regions (lower tail / central / upper tail) are evaluated
/// unconditionally and combined with `0.0/1.0` masks derived from the comparisons,
/// so there is no data-dependent control flow. `ln` uses the branchless [`ct_ln`]
/// and `sqrt` is a single hardware instruction, so the whole function is constant-time.
#[inline]
#[allow(clippy::excessive_precision)] // keep Acklam's published constants verbatim
fn inv_norm_cdf(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
        1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
        6.680131188771972e+01, -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
        -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const P_LOW: f64 = 0.02425;
    const P_HIGH: f64 = 1.0 - P_LOW;

    // Region selectors as 0.0/1.0 (comparisons compile to setcc, not a branch).
    let is_low = (p < P_LOW) as u8 as f64;
    let is_high = (p > P_HIGH) as u8 as f64;
    let is_cen = 1.0 - is_low - is_high;

    // Central region: rational in r = (p-0.5)².
    let q = p - 0.5;
    let r = q * q;
    let cen_num = (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q;
    let cen_den = ((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0;
    let cen = cen_num / cen_den;

    // Tail region: q = sqrt(-2·ln(pp)), pp = p (lower) or 1-p (upper / central filler).
    let pp = is_low * p + (1.0 - is_low) * (1.0 - p);
    let qt = (-2.0 * ct_ln(pp)).sqrt();
    let tail_num = ((((C[0] * qt + C[1]) * qt + C[2]) * qt + C[3]) * qt + C[4]) * qt + C[5];
    let tail_den = (((D[0] * qt + D[1]) * qt + D[2]) * qt + D[3]) * qt + 1.0;
    let tail = tail_num / tail_den;
    // Lower tail uses +tail (negative), upper tail uses −tail (positive).
    let tail_signed = tail * (is_low - is_high);

    is_cen * cen + tail_signed
}

/// Sample a polynomial with coefficients from a rounded continuous Gaussian, using
/// **constant-time inverse-transform sampling** (see module docs and
/// [`inv_norm_cdf`]).
///
/// Used for the flood width only. Unlike a rejection-based (Box–Muller/ziggurat)
/// sampler, every coefficient costs exactly one uniform draw, one probit evaluation
/// and one rounding, with no data-dependent branch or memory indexing.
pub fn sample_gaussian<R: Rng + CryptoRng>(rng: &mut R, sigma: f64) -> [i64; N] {
    // 2⁻⁵³, for mapping a 53-bit mantissa into (0,1).
    const INV_2_53: f64 = 1.0 / 9_007_199_254_740_992.0;
    let mut coeffs = [0i64; N];
    for c in coeffs.iter_mut() {
        // 53-bit uniform in (0,1); the +0.5 keeps p strictly inside the open interval.
        let m = (rng.gen::<u64>() >> 11) as f64;
        let p = (m + 0.5) * INV_2_53;
        *c = (inv_norm_cdf(p) * sigma).round() as i64;
    }
    coeffs
}

// ---------------------------------------------------------------------------
// Exact discrete Gaussian: Karney's algorithm
// ---------------------------------------------------------------------------

/// Greatest common divisor (Euclid), usable in `const` context.
const fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

/// Tail cut-off of the rejection sampler, in units of σ. The integer part `k`
/// of a normal deviate exceeds `K_MAX` with probability about `exp(-K_MAX²/2)`,
/// `10^-543` for `K_MAX = 50`: the algorithm never gets there, but every
/// intermediate value is sized so that it could without overflowing.
const K_MAX: u64 = 50;

/// Parameters of the exact discrete-Gaussian flood sampler
/// ([`sample_gaussian_karney`]): `σ = num/den` in lowest terms, validated at
/// construction so that sampling itself cannot fail.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KarneyParams {
    /// Numerator of σ (lowest terms).
    sig: u64,
    /// Denominator of σ (lowest terms).
    den: u64,
    /// `⌈σ⌉`.
    isig: u64,
}

impl KarneyParams {
    /// The scheme's flood width, `σ_f = SIGMA_FLOOD_NUM / SIGMA_FLOOD_DEN`
    /// (`= 93 733.236` exactly, see [`crate::params`]).
    pub const FLOOD: KarneyParams = match KarneyParams::new(SIGMA_FLOOD_NUM, SIGMA_FLOOD_DEN) {
        Ok(p) => p,
        Err(_) => panic!("SIGMA_FLOOD_NUM / SIGMA_FLOOD_DEN is not a valid Karney parameter"),
    };

    /// Parameters for `σ = sigma_num / sigma_den`.
    ///
    /// Both must be nonzero, and after reduction to lowest terms the numerator
    /// must be below `2^57`, so that every intermediate value of the algorithm
    /// (which sizes the tail at `K_MAX = 50` standard deviations) fits in 64
    /// bits. The scheme's flood width has a 25-bit numerator.
    pub const fn new(sigma_num: u64, sigma_den: u64) -> Result<Self, &'static str> {
        if sigma_num == 0 || sigma_den == 0 {
            return Err("KarneyParams: sigma must be positive");
        }
        let g = gcd_u64(sigma_num, sigma_den);
        let sig = sigma_num / g;
        let den = sigma_den / g;
        if sig >= (1u64 << 57) {
            return Err("KarneyParams: sigma numerator too large (>= 2^57 in lowest terms)");
        }
        // ⌈σ⌉. The largest result is about (K_MAX + 1)·⌈σ⌉ and the largest
        // intermediate product is i·den ≤ (K_MAX + 1)·sig + den; both must fit
        // (they do whenever sig < 2^57, but check explicitly).
        let isig = sig.div_ceil(den);
        let max_i = match isig.checked_mul(K_MAX + 1) {
            Some(v) => v,
            None => return Err("KarneyParams: sigma too large"),
        };
        if max_i.checked_mul(den).is_none() || sig.checked_mul(K_MAX + 1).is_none() {
            return Err("KarneyParams: sigma too large");
        }
        Ok(KarneyParams { sig, den, isig })
    }

    /// `σ` as a floating-point number (for reporting; the sampler never uses it).
    pub fn sigma(&self) -> f64 {
        self.sig as f64 / self.den as f64
    }

    /// `σ` as the reduced fraction `(numerator, denominator)`.
    pub fn fraction(&self) -> (u64, u64) {
        (self.sig, self.den)
    }
}

/// Number of entries of the precomputed `⌈k·σ⌉` table (k is the integer part
/// of the deviate in units of σ; `k ≥ 64` has probability below `10^-889`).
const CEIL_TABLE: usize = 64;

/// Random bytes fetched from the RNG at a time (eight ChaCha blocks for the
/// `rand` `StdRng`); one refill branch per 512 digits.
const RAND_BUF: usize = 512;

/// Streaming state of Karney's Algorithm D for one polynomial: the random
/// source, a buffer of random bytes, and the lazily extended base-256 digits
/// of the current uniform deviate (Karney's "u-rand").
///
/// Port of the reference implementation (`RandomLib::DiscreteNormal`,
/// Copyright (c) Charles Karney 2013, MIT/X11 license,
/// <https://randomlib.sourceforge.net/>) with changes that preserve
/// exactness: 8-bit digits read from a 512-byte random buffer (a tie, which
/// needs a further digit, has probability `2^-8` per comparison, and the
/// refill branch is taken once per 512 digits); the Bernoulli(e^{-1/2})
/// trial evaluated without data-dependent branches on a window of four
/// digits, with the digit-by-digit procedure as its exact continuation
/// ([`Self::exp_prob_h`]); the offset `j` and the sign drawn together from
/// one 32-bit chunk; a precomputed table of `⌈kσ⌉`; and the digits of `p/σ`
/// computed with a floating-point estimate that an exact integer step
/// corrects. Every change is checked against the reference formulation or
/// the exact probability mass function in the tests.
struct Karney<'a, R: Rng + CryptoRng> {
    p: KarneyParams,
    rng: &'a mut R,
    buf: [u8; RAND_BUF],
    pos: usize,
    /// Digits of the current uniform deviate `V`, most significant first
    /// (`v[..m]`; the vector only grows on ties, probability `2^-8` each).
    v: Vec<u8>,
    /// Number of digits of `V` sampled so far.
    m: usize,
    /// `⌈k·σ⌉` for `k < CEIL_TABLE`.
    ceil_k: [u64; CEIL_TABLE],
    /// `256 / sig` as a float, for the quotient estimate in `leading_digit`.
    inv_sig_256: f64,
}

impl<'a, R: Rng + CryptoRng> Karney<'a, R> {
    fn new(rng: &'a mut R, p: KarneyParams) -> Self {
        let mut ceil_k = [0u64; CEIL_TABLE];
        for (k, c) in ceil_k.iter_mut().enumerate() {
            *c = Self::ceil_div(p.sig * k as u64, p.den);
        }
        let mut st = Karney {
            p,
            rng,
            buf: [0u8; RAND_BUF],
            pos: RAND_BUF,
            v: Vec::with_capacity(16),
            m: 0,
            ceil_k,
            inv_sig_256: 256.0 / p.sig as f64,
        };
        st.refill();
        st
    }

    #[inline]
    fn ceil_div(n: u64, d: u64) -> u64 {
        n.div_ceil(d)
    }

    #[inline(never)]
    fn refill(&mut self) {
        self.rng.fill_bytes(&mut self.buf);
        self.pos = 0;
    }

    /// `K` consecutive uniform random bytes (`K ≤ RAND_BUF`).
    #[inline]
    fn bytes<const K: usize>(&mut self) -> [u8; K] {
        if self.pos + K > RAND_BUF {
            self.refill();
        }
        let mut out = [0u8; K];
        out.copy_from_slice(&self.buf[self.pos..self.pos + K]);
        self.pos += K;
        out
    }

    /// One uniform random digit (8 bits).
    #[inline]
    fn digit(&mut self) -> u8 {
        if self.pos == RAND_BUF {
            self.refill();
        }
        let d = self.buf[self.pos];
        self.pos += 1;
        d
    }

    /// One uniform 32-bit chunk.
    #[inline]
    fn chunk32(&mut self) -> u32 {
        u32::from_le_bytes(self.bytes::<4>())
    }

    /// Uniform integer in `[0, n)`, `0 < n ≤ 256`, from one digit: Lemire's
    /// multiply-shift with exact rejection of the biased residue classes.
    #[inline]
    fn uniform_below_u8(&mut self, n: u16) -> u8 {
        let mut m = u16::from(self.digit()) * n;
        let mut l = m as u8;
        if u16::from(l) < n {
            let t = (n.wrapping_neg() % n) as u8; // 256 mod n
            while l < t {
                m = u16::from(self.digit()) * n;
                l = m as u8;
            }
        }
        (m >> 8) as u8
    }

    /// Uniform integer in `[0, n)`, `0 < n ≤ 2^32`, from 32-bit chunks.
    #[inline]
    fn uniform_below_u32(&mut self, n: u32) -> u32 {
        let n64 = u64::from(n);
        let mut m = u64::from(self.chunk32()) * n64;
        let mut l = m as u32;
        if l < n {
            let t = n.wrapping_neg() % n; // 2^32 mod n
            while l < t {
                m = u64::from(self.chunk32()) * n64;
                l = m as u32;
            }
        }
        (m >> 32) as u32
    }

    /// Uniform integer in `[0, n)` for any `n > 0`.
    #[inline]
    fn uniform_below_u64(&mut self, n: u64) -> u64 {
        if n <= u64::from(u32::MAX) {
            return u64::from(self.uniform_below_u32(n as u32));
        }
        let n128 = u128::from(n);
        let mut m = u128::from(self.rng.next_u64()) * n128;
        let mut l = m as u64;
        if l < n {
            let t = n.wrapping_neg() % n; // 2^64 mod n
            while l < t {
                m = u128::from(self.rng.next_u64()) * n128;
                l = m as u64;
            }
        }
        (m >> 64) as u64
    }

    /// Append a digit to `V`.
    #[inline]
    fn push(&mut self, d: u8) {
        if self.m < self.v.len() {
            self.v[self.m] = d;
        } else {
            self.v.push(d);
        }
        self.m += 1;
    }

    /// Start a fresh deviate `V` with first digit `d`.
    #[inline]
    fn reset_v(&mut self, d: u8) {
        self.m = 0;
        self.push(d);
    }

    /// Next base-256 digit of `p/σ` (for `0 ≤ p < sig`), replacing `p` by the
    /// remainder: `p/sig = (digit + p'/sig) / 256`.
    ///
    /// For `sig < 2^52` the quotient `⌊256p/sig⌋` is estimated in floating
    /// point and corrected by at most one in exact integer arithmetic (the
    /// estimate is within one of the truth since `p < 2^52` is exact in an
    /// `f64` and the product's relative error is below `2^-50`), which avoids
    /// a 64-bit division on the hot path; the result is exact either way.
    #[inline]
    fn leading_digit(&self, p: &mut u64) -> u8 {
        let sig = self.p.sig;
        if sig < (1u64 << 52) {
            let num = *p << 8;
            let mut d = (*p as f64 * self.inv_sig_256) as u64;
            let mut r = num.wrapping_sub(d.wrapping_mul(sig));
            if r >= sig {
                // Estimate too low (r ≥ sig) or too high (r wrapped negative).
                if (r as i64) < 0 {
                    d -= 1;
                    r = r.wrapping_add(sig);
                } else {
                    d += 1;
                    r -= sig;
                }
            }
            *p = r;
            d as u8
        } else {
            let num = u128::from(*p) << 8;
            let d = num / u128::from(sig);
            *p = (num - d * u128::from(sig)) as u64;
            d as u8
        }
    }

    // --- u-rand comparisons ---------------------------------------------

    /// Draw a fresh uniform `V'`; return `V' < V`, and if so make `V'` the
    /// current deviate (the comparison extends either deviate lazily as needed).
    #[inline]
    fn less_than(&mut self) -> bool {
        self.less_than_from(0)
    }

    /// [`Self::less_than`] when the first `j` digits of `V'` are already known
    /// to equal those of `V` (used to resume a comparison after a tie).
    fn less_than_from(&mut self, j: usize) -> bool {
        let mut j = j;
        loop {
            if j == self.m {
                let d = self.digit();
                self.push(d);
            }
            let w = self.digit();
            let vj = self.v[j];
            if w > vj {
                return false;
            }
            if w < vj {
                self.v[j] = w;
                self.m = j + 1;
                return true;
            }
            j += 1;
        }
    }

    /// `V < x` for the rational `x = (x0 + p/sig) / 256`, drawing `V`'s digits
    /// lazily (its first digit if none has been drawn yet).
    fn less_than_x(&mut self, x0: u8, p: u64) -> bool {
        if self.m == 0 {
            let d = self.digit();
            self.push(d);
        }
        let v0 = self.v[0];
        if v0 != x0 {
            return v0 < x0;
        }
        self.less_than_x_tail(p)
    }

    /// [`Self::less_than_x`] after a tie on the first digit: compare the
    /// remaining digits of `V` (drawn lazily) with those of `p/sig`.
    fn less_than_x_tail(&mut self, mut p: u64) -> bool {
        let mut j = 1;
        loop {
            if p == 0 {
                return false;
            }
            if j == self.m {
                let d = self.digit();
                self.push(d);
            }
            let x = self.leading_digit(&mut p);
            let vj = self.v[j];
            if vj != x {
                return vj < x;
            }
            j += 1;
        }
    }

    /// Bernoulli trial: true with probability `x = (x0 + p/sig) / 256`.
    fn bernoulli(&mut self, x0: u8, p: u64) -> bool {
        let w = self.digit();
        if w != x0 {
            return w < x0;
        }
        self.bernoulli_tail(p)
    }

    /// [`Self::bernoulli`] after a tie on the first digit.
    fn bernoulli_tail(&mut self, mut p: u64) -> bool {
        loop {
            if p == 0 {
                return false;
            }
            let x = self.leading_digit(&mut p);
            let w = self.digit();
            if w != x {
                return w < x;
            }
        }
    }

    // --- Algorithm H: Bernoulli(exp(-1/2)) -------------------------------
    //
    // With uniform deviates U₁, U₂, …, let n be the length of the longest run
    // 1/2 > U₁ > U₂ > … > Uₙ; the trial is true iff n is even. The deviates
    // are compared digit by digit (base 256); a tie on a digit is resolved by
    // the next digits of both deviates.

    /// Reference, digit-by-digit evaluation of one trial (the formulation of
    /// the reference implementation; the differential test checks
    /// [`Self::exp_prob_h`] against it on identical random streams).
    #[cfg(test)]
    fn exp_prob_h_ref(&mut self) -> bool {
        let w = self.digit();
        self.reset_v(w);
        if w & 0x80 != 0 {
            return true; // U₁ ≥ 1/2: n = 0
        }
        let mut odd = false;
        loop {
            if !self.less_than() {
                return odd; // n - 1 successes; n even ⟺ successes odd
            }
            odd = !odd;
        }
    }

    /// Branch-free evaluation of one trial on a window of four leading digits.
    /// Returns `(result, needs_slow_path, successes, tie)`: the result is valid
    /// only if `needs_slow_path` is false, i.e. if the run ended within the
    /// window without a tie between two compared digits.
    #[inline]
    fn h_window(c: [u8; 4]) -> (bool, bool, u32, bool) {
        let (w0, w1, w2, w3) = (c[0], c[1], c[2], c[3]);
        let hi = w0 & 0x80 != 0;
        let lt1 = u32::from(w1 < w0);
        let lt2 = u32::from(w2 < w1);
        let lt3 = u32::from(w3 < w2);
        let eq1 = u32::from(w1 == w0);
        let eq2 = u32::from(w2 == w1);
        let eq3 = u32::from(w3 == w2);
        let r2 = lt1; // second comparison reached
        let r3 = lt1 & lt2; // third comparison reached
        let successes = lt1 + r2 * lt2 + r3 * lt3;
        let tie = (eq1 | (r2 & eq2) | (r3 & eq3)) != 0;
        let exhausted = (r3 & lt3) != 0;
        let slow = !hi & (tie | exhausted);
        (hi | (successes & 1 == 1), slow, successes, tie)
    }

    /// Exact continuation of a trial whose first four digits `c` were read
    /// by [`Self::h_window`]: `successes` comparisons succeeded before either
    /// a tie (`tie`) or the end of the window.
    #[cold]
    #[inline(never)]
    fn exp_prob_h_slow(&mut self, c: [u8; 4], successes: u32, tie: bool) -> bool {
        // The deviate at the end of the run has first digit c[successes]; on a
        // tie its challenger has the same first digit.
        let mut odd = successes & 1 == 1;
        self.reset_v(c[successes as usize]);
        if tie {
            if !self.less_than_from(1) {
                return odd;
            }
            odd = !odd;
        }
        loop {
            if !self.less_than() {
                return odd;
            }
            odd = !odd;
        }
    }

    /// One trial: the four-digit window, with the exact continuation on the
    /// rare occasions it does not decide. The continuation first gives back
    /// the window digits the reference procedure would not have consumed yet
    /// (a tie at comparison `successes + 1` has used `successes + 2` digits),
    /// so the trial is the same function of the digit stream as
    /// [`Self::exp_prob_h_ref`], not merely the same distribution.
    #[inline]
    fn exp_prob_h(&mut self) -> bool {
        let c = self.bytes::<4>();
        let (h, slow, successes, tie) = Self::h_window(c);
        if slow {
            let used = if tie { successes as usize + 2 } else { 4 };
            self.pos -= 4 - used; // bytes::<4> left pos ≥ 4, so this cannot underflow
            return self.exp_prob_h_slow(c, successes, tie);
        }
        h
    }

    /// Integer part `k ≥ 0` of the deviate, with probability
    /// `exp(-k/2)(1 - exp(-1/2))`: the number of consecutive true trials.
    #[inline]
    fn sample_k(&mut self) -> u32 {
        let mut k = 0u32;
        while self.exp_prob_h() {
            k += 1;
        }
        k
    }

    // --- Algorithm B: acceptance of the fractional part ------------------

    /// Algorithm C: `-1` with probability `1/m`, `0` with probability `1/m`,
    /// `+1` otherwise. `m = 2k + 2 ≤ 2·K_MAX + 2` in any realizable run; the
    /// wide path keeps the sampler exact even for an absurd `k`.
    #[inline]
    fn choose(&mut self, m: u32) -> i32 {
        let u = if m <= 256 {
            u32::from(self.uniform_below_u8(m as u16))
        } else {
            self.uniform_below_u32(m)
        };
        match u {
            0 => 0,
            1 => -1,
            _ => 1,
        }
    }

    /// Algorithm B: Bernoulli trial with probability
    /// `exp(-x (2k + x) / (2k + 2))`, `x = (x0 + p/sig) / 256`.
    fn b_trial(&mut self, k: u32, x0: u8, p: u64) -> bool {
        let h = 2 * k + 2;
        let mut n: u32 = 0;
        self.m = 0;
        loop {
            // For k = 0 the 3-way selector is a coin (−1 or 0, never +1) and
            // goes first, failing outright half of the time; for k > 0 it
            // goes after the run test.
            let mut f = if k == 0 {
                if self.digit() & 0x80 != 0 { -1 } else { 0 }
            } else {
                0
            };
            if f < 0 {
                break;
            }
            let lt = if n == 0 {
                self.less_than_x(x0, p)
            } else {
                self.less_than()
            };
            if !lt {
                break;
            }
            if k != 0 {
                f = self.choose(h);
                if f < 0 {
                    break;
                }
            }
            if f == 0 && !self.bernoulli(x0, p) {
                break;
            }
            n += 1;
        }
        n % 2 == 0
    }

    // --- Algorithm D -----------------------------------------------------

    /// One sample from `D_{Z,σ}`.
    fn sample(&mut self) -> i64 {
        let sig = self.p.sig;
        let den = self.p.den;
        loop {
            // Integer part k ≥ 0 with probability exp(-k/2)(1 - exp(-1/2)) ...
            let k = self.sample_k();
            // ... accepted with probability exp(-k(k-1)/2).
            let mut t = k * k.saturating_sub(1);
            let mut ok = true;
            while t > 0 {
                if !self.exp_prob_h() {
                    ok = false;
                    break;
                }
                t -= 1;
            }
            if !ok {
                continue;
            }
            // Sign and the offset j ∈ [0, ⌈σ⌉) from one draw.
            let u = self.uniform_below_u64(2 * self.p.isig);
            let neg = u & 1 == 1;
            let j = u >> 1;
            // x = (i - kσ)/σ with i = ⌈kσ⌉ + j, i.e. x = xn/sig with
            // xn = i·den - k·sig ≥ 0.
            let k_sig = sig * u64::from(k);
            let i0 = if (k as usize) < CEIL_TABLE {
                self.ceil_k[k as usize]
            } else {
                Self::ceil_div(k_sig, den)
            };
            let i = i0 + j;
            let xn = i * den - k_sig;
            if xn >= sig || (k == 0 && neg && xn == 0) {
                // x ≥ 1 is counted with the next k; −0 is not counted twice.
                continue;
            }
            if xn > 0 {
                // Accept x with probability exp(-x(2k+x)/2), as k + 1 trials
                // with probability exp(-x(2k+x)/(2k+2)) each.
                let mut p = xn;
                let x0 = self.leading_digit(&mut p);
                let mut h = k + 1;
                let mut acc = true;
                while h > 0 {
                    if !self.b_trial(k, x0, p) {
                        acc = false;
                        break;
                    }
                    h -= 1;
                }
                if !acc {
                    continue;
                }
            }
            let i = i as i64;
            return if neg { -i } else { i };
        }
    }
}

/// Sample a polynomial with coefficients from the **exact** discrete Gaussian
/// `D_{Z,σ}`, `Pr[i] ∝ exp(-i²/(2σ²))`, `σ = num/den` given by `params`, using
/// Karney's Algorithm D (see the module documentation).
///
/// Integer arithmetic and uniform random digits only; the output law is
/// exactly `D_{Z,σ}` for an ideal random source. **Variable-time** (rejection
/// sampling): use it only where the sampled values are not secret, as for the
/// flood packet (an encryption of zero).
pub fn sample_gaussian_karney<R: Rng + CryptoRng>(rng: &mut R, params: &KarneyParams) -> [i64; N] {
    let mut st = Karney::new(rng, *params);
    let mut coeffs = [0i64; N];
    for c in coeffs.iter_mut() {
        *c = st.sample();
    }
    coeffs
}

// ---------------------------------------------------------------------------
// Flood sampler selection
// ---------------------------------------------------------------------------

/// How the flood packet of an encryption is sampled, and with which width.
///
/// One value fully specifies the flood packet, so a width and a sampler can
/// never disagree. The published implementation is [`FloodSampler::INVERSE_CDF`]
/// (the default). [`FloodSampler::KARNEY`] samples the exact discrete Gaussian
/// at the same width; see the [module documentation](self) for the trade-off.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FloodSampler {
    /// No flood packet (`κ_f = 0`): the base encryption only.
    None,
    /// Constant-time inverse-CDF sampler ([`sample_gaussian`]): a rounded
    /// continuous Gaussian of standard deviation `sigma`, up to the probit
    /// approximation error.
    InverseCdf {
        /// Standard deviation of the flood packet.
        sigma: f64,
    },
    /// Exact discrete Gaussian `D_{Z,σ}` ([`sample_gaussian_karney`]),
    /// variable-time.
    Karney(KarneyParams),
}

impl FloodSampler {
    /// The published configuration: inverse-CDF sampler at `σ_f = SIGMA_FLOOD`.
    pub const INVERSE_CDF: FloodSampler = FloodSampler::InverseCdf { sigma: SIGMA_FLOOD };
    /// The exact sampler at `σ_f = SIGMA_FLOOD_NUM / SIGMA_FLOOD_DEN`.
    pub const KARNEY: FloodSampler = FloodSampler::Karney(KarneyParams::FLOOD);

    /// The flood packet's standard deviation (`0.0` for [`FloodSampler::None`]).
    pub fn sigma(&self) -> f64 {
        match self {
            FloodSampler::None => 0.0,
            FloodSampler::InverseCdf { sigma } => *sigma,
            FloodSampler::Karney(p) => p.sigma(),
        }
    }

    /// Sample one flood polynomial (all zeros for [`FloodSampler::None`]).
    pub fn sample<R: Rng + CryptoRng>(&self, rng: &mut R) -> [i64; N] {
        match self {
            FloodSampler::None => [0i64; N],
            FloodSampler::InverseCdf { sigma } => sample_gaussian(rng, *sigma),
            FloodSampler::Karney(p) => sample_gaussian_karney(rng, p),
        }
    }
}

impl Default for FloodSampler {
    fn default() -> Self {
        FloodSampler::INVERSE_CDF
    }
}

/// Sample a polynomial uniformly in [0, modulus).
///
/// Uses rejection sampling (`gen_range`); its timing depends only on the public
/// modulus and on fresh randomness, never on secret data (the sampled values are
/// the *public* polynomial `a`).
pub fn sample_uniform<R: Rng + CryptoRng>(rng: &mut R, modulus: u32) -> [u32; N] {
    let mut coeffs = [0u32; N];
    for c in coeffs.iter_mut() {
        *c = rng.gen_range(0..modulus);
    }
    coeffs
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn cdt_thresholds_monotonic() {
        let cdt = cdt_table();
        for w in cdt.windows(2) {
            assert!(w[1] >= w[0], "CDT thresholds must be non-decreasing");
        }
        // Last threshold should be essentially the full mass (≈ 2^63).
        assert!(cdt[CDT_LEN - 1] > (1u64 << 62), "CDT does not cover the full mass");
    }

    #[test]
    fn ct_ln_matches_std() {
        // Across the full input range the sampler can feed (p ∈ [2⁻⁵⁴, 1)).
        let xs = [
            5.6e-17, 1e-12, 1e-6, 0.01, 0.02425, 0.1, 0.3, 0.5, 0.7, 0.97575, 0.999, 0.999_999,
        ];
        for &x in &xs {
            let err = (ct_ln(x) - x.ln()).abs();
            assert!(err < 1e-6, "ct_ln({x}) err = {err}");
        }
    }

    #[test]
    fn inv_norm_cdf_known_quantiles() {
        // (p, expected Φ⁻¹(p)) — spans central and both tails.
        let cases = [
            (0.5, 0.0),
            (0.975, 1.959963985),
            (0.025, -1.959963985),
            (0.8413447461, 1.0),
            (0.1586552539, -1.0),
            (0.99, 2.326347874),
            (0.001, -3.090232306),
            (0.9999, 3.719016485),
        ];
        for (p, want) in cases {
            let got = inv_norm_cdf(p);
            assert!((got - want).abs() < 1e-4, "Φ⁻¹({p}) = {got}, want {want}");
        }
    }

    #[test]
    fn flood_sampler_moments() {
        let mut rng = StdRng::seed_from_u64(11);
        let sigma = 1000.0;
        let mut count = 0u64;
        let mut sum = 0.0f64;
        let mut sumsq = 0.0f64;
        for _ in 0..200 {
            for &x in sample_gaussian(&mut rng, sigma).iter() {
                sum += x as f64;
                sumsq += (x as f64) * (x as f64);
                count += 1;
            }
        }
        let mean = sum / count as f64;
        let var = sumsq / count as f64 - mean * mean;
        assert!(mean.abs() < 0.05 * sigma, "flood mean too far from 0: {mean}");
        assert!((var / (sigma * sigma) - 1.0).abs() < 0.05, "flood variance off: {var}");
    }

    /// Exact pmf of D_{Z,σ} on [-support, support] (σ as f64 is fine here:
    /// the test compares against the mathematical target, not the sampler).
    fn exact_pmf(sigma: f64, support: i64) -> Vec<f64> {
        let w: Vec<f64> = (-support..=support)
            .map(|i| (-(i as f64) * (i as f64) / (2.0 * sigma * sigma)).exp())
            .collect();
        let z: f64 = w.iter().sum();
        w.into_iter().map(|x| x / z).collect()
    }

    /// Chi-square goodness of fit of the Karney sampler against the exact
    /// discrete Gaussian, pooling tail bins with expected count below 20.
    /// Returns (statistic, degrees of freedom).
    fn karney_chi_square(num: u64, den: u64, samples: usize, seed: u64) -> (f64, usize) {
        let params = KarneyParams::new(num, den).unwrap();
        let sigma = params.sigma();
        let support = (12.0 * sigma).ceil() as i64 + 2;
        let pmf = exact_pmf(sigma, support);
        let mut counts = vec![0u64; pmf.len()];
        let mut rng = StdRng::seed_from_u64(seed);
        let mut st = Karney::new(&mut rng, params);
        for _ in 0..samples {
            let x = st.sample();
            assert!(x.abs() <= support, "sample {x} outside the 12σ support (σ={sigma})");
            counts[(x + support) as usize] += 1;
        }
        // Pool bins into cells with expected count ≥ 20, walking outward from
        // the centre so that pooling only merges tail bins.
        let n = samples as f64;
        let mut chi = 0.0f64;
        let mut cells = 0usize;
        let (mut acc_o, mut acc_e) = (0.0f64, 0.0f64);
        for (c, p) in counts.iter().zip(pmf.iter()) {
            acc_o += *c as f64;
            acc_e += p * n;
            if acc_e >= 20.0 {
                chi += (acc_o - acc_e).powi(2) / acc_e;
                cells += 1;
                acc_o = 0.0;
                acc_e = 0.0;
            }
        }
        if acc_e > 0.0 {
            chi += (acc_o - acc_e).powi(2) / acc_e;
            cells += 1;
        }
        (chi, cells - 1)
    }

    /// The chi-square statistic has mean df and standard deviation sqrt(2 df);
    /// df + 6 sqrt(2 df) is a p ~ 1e-9 threshold, far above any seed-to-seed
    /// variation yet sensitive to a wrong distribution (a 1% relative error on
    /// the central bins at 2e6 samples moves the statistic by hundreds).
    fn chi_threshold(df: usize) -> f64 {
        df as f64 + 6.0 * (2.0 * df as f64).sqrt()
    }

    /// An RNG that replays a fixed byte stream (test only).
    struct Scripted {
        bytes: Vec<u8>,
        pos: usize,
    }
    impl rand::RngCore for Scripted {
        fn next_u32(&mut self) -> u32 {
            let mut b = [0u8; 4];
            self.fill_bytes(&mut b);
            u32::from_le_bytes(b)
        }
        fn next_u64(&mut self) -> u64 {
            let mut b = [0u8; 8];
            self.fill_bytes(&mut b);
            u64::from_le_bytes(b)
        }
        fn fill_bytes(&mut self, dest: &mut [u8]) {
            for d in dest.iter_mut() {
                *d = self.bytes[self.pos];
                self.pos += 1;
            }
        }
        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
            self.fill_bytes(dest);
            Ok(())
        }
    }
    impl rand::CryptoRng for Scripted {}

    #[test]
    fn karney_window_trial_matches_reference_on_identical_streams() {
        // Both evaluations of the Bernoulli(exp(-1/2)) trial are deterministic
        // functions of the digit stream, so they must agree stream by stream,
        // including on ties and on runs longer than the window (together about
        // 0.5% of trials), which take the exact continuation.
        use rand::RngCore;
        let params = KarneyParams::FLOOD;
        let mut rng = StdRng::seed_from_u64(99);
        let (mut ties_or_long, mut trues) = (0u32, 0u32);
        const TRIALS: u32 = 400_000;
        for _ in 0..TRIALS {
            let mut stream = vec![0u8; RAND_BUF + 64];
            rng.fill_bytes(&mut stream);
            let mut a = Scripted { bytes: stream.clone(), pos: 0 };
            let mut b = Scripted { bytes: stream, pos: 0 };
            let mut fast = Karney::new(&mut a, params);
            let mut reference = Karney::new(&mut b, params);
            let c = [fast.buf[0], fast.buf[1], fast.buf[2], fast.buf[3]];
            let (_, slow, _, _) = Karney::<Scripted>::h_window(c);
            ties_or_long += u32::from(slow);
            let hf = fast.exp_prob_h();
            let hr = reference.exp_prob_h_ref();
            assert_eq!(hf, hr, "window trial disagrees with the reference on {c:?}");
            trues += u32::from(hf);
        }
        assert!(ties_or_long > 1_000, "slow path exercised {ties_or_long} times only");
        // Pr[true] = exp(-1/2) = 0.6065; 400k trials give a standard error of 0.0008.
        let p = trues as f64 / TRIALS as f64;
        assert!((p - 0.606_53).abs() < 0.005, "Bernoulli(exp(-1/2)) frequency {p}");
    }

    #[test]
    fn karney_params_validation() {
        assert!(KarneyParams::new(0, 1).is_err());
        assert!(KarneyParams::new(1, 0).is_err());
        assert!(KarneyParams::new(1 << 57, 1).is_err());
        let p = KarneyParams::new(14, 6).unwrap();
        assert_eq!(p.fraction(), (7, 3), "fraction must be reduced");
        assert_eq!(p.isig, 3, "⌈7/3⌉ = 3");
        assert!((p.sigma() - 7.0 / 3.0).abs() < 1e-12);
        let f = KarneyParams::FLOOD;
        assert_eq!(f.fraction(), (23_433_309, 250));
        assert_eq!(f.sigma(), SIGMA_FLOOD);
        assert_eq!(f.isig, 93_734);
    }

    #[test]
    fn karney_leading_digits_expand_fraction() {
        // The base-256 digits of p/sig must reconstruct p/sig.
        let params = KarneyParams::new(23_433_309, 250).unwrap();
        let mut rng = StdRng::seed_from_u64(3);
        let st = Karney::new(&mut rng, params);
        for &p0 in &[1u64, 2, 1_000_000, 23_433_308] {
            let mut p = p0;
            let mut acc = 0.0f64;
            let mut scale = 1.0f64 / 256.0;
            for _ in 0..8 {
                let d = st.leading_digit(&mut p);
                acc += d as f64 * scale;
                scale /= 256.0;
            }
            let target = p0 as f64 / 23_433_309.0;
            assert!((acc - target).abs() < 1e-15, "digits of {p0}/sig: {acc} vs {target}");
        }
    }

    #[test]
    fn karney_matches_exact_pmf_small_sigma() {
        // σ = 1 (heavy rejection), 7/3 (non-integer rational), 10 and 21/2:
        // the regime where every branch of the algorithm (k up to ~6, x = 0,
        // x ≥ 1, the −0 exclusion) is exercised many times.
        for (num, den, seed) in [(1u64, 1u64, 1u64), (7, 3, 2), (10, 1, 3), (21, 2, 4)] {
            let (chi, df) = karney_chi_square(num, den, 2_000_000, seed);
            assert!(
                chi < chi_threshold(df),
                "Karney sampler off the exact pmf for σ={num}/{den}: chi²={chi:.1}, df={df}"
            );
        }
    }

    #[test]
    fn karney_flood_width_moments_and_symmetry() {
        // At σ_f ≈ 93 733 the pmf is too spread for a per-value test; check the
        // first two moments, symmetry, and the tail, over 200 polynomials.
        let mut rng = StdRng::seed_from_u64(11);
        let params = KarneyParams::FLOOD;
        let sigma = params.sigma();
        let (mut n, mut sum, mut sumsq, mut pos) = (0u64, 0.0f64, 0.0f64, 0u64);
        let mut max_abs = 0i64;
        for _ in 0..200 {
            for &x in sample_gaussian_karney(&mut rng, &params).iter() {
                n += 1;
                sum += x as f64;
                sumsq += (x as f64) * (x as f64);
                pos += u64::from(x > 0);
                max_abs = max_abs.max(x.abs());
            }
        }
        let nf = n as f64;
        let mean = sum / nf;
        let var = sumsq / nf - mean * mean;
        // Standard errors: mean σ/√n ≈ 104, variance σ²√(2/n) ≈ 0.16%.
        assert!(mean.abs() < 5.0 * sigma / nf.sqrt(), "flood mean too far from 0: {mean}");
        assert!((var / (sigma * sigma) - 1.0).abs() < 0.01, "flood variance off: {var}");
        let frac_pos = pos as f64 / nf;
        assert!((frac_pos - 0.5).abs() < 0.005, "flood sign asymmetry: {frac_pos}");
        assert!(max_abs < (8.0 * sigma) as i64, "implausible tail value {max_abs}");
    }

    #[test]
    fn flood_sampler_dispatch() {
        let mut rng = StdRng::seed_from_u64(5);
        assert!(FloodSampler::None.sample(&mut rng).iter().all(|&x| x == 0));
        let a = FloodSampler::INVERSE_CDF.sample(&mut rng);
        let b = FloodSampler::KARNEY.sample(&mut rng);
        let var = |v: &[i64; N]| v.iter().map(|&x| (x as f64).powi(2)).sum::<f64>() / N as f64;
        for v in [&a, &b] {
            assert!((var(v) / (SIGMA_FLOOD * SIGMA_FLOOD) - 1.0).abs() < 0.1);
        }
    }

    #[test]
    fn cdt_sampler_moments() {
        let mut rng = StdRng::seed_from_u64(7);
        let mut count = 0u64;
        let mut sum = 0.0f64;
        let mut sumsq = 0.0f64;
        let mut max_abs = 0i64;
        for _ in 0..200 {
            for &x in sample_gaussian_cdt(&mut rng).iter() {
                sum += x as f64;
                sumsq += (x * x) as f64;
                max_abs = max_abs.max(x.abs());
                count += 1;
            }
        }
        let mean = sum / count as f64;
        let var = sumsq / count as f64 - mean * mean;
        // Mean ≈ 0, variance ≈ σ² = 10.24 for the discrete Gaussian.
        assert!(mean.abs() < 0.1, "mean too far from 0: {mean}");
        assert!((var - SIGMA * SIGMA).abs() < 1.0, "variance off: {var}");
        assert!(max_abs < CDT_LEN as i64, "magnitude exceeded table support");
    }
}
