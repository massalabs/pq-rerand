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
//!   (>10^5 entries), so [`sample_gaussian`] uses **constant-time inverse-transform
//!   sampling**: a 53-bit uniform mantissa is mapped through a branchless rational
//!   approximation of the probit ([`inv_norm_cdf`]) and rounded. There is no
//!   rejection loop, no data-dependent branch and no table indexing, so timing and
//!   memory-access patterns are independent of the sampled value. (This packet only
//!   ever processes public-domain randomness — an encryption of zero — and never
//!   secret-dependent data.)

use crate::params::{N, SIGMA};
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
