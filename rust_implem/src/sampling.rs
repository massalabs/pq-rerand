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
//!   (>10^5 entries), so [`add_flood`] uses **constant-time inverse-transform
//!   sampling**: a 53-bit uniform (a sign and a 52-bit mantissa) is mapped through a
//!   branchless rational approximation of the probit (`inv_norm_cdf`) and rounded.
//!   There is no rejection loop, no data-dependent branch and no table indexing, so
//!   control flow and memory-access patterns are independent of the sampled value.
//!   (The flood triple is part of the encryption randomness: it is as secret as the
//!   base triple and reveals the plaintext just the same.)
//!
//!   The output is a rounded Gaussian with 53 bits of resolution, not an exact
//!   discrete Gaussian: it is symmetric and bounded by `8.3 * sigma_f`, and its far
//!   tail is quantized by the 53-bit input (beyond about `7 * sigma_f` an output
//!   value is produced by at most a few inputs). Security does not rely on its
//!   shape: flooding preserves pseudorandomness for any flood distribution (paper
//!   Lemma 6.1).
//!
//! # Use of the generator
//!
//! The samplers read the generator only through `next_u32` and `next_u64`, in the
//! fixed order documented on each function. Their output for a given generator
//! stream is therefore determined by this crate alone: it does not vary with the
//! platform or with the version of the `rand` crates. This is what makes a key pair
//! derived from a seeded generator reproducible.

use crate::params::{N, SIGMA_FLOOD};
use core::hint::select_unpredictable;
use rand_core::{CryptoRng, Rng};
use zeroize::Zeroizing;

/// Run `draw`, which reads a generator, then overwrite the stack that it used,
/// also when `draw` unwinds.
///
/// A generator computes on the stack: a ChaCha implementation, for one, spills
/// words of its key there, and that key reproduces everything the generator
/// returned, long after the generator itself is gone. `draw` therefore runs in a
/// stack frame of its own, which is then overwritten together with those of the
/// functions that it called. This is best effort: it overwrites the 8 KiB below
/// the frame of the caller, many times what the samplers use with a ChaCha
/// generator.
pub(crate) fn wiped<T>(draw: impl FnOnce() -> T) -> T {
    #[inline(never)]
    fn framed<T>(draw: impl FnOnce() -> T) -> T {
        draw()
    }
    #[inline(never)]
    fn wipe() {
        core::hint::black_box(&mut [0u8; 8 * 1024]);
    }
    struct Wipe;
    impl Drop for Wipe {
        fn drop(&mut self) {
            wipe();
        }
    }
    let _wipe = Wipe;
    framed(draw)
}

/// The generator is read a few words at a time, ahead of the arithmetic on them:
/// with the generator out of the way, the compiler can vectorize that arithmetic.
const CHUNK: usize = 16;
const _: () = assert!(N.is_multiple_of(CHUNK));

/// Replace every coefficient `c` of `out` by `sample(word, c)`, for consecutive
/// `next_u64` words of `rng`.
#[inline(always)]
fn sample_with<R: Rng + CryptoRng>(rng: &mut R, out: &mut [i32; N], sample: impl Fn(u64, i32) -> i32) {
    let mut draws = Zeroizing::new([0u64; CHUNK]);
    for chunk in out.as_chunks_mut::<CHUNK>().0 {
        draws.iter_mut().for_each(|draw| *draw = rng.next_u64());
        for (c, &draw) in chunk.iter_mut().zip(draws.iter()) {
            *c = sample(draw, *c);
        }
    }
}

/// Cumulative distribution table for the centered discrete Gaussian `D_{Z,sigma}`
/// at `sigma = SIGMA` (16/5 exactly), folded onto non-negative magnitudes: the
/// magnitude `k` has probability `(CDT[k] - CDT[k-1]) / 2^63`, with `CDT[-1] = 0`
/// and `CDT[29] = 2^63`.
///
/// Every magnitude `k >= 1` has its probability in `D_{Z,sigma}` rounded to the
/// nearest multiple of 2^-63, and zero takes what remains. The probability of 30
/// rounds to zero (the tail beyond 29 has mass 2^-65.3), so 29 is the largest
/// magnitude. Rounding the cumulative values instead would leave two rounding
/// errors in the probability of every nonzero magnitude, and a distribution that
/// is further from `D_{Z,sigma}`. No entry is below the exact cumulative value:
/// no tail `|x| > k` is more likely than in `D_{Z,sigma}`, and a sample is
/// sub-Gaussian with parameter `sigma`, like an exact one.
///
/// The table is a constant computed in exact arithmetic: it involves no floating
/// point and is the same on every platform. The test
/// `cdt_matches_exact_distribution` recomputes it. Per sample, the sampled
/// distribution is within statistical distance 2^-60.96 of `D_{Z,sigma}`; its
/// Renyi divergence from it is below 1 + 2.65e-20 at order 2 and below
/// 1 + 1.19e-5 at order infinity.
const CDT: [u64; 29] = [
    1_149_872_835_429_266_009, 3_340_023_666_152_832_879, 5_231_742_854_224_525_756,
    6_713_673_034_491_318_535, 7_766_573_326_200_196_560, 8_445_050_402_542_556_635,
    8_841_576_285_654_612_685, 9_051_758_678_878_186_098, 9_152_802_451_769_415_981,
    9_196_859_074_767_746_707, 9_214_281_206_174_004_122, 9_220_529_764_022_708_443,
    9_222_562_339_745_873_208, 9_223_161_995_634_596_966, 9_223_322_447_917_711_090,
    9_223_361_386_320_111_734, 9_223_369_956_674_611_013, 9_223_371_667_508_612_692,
    9_223_371_977_254_386_297, 9_223_372_028_116_140_534, 9_223_372_035_690_845_300,
    9_223_372_036_713_969_872, 9_223_372_036_839_307_003, 9_223_372_036_853_232_779,
    9_223_372_036_854_636_069, 9_223_372_036_854_764_321, 9_223_372_036_854_774_952,
    9_223_372_036_854_775_751, 9_223_372_036_854_775_805,
];

/// Largest magnitude of a base-width sample (about 9 sigma).
pub(crate) const MAX_BASE_SAMPLE: u32 = CDT.len() as u32;

/// The base-width sample selected by one 64-bit draw: its top 63 bits are a
/// uniform `u` and its low bit is the sign.
#[inline(always)]
fn cdt_sample(draw: u64) -> i32 {
    let u = draw >> 1; // 63-bit uniform
    let mut mag = 0i32;
    for &threshold in CDT.iter() {
        // Unsigned comparison of two values below 2^63: nothing can overflow.
        mag += (u >= threshold) as i32;
    }
    // Uniform sign applied unconditionally (negating 0 is harmless, so no bias).
    select_unpredictable(draw & 1 == 1, -mag, mag)
}

/// Constant-time discrete-Gaussian sampler for the base width `σ = SIGMA`.
///
/// For each coefficient it draws a 63-bit uniform `u`, computes the magnitude as the
/// number of table thresholds `u` meets or exceeds (a fixed-length, branchless scan),
/// and applies an independent uniform sign. The control flow, loop count and table
/// access pattern do not depend on the sampled value.
///
/// Each coefficient consumes one `next_u64` (see `cdt_sample` for its bits).
pub fn sample_gaussian_cdt<R: Rng + CryptoRng>(rng: &mut R, out: &mut [i32; N]) {
    wiped(|| sample_with(rng, out, |draw, _| cdt_sample(draw)));
}

/// Smallest probability fed to the probit: `2^-54`, the midpoint of the first of
/// the `2^53` cells that partition `(0, 1)`.
const MIN_P: f64 = 1.0 / 18_014_398_509_481_984.0;

/// Branchless, constant-time natural logarithm for a normal positive `f64`.
///
/// Decomposes `x = m*2^e` by bit manipulation (`e` from the exponent field, `m` in
/// `[1,2)` from the mantissa) and evaluates `ln(m) = 2*atanh((m-1)/(m+1))` with a
/// fixed-degree odd series (`t <= 1/3`, truncation error `< 1.1e-13`). No branches
/// and no table lookups, so it does not inherit libm's data-dependent timing. Valid
/// for finite normal `x > 0`; this module only feeds it `p` in `[2^-54, 1/2)`.
#[inline(always)]
fn ct_ln(x: f64) -> f64 {
    const LN2: f64 = core::f64::consts::LN_2;
    debug_assert!((MIN_P..1.0).contains(&x), "ct_ln: input outside [2^-54, 1)");
    let bits = x.to_bits();
    // The exponent goes through i32: a 64-bit integer-to-float conversion is a
    // variable-time library call on 32-bit targets.
    let e = (((bits >> 52) & 0x7ff) as i32 - 1023) as f64;
    let m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let t = (m - 1.0) / (m + 1.0);
    let t2 = t * t;
    let t4 = t2 * t2;
    // 2t*(1 + t^2/3 + t^4/5 + ... + t^22/23): the terms of even and of odd rank
    // form two chains in t^4, half as long as one chain in t^2.
    let even = ((((1.0 / 21.0 * t4 + 1.0 / 17.0) * t4 + 1.0 / 13.0) * t4 + 1.0 / 9.0) * t4 + 1.0 / 5.0) * t4 + 1.0;
    let odd = ((((1.0 / 23.0 * t4 + 1.0 / 19.0) * t4 + 1.0 / 15.0) * t4 + 1.0 / 11.0) * t4 + 1.0 / 7.0) * t4 + 1.0 / 3.0;
    e * LN2 + 2.0 * t * (even + t2 * odd)
}

/// Branchless approximation of the standard-normal quantile `Phi^-1(p)` on the
/// lower half `p` in `(0, 1/2)`, where it is negative, via Acklam's rational
/// approximation with mask-selected regions (relative error `< 1.15e-9`).
///
/// The sampler obtains the upper half by symmetry, so `1 - p` is never formed and
/// `p` stays strictly inside the open interval. Both Acklam regions (lower tail /
/// central) are evaluated unconditionally and combined with a `0.0/1.0` mask
/// derived from the comparison, so there is no data-dependent control flow. `ln`
/// uses the branchless [`ct_ln`]. The square root and the three divisions are
/// single hardware instructions; whether their latency depends on the operands is
/// a property of the CPU, which `examples/leakage.rs` measures.
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

    // Region selector as 0.0/1.0 (the comparison compiles to setcc, not a branch).
    let is_low = (p < P_LOW) as u8 as f64;

    // Central region: rational in r = (p-0.5)².
    let q = p - 0.5;
    let r = q * q;
    let cen_num = (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q;
    let cen_den = ((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0;
    let cen = cen_num / cen_den;

    // Lower tail: rational in sqrt(-2 ln p).
    let qt = (-2.0 * ct_ln(p)).sqrt();
    let tail_num = ((((C[0] * qt + C[1]) * qt + C[2]) * qt + C[3]) * qt + C[4]) * qt + C[5];
    let tail_den = (((D[0] * qt + D[1]) * qt + D[2]) * qt + D[3]) * qt + 1.0;
    let tail = tail_num / tail_den;

    is_low * tail + (1.0 - is_low) * cen
}

/// The flood sample selected by one 64-bit draw: bits 12 to 63 are a 52-bit
/// mantissa `m`, bit 11 is the sign and the low 11 bits are unused.
///
/// The magnitude is `round(-SIGMA_FLOOD * Phi^-1(p))` with `p = (m + 1/2) / 2^53`,
/// which lies in `[2^-54, 1/2 - 2^-54]`: the output is exactly symmetric and no
/// sample exceeds `8.3 * SIGMA_FLOOD` in magnitude.
#[inline(always)]
fn flood_sample(draw: u64) -> i32 {
    // 2^-53, half the spacing of the doubles in [1, 2).
    const HALF_ULP: f64 = 1.0 / 9_007_199_254_740_992.0;
    // p from the bits of 1 + m/2^52. Every step is exact and there is no
    // integer-to-float conversion (a library call on some targets).
    let m = f64::from_bits(0x3ff0_0000_0000_0000 | (draw >> 12));
    let p = (m - 1.0 + HALF_ULP) * 0.5;
    // Phi^-1(p) < 0. Round the magnitude to nearest by adding 1/2 and truncating
    // through i32 (one instruction, where the library `round` and a 64-bit
    // conversion may branch on the value), then apply the sign.
    let mag = (0.5 - inv_norm_cdf(p) * SIGMA_FLOOD) as i32;
    select_unpredictable((draw >> 11) & 1 == 1, -mag, mag)
}

/// Add an independent flood-width sample to every coefficient of `acc`: a rounded
/// continuous Gaussian of standard deviation [`SIGMA_FLOOD`], by **constant-time
/// inverse-transform sampling** (see module docs and `inv_norm_cdf`).
///
/// Unlike a rejection-based (Box-Muller/ziggurat) sampler, every coefficient costs
/// exactly one uniform draw, one probit evaluation and one rounding, with no
/// data-dependent branch or memory indexing.
///
/// Each coefficient consumes one `next_u64` (see `flood_sample` for its bits). The
/// addition wraps, which the sums formed in this crate (below 2^20) never do.
pub fn add_flood<R: Rng + CryptoRng>(rng: &mut R, acc: &mut [i32; N]) {
    wiped(|| sample_with(rng, acc, |draw, c| c.wrapping_add(flood_sample(draw))));
}

/// Sample a polynomial uniformly in [0, modulus).
///
/// Rejection sampling: each coefficient is the first `next_u32` draw below
/// `modulus`. Its timing depends only on the public modulus and on the sampled
/// values, which are public themselves (they form the polynomial `a`).
pub(crate) fn sample_uniform<R: Rng + CryptoRng>(rng: &mut R, modulus: u32, out: &mut [u32; N]) {
    wiped(|| {
        for c in out.iter_mut() {
            *c = loop {
                let draw = rng.next_u32();
                if draw < modulus {
                    break draw;
                }
            };
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{Q2, SIGMA, T};
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use rand_core::{Infallible, TryCryptoRng, TryRng};

    /// Generator that replays a script of words in a loop, to drive a sampler at
    /// chosen points of its input space.
    struct Script(Vec<u64>, usize);

    impl TryRng for Script {
        type Error = Infallible;
        fn try_next_u32(&mut self) -> Result<u32, Infallible> {
            self.try_next_u64().map(|x| x as u32)
        }
        fn try_next_u64(&mut self) -> Result<u64, Infallible> {
            self.1 += 1;
            Ok(self.0[(self.1 - 1) % self.0.len()])
        }
        fn try_fill_bytes(&mut self, _: &mut [u8]) -> Result<(), Infallible> {
            unimplemented!("the samplers only draw whole words")
        }
    }

    impl TryCryptoRng for Script {}

    #[test]
    fn cdt_thresholds_monotonic() {
        for w in CDT.windows(2) {
            assert!(w[1] > w[0], "CDT thresholds must be increasing");
        }
        assert!(CDT[CDT.len() - 1] < 1 << 63, "a threshold of 2^63 can never be met");
    }

    /// Recompute the table from the definition of the discrete Gaussian, in
    /// fixed-point integer arithmetic with 100 fractional bits. The accumulated
    /// truncation error is below 2^-90, far less than the distance (0.0065 units of
    /// 2^-63) between any exact probability and the nearest rounding boundary.
    #[test]
    fn cdt_matches_exact_distribution() {
        const ONE: u128 = 1 << 100;
        // floor(a * b / 2^100), for a, b <= 2^100.
        fn mul(a: u128, b: u128) -> u128 {
            const LOW: u128 = (1 << 64) - 1;
            let (ah, al, bh, bl) = (a >> 64, a & LOW, b >> 64, b & LOW);
            ((ah * bh) << 28) + ((ah * bl + al * bh + ((al * bl) >> 64)) >> 36)
        }
        // rho(k) = exp(-k^2 / (2 sigma^2)) = q^(k^2), with q = exp(-25/512).
        assert!((2.0 * SIGMA * SIGMA - 512.0 / 25.0).abs() < 1e-12);
        let (mut term, mut q) = (ONE, ONE);
        for n in 1..40u128 {
            term = term * 25 / (512 * n);
            q = if n % 2 == 1 { q - term } else { q + term };
        }
        // Folded weights: w(0) = rho(0), w(k) = 2 rho(k).
        // rho(k) = rho(k-1) * q^(2k-1). Beyond k = 47 the weights are below 2^-160.
        let q2 = mul(q, q);
        let (mut rho, mut step) = (ONE, q);
        let mut weights = vec![ONE];
        for _ in 1..48 {
            rho = mul(rho, step);
            step = mul(step, q2);
            weights.push(2 * rho);
        }
        let norm: u128 = weights.iter().sum();
        // round(2^63 * w / norm), from the quotient with one more bit.
        let rounded = |w: u128| {
            let (mut rem, mut quo) = (w, 0u128);
            for _ in 0..4 {
                rem <<= 16;
                quo = (quo << 16) | (rem / norm);
                rem %= norm;
            }
            (quo + 1) >> 1
        };
        // Probabilities of the magnitudes 1, 2, ...: zero takes what remains.
        let probabilities: Vec<u128> = weights[1..].iter().map(|&w| rounded(w)).collect();
        let mut threshold = (1 << 63) - probabilities.iter().sum::<u128>();
        for (k, &entry) in CDT.iter().enumerate() {
            assert_eq!(entry as u128, threshold, "CDT[{k}]");
            threshold += probabilities[k];
        }
        // The table is complete: beyond its largest magnitude the probabilities
        // round to zero.
        assert_eq!(threshold, 1 << 63);
    }

    /// The scan returns the number of thresholds at or below the draw, at every
    /// threshold boundary, and the low bit of the draw selects the sign.
    #[test]
    fn cdt_sampler_boundaries() {
        assert_eq!(cdt_sample(0), 0);
        assert_eq!(cdt_sample(1), 0);
        // An all-ones word gives u = 2^63 - 1, the largest magnitude, and a minus sign.
        assert_eq!(cdt_sample(u64::MAX), -(MAX_BASE_SAMPLE as i32));
        assert_eq!(cdt_sample(u64::MAX - 1), MAX_BASE_SAMPLE as i32);
        for (k, &threshold) in CDT.iter().enumerate() {
            // u = CDT[k]: magnitude k + 1, plus sign (low bit clear).
            assert_eq!(cdt_sample(threshold << 1), k as i32 + 1, "u = CDT[{k}]");
            // u = CDT[k] - 1: magnitude k, minus sign (low bit set).
            assert_eq!(cdt_sample(((threshold - 1) << 1) | 1), -(k as i32), "u = CDT[{k}] - 1");
        }
        // The sampler maps each generator word to its sample, in order.
        let draws: Vec<u64> = (0..5 * CHUNK as u64 / 2).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        let mut out = [0i32; N];
        sample_gaussian_cdt(&mut Script(draws.clone(), 0), &mut out);
        for (i, &x) in out.iter().enumerate() {
            assert_eq!(x, cdt_sample(draws[i % draws.len()]), "coefficient {i}");
        }
    }

    #[test]
    fn ct_ln_matches_std() {
        // Across the documented domain [2^-54, 1), including both ends of a binade.
        let xs = [
            MIN_P, 5.6e-17, 1e-12, 1e-6, 0.01, 0.02425, 0.1, 0.25, 0.3, 0.499_999_999_999,
            0.5, 0.7, 0.97575, 0.999, 0.999_999,
        ];
        for &x in &xs {
            let err = (ct_ln(x) - x.ln()).abs();
            assert!(err < 1e-12, "ct_ln({x}) err = {err}");
        }
    }

    #[test]
    fn inv_norm_cdf_known_quantiles() {
        // (p, Phi^-1(p)) over the central region and the lower tail.
        let cases = [
            (0.25, -0.674_489_750_196_081_7),
            (0.1, -1.281_551_565_544_600_5),
            (0.025, -1.959_963_984_540_054_2),
            (0.01, -2.326_347_874_040_841),
            (0.001, -3.090_232_306_167_813_5),
            (1e-6, -4.753_424_308_822_899),
            (1e-10, -6.361_340_902_404_056),
            (MIN_P, -8.292_361_075_813_595),
        ];
        for (p, want) in cases {
            let got = inv_norm_cdf(p);
            assert!((got / want - 1.0).abs() < 2e-9, "quantile({p}) = {got}, want {want}");
        }
    }

    /// Flood sample for the 53-bit input `v`: its low bit is the sign, the other 52
    /// bits are the mantissa.
    fn flood_at(v: u64) -> i32 {
        flood_sample(v << 11)
    }

    /// The sampler at the ends and in the middle of its input space, against the
    /// rounded quantiles of the exact Gaussian. No input reaches p = 0 or p = 1.
    #[test]
    fn flood_sampler_edge_inputs() {
        const TOP: u64 = (1 << 53) - 1;
        let cases = [
            (0, 777_270),
            (1, -777_270),
            (2, 764_930),
            (4, 759_125),
            (6, 755_278),
            ((1 << 52) - 1, -63_222),
            (1 << 52, 63_222),
            (TOP - 1, 0),
            (TOP, 0),
        ];
        for (v, want) in cases {
            assert_eq!(flood_at(v), want, "input {v}");
        }
        // The unused low bits do not matter, and an all-ones word gives the
        // smallest magnitude.
        assert_eq!(flood_sample(0x7ff), 777_270);
        assert_eq!(flood_sample(u64::MAX), 0);
        // The sampler adds the sample of each generator word, in order, to what
        // the buffer holds.
        let draws: Vec<u64> = (0..5 * CHUNK as u64 / 2).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        let mut acc: [i32; N] = std::array::from_fn(|i| i as i32);
        add_flood(&mut Script(draws.clone(), 0), &mut acc);
        for (i, &x) in acc.iter().enumerate() {
            assert_eq!(x, i as i32 + flood_sample(draws[i % draws.len()]), "coefficient {i}");
        }
        let mut acc = [5i32; N];
        add_flood(&mut Script(vec![0, 1 << 11], 0), &mut acc);
        assert_eq!(acc[..4], [777_275, -777_265, 777_275, -777_265]);
    }

    /// The magnitude is non-increasing in the mantissa and the sign bit mirrors it:
    /// in the far tail, across the boundary between the two regions of the
    /// approximation (p = 0.02425), over the whole range and around the center.
    #[test]
    fn flood_sampler_is_monotonic() {
        let boundary = (0.02425 * (1u64 << 53) as f64) as u64;
        let mut mantissas: Vec<u64> = (0..100).collect();
        mantissas.extend((0..200).map(|j| boundary - (100 << 28) + (j << 28)));
        mantissas.extend((0..52).flat_map(|i| [1u64 << i, 3u64 << i >> 1]));
        mantissas.extend((0..200).map(|j| (1u64 << 52) - 1 - (j << 33)));
        mantissas.sort_unstable();
        let mut previous = i32::MAX;
        for m in mantissas {
            let x = flood_at(m << 1);
            assert_eq!(flood_at((m << 1) | 1), -x, "mantissa {m}: not symmetric");
            assert!((0..=previous).contains(&x), "mantissa {m}: {x} after {previous}");
            previous = x;
        }
    }

    /// Draws at or above the modulus are rejected, the others are kept in order.
    #[test]
    fn uniform_sampler_keeps_the_draws_below_the_modulus() {
        for p in [T, Q2] {
            let draws = [p, u32::MAX, p - 1, 0, 7, p + 1, 5].map(u64::from).to_vec();
            let mut out = [0u32; N];
            sample_uniform(&mut Script(draws, 0), p, &mut out);
            assert_eq!(out[..8], [p - 1, 0, 7, 5, p - 1, 0, 7, 5]);
        }
    }

    #[test]
    fn uniform_sampler_distribution() {
        let mut rng = StdRng::seed_from_u64(3);
        for p in [T, Q2] {
            let mut out = [0u32; N];
            sample_uniform(&mut rng, p, &mut out);
            assert!(out.iter().all(|&x| x < p));
            // Each of 16 equal buckets expects N/16 = 256 samples (standard deviation 15.5).
            let mut buckets = [0i32; 16];
            for &x in out.iter() {
                buckets[(x as u64 * 16 / p as u64) as usize] += 1;
            }
            assert!(buckets.iter().all(|&b| (b - 256).abs() < 80), "{buckets:?}");
        }
    }

    #[test]
    fn flood_sampler_moments() {
        let mut rng = StdRng::seed_from_u64(11);
        let mut count = 0u64;
        let mut sum = 0.0f64;
        let mut sumsq = 0.0f64;
        for _ in 0..200 {
            let mut samples = [0i32; N];
            add_flood(&mut rng, &mut samples);
            for &x in samples.iter() {
                sum += x as f64;
                sumsq += (x as f64) * (x as f64);
                count += 1;
            }
        }
        let mean = sum / count as f64;
        let var = sumsq / count as f64 - mean * mean;
        assert!(mean.abs() < 0.01 * SIGMA_FLOOD, "flood mean too far from 0: {mean}");
        assert!((var / (SIGMA_FLOOD * SIGMA_FLOOD) - 1.0).abs() < 0.01, "flood variance off: {var}");
    }

    #[test]
    fn cdt_sampler_moments() {
        let mut rng = StdRng::seed_from_u64(7);
        let mut count = 0u64;
        let mut sum = 0.0f64;
        let mut sumsq = 0.0f64;
        let mut max_abs = 0i32;
        for _ in 0..200 {
            let mut samples = [0i32; N];
            sample_gaussian_cdt(&mut rng, &mut samples);
            for &x in samples.iter() {
                sum += x as f64;
                sumsq += (x * x) as f64;
                max_abs = max_abs.max(x.abs());
                count += 1;
            }
        }
        let mean = sum / count as f64;
        let var = sumsq / count as f64 - mean * mean;
        // Mean ≈ 0, variance ≈ σ² = 10.24 for the discrete Gaussian.
        assert!(mean.abs() < 0.02, "mean too far from 0: {mean}");
        assert!((var - SIGMA * SIGMA).abs() < 0.1, "variance off: {var}");
        assert!(max_abs <= MAX_BASE_SAMPLE as i32, "magnitude exceeded table support");
    }
}
