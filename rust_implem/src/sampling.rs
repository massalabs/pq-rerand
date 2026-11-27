//! Gaussian and uniform sampling.
//!
//! Two error widths are used by the scheme:
//!
//! * **Base width** (`σ = SIGMA`, used by keygen, base encryption and
//!   re-randomization): sampled with a **constant-time** cumulative-distribution-table
//!   (CDT) discrete Gaussian sampler ([`sample_gaussian_cdt`]). The per-coefficient
//!   lookup scans the whole table unconditionally with branchless comparisons, so its
//!   running time and memory-access pattern are independent of the secret value.
//! * **Flood width** (`sigma_f`, about 93,733): a table is impractical at this
//!   width (more than 10^5 entries), so [`add_flood`] assembles each sample from
//!   eight samples of a narrower table, as the digits of a number in base 4:
//!   `x_0 + 4 x_1 + ... + 4^7 x_7`. Like the base sampler, it only compares,
//!   selects and adds integers: there is no floating point, no rejection loop and
//!   no branch or memory address that depends on the data, so control flow and
//!   memory-access patterns are independent of the sampled value. (The flood
//!   triple is part of the encryption randomness: it is as secret as the base
//!   triple and reveals the plaintext just the same.)
//!
//!   A sum of independent discrete Gaussians is a discrete Gaussian again, up to
//!   a negligible factor, when each term is wide enough next to the step of the
//!   following one, as the digits are (see `flood_sample`). What the table of the
//!   digits rounds away remains: no value is more likely than under
//!   `D_{Z,sigma_f}` by a factor above 1 + 6.5e-8, values far in the tail are less
//!   likely, none exceeds `7.93 * sigma_f`, and a sample is sub-Gaussian with
//!   parameter `sigma_f`, like an exact one. Security does not rely on the shape:
//!   flooding preserves pseudorandomness for any flood distribution (paper
//!   Lemma 6.1).
//!
//! # Use of the generator
//!
//! The samplers read the generator only through `next_u32` and `next_u64`, in the
//! fixed order documented on each function. Their output for a given generator
//! stream is therefore determined by this crate alone: it does not vary with the
//! platform or with the version of the `rand` crates. This is what makes a key pair
//! derived from a seeded generator reproducible.

use crate::params::N;
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

/// The generator is read for a few coefficients at a time, ahead of the arithmetic
/// on its words: with the generator out of the way, the compiler can vectorize
/// that arithmetic.
const CHUNK: usize = 16;
const _: () = assert!(N.is_multiple_of(CHUNK));

/// Replace every coefficient `c` of `out` by `sample(draws, c)`, for consecutive
/// groups of `W` `next_u64` words of `rng`.
#[inline(always)]
fn sample_with<const W: usize, R: Rng + CryptoRng>(
    rng: &mut R,
    out: &mut [i32; N],
    sample: impl Fn(&[u64; W], i32) -> i32,
) {
    let mut draws = Zeroizing::new([[0u64; W]; CHUNK]);
    for chunk in out.as_chunks_mut::<CHUNK>().0 {
        draws.as_flattened_mut().iter_mut().for_each(|draw| *draw = rng.next_u64());
        for (c, draws) in chunk.iter_mut().zip(draws.iter()) {
            *c = sample(draws, *c);
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
/// `cdt_matches_exact_distribution` recomputes it, and so does
/// `tools/lean/cdt_table.lean`, which also checks the tails and measures the
/// distances that follow. Per sample, the sampled distribution is within
/// statistical distance 2^-60.96 of `D_{Z,sigma}`; its Renyi divergence from it is
/// below 1 + 2.65e-20 at order 2 and below 1 + 1.19e-5 at order infinity.
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
    wiped(|| sample_with(rng, out, |&[draw], _| cdt_sample(draw)));
}

/// Cumulative distribution table of one digit of a flood sample, laid out like
/// [`CDT`] with 31-bit thresholds: the magnitude `k` has probability
/// `(FLOOD_CDT[k] - FLOOD_CDT[k-1]) / 2^31`, with `FLOOD_CDT[-1] = 0` and
/// `FLOOD_CDT[34] = 2^31`. The digits follow the centered discrete Gaussian of
/// parameter `sigma_d = SIGMA_FLOOD / sqrt(S)`, about 5.54, where
/// `S = 1 + 16 + ... + 16^7 = 286,331,153` is the sum of the squares of their
/// weights `4^i`: the weighted sum then has the parameter `SIGMA_FLOOD`.
///
/// Every magnitude `k >= 1` has its probability in that distribution rounded
/// down to a multiple of 2^-31, and zero takes what remains: 34 is the largest
/// magnitude. No digit other than zero is more likely than in the discrete
/// Gaussian, zero is by a factor below 1 + 1.14e-7, and no tail `|x| > k` is
/// more likely: a digit is sub-Gaussian with parameter `sigma_d`, and a flood
/// sample with parameter `SIGMA_FLOOD`.
const FLOOD_CDT: [u32; 34] = [
    154_660_945, 458_983_279, 748_786_686, 1_015_914_794, 1_254_246_669, 1_460_068_347,
    1_632_115_122, 1_771_318_060, 1_880_335_682, 1_962_975_899, 2_023_612_191, 2_066_676_776,
    2_096_281_044, 2_115_979_620, 2_128_666_703, 2_136_575_950, 2_141_348_549, 2_144_136_090,
    2_145_712_009, 2_146_574_378, 2_147_031_149, 2_147_265_329, 2_147_381_540, 2_147_437_360,
    2_147_463_312, 2_147_474_991, 2_147_480_078, 2_147_482_223, 2_147_483_098, 2_147_483_443,
    2_147_483_575, 2_147_483_623, 2_147_483_640, 2_147_483_646,
];

/// The flood sample selected by four 64-bit draws: `x_0 + 4 x_1 + ... + 4^7 x_7`,
/// where the low half of draw `i` selects the digit `x_2i` and its high half
/// `x_(2i+1)`, each as a draw selects a sample in `cdt_sample`: the top 31 bits of
/// the half are a uniform and its low bit is the sign.
///
/// For independent discrete Gaussians `y`, `y'` of parameter `s` and an integer
/// `k` with `1 + k^2 <= s^2`, the sum `y + k y'` is a discrete Gaussian of
/// parameter `s sqrt(1 + k^2)` up to a factor that differs from 1 by less than
/// `2.1 exp(-2 pi^2 s^2 / (1 + k^2))`. The digits are combined in three such
/// steps, with `k` = 4, 16 and 256, and `s^2 / (1 + k^2)` is at least 1.8 in
/// each: over the three steps (the first applies to four pairs of digits, the
/// second to two pairs of sums) these factors stay within 3e-15 of 1. What
/// separates a flood sample from the discrete Gaussian is therefore the rounding
/// of the table, whose effect the test `flood_sample_distribution` computes from
/// the exact distribution of the sample.
#[inline(always)]
fn flood_sample(draws: &[u64; 4]) -> i32 {
    let halves: [u32; 8] = core::array::from_fn(|i| (draws[i / 2] >> (32 * (i % 2))) as u32);
    // The eight scans advance together, threshold by threshold. Both sides of a
    // comparison are below 2^31: comparing them as signed integers gives the same
    // result, and more targets have an instruction for it.
    let mut digits = [0i32; 8];
    for &threshold in FLOOD_CDT.iter() {
        for (digit, &half) in digits.iter_mut().zip(halves.iter()) {
            *digit += ((half >> 1) as i32 >= threshold as i32) as i32;
        }
    }
    // Horner's rule from the last digit down, each with its sign: the result is
    // at most 34 * (4^8 - 1) / 3 = 742,730 in magnitude.
    let signed = |i: usize| select_unpredictable(halves[i] & 1 == 1, -digits[i], digits[i]);
    (0..8).rev().fold(0, |x, i| 4 * x + signed(i))
}

/// Add an independent flood-width sample to every coefficient of `acc` (see the
/// module docs and `flood_sample`). Every coefficient costs the same four draws,
/// comparisons and additions, with no data-dependent branch or memory indexing.
///
/// Each coefficient consumes four `next_u64` (see `flood_sample` for their bits).
/// The addition wraps, which the sums formed in this crate (below 2^20) never do.
pub fn add_flood<R: Rng + CryptoRng>(rng: &mut R, acc: &mut [i32; N]) {
    wiped(|| sample_with(rng, acc, |draws, c| c.wrapping_add(flood_sample(draws))));
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
    use crate::params::{Q2, SIGMA, SIGMA_FLOOD, T};
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
    ///
    /// The weights `exp(-k^2 / (2 sigma^2))` are powers of `exp(-1 / (2 sigma^2))`,
    /// which is evaluated by its Taylor series. `tools/lean/cdt_table.lean` builds the
    /// same table in another way (512 fractional bits, one exponential per magnitude),
    /// then checks its tails and measures its distance from the discrete Gaussian.
    #[test]
    fn cdt_matches_exact_distribution() {
        const ONE: u128 = 1 << 100;
        // floor(a * b / 2^100), for a, b <= 2^100. The product may not fit in a u128:
        // from the 64-bit halves, it is ah*bh*2^128 + (ah*bl + al*bh)*2^64 + al*bl.
        fn mul(a: u128, b: u128) -> u128 {
            const LOW: u128 = (1 << 64) - 1;
            let (ah, al, bh, bl) = (a >> 64, a & LOW, b >> 64, b & LOW);
            ((ah * bh) << 28) + ((ah * bl + al * bh + ((al * bl) >> 64)) >> 36)
        }
        // rho(k) = exp(-k^2 / (2 sigma^2)) = q^(k^2), with q = exp(-25/512).
        assert!((2.0 * SIGMA * SIGMA - 512.0 / 25.0).abs() < 1e-12);
        // Taylor series: q = 1 - x + x^2/2! - x^3/3! + ... at x = 25/512, where `term`
        // is x^n / n!. It falls below 2^-100 at n = 15: the later terms add nothing.
        let (mut term, mut q) = (ONE, ONE);
        for n in 1..40u128 {
            term = term * 25 / (512 * n);
            q = if n % 2 == 1 { q - term } else { q + term };
        }
        // Folded weights: w(0) = rho(0), w(k) = 2 rho(k) for k and -k.
        // rho(k) = rho(k-1) * q^(2k-1), since k^2 - (k-1)^2 = 2k - 1: `step` is that
        // power of q. Beyond k = 47 the weights are below 2^-160.
        let q2 = mul(q, q);
        let (mut rho, mut step) = (ONE, q);
        let mut weights = vec![ONE];
        for _ in 1..48 {
            rho = mul(rho, step);
            step = mul(step, q2);
            weights.push(2 * rho);
        }
        // The probability of a magnitude is its weight over the sum of the weights.
        let norm: u128 = weights.iter().sum();
        // round(2^63 * w / norm), from the quotient with one more bit: the long
        // division gives floor(2^64 * w / norm), 16 bits at a time.
        let rounded = |w: u128| {
            let (mut rem, mut quo) = (w, 0u128);
            for _ in 0..4 {
                rem <<= 16;
                quo = (quo << 16) | (rem / norm);
                rem %= norm;
            }
            (quo + 1) >> 1
        };
        // Probabilities of the magnitudes 1, 2, ..., in units of 2^-63: zero takes
        // what remains.
        let probabilities: Vec<u128> = weights[1..].iter().map(|&w| rounded(w)).collect();
        // The thresholds are the cumulative sums: the first is the probability of zero.
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

    /// Recompute the table from the definition of the discrete Gaussian. Double
    /// precision is enough here: its error is below 1e-6 units of 2^-31, and no
    /// probability is within 5e-3 of a whole number of units.
    #[test]
    fn flood_cdt_matches_exact_distribution() {
        // rho(k) = exp(-k^2 / (2 sigma_d^2)), with sigma_d^2 = SIGMA_FLOOD^2 / S.
        let weights: f64 = (0..8).map(|i| 16f64.powi(i)).sum();
        let rho = |k: usize| (-((k * k) as f64) * weights / (2.0 * SIGMA_FLOOD * SIGMA_FLOOD)).exp();
        let norm = 1.0 + 2.0 * (1..200).map(rho).sum::<f64>();
        // Probabilities of the magnitudes 1, 2, ..., in units of 2^-31, rounded down
        // until one vanishes: zero takes what remains.
        let exact = |k: usize| 2.0 * rho(k) / norm * (1u64 << 31) as f64;
        let probabilities: Vec<u32> = (1..).map(|k| exact(k) as u32).take_while(|&p| p > 0).collect();
        assert!((1..=probabilities.len() + 1).all(|k| (exact(k) - exact(k).round()).abs() > 5e-3));
        let mut threshold = (1 << 31) - probabilities.iter().sum::<u32>();
        for (k, &entry) in FLOOD_CDT.iter().enumerate() {
            assert_eq!(entry, threshold, "FLOOD_CDT[{k}]");
            threshold += probabilities[k];
        }
        assert_eq!(threshold, 1 << 31);
        // Zero is more likely than in the discrete Gaussian, by less than 1.14e-7.
        let excess = FLOOD_CDT[0] as f64 / (1u64 << 31) as f64 * norm - 1.0;
        assert!((0.0..1.14e-7).contains(&excess), "{excess}");
    }

    /// The exact distribution of a flood sample, computed from the table, against
    /// the discrete Gaussian of parameter `SIGMA_FLOOD`. In double precision: the
    /// sums below only add positive terms, and are exact to 1e-11.
    #[test]
    fn flood_sample_distribution() {
        // One digit: the value d is at index d + 34.
        let top = FLOOD_CDT.len();
        let mut dist = vec![0.0f64; 2 * top + 1];
        let mut below = 0;
        for (k, threshold) in FLOOD_CDT.iter().copied().chain([1 << 31]).enumerate() {
            // Half of the probability of the magnitude k on each sign (both on 0).
            let half = (threshold - below) as f64 / (1u64 << 32) as f64;
            dist[top + k] += half;
            dist[top - k] += half;
            below = threshold;
        }
        // y + k y' for two independent copies of what precedes: digits by pairs,
        // then pairs by pairs, then the two halves of the sample.
        for k in [4, 16, 256] {
            let mut sum = vec![0.0; (k + 1) * (dist.len() - 1) + 1];
            for (j, &high) in dist.iter().enumerate() {
                for (to, &low) in sum[k * j..].iter_mut().zip(dist.iter()) {
                    *to += low * high;
                }
            }
            dist = sum;
        }
        let reach = (dist.len() / 2) as i32;
        assert_eq!(reach, 742_730);
        assert!((dist.iter().sum::<f64>() - 1.0).abs() < 1e-11);

        // The discrete Gaussian: at this width its normalization is that of the
        // continuous one, to more than 300 digits.
        let variance = SIGMA_FLOOD * SIGMA_FLOOD;
        let ideal = |x: i32| (-(x as f64).powi(2) / (2.0 * variance)).exp() / (2.0 * core::f64::consts::PI * variance).sqrt();
        let ratios = || (-reach..=reach).zip(dist.iter()).map(|(x, &p)| (x, p / ideal(x)));
        // No value is more likely than in the discrete Gaussian by a factor above
        // 1 + 6.5e-8, reached at zero. This bounds the Renyi divergence of order
        // infinity of the sample from the discrete Gaussian, which is what the
        // analysis of the flood packet uses of its condition on the two
        // distributions (paper section 6: a log-ratio within 7.2e-5).
        let (at, largest) = ratios().fold((0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
        assert!(at == 0 && largest < 1.0 + 6.5e-8, "ratio {largest} at {at}");
        // Within 4 sigma_f the two distributions agree to 1e-5 either way. Further
        // out the rounding of the table makes a value less likely, never more.
        let bulk = (4.0 * SIGMA_FLOOD) as i32;
        assert!(ratios().filter(|&(x, _)| x.abs() <= bulk).all(|(_, r)| (r - 1.0).abs() < 1e-5));
        // The variance is that of the discrete Gaussian to 2e-7: the rounding of the
        // table lowers it.
        let second: f64 = (-reach..=reach).zip(dist.iter()).map(|(x, &p)| (x as f64).powi(2) * p).sum();
        assert!((second / variance - 1.0).abs() < 2e-7, "{second}");
    }

    /// Each half of a draw selects its digit like `cdt_sample`, and the digits have
    /// the weights 1, 4, ..., 4^7 in the order of the halves.
    #[test]
    fn flood_sampler_digits() {
        assert_eq!(flood_sample(&[0; 4]), 0);
        assert_eq!(flood_sample(&[1 | 1 << 32; 4]), 0);
        // All-ones words: every digit is the largest magnitude, with a minus sign.
        let top = FLOOD_CDT.len() as i32;
        assert_eq!(flood_sample(&[u64::MAX; 4]), -top * 21_845);
        assert_eq!(flood_sample(&[u64::MAX - 1 - (1 << 32); 4]), top * 21_845);
        for position in 0..8 {
            let weight = 4i32.pow(position);
            let alone = |half: u32| {
                let mut draws = [0u64; 4];
                draws[position as usize / 2] = (half as u64) << (32 * (position % 2));
                flood_sample(&draws)
            };
            for (k, &threshold) in FLOOD_CDT.iter().enumerate() {
                // u = FLOOD_CDT[k]: magnitude k + 1, plus sign (low bit clear).
                assert_eq!(alone(threshold << 1), (k as i32 + 1) * weight, "digit {position}, u = FLOOD_CDT[{k}]");
                // u = FLOOD_CDT[k] - 1: magnitude k, minus sign (low bit set).
                assert_eq!(alone(((threshold - 1) << 1) | 1), -(k as i32) * weight, "digit {position}, u = FLOOD_CDT[{k}] - 1");
            }
        }
        // All eight digits at once, with different magnitudes and signs.
        let digits = [1, -6, -34, 13, 21, 2, -29, -10];
        let half = |d: i32| ((FLOOD_CDT[d.unsigned_abs() as usize - 1] << 1) | (d < 0) as u32) as u64;
        let draws: [u64; 4] = std::array::from_fn(|i| half(digits[2 * i]) | half(digits[2 * i + 1]) << 32);
        assert_eq!(flood_sample(&draws), (0..8).map(|i| digits[i] * 4i32.pow(i as u32)).sum::<i32>());
        // The sampler adds the sample of each group of four generator words, in
        // order, to what the buffer holds.
        let draws: Vec<u64> = (0..5 * 4 * CHUNK as u64 / 2).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        let mut acc: [i32; N] = std::array::from_fn(|i| i as i32);
        add_flood(&mut Script(draws.clone(), 0), &mut acc);
        let groups = draws.as_chunks::<4>().0;
        for (i, &x) in acc.iter().enumerate() {
            assert_eq!(x, i as i32 + flood_sample(&groups[i % groups.len()]), "coefficient {i}");
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
