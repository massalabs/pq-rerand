//! Gaussian and uniform sampling.
//!
//! ⚠️ NOT CONSTANT-TIME. Uses floating-point Box-Muller.
//! A production implementation needs a proper discrete Gaussian sampler.

use crate::params::N;
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Sample a polynomial with coefficients from a discrete Gaussian
/// approximated by rounding continuous Gaussian samples.
/// Returns signed integer coefficients.
pub fn sample_gaussian<R: Rng>(rng: &mut R, sigma: f64) -> [i64; N] {
    let normal = Normal::new(0.0, sigma).unwrap();
    let mut coeffs = [0i64; N];
    for i in 0..N {
        coeffs[i] = normal.sample(rng).round() as i64;
    }
    coeffs
}

/// Sample a polynomial uniformly in [0, modulus).
pub fn sample_uniform<R: Rng>(rng: &mut R, modulus: u32) -> [u32; N] {
    let mut coeffs = [0u32; N];
    for i in 0..N {
        coeffs[i] = rng.gen_range(0..modulus);
    }
    coeffs
}
