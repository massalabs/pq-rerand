//! Scheme parameters and precomputed constants.

/// Ring dimension (power of 2).
pub const N: usize = 4096;

/// Plaintext modulus. Prime, ≡ 1 (mod 8192). Fits in u32.
/// t = 262143 * 8192 + 1 = 2_147_565_569
pub const T: u32 = 2_147_565_569;

/// Second CRT prime (noise-limb modulus). Prime, ≡ 1 (mod 8192). Fits in u32.
/// q2 = 524413 * 8192 + 1 = 4_294_828_033
pub const Q2: u32 = 4_294_828_033;

/// Scaling factor residue: Δ_t = q₂ mod t. The message embedding in the t-limb
/// uses Δ_t · M (mod t). Precomputed.
pub const DELTA_T: u32 = {
    // q2 mod t = 4_294_828_033 mod 2_147_565_569 = 2_147_262_464
    (Q2 as u64 % T as u64) as u32
};

/// Modular inverse of DELTA_T mod t, computed via Fermat's little theorem.
/// INV_DELTA_T = DELTA_T^(t-2) mod t. Used in decryption.
pub const INV_DELTA_T: u32 = 1_987_184_532;

/// Narrow Gaussian standard deviation (for keygen, rerand, base encrypt).
pub const SIGMA: f64 = 3.2;

/// Default re-randomization budget.
pub const K_MAX: u64 = 858_000_000;

/// Flood Gaussian standard deviation: σ · √k_max.
pub const SIGMA_FLOOD: f64 = SIGMA * 29_292.4; // √858_000_000 ≈ 29292.4, so σ_flood ≈ 93735.7

/// Number of slots for 15.5 MiB plaintext.
pub const NUM_SLOTS: usize = 1024;

/// Bytes of plaintext per slot (4096 coefficients × 31 bits / 8).
pub const SLOT_BYTES: usize = N * 31 / 8; // 15_872

/// Total plaintext capacity in bytes.
pub const PLAINTEXT_BYTES: usize = NUM_SLOTS * SLOT_BYTES; // 16_252_928

/// Bits per plaintext coefficient.
pub const BITS_PER_COEFF: usize = 31;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_t() {
        assert_eq!(DELTA_T, (Q2 as u64 % T as u64) as u32);
        assert_eq!(DELTA_T, 2_147_262_464);
    }

    #[test]
    fn test_inv_delta_t() {
        let product = (DELTA_T as u64 * INV_DELTA_T as u64) % T as u64;
        assert_eq!(product, 1, "INV_DELTA_T is not the inverse of DELTA_T mod T");
    }

    #[test]
    fn test_primes_ntt_friendly() {
        assert_eq!(T as u64 % (2 * N as u64), 1, "T is not ≡ 1 (mod 2N)");
        assert_eq!(Q2 as u64 % (2 * N as u64), 1, "Q2 is not ≡ 1 (mod 2N)");
    }

    #[test]
    fn test_slot_bytes() {
        assert_eq!(SLOT_BYTES, 15_872);
        assert_eq!(PLAINTEXT_BYTES, 16_252_928);
    }

    #[test]
    fn test_sigma_flood() {
        let expected = SIGMA * (K_MAX as f64).sqrt();
        let diff = (SIGMA_FLOOD - expected).abs();
        assert!(diff < 10.0, "SIGMA_FLOOD is off: {} vs {}", SIGMA_FLOOD, expected);
    }
}
