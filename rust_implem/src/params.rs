//! Scheme parameters and precomputed constants.

/// Ring dimension (power of 2).
pub const N: usize = 4096;

/// Plaintext modulus. Prime, ≡ 1 (mod 8192). Fits in u32.
/// t = 262154 * 8192 + 1 = 2_147_565_569
pub const T: u32 = 2_147_565_569;

/// Second CRT prime (noise-limb modulus). Prime, ≡ 1 (mod 8192). Fits in u32.
/// q2 = 524271 * 8192 + 1 = 4_294_828_033
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
/// Matches the paper's σ = 3.2.
pub const SIGMA: f64 = 3.2;

/// Flood "virtual re-randomization count" κ_f (paper §6). The flood packet is
/// equivalent, for the noise budget, to κ_f extra re-randomizations.
pub const KAPPA_F: u64 = 858_000_000;

/// Flood Gaussian standard deviation σ_f = σ · √κ_f ≈ 93,733 (paper §6).
/// 3.2 * sqrt(858_000_000) = 3.2 * 29291.6370... = 93733.2385...
pub const SIGMA_FLOOD: f64 = 93_733.238_501_611_59;

/// Noise-budget threshold: a coefficient decrypts correctly iff |ν| < q₂/2
/// (paper Theorem 4.3). With q₂ ≈ 2³², q₂/2 ≈ 2.147 × 10⁹.
pub const NOISE_THRESHOLD: u32 = Q2 / 2;

/// Conservative limit on the number of *additional* public re-randomizations of
/// one ciphertext on top of the flood baseline. The fixed-key batch bound (paper
/// Theorem 4.6, Remark 4.7 and Table 7; p_batch = 2^-106, B = 1024, good key)
/// gives k_max = 14,151,421,672 with flooding and 15,009,421,672 without; this
/// constant rounds the former down. Nothing in this crate counts
/// re-randomizations: see [`rerandomize_slot`](crate::rerandomize::rerandomize_slot).
pub const K_MAX_RERAND: u64 = 14_000_000_000;

/// Number of slots for 15.5 MiB plaintext.
pub const NUM_SLOTS: usize = 1024;

/// Bytes of plaintext per slot (N coefficients of BITS_PER_COEFF bits each).
pub const SLOT_BYTES: usize = N * BITS_PER_COEFF / 8; // 15_872
const _: () = assert!((N * BITS_PER_COEFF).is_multiple_of(8), "a slot must hold whole bytes");

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
        // σ_f must equal σ·√κ_f to match the paper's flood calibration.
        let expected = SIGMA * (KAPPA_F as f64).sqrt();
        let diff = (SIGMA_FLOOD - expected).abs();
        assert!(diff < 1e-6, "SIGMA_FLOOD is off: {} vs {}", SIGMA_FLOOD, expected);
    }

    #[test]
    fn test_noise_threshold() {
        assert_eq!(NOISE_THRESHOLD, Q2 / 2);
    }
}
