//! Decryption using the noise-limb trick.
//!
//! Because Δ·M ≡ 0 (mod q₂), the q₂-limb of v = c0 - c1·s is pure noise.
//! Center-lift it to recover the signed noise, then subtract from the t-limb
//! and multiply by INV_DELTA_T to recover M.

use crate::params::{N, T, Q2, INV_DELTA_T, NOISE_THRESHOLD};
use crate::ntt::{barrett_mu, barrett_reduce, csub, pointwise_mul, submod};
use crate::poly::{NttContext, secret};
use crate::encoding::Plaintext;
use crate::keygen::SecretKey;
use crate::encrypt::CiphertextSlot;
use core::hint::select_unpredictable;

/// Decrypt a single ciphertext slot, returning the message polynomial
/// coefficients in [0, t).
///
/// The result is only the plaintext that was encrypted if the ciphertext is
/// authentic. For a ciphertext that someone crafted or modified it can be a
/// function of the secret key (for instance its negation): treat it as key
/// material until the payload it carries has been authenticated, and never
/// reveal it. It is returned in heap storage that is zeroized on drop, as are the
/// intermediate values. See the crate documentation on decryption outcomes.
pub fn decrypt_slot(
    ctx: &NttContext,
    sk: &SecretKey,
    ct: &CiphertextSlot,
) -> Plaintext {
    // Compute c1 * s in each limb via NTT. c1 is in coefficient domain, so it is
    // forward-transformed here; s is already cached in the NTT domain (see
    // `SecretKey`), saving two forward transforms per decryption.
    let (mut cs_t, mut cs_q2) = (secret(), secret());
    cs_t.copy_from_slice(&*ct.c1_t);
    cs_q2.copy_from_slice(&*ct.c1_q2);
    ctx.tables_t.forward(&mut cs_t);
    ctx.tables_q2.forward(&mut cs_q2);

    // c1*s in NTT domain
    pointwise_mul(&mut cs_t, &sk.s_ntt_t, T);
    pointwise_mul(&mut cs_q2, &sk.s_ntt_q2, Q2);

    // INTT to get c1*s in coefficient domain
    ctx.tables_t.inverse(&mut cs_t);
    ctx.tables_q2.inverse(&mut cs_q2);

    // Noise-limb trick: v_q2 is pure noise. The center-lift and final reduction
    // run on secret-dependent values, so they are done branchlessly (constant-time):
    // no data-dependent control flow and no data-indexed memory access.
    let mu_t = barrett_mu(T);
    let t_u = T as u64;
    let mut message = Plaintext(secret());
    for i in 0..N {
        // v = c0 - c1*s in each limb
        let v_t = submod(ct.c0_t[i], cs_t[i], T);
        let v = submod(ct.c0_q2[i], cs_q2[i], Q2);
        // The center-lift maps v to v - q2 iff v > q2/2.
        let noise = v as i64 - select_unpredictable(v > NOISE_THRESHOLD, Q2 as i64, 0); // in [-q2/2, q2/2]

        // M = (v_t - noise) * Delta_t^-1 mod t. The value (v_t - noise) is in (-T, 2T);
        // shift by +T into (0, 3T) and reduce with two branchless subtractions.
        let u = ((v_t as i64 - noise) + T as i64) as u64;
        let u = csub(csub(u, t_u), t_u);
        message.0[i] = barrett_reduce(u * INV_DELTA_T as u64, T, mu_t);
    }
    message
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::{decode, encode};
    use crate::encrypt::{Flooding, encrypt_slot};
    use crate::keygen::keygen;
    use crate::params::{DELTA_T, SIGMA, SIGMA_FLOOD, SLOT_BYTES};
    use crate::poly::NTT_CONTEXT;
    use crate::rerandomize::rerandomize_slot;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    /// The decryption noise of `ct`: the centered value of `c0 - c1*s` in the q2 limb.
    fn noise(sk: &SecretKey, ct: &CiphertextSlot) -> Vec<i64> {
        let mut cs = *ct.c1_q2;
        NTT_CONTEXT.tables_q2.forward(&mut cs);
        pointwise_mul(&mut cs, &sk.s_ntt_q2, Q2);
        NTT_CONTEXT.tables_q2.inverse(&mut cs);
        (0..N)
            .map(|i| {
                let v = submod(ct.c0_q2[i], cs[i], Q2) as i64;
                if v > NOISE_THRESHOLD as i64 { v - Q2 as i64 } else { v }
            })
            .collect()
    }

    fn std_dev(xs: &[i64]) -> f64 {
        (xs.iter().map(|&x| (x as f64) * (x as f64)).sum::<f64>() / xs.len() as f64).sqrt()
    }

    /// Decryption is exact while the noise magnitude is at most floor(q2/2), and
    /// off by one unit of the plaintext just beyond, on either side.
    #[test]
    fn noise_threshold_is_exact() {
        let (sk, _pk) = keygen(&mut StdRng::seed_from_u64(1), &NTT_CONTEXT);
        let data: [u8; SLOT_BYTES] = std::array::from_fn(|i| (i * 29 + 5) as u8);
        let message = encode(&data);
        let half = NOISE_THRESHOLD as i64;
        // (noise, error it causes on the plaintext coefficient)
        let cases = [
            (0, 0), (1, 0), (-1, 0),
            (half - 1, 0), (half, 0), (half + 1, 1),
            (-(half - 1), 0), (-half, 0), (-(half + 1), -1),
        ];
        // With c1 = 0, c0 is what decryption works on: Delta_t*M + noise in the t
        // limb and the noise alone in the q2 limb.
        let mut ct = CiphertextSlot::zeroed();
        for i in 0..N {
            let nu = cases[i % cases.len()].0;
            let embedded = (DELTA_T as u64 * message.0[i] as u64 % T as u64) as i64;
            ct.c0_t[i] = (embedded + nu).rem_euclid(T as i64) as u32;
            ct.c0_q2[i] = nu.rem_euclid(Q2 as i64) as u32;
        }
        let out = decrypt_slot(&NTT_CONTEXT, &sk, &ct);
        for i in 0..N {
            let want = (message.0[i] as i64 + cases[i % cases.len()].1).rem_euclid(T as i64);
            assert_eq!(out.0[i] as i64, want, "coefficient {i}, noise {}", cases[i % cases.len()].0);
        }
    }

    /// A fresh ciphertext carries noise of the base width without flooding, and of
    /// the flood width with it: sigma * sqrt(|e|^2 + |s|^2 + 1), about sigma * sqrt(2 N) * sigma.
    #[test]
    fn flooding_sets_the_noise_width() {
        let mut rng = StdRng::seed_from_u64(2);
        let (sk, pk) = keygen(&mut rng, &NTT_CONTEXT);
        let message = encode(&[0x5a; SLOT_BYTES]);
        let key_norm = SIGMA * (2.0 * N as f64).sqrt();
        for (flood, width) in [(Flooding::Off, SIGMA), (Flooding::Standard, SIGMA_FLOOD)] {
            let ct = encrypt_slot(&mut rng, &NTT_CONTEXT, &pk, &message, flood);
            let ratio = std_dev(&noise(&sk, &ct)) / (width * key_norm);
            assert!((0.9..1.1).contains(&ratio), "{flood:?}: noise is {ratio} times the expected width");
            let decoded = decode(&decrypt_slot(&NTT_CONTEXT, &sk, &ct)).expect("valid slot");
            assert!(decoded.iter().all(|&b| b == 0x5a));
        }
    }

    /// Re-randomization changes every component of the ciphertext beyond
    /// recognition, yet only adds noise of the base width.
    #[test]
    fn rerandomization_adds_narrow_noise() {
        let mut rng = StdRng::seed_from_u64(3);
        let (sk, pk) = keygen(&mut rng, &NTT_CONTEXT);
        let message = encode(&[0xa5; SLOT_BYTES]);
        let ct = encrypt_slot(&mut rng, &NTT_CONTEXT, &pk, &message, Flooding::Standard);
        let ct2 = rerandomize_slot(&mut rng, &NTT_CONTEXT, &pk, &ct);

        let added: Vec<i64> = noise(&sk, &ct2).iter().zip(noise(&sk, &ct)).map(|(a, b)| a - b).collect();
        let ratio = std_dev(&added) / (SIGMA * SIGMA * (2.0 * N as f64).sqrt());
        assert!((0.9..1.1).contains(&ratio), "added noise is {ratio} times the expected width");

        for (old, new, p) in [
            (&ct.c0_t, &ct2.c0_t, T), (&ct.c0_q2, &ct2.c0_q2, Q2),
            (&ct.c1_t, &ct2.c1_t, T), (&ct.c1_q2, &ct2.c1_q2, Q2),
        ] {
            // The difference is uniform modulo p: about half its values exceed p/2.
            let large = (0..N).filter(|&i| submod(new[i], old[i], p) > p / 2).count();
            assert!((N / 2 - 200..N / 2 + 200).contains(&large), "{large} large differences");
        }
    }
}
