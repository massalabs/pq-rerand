//! Decryption using the noise-limb trick.
//!
//! Because Δ·M ≡ 0 (mod q₂), the q₂-limb of v = c0 - c1·s is pure noise.
//! Center-lift it to recover the signed noise, then subtract from the t-limb
//! and multiply by INV_DELTA_T to recover M.

use crate::params::{N, T, Q2, INV_DELTA_T, NOISE_THRESHOLD, SLOT_BYTES};
use crate::ntt::{barrett_mu, barrett_reduce, csub, pointwise_mul, submod};
use crate::poly::{NTT_CONTEXT, secret};
use crate::encoding::{Plaintext, decode};
use crate::keygen::SecretKey;
use crate::encrypt::CiphertextSlot;
use core::hint::{black_box, select_unpredictable};
use zeroize::Zeroize;

/// Decrypt a single ciphertext slot and pass its payload to `open`, which
/// authenticates it. Returns what `open` returns, or `None` if the slot did not
/// decrypt to a payload.
///
/// `open` is the authenticated decryption of the application (see the crate
/// documentation): it must return `None` unless the payload authenticates, and
/// what it returns is all that leaves this function. The payload that it
/// receives is not authentic yet: for a ciphertext that someone crafted or
/// modified it can be a function of the secret key (for instance a constant
/// minus the key). `open` may work in the buffer that it is given, as an AEAD
/// that decrypts in place does: like the intermediate values, that buffer only
/// exists in heap storage that is zeroized before returning.
///
/// A slot that did not decrypt to a payload (wrong key, modified ciphertext,
/// exhausted noise budget) is not reported apart from one that `open` rejects:
/// `open` still runs, on a slot of zeros, and its result is discarded. The two
/// cases give the same `None` after the same work, provided that the time `open`
/// takes to reject does not depend on its input. That a call failed at all is
/// still an outcome to hide: see the crate documentation.
pub fn decrypt_slot<T>(
    sk: &SecretKey,
    ct: &CiphertextSlot,
    open: impl FnOnce(&mut [u8; SLOT_BYTES]) -> Option<T>,
) -> Option<T> {
    let mut payload = Payloads::new(1);
    let payload = &mut payload.slots()[0];
    let valid = decrypt_into(sk, ct, payload);
    // The verdict is read inside the closure, through `black_box`: only once
    // `open` has accepted. A compiler is otherwise free to test it before the
    // answer of `open`, in a decryption that fails.
    open(payload).filter(|_| black_box(valid))
}

/// Heap storage for the payloads of decrypted slots, zeroized on drop.
pub(crate) struct Payloads(Vec<u8>);

impl Payloads {
    /// Zeroed storage for `slots` payloads. It is allocated as bytes: a vector of
    /// slots can be built from a slot on the stack.
    pub(crate) fn new(slots: usize) -> Self {
        Payloads(vec![0; slots.checked_mul(SLOT_BYTES).expect("capacity overflow")])
    }

    /// The payloads, slot by slot.
    pub(crate) fn slots(&mut self) -> &mut [[u8; SLOT_BYTES]] {
        self.0.as_chunks_mut().0
    }
}

impl Drop for Payloads {
    fn drop(&mut self) {
        // As a slice, the form that `zeroize` clears fastest. The vector has no
        // spare capacity to clear.
        self.0.as_mut_slice().zeroize();
    }
}

/// Decrypt `ct` into `payload`. Returns whether it decrypted to a payload: if
/// not, `payload` is all zero (see [`decode`]).
pub(crate) fn decrypt_into(sk: &SecretKey, ct: &CiphertextSlot, payload: &mut [u8; SLOT_BYTES]) -> bool {
    decode(&decrypt_message(sk, ct), payload)
}

/// The message polynomial that `ct` decrypts to, with coefficients in [0, t).
/// It is returned in heap storage that is zeroized on drop, as are the
/// intermediate values.
fn decrypt_message(sk: &SecretKey, ct: &CiphertextSlot) -> Plaintext {
    let ctx = &NTT_CONTEXT;
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
    use crate::encoding::encode;
    use crate::encrypt::{Flooding, encrypt_slot};
    use crate::keygen::keygen;
    use crate::ntt::addmod;
    use crate::params::{DELTA_T, SIGMA, SIGMA_FLOOD};
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
        let (sk, _pk) = keygen(&mut StdRng::seed_from_u64(1));
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
        let out = decrypt_message(&sk, &ct);
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
        let (sk, pk) = keygen(&mut rng);
        let data = [0x5a; SLOT_BYTES];
        let key_norm = SIGMA * (2.0 * N as f64).sqrt();
        for (flood, width) in [(Flooding::Off, SIGMA), (Flooding::Standard, SIGMA_FLOOD)] {
            let ct = encrypt_slot(&mut rng, &pk, &data, flood);
            let ratio = std_dev(&noise(&sk, &ct)) / (width * key_norm);
            assert!((0.9..1.1).contains(&ratio), "{flood:?}: noise is {ratio} times the expected width");
            assert!(decrypt_slot(&sk, &ct, |payload| Some(*payload)) == Some(data));
        }
    }

    /// One failure: `open` rejecting the payload, or a slot that did not decrypt
    /// to one, which `open` then sees as zeros and which fails whatever `open`
    /// answers.
    #[test]
    fn decryption_has_one_failure() {
        let mut rng = StdRng::seed_from_u64(4);
        let (sk, pk) = keygen(&mut rng);
        let data: [u8; SLOT_BYTES] = std::array::from_fn(|i| (i * 31 + 7) as u8);
        let mut ct = encrypt_slot(&mut rng, &pk, &data, Flooding::Standard);
        // What `open` is given, and the outcome when it accepts.
        let accepted = |ct: &CiphertextSlot| {
            let mut given = [0xaa; SLOT_BYTES];
            let outcome = decrypt_slot(&sk, ct, |payload| {
                given = *payload;
                Some(())
            });
            (outcome, given)
        };
        assert!(accepted(&ct) == (Some(()), data));
        assert_eq!(decrypt_slot(&sk, &ct, |_| None::<()>), None);

        // The last coefficient of the message, m, is moved to t - 1, which is not
        // below 2^31, by adding Delta_t * (t - 1 - m) to c0 in the t limb.
        let m = encode(&data).0[N - 1];
        let shift = (DELTA_T as u64 * (T - 1 - m) as u64 % T as u64) as u32;
        ct.c0_t[N - 1] = addmod(ct.c0_t[N - 1], shift, T);
        assert!(accepted(&ct) == (None, [0; SLOT_BYTES]));
    }

    /// Re-randomization changes every component of the ciphertext beyond
    /// recognition, yet only adds noise of the base width.
    #[test]
    fn rerandomization_adds_narrow_noise() {
        let mut rng = StdRng::seed_from_u64(3);
        let (sk, pk) = keygen(&mut rng);
        let ct = encrypt_slot(&mut rng, &pk, &[0xa5; SLOT_BYTES], Flooding::Standard);
        let ct2 = rerandomize_slot(&mut rng, &pk, &ct);

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
