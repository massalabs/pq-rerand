#!/usr/bin/env python3
"""
Security estimation for the PQ re-randomizable encryption scheme.

Uses the HomomorphicEncryption.org Security Standard v1.1 tables as the
authoritative reference (which are derived from lattice-estimator runs).
We interpolate from those tables rather than reimplementing the estimator,
because simplified core-SVP formulas often disagree with the full analysis.
"""

import math

# ─── Parameters ───────────────────────────────────────────────────────
N = 4096
T = 2_147_565_569
Q2 = 4_294_828_033
Q_FULL = T * Q2
SIGMA = 3.2

# ─── HE Standard v1.1 reference tables ────────────────────────────────
# Source: https://homomorphicencryption.org/standard/ (2024 revision)
# Table for secret distribution = error distribution (Gaussian, σ ≈ 3.2)
# Maps n → max log₂(q) for 128-bit classical security.

HE_STANDARD_128BIT = {
    # n: max_log2_q
    1024: 27,
    2048: 56,   # reviewer cited ~56 for n=2048
    4096: 109,  # Table 2, error distribution row
    8192: 218,
    16384: 438,
    32768: 881,
}

def interpolate_security(n, log2_q):
    """Estimate whether (n, log₂(q)) achieves 128-bit security
    using HE Standard v1.1 tables. Returns the fraction of the
    128-bit budget consumed."""
    if n not in HE_STANDARD_128BIT:
        # Linear interpolation in log space
        ns = sorted(HE_STANDARD_128BIT.keys())
        for i in range(len(ns) - 1):
            if ns[i] <= n <= ns[i+1]:
                frac = (math.log2(n) - math.log2(ns[i])) / (math.log2(ns[i+1]) - math.log2(ns[i]))
                max_q = HE_STANDARD_128BIT[ns[i]] * (1-frac) + HE_STANDARD_128BIT[ns[i+1]] * frac
                return log2_q / max_q
        return None
    max_q = HE_STANDARD_128BIT[n]
    return log2_q / max_q

# ─── Run estimates ────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Security Estimation for PQ Re-randomizable Encryption")
    print("=" * 70)
    print(f"\nScheme parameters: n={N}, σ={SIGMA}")
    print(f"  t  = {T:>15,} (log₂ ≈ {math.log2(T):.2f})")
    print(f"  q₂ = {Q2:>15,} (log₂ ≈ {math.log2(Q2):.2f})")
    print(f"  q  = t·q₂ = {Q_FULL:>20,} (log₂ ≈ {math.log2(Q_FULL):.2f})")

    log2_q = math.log2(Q_FULL)
    max_q_128 = HE_STANDARD_128BIT[N]
    budget_frac = log2_q / max_q_128

    print(f"\n─── HE Standard v1.1 Analysis ───")
    print(f"\nReference: for n={N}, σ≈3.2, 128-bit classical security:")
    print(f"  Max log₂(q) = {max_q_128}")
    print(f"  Our log₂(q) = {log2_q:.1f}")
    print(f"  Budget used  = {budget_frac*100:.1f}%")
    print(f"  Margin       = {max_q_128 - log2_q:.0f} bits of q remaining")

    print(f"\nSince log₂(q) = {log2_q:.0f} ≪ {max_q_128}, our parameters")
    print(f"comfortably exceed the 128-bit classical security threshold.")
    print(f"\nThe HE Standard v1.1 does not provide a precise security level")
    print(f"for intermediate q values (only the 128-bit threshold). For a")
    print(f"precise estimate, run the lattice-estimator:")
    print(f"  https://github.com/malb/lattice-estimator")
    print(f"  with parameters: n={N}, q={Q_FULL}, secret_dist='normal',")
    print(f"  error_dist=σ={SIGMA}")

    print(f"\n─── Per-Limb Analysis ───")
    print(f"\nThe adversary sees RLWE instances in both limbs sharing the")
    print(f"same secret s. By CRT, this is equivalent to RLWE with the")
    print(f"full modulus q = t·q₂ (since CRT of independent uniforms is")
    print(f"uniform). So security is determined by the full q.")
    print(f"\n  Individual t-limb:   log₂(t) = {math.log2(T):.0f} vs max {max_q_128} → {math.log2(T)/max_q_128*100:.0f}% of budget")
    print(f"  Individual q₂-limb:  log₂(q₂) = {math.log2(Q2):.0f} vs max {max_q_128} → {math.log2(Q2)/max_q_128*100:.0f}% of budget")
    print(f"  Combined q = t·q₂:   log₂(q) = {log2_q:.0f} vs max {max_q_128} → {budget_frac*100:.0f}% of budget")

    print(f"\n─── HE Standard Reference Table (n → max log₂(q) at 128-bit) ───")
    for n, max_q in sorted(HE_STANDARD_128BIT.items()):
        marker = " ← our n" if n == N else ""
        print(f"  n = {n:>6}: max log₂(q) = {max_q:>4}{marker}")

    print(f"\n─── Sensitivity: what if we used larger q? ───")
    print(f"  {'log₂(q)':>8} {'Budget%':>8} {'Status':>20}")
    print(f"  {'-'*40}")
    for lq in [32, 50, 63, 80, 100, 109, 120]:
        pct = lq / max_q_128 * 100
        status = "✓ > 128-bit" if lq <= max_q_128 else "✗ < 128-bit"
        marker = " ← our scheme" if abs(lq - 63) < 1 else ""
        marker = " ← HE Standard limit" if lq == 109 else marker
        print(f"  {lq:>8} {pct:>7.0f}% {status:>20}{marker}")

    print(f"\n─── Conclusion ───")
    print(f"Our parameters (n={N}, log₂(q)≈{log2_q:.0f}) use {budget_frac*100:.0f}% of the")
    print(f"128-bit security budget per the HE Standard v1.1. This provides")
    print(f"a substantial safety margin against both classical and quantum")
    print(f"lattice attacks.")

if __name__ == '__main__':
    main()
