#!/usr/bin/env sage
"""
Concrete security estimation for the two-limb CRT RLWE scheme.

Requires SageMath + lattice-estimator (https://github.com/malb/lattice-estimator).

Install:
  git clone https://github.com/malb/lattice-estimator.git
  cd lattice-estimator
  sage -python -m pytest  # verify installation

Run:
  sage lattice_estimate_sage.py
"""

import sys
sys.path.insert(0, ".")  # if run from lattice-estimator directory

from estimator import *

# === Scheme parameters ===
n = 4096                    # Ring dimension (x^n + 1)
t = 2_147_565_569           # Plaintext modulus (32-bit NTT prime)
q2 = 4_294_828_033          # Noise modulus (32-bit NTT prime)
q = t * q2                  # Full modulus  (approx 2^63)
sigma = 3.2                 # Gaussian std dev

print(f"=== Two-Limb CRT RLWE Parameters ===")
print(f"  n     = {n}")
print(f"  t     = {t}  (log2 = {t.bit_length()})")
print(f"  q2    = {q2} (log2 = {q2.bit_length()})")
print(f"  q     = {q}  (log2 = {q.bit_length()})")
print(f"  sigma = {sigma}")
print()

# Define the LWE parameters
# Secret distribution: discrete Gaussian with sigma=3.2
# Error distribution: discrete Gaussian with sigma=3.2
params = LWE.Parameters(
    n=n,
    q=q,
    Xs=ND.DiscreteGaussian(sigma),
    Xe=ND.DiscreteGaussian(sigma),
)

print(f"LWE Parameters: {params}")
print()

# Run all standard attacks
print("=== Running lattice estimator (this may take a few minutes) ===")
print()
results = LWE.estimate(params)

print()
print("=== Summary ===")
for attack_name, result in sorted(results.items(), key=lambda x: x[1].get("rop", float("inf"))):
    rop = result.get("rop", "?")
    if isinstance(rop, (int, float)):
        import math
        log2_rop = math.log2(rop) if rop > 0 else "?"
        print(f"  {attack_name:30s}: log2(rop) = {log2_rop:.1f}")
    else:
        print(f"  {attack_name:30s}: {rop}")

print()
print("=== HE Standard v1.1 Cross-reference ===")
print(f"  For n={n}, sigma={sigma} (Gaussian secret):")
print(f"  128-bit classical budget: log2(q) <= 109")
print(f"  128-bit quantum budget:   log2(q) <= 101")
print(f"  Our log2(q) = {q.bit_length()} ({q.bit_length()/109*100:.0f}% of classical budget)")
