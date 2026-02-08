#!/usr/bin/env python3
"""
Generate publication-quality figures for the paper.
Run from the pq-rerand/ directory: python3 tools/make_figures.py

Fixed-key simulation: 20 keys x 10 trials per key, matching Theorem 4.6's
operational setting (key fixed, randomness over encryption/rerand only).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'legend.fontsize': 9.5,
    'figure.dpi': 150,
})
import matplotlib.pyplot as plt
import os, math

# Parameters
N = 4096
SIGMA = 3.2
K_MAX = 858_000_000
Q2 = 4_294_828_033
THRESHOLD = Q2 / 2
# Fixed-key sub-Gaussian bound: sigma_k * sqrt(2*(4*N*SIGMA^2 + 1)*ln(2/delta))
# with delta = 2^{-128}
SG_CONST_SQ = 2 * (4 * N * SIGMA**2 + 1) * (math.log(2) * 129)
SIGMA_FLOOD = SIGMA * math.sqrt(K_MAX)

NUM_KEYS = 10
NUM_TRIALS_PER_KEY = 5

os.makedirs('paper/figures', exist_ok=True)

# ─── Negacyclic multiplication via FFT ────────────────────────────────

def negacyclic_mul(a, b):
    n = len(a)
    fa = np.fft.rfft(a.astype(np.float64), 2*n)
    fb = np.fft.rfft(b.astype(np.float64), 2*n)
    c = np.round(np.fft.irfft(fa * fb, 2*n)).astype(np.int64)
    return c[:n] - c[n:]

# ─── Fixed-key simulation ─────────────────────────────────────────────

def simulate_noise_fixed_key(s, e, k, rng):
    """Simulate noise for a FIXED key (s, e) with fresh enc/rerand randomness."""
    # Fresh encryption
    r = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
    e1 = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
    e2 = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
    noise = negacyclic_mul(e, r) + e2 - negacyclic_mul(e1, s)
    # Flood (one-time smudge)
    rf = np.round(rng.normal(0, SIGMA_FLOOD, N)).astype(np.int64)
    e1f = np.round(rng.normal(0, SIGMA_FLOOD, N)).astype(np.int64)
    e2f = np.round(rng.normal(0, SIGMA_FLOOD, N)).astype(np.int64)
    noise += negacyclic_mul(e, rf) + e2f - negacyclic_mul(e1f, s)
    # k re-randomizations
    for _ in range(k):
        rr = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
        er1 = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
        er2 = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
        noise += negacyclic_mul(e, rr) + er2 - negacyclic_mul(er1, s)
    return noise

print(f"Fixed-key simulation: {NUM_KEYS} keys x {NUM_TRIALS_PER_KEY} trials")
print("Generating keys...")
rng = np.random.default_rng(seed=2024)

# Pre-generate fixed keys
keys = []
for ki in range(NUM_KEYS):
    s = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
    e = np.round(rng.normal(0, SIGMA, N)).astype(np.int64)
    keys.append((s, e))
    norm_s2 = float(np.sum(s**2))
    norm_e2 = float(np.sum(e**2))
    expected = N * SIGMA**2
    print(f"  Key {ki:2d}: ||s||^2={norm_s2:.0f}, ||e||^2={norm_e2:.0f} "
          f"(expected ~{expected:.0f}, good-key bound {2*expected:.0f})")

SIM_K = [0, 1, 5, 10, 50, 100, 500, 1000]

print("\nSimulating noise (this takes a few minutes)...")

# For each k: for each key, run NUM_TRIALS_PER_KEY trials
# Store: per-key means and the cross-key mean/std
sim_results = {}       # k -> (cross_key_mean, cross_key_std)
sim_per_key = {}       # k -> list of NUM_KEYS per-key means

for k in SIM_K:
    per_key_means = []
    for key_idx, (s, e) in enumerate(keys):
        trial_maxes = []
        for _ in range(NUM_TRIALS_PER_KEY):
            noise = simulate_noise_fixed_key(s, e, k, rng)
            trial_maxes.append(np.max(np.abs(noise)))
        per_key_means.append(np.mean(trial_maxes))
    per_key_means = np.array(per_key_means)
    cross_mean = np.mean(per_key_means)
    cross_std = np.std(per_key_means)
    sim_results[k] = (cross_mean, cross_std)
    sim_per_key[k] = per_key_means
    print(f"  k={k:>6}: cross-key mean |ν|_∞ = {cross_mean:.3e} ± {cross_std:.3e} "
          f"(min={np.min(per_key_means):.3e}, max={np.max(per_key_means):.3e})")

# Fit A from simulated data: noise = A * sqrt(k_max + 1 + k)
A_vals = [sim_results[k][0] / np.sqrt(K_MAX + 1 + k) for k in SIM_K]
A = np.mean(A_vals)
print(f"  Fitted A = {A:.2f}")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1: Noise growth — full extrapolation to k_max (log-x scale)
# ═══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 4.5))

# Extrapolated model
k_extrap = np.logspace(0, np.log10(2 * K_MAX), 500)
noise_model = A * np.sqrt(K_MAX + 1 + k_extrap)
# Sub-Gaussian bound: B_coeff(k) = sigma_k * sqrt(SG_CONST_SQ)
# where sigma_k = sigma * sqrt(k + 1 + K_MAX)
theoretical = SIGMA * np.sqrt(k_extrap + 1 + K_MAX) * np.sqrt(SG_CONST_SQ)

ax.semilogx(k_extrap, noise_model / 1e8, 'b-', linewidth=1.5,
            label='Empirical fit: $A\\sqrt{\\kappa_f+1+k}$')
ax.semilogx(k_extrap, theoretical / 1e8, 'r--', linewidth=1.2,
            label='Sub-Gaussian bound ($\\delta_{\\mathrm{coeff}}=2^{-128}$)')
ax.axhline(THRESHOLD / 1e8, color='black', linestyle=':', linewidth=1.5,
           label=f'Threshold $q_2/2$')

# Simulated data points — show cross-key spread
k_pts = np.array(SIM_K[1:])  # skip k=0 for log scale
means = np.array([sim_results[k][0] for k in SIM_K[1:]])
stds = np.array([sim_results[k][1] for k in SIM_K[1:]])
ax.errorbar(k_pts, means / 1e8, yerr=stds / 1e8, fmt='ko', markersize=4,
            capsize=3, label=f'Fixed-key sim ({NUM_KEYS} keys, {NUM_TRIALS_PER_KEY} trials/key)', zorder=5)

ax.axvline(K_MAX, color='gray', linestyle='-.', alpha=0.7, linewidth=1)
ax.text(K_MAX * 0.6, 0.3, '$k_{\\max}$', color='gray', fontsize=10, ha='right')

ax.set_xlabel('Number of additional re-randomizations ($k$)')
ax.set_ylabel('Max $|\\nu|$ per coefficient  ($\\times 10^8$)')
ax.set_title('Noise Growth vs Re-randomizations (Fixed-Key)')
ax.legend(loc='upper left', framealpha=0.9)
ax.set_ylim(0, THRESHOLD / 1e8 * 1.15)
ax.set_xlim(1, 2 * K_MAX)
ax.grid(True, alpha=0.2, which='both')

fig.tight_layout()
fig.savefig('paper/figures/Fig1.png', dpi=200)
print("Saved paper/figures/Fig1.png  (noise growth)")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2: Headroom vs k (log-x)
# ═══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 4))

headroom_emp = THRESHOLD / (A * np.sqrt(K_MAX + 1 + k_extrap))
headroom_theory = THRESHOLD / (SIGMA * np.sqrt(k_extrap + 1 + K_MAX) * np.sqrt(SG_CONST_SQ))

ax.semilogx(k_extrap, headroom_emp, 'b-', linewidth=1.5,
            label='Empirical headroom')
ax.semilogx(k_extrap, headroom_theory, 'r--', linewidth=1.2,
            label='Sub-Gaussian headroom (worst-case)')
ax.axhline(1.0, color='red', linestyle=':', linewidth=1, alpha=0.7,
           label='Failure boundary')
ax.axvline(K_MAX, color='gray', linestyle='-.', alpha=0.7, linewidth=1)
ax.text(K_MAX * 0.6, 8, '$k_{\\max}$', color='gray', fontsize=10, ha='right')

ax.set_xlabel('Number of additional re-randomizations ($k$)')
ax.set_ylabel('Headroom ($q_2/2 \\;/\\; \\|\\nu\\|_\\infty$)')
ax.set_title('Decryption Headroom vs Re-randomizations (Fixed-Key)')
ax.legend(loc='upper right', framealpha=0.9)
ax.set_ylim(0, 25)
ax.set_xlim(1, 2 * K_MAX)
ax.grid(True, alpha=0.2, which='both')

fig.tight_layout()
fig.savefig('paper/figures/Fig2.png', dpi=200)
print("Saved paper/figures/Fig2.png  (headroom)")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3: Noise histogram at k=5000 (single fixed key)
# ═══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 4))

# Use first key for histogram
s0, e0 = keys[0]
noise_sample = simulate_noise_fixed_key(s0, e0, 1000, rng)
sigma_emp = np.std(noise_sample)

ax.hist(noise_sample / 1e6, bins=120, density=True, alpha=0.7, color='steelblue',
        edgecolor='none', label=f'Noise after $k=1000$ rerands (1 fixed key)')

# Gaussian overlay
from scipy.stats import norm
x = np.linspace(-4 * sigma_emp, 4 * sigma_emp, 500) / 1e6
pdf = norm.pdf(x, 0, sigma_emp / 1e6)
ax.plot(x, pdf, 'r-', linewidth=1.5, label=f'Gaussian fit ($\\sigma={sigma_emp/1e6:.2f} \\times 10^6$)')

ax.set_xlabel('Noise value ($\\times 10^6$)')
ax.set_ylabel('Density')
ax.set_title('Per-Coefficient Noise Distribution After 1,000 Re-randomizations')
ax.legend(framealpha=0.9)

fig.tight_layout()
fig.savefig('paper/figures/Fig3.png', dpi=200)
print("Saved paper/figures/Fig3.png  (noise histogram)")

# ═══════════════════════════════════════════════════════════════════════
# FIGURE 4: Inter-key variability at k=1000
# ═══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 3.5))

k_show = 1000
if k_show in sim_per_key:
    per_key = sim_per_key[k_show] / 1e8
    ax.bar(range(NUM_KEYS), per_key, color='steelblue', alpha=0.8,
           edgecolor='navy', linewidth=0.5)
    ax.axhline(np.mean(per_key), color='red', linestyle='--', linewidth=1.2,
               label=f'Cross-key mean = {np.mean(per_key):.3f}')
    ax.set_xlabel('Key index')
    ax.set_ylabel('Mean max $|\\nu|$ ($\\times 10^8$)')
    ax.set_title(f'Inter-Key Variability at $k={k_show}$ ({NUM_TRIALS_PER_KEY} trials/key)')
    ax.legend(framealpha=0.9)
    ax.set_xticks(range(0, NUM_KEYS, 2))
    ax.grid(True, alpha=0.2, axis='y')

fig.tight_layout()
fig.savefig('paper/figures/Fig4.png', dpi=200)
print("Saved paper/figures/Fig4.png  (key variability)")

print("\nDone. All figures in paper/figures/")
