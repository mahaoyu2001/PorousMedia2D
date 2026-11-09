R5Q9: leading hard-core pair-correlation prefactor

Analytical result: C_i,0 = -1 for the zero-extended circular-pore Kirsch
disturbance under biaxial normal loading, for isotropic step pair correlation
g_2,0(r)=0 inside an exclusion radius D>=2a and 1 outside.
The derivation is in Section 2 of Supplementary Material 3.

Run verify_pair_prefactor in MATLAB (no toolbox required).
It writes Pair_Prefactor_FFT_Convergence.csv beside the script.
The numerical check uses a=1, Sigma1=0, Sigma2=1. Components 1,2,3 mean
11,22,12. An overall mean-field stress factor cancels in the coefficient.
The square has half-width R; dx is the midpoint-grid spacing. The kernel is
zero inside the pore. FFT autocorrelation uses periodic continuation of the
truncated square. Increasing R checks its boundary/truncation effects;
decreasing dx checks kernel sampling, including the jump at r=a.
The disk integration uses 8-by-8 subcell weights at the exclusion boundary.
D/a=2 and 2.04 are both checked (the latter includes the generator's 1.02
spacing guard). Both have the same infinite-domain analytical coefficient.
The supplied finite-box values approach -1; they are numerical checks of the
identity, not improved estimates to substitute for its exact value.

Scope: this is the leading low-density pair contribution to an unconditioned,
prescribed superposition at equal center intensity. It is not an empirical
finite-density g_2 estimate, a matrix-conditioned variance correction, or a
total FEM error estimate. No existing statistics/figures should be adjusted
using this coefficient. No new finite-element calculations are performed.
