# Exact finite-window VAR(p) population spectra and conditional canonical coherence

**SciFort extension, not a `scipy.signal` function.** This extends the
[VAR(1) spectral target](signal_var1_smoothed_coverage.md) to stable, real
multivariate VAR(p) models and evaluates frequency-band integrals **without
approximating an integral on a sampled FFT grid**. The output can be used to
specify the appropriate population estimand for conditional canonical
coherence and aggregate group coherence.

The companion [fixed-weight multitaper extension](signal_varp_multitaper_targets.md) supports deterministic DPSS and arbitrary taper banks, with exact band integration and conditional canonical coherence of the combined spectral matrix.

## Public Fortran API

```fortran
use scifort_signal, only : dp, varp_taper_spectrum_result, &
    varp_expected_taper_spectrum, varp_expected_band_taper_spectrum, &
    varp_expected_canonical_coherence, varp_expected_band_canonical_coherence, &
    canonical_coherence_result, varp_spectrum_success

type(varp_taper_spectrum_result) :: spectrum, integrated
type(canonical_coherence_result) :: modes, band_modes
real(dp) :: a(5,5,2), sigma(5,5), taper(24), bands(2,2)

! Initialize stable VAR(2) coefficients a(:,:,1:2), positive-definite
! innovation covariance sigma, nonzero taper, and valid band boundaries.
spectrum = varp_expected_taper_spectrum(a,sigma,taper, &
    [0.0_dp,0.125_dp,0.25_dp], normalized=.true.)

bands(:,1) = [0.0_dp,0.15_dp]
bands(:,2) = [0.15_dp,0.4_dp]
integrated = varp_expected_band_taper_spectrum(a,sigma,taper,bands, &
    normalized=.true.)

modes = varp_expected_canonical_coherence(a,sigma,taper, &
    [0.0_dp,0.125_dp,0.25_dp],[1,2],[3,4],conditioned_on=[5])
band_modes = varp_expected_band_canonical_coherence(a,sigma,taper, &
    bands,[1,2],[3,4],conditioned_on=[5],average=.true.)

if (integrated%status == varp_spectrum_success) then
    print *, real(integrated%power(1,1,1),dp)
end if
```

All routines are available both in their own modules and through the
`scifort_signal` umbrella module. No C ABI or additional Fortran dependencies
are introduced.

- `a(n,n,p)` contains coefficients from lag one through lag p, with `p>=1`.
  `innovation_covariance(n,n)` must be finite, symmetric, and **strictly**
  positive definite. The process must be stationary and causal.
- `taper(length)` must be finite and have nonzero energy. Fourier
  `frequencies` and `bands(2,n_band)` are in **cycles per sample**, not Hz.
  Frequencies may be negative and do not need to be sorted. Bands may
  overlap; each lower boundary must be strictly less than its upper boundary.
  The spectrum is periodic with period one cycle/sample.
- `normalized=.true.` divides by `sum(taper**2)`; the default is the
  **raw** expected Fourier outer product. `average=.true.` divides an exact
  band integral by that band's width, giving mean spectral density rather
  than integrated power. Both options preserve coherence without fixed
  additive regularization; units still matter for comparing spectral power.
- `result%power(f,i,j)` follows the convention
  `E[Z_i(f)*conjg(Z_j(f))]`. Band results store the integrated cross-spectrum
  in `%power(b,i,j)`; `%bands` stores the original limits and `%frequencies`
  contains band midpoints as labels (not sampled spectral frequencies).
  `%lag0_covariance` is the stationary covariance of the observed channels.
- Status codes: `varp_spectrum_success=0`,
  `varp_spectrum_invalid_input=1`, `varp_spectrum_nonstationary=2`.
  Invalid or nonconverged calls do not allocate the cross-spectral arrays.
- Canonical wrappers accept arbitrary nonempty disjoint groups and an
  optional disjoint `conditioned_on` set. They delegate to the existing
  `canonical_coherence_matrix` solver and can accept fixed `ridge` and
  `shrinkage` parameters and `return_vectors=.true.`. Invalid groups and
  singular cells are reported using its existing statuses and validity masks.
  The wrappers do not fit a VAR model or estimate confidence intervals.

## Mathematical definition

For a zero-mean stationary real VAR(p),

\[
x_t=\sum_{j=1}^p A_jx_{t-j}+\varepsilon_t,\qquad
E[\varepsilon_t\varepsilon_t^T]=\Sigma,
\]

form the state vector `y(t)=(x(t),...,x(t-p+1))` and companion transition
matrix `F`. Its innovation covariance is `Q=diag(Sigma,0,...,0)` (singular
for p>1). Solve the stationary discrete Lyapunov equation

\[
G=FGF^T+Q
\]

by **Smith doubling**, which allows the singular state innovation `Q`. The
original-channel lag covariances are the leading `n x n` blocks of `F^h G`
for `h>=0`. The algorithm propagates `F^h G(:,1:n)` rather than repeatedly
forming full matrix powers.

With taper weights `w_t` and autocorrelations
`c_h=sum_{t=0}^{L-1-h}w_{t+h}w_t`, the exact expected finite-window
cross-spectrum is

\[
S_L(f)=c_0\Gamma_0+\sum_{h=1}^{L-1}c_h\left[
  \Gamma_h e^{-2\pi i fh}+\Gamma_h^T e^{2\pi i fh}\right].
\]

For a band `[a,b]` (cycles/sample), integration of each exponential is
**analytic**:

\[
\begin{split}
I_0(a,b)&=b-a,\\
I_h(a,b)&=\int_a^b e^{-2\pi ifh}\,df\\
&=(b-a)\,\operatorname{sinc}(h(b-a))
  e^{-\pi i h(a+b)},\quad h\ge1,
\end{split}
\]

where `sinc(x)=sin(pi*x)/(pi*x)`. For very narrow bands a Taylor series
avoids cancellation in the sinc ratio. Thus the *exact* band-integrated
population matrix is

\[
S_{[a,b]}=c_0(b-a)\Gamma_0+
\sum_{h=1}^{L-1}c_h\left[I_h\Gamma_h+
\overline{I_h}\Gamma_h^T\right].
\]

Integration and canonical coherence are **not interchangeable**:
coherence of the integrated matrix generally differs from a mean of
per-frequency coherence values. For group A and B conditioned on C, the
canonical wrappers remove C by spectral Schur complements, whiten the
remaining A and B spectral covariance blocks, and compute the squared
singular values in descending order. The aggregate group coherence is
`1-prod(1-squared_modes)`.

These are population *plug-in* functionals, such as
`canonical(E[tapered cross spectrum])`. They are not the expectation of
sample canonical correlations, and do not automatically account for
estimator-specific adaptive DPSS weighting, detrending, ridge-selection,
shrinkage-selection, overlapping windows, or bootstrap construction. In
particular, finite-window exactness **does not establish 95% coverage** of
any bootstrap or jackknife interval.

## Verification

```bash
python tools/run_signal_varp_population.py
python tools/check_signal_varp_population_public.py
```

The first script builds with checked and optimized GNU Fortran, executes
native analytic and invalid-input tests, and compares 80 random VAR(p)
processes against an independent SciPy `solve_discrete_lyapunov` plus
**direct sample-pair** NumPy sum and analytically integrated sample-pair
band kernels. It also simulates 20,000 stationary-approximated VAR(2)
windows to verify the finite-taper covariance prediction empirically.

The second script compiles the complete `scifort_signal` public dependency
closure and all existing public tests, then compares conditioned canonical
coherence for 50 random VAR(p) processes to independent SciPy Cholesky and
singular-value calculations, at multiple frequencies and bands.

Both scripts require GNU Fortran, NumPy and SciPy for tests; the compiled
Fortran core itself requires no SciPy, Python, BLAS, or LAPACK dependency.

## Limitations and next steps

This direct covariance formulation is intended for **moderate n, p, and
window lengths**: Smith doubling forms dense `(n*p) x (n*p)` products.
It detects divergence/nonconvergence but does not provide a numerical
certificate of stability for every near-unit-root model. Input data must
be finite and the innovation covariance strictly positive definite;
complex-valued VAR coefficients, exogenous terms, and nonstationary models
are not supported. Fourier output is two-sided and does not include a
one-sided PSD scaling or automatic conversion to Hz.

Useful follow-up studies would simulate conditional canonical bootstrap
intervals under VAR(2)/VAR(3), vary smoother/taper and block lengths,
and distinguish finite-window population targets from finite-sample
expectations of the *nonlinear* coherence estimator. Rank-tied modes
require subspace-level rather than individual-mode inference.
