# Exact finite-window VAR(1) spectral targets and bootstrap coverage

**SciFort extension, not a `scipy.signal` API.**
`var1_expected_taper_spectrum` computes the exact expected *raw* tapered
Fourier outer product of a stationary, real, zero-mean VAR(1) process. It
allows numerical tests of spectral-coherence confidence intervals against the
population functional of the *same finite-window spectral estimator* rather
than against an unrelated population frequency or spectrum.

## Fortran interface

```fortran
use scifort_signal, only : dp, var1_taper_spectrum_result, &
    var1_expected_taper_spectrum, var1_spectrum_success

type(var1_taper_spectrum_result) :: target
real(dp) :: a(3,3), sigma(3,3), taper(12)

! Set a, sigma and taper. Frequencies are in cycles/sample.
target = var1_expected_taper_spectrum(a,sigma,taper,[0.0_dp,0.125_dp])
if (target%status == var1_spectrum_success) then
    ! power(f,i,j) = E[Z_i(f) conjg(Z_j(f))]
    print *, target%power(1,1,2)
    print *, target%lag0_covariance(1,1)
end if
```

Inputs are a square real VAR(1) matrix `a`, its symmetric strictly
positive-definite innovation covariance `innovation_covariance`, nonzero
finite taper `taper`, and finite frequencies in cycles per sample. The
optional `normalized=.true.` divides the returned matrices by the sum of
squared taper weights; the default returns the *raw* spectral outer product
and matches the simulated driver directly. Values at frequencies differing
by an integer are mathematically equivalent. There is no factor of two for
one-sided real spectra, no sampling-rate conversion, and no detrending.

`result%power(f,i,j)` is allocated on success, along with
`result%frequencies` and `result%lag0_covariance`. Status constants are
`var1_spectrum_success=0`, `var1_spectrum_invalid_input=1`, and
`var1_spectrum_nonstationary=2`. Failed calls do not allocate spectral
arrays. The routine uses the Smith-doubling discrete Lyapunov iteration
and returns status 2 on divergence or nonconvergence. It requires a stable
VAR(1), finite input and nonsingular positive-definite innovation covariance.
No BLAS/LAPACK dependency is introduced. It is intended for moderate channel
counts and taper lengths, not high-order VAR models.

## Population quantity

For `x(t) = A*x(t-1) + e(t)`, with `E[e(t)e(t)^T]=Sigma`, the stationary
lag-zero covariance satisfies

\[
\Gamma_0=A\Gamma_0 A^T+\Sigma,
\qquad \Gamma_h=A^h\Gamma_0\quad(h\geq 0).
\]

For a taper of length `L` with weights `w_t`, define
`Z(f)=sum_t w_t*x(t)*exp(-2*pi*i*f*t)`. With
`c_h=sum_{t=0}^{L-1-h} w_{t+h}*w_t`, the target is

\[
Q(f)=c_0\Gamma_0+\sum_{h=1}^{L-1}c_h\left[
\Gamma_h e^{-2\pi ifh}+\Gamma_h^T e^{2\pi ifh}\right].
\]

This formula retains the finite Hann-window spectral smoothing and complex
cross-channel phase, even when the underlying time series are real. In
particular, squared coherence of channels 1 and 2 is

\[
C(Q)=\frac{|Q_{12}|^2}{Q_{11}Q_{22}}.
\]

The relevant population *plug-in functional* is `C(E[Z Z*])`, which in
general is **not** `C(S(0))` (unsmoothed zero-frequency coherence),
`E[C(Z Z*)]`, or `C(E[sample ratio])`. A nonzero `ridge`, channel-dependent
shrinkage, detrending, different DPSS weights, frequency-band integration,
or overlapping windows changes the estimand and must be modeled separately.
The simulation below uses zero ridge/shrinkage and nonoverlapping Hann
windows; its target corresponds exactly to those choices.

## Independent numerical checks

Run:

```bash
python tools/run_signal_var1_smoothed_coverage.py
python tools/check_signal_var1_smoothed_public.py
```

The first script builds checked and optimized GNU Fortran versions, runs the
native analytic/error-handling tests, compares 80 seeded randomized cases
(240 spectral matrices) to *SciPy* `solve_discrete_lyapunov` plus independent
sample-pair covariance summation in NumPy, and checks empirical moments from
60,000 simulated windows started in the stationary distribution. The public
script compiles the full `scifort_signal` dependency closure and runs its
previous public regressions plus the new interface test.

Run the four-model coverage study with:

```bash
python tools/run_signal_var1_smoothed_coverage.py --checked-only --study --trials 64
```

The study generates 64 independent datasets per model, with 48 complete
nonoverlapping 12-sample periodic-Hann windows, burn-in of 400 observations,
199 fixed-seed circular moving-block draws of length four, and the five
existing `canonical_moving_block_bootstrap` interval methods. All five use
the same data and bootstrap seeds. The four VAR(1) models include a
persistent target with white nuisance channel, a moderate target with
persistent nuisance, a lag-induced target with independent innovations,
and an opposite-sign autoregressive pair. The third channel is *not*
conditioned upon, so the bootstrap group statistic is the unconditional
one-mode squared coherence of channels 1 and 2.

## Reproducible 64-dataset findings (seed 20261030)

| Model | Finite-window target | Unsmoothed `C(S(0))` | Percentile valid/covered | Basic valid/covered | BC valid/covered | BCa valid/covered | Studentized valid/covered |
|---|---:|---:|---:|---:|---:|---:|---:|
| Persistent target | 0.42978 | 0.61011 | 64 / 59 | 64 / 52 | 64 / 59 | 64 / 59 | 63 / 59 |
| Moderate target | 0.35079 | 0.37985 | 64 / 60 | 64 / 60 | 64 / 60 | 64 / 60 | 60 / 57 |
| Lag-induced target | 0.14067 | 0.27468 | 64 / 57 | 64 / 45 | 64 / 55 | 63 / 56 | 33 / 33 |
| Opposite-sign AR | 0.08706 | 0.17672 | 64 / 64 | 64 / 54 | 64 / 61 | 64 / 60 | 30 / 29 |

`valid / covered` expresses the number of intervals for which bounds are
identified and the number of those that contain the smoothed population
target. For example, the persistent-target percentile estimate covers
`59/64 = 92.2%`, compared with the much lower apparent coverage when
comparing to the *different* unsmoothed zero-frequency target.
The study prints Wilson intervals for Monte Carlo uncertainty, average
interval width, and valid fractions for every method. Studentized intervals
become unidentifiable relatively often near weak coherence, so reporting
coverage conditional on their validity alone could be misleading.

**These are small exploratory simulations, not a demonstration of calibrated
95% confidence intervals.** Each estimated coverage proportion has
substantial Monte Carlo uncertainty. The fixed block size, bootstrap count,
window length, model coefficients, and parameter boundaries affect the
results. Future studies should vary those inputs, use substantially more
replications, evaluate bias in the estimated coherence itself, and explore
conditional and multi-mode canonical coherence. The analytic population
calculation is exact under the stipulated real stationary VAR(1) model,
up to floating-point error; the *coverage estimates* are not exact.
