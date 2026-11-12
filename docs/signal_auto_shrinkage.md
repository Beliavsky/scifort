# Automatic diagonal spectral shrinkage (SciFort extension)

SciFort's `partial_coherence_matrix`, `partial_coherence_tf_matrix`,
`multitaper_partial_coherence`, and `multitaper_tf_partial_coherence`
support **opt-in frequency-local diagonal-target shrinkage**. A separate
predictive leave-one-out option, `cv_shrinkage=.true.`, supports adaptive taper
weights; see [cross-validated shrinkage](signal_cv_shrinkage.md). This estimates
how much to suppress off-diagonal cross-spectral terms before computing
complex partial coherency and squared partial coherence. It is an original
SciFort numerical extension, **not** a public `scipy.signal` API.

## Quick start

```fortran
use scifort_signal, only : dp, multitaper_partial_coherence, &
    partial_coherence_result

real(dp) :: x(512,4)  ! (sample, channel)
type(partial_coherence_result) :: result

! Fill x with finite observations.
result = multitaper_partial_coherence(x, 3.5_dp, kmax=6, &
    weighting='equal', auto_shrinkage=.true., &
    ridge=0.01_dp, return_coherency=.true.)

! Per-frequency shrinkage coefficients, between 0 and 1:
print *, result%shrinkage_used
! Phase-preserving partial coherency of channels 1 and 3:
print *, result%coherency(:,1,3)
```

`result%shrinkage_used(f)` reports the estimated coefficient for each
frequency bin, or NaN if its input spectrum is entirely zero. For time-frequency
results, `result%shrinkage_used(f,t)` is chosen independently for each
frequency and time window. With explicit `shrinkage=alpha`, these fields
report `alpha` instead; without either option they report zero on nonzero
spectral bins. The new arrays are allocated even when automatic shrinkage is
disabled, but cost only one real number per frequency/time cell.

For an already computed cross-spectrum matrix, **supply the effective number
of independent taper spectra explicitly**:

```fortran
spectral = multitaper_spectral_matrix(x, 3.5_dp, kmax=6, weighting='equal')
result = partial_coherence_matrix(spectral, auto_shrinkage=.true., &
    effective_tapers=6.0_dp, return_coherency=.true.)
```

For an existing 4-D multitaper spectral matrix spectrogram:

```fortran
result_tf = partial_coherence_tf_matrix(spectral_tf, &
    auto_shrinkage=.true., effective_tapers=6.0_dp)
```

The signal-level calls infer the number of tapers from `kmax` or the existing
DPSS default `min(nsamples,max(1,int(2*nw)-1))`, using `nperseg` for moving
windows. Automatic shrinkage accepts both `weighting='equal'` (default) and
`weighting='eigen'`. Equal weights use the taper count; deterministic DPSS
eigenvalue weights use the effective count below. `weighting='adaptive'`
remains excluded because its weights depend on the observations. All modes
remain available with fixed shrinkage.

### Eigenvalue-weighted example

```fortran
use scifort_signal, only : dp, multitaper_partial_coherence, &
    multitaper_spectral_matrix, partial_coherence_matrix, &
    multitaper_matrix_result, partial_coherence_result

real(dp) :: x(512,4) ! populate with signal observations
type(multitaper_matrix_result) :: spectral
type(partial_coherence_result) :: conditional

conditional = multitaper_partial_coherence(x,3.5_dp,kmax=6, &
    weighting='eigen',auto_shrinkage=.true.,ridge=0.02_dp)
! Or reuse a previously estimated spectral matrix:
spectral = multitaper_spectral_matrix(x,3.5_dp,kmax=6,weighting='eigen')
conditional = partial_coherence_matrix(spectral,auto_shrinkage=.true., &
    effective_tapers=spectral%effective_tapers,ridge=0.02_dp)
```

The static and moving-window spectral matrix result types report
`%effective_tapers` for equal or eigen weighting and zero for invalid or
adaptive-weighted outputs. The public helper
`multitaper_effective_tapers(weights)` computes the equivalent count for
any finite nonnegative fixed weights, returning zero for empty, nonfinite,
negative or all-zero inputs. Matrix-level automatic shrinkage retains its
requirement for an explicit `effective_tapers=` argument.


## Mathematical model

At one frequency/time bin, write the cross-spectrum estimate as a Hermitian
matrix `S` of size `p` channels by `p` channels. The shrinkage target is the
matrix of individual channel powers:

\[
  T=\operatorname{diag}(S),\qquad S_\alpha=(1-\alpha)S+\alpha T.
\]

Assume `S` estimates the covariance of `n` **independent, zero-mean, proper
complex Gaussian** Fourier-coefficient vectors. Under the complex Wishart
model, for distinct channels `i,j`,

\[
 E[S_{ii}S_{jj}]=\Sigma_{ii}\Sigma_{jj}+|\Sigma_{ij}|^2/n,
 \qquad
 E[|S_{ij}|^2]=|\Sigma_{ij}|^2+\Sigma_{ii}\Sigma_{jj}/n.
\]

A moment plug-in estimate for the diagonally targeted minimum-risk
shrinkage factor is

\[
 D=\sum_{i\ne j} S_{ii}S_{jj},\qquad
 O=\sum_{i\ne j}|S_{ij}|^2,
 \qquad
 \hat\alpha=\operatorname{clip}_{[0,1]}
 \left(\frac{nD-O}{(n^2-1)O}\right).
\]

For **fixed deterministic, nonnegative weights**, independent proper complex
Gaussian Fourier vectors also obey the moment model with the effective
sample count

\[
 n_{\mathrm{eff}}=\frac{(\sum_k w_k)^2}{\sum_k w_k^2}.
\]

This reduces to the number of tapers for equal weights and can be fractional
for DPSS eigenvalue weights. It is derived from the variance of a weighted
sample covariance of independent observations, and does not automatically
account for correlated tapered spectra. Data-dependent adaptive weights
invalidate the independence-of-weights assumption.

The code evaluates `D` and `O` after dividing by the largest spectral
matrix magnitude to avoid over/underflow from spectral units. For `n=1`,
it uses full diagonal shrinkage (`alpha=1`); when `O=0`, shrinkage is
unnecessary and it uses `alpha=0`. For a zero spectral matrix, the bin is
singular and its reported coefficient is NaN. The result is then inverted
by the **existing Hermitian Cholesky** algorithm, after any separately
requested trace-scaled `ridge` loading.

This is a **Gaussian-Wishart moment-based plug-in estimate**, not an exact
oracle, a full Ledoit-Wolf estimator, or a calibrated significance test.
Automatic shrinkage may stabilize inversions when there are fewer independent
spectral vectors than channels, but it **does not guarantee** an invertible
matrix. If a bin remains singular, `valid(f)` or `valid(f,t)` is false and
`status=partial_coherence_singular`; its precision/coherence values are NaN.

## Statistical limitations and validation

- Orthogonal DPSS tapers are not automatically **statistically independent**
  for colored, nonstationary, strongly line-dominated, or non-Gaussian
  processes. Treat the fixed-weight effective count as an approximation to the
  independent sample count; choose manual `effective_tapers` for input
  matrices only when it has a defensible interpretation.
- Proper-complex Gaussian assumptions may fail at DC and real-valued Nyquist
  bins. The formula remains a defined numerical regularizer there, but its
  nominal risk calculation may not apply.
- Automatic shrinkage makes no claims about conditional independence,
  directed connectivity, p-values, confidence intervals or causal effects.
- `auto_shrinkage=.true.` cannot be combined with explicit `shrinkage=`.
  Matrix-level calls require `effective_tapers >= 1` and finite. Supplying
  `effective_tapers` without enabling automatic shrinkage is rejected.
- Existing `shrinkage=` and `ridge=` semantics, positional APIs, `valid`
  masks, and `return_coherency` behavior remain unchanged. Float estimates
  are computed separately for every valid frequency/time cell.

The test suite builds checked and targeted `-O3` versions of the module,
checks the legacy native programs, and independently compares 196 randomized
Hermitian real/complex spectral cases (949 invertible and eight zero/singular
cells, 85,545 precision/coherency values per build) against NumPy. The maximum
scaled difference observed was `4.27e-15` for the original matrix-only regression.

```bash
python tools/run_signal_partial_shrinkage_auto.py
python tools/run_signal_partial_shrinkage_auto.py --legacy-references
python tools/check_signal_partial_shrinkage_auto_public.py
```

The `--legacy-references` option also checks the previous 288-, 252-, and
120-case independent NumPy/SciPy reference suites in checked mode. The
public-module checker builds the full `scifort_signal` dependency closure
and verifies both old and new public interfaces. NumPy/SciPy are only test
references; the Fortran implementation adds no runtime dependencies.
