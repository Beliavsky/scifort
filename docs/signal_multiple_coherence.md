# Multichannel spectral multiple coherence (SciFort extension)

SciFort now estimates **squared multiple coherence** between any one channel
and the *joint linear predictor* formed from every other channel. It also
returns the residual (conditional) spectral power and power explained by that
joint predictor. This is a **SciFort extension**: it is not a routine in
`scipy.signal`. It complements pairwise ordinary coherence and pairwise
partial coherence, which answer different questions.

## Mathematical definition

For a positive-definite Hermitian channel cross-spectral matrix `S(f)` of
size `nchannel` by `nchannel`, define `P(f) = inverse(S(f))`. Then the
conditional (unexplained) spectral power for channel `i`, given all the
other channels, is the Schur complement

\[
S_{i\mid-i}(f) = S_{ii}(f) - S_{i,-i}(f)\,S_{-i,-i}(f)^{-1}\,S_{-i,i}(f)
             = \frac{1}{P_{ii}(f)}.
\]

Its **squared multiple coherence** is

\[
R_i^2(f) = 1-\frac{S_{i\mid-i}(f)}{S_{ii}(f)}
= 1-\frac{1}{S_{ii}(f)P_{ii}(f)} \quad \in[0,1].
\]

This is a dimensionless, frequency-local **linear association measure**, not
a causal measure or a statistical significance probability. For two channels,
it is identical to ordinary magnitude-squared coherence; for three or more
channels, it differs from any pairwise or partial-coherence entry.

The implementation reuses SciFort's existing Hermitian Cholesky precision
solver (no new LAPACK dependency). It obtains each residual spectrum from the
precision diagonal and computes the ratio as `(1/S_ii)/P_ii` to avoid overflow
in the product `S_ii*P_ii`.

## API examples

Compute from previously estimated cross-spectral matrices:

```fortran
use scifort_signal, only : dp, multitaper_spectral_matrix, &
    multiple_coherence_matrix, multiple_coherence_result, multitaper_matrix_result

real(dp) :: x(512,4)        ! data: (sample,channel)
type(multitaper_matrix_result) :: spectrum
type(multiple_coherence_result) :: out

! Fill x with observed data.
spectrum = multitaper_spectral_matrix(x,3.5_dp,kmax=6, &
    weighting='adaptive',cv_shrinkage=.true.)
out = multiple_coherence_matrix(spectrum,cv_shrinkage=.true.,ridge=0.02_dp)
print *, shape(out%coherence)          ! (frequency,channel)
print *, shape(out%conditional_power) ! (frequency,channel)
print *, out%valid                     ! frequency-by-frequency inversion status
```

Or estimate spectral matrices and multiple coherence together:

```fortran
use scifort_signal, only : dp, multitaper_multiple_coherence, &
    multitaper_tf_multiple_coherence, multiple_coherence_result, &
    multiple_coherence_tf_result

real(dp) :: x(512,4)
type(multiple_coherence_result) :: static
type(multiple_coherence_tf_result) :: moving

static = multitaper_multiple_coherence(x,3.5_dp,kmax=6, &
    weighting='eigen',auto_shrinkage=.true.,ridge=0.02_dp)

moving = multitaper_tf_multiple_coherence(x,3.5_dp,nperseg=128, &
    noverlap=64,kmax=6,weighting='adaptive', &
    cv_shrinkage=.true.,ridge=0.02_dp)
print *, shape(moving%coherence) ! (frequency,time,channel)
```

The same entry points accept `complex(dp)` input. `axis=2` selects samples
on dimension 2 instead of 1. The direct multitaper functions support `fs`,
`nfft`, `detrend`, `detrend_mode`, and equal/eigen/adaptive DPSS weighting.
The time-frequency version accepts `noverlap` and `nperseg` as in the existing
multitaper spectral matrix spectrogram.

### Output fields

- `coherence`: squared multiple coherence, in `[0,1]`, with dimensions
  `(frequency,channel)` or `(frequency,time,channel)`.
- `conditional_power`: unmodeled power, `1/P_ii`, in the same physical power
  units as the spectral matrix.
- `explained_power`: regularized diagonal PSD minus conditional power.
- `frequencies` and optional `times`: copied frequency/segment-center grids.
- `shrinkage_used`: the actual coefficient at each frequency (and time cell).
- `valid`: per-frequency or per-frequency/time inversion validity.
- `status`: `partial_coherence_success`, `partial_coherence_invalid_input`,
  or `partial_coherence_singular` from the shared precision kernel.

Invalid inputs return zero-size arrays and the invalid-input status. If a
frequency/time cell is singular, its multiple coherence and residual/explained
spectra are NaN; other valid cells are retained. A zero spectral matrix is
singular even when relative ridge loading is requested, because its trace is
zero.

## Regularization and interpretation

All outputs are computed for the **conditioned** spectral matrix

\[
\widetilde S_{ij} = \begin{cases}
S_{ij} + \lambda\bar S, & i=j,\\
(1-\alpha)S_{ij}, & i\ne j,
\end{cases}
\quad \bar S=\frac1p\sum_i S_{ii},
\]

with relative `ridge=lambda >= 0` and diagonal-target
`shrinkage=alpha in [0,1]`. Shrinkage is applied to the off-diagonals
before trace ridge. If neither option is provided, the method uses the
unregularized original spectrum.

`auto_shrinkage=.true.` uses the existing Gaussian-Wishart effective-sample
plug-in, supported for equal and eigen tapers; at matrix level it needs
`effective_tapers=`. `cv_shrinkage=.true.` reuses the per-frequency
leave-one-taper-out coefficients stored in the spectral matrix, including
adaptive taper weights. These flags are mutually exclusive and cannot be
combined with explicit `shrinkage=`; all existing validation rules apply.

**Important:** `conditional_power` and `explained_power` refer to the
**regularized** matrix, not to unbiased conditional prediction error of the
unknown true spectrum. Trace ridge increases the target channel's diagonal,
so it can change these power quantities even if its cross spectrum is zero.
Shrinkage will generally reduce multiple coherence. For a rank-deficient
spectral matrix, use sufficient regularization if an inverse is desired.

No confidence limits, p-values, or calibrated significance inferences are
reported. Multitaper observations may remain dependent under colored noise;
frequency-local CV chooses predictive shrinkage rather than guaranteeing an
unbiased covariance-risk optimum. See [CV shrinkage](signal_cv_shrinkage.md).

## Reproducible validation

```bash
python tools/run_signal_multiple_coherence.py
python tools/check_signal_multiple_coherence_public.py
```

The native checked and optimized regressions cover known two-channel
examples, complex-phase cross spectra, singular cells, input validation,
real/complex static and moving-window functions, automatic fixed-weight
shrinkage and adaptive CV. The independent NumPy reference compares
**700 frequency/time/channel/configuration cases** against Schur complements
obtained by independently solving a `(nchannel-1)` regression system; the
maximum normalized difference observed in the development check was below
`3e-15`. Public-API compilation covers the complete `scifort_signal` module
closure. NumPy and SciPy are development/test tools only; no new runtime
Fortran dependencies are required.
