# Band-integrated spectral matrices, partial coherence and multiple coherence

This is an **original SciFort extension**, not a `scipy.signal` function.
It aggregates an estimated multichannel cross-spectral density matrix **before**
computing band-level partial or multiple coherence. Unlike averaging squared
coherence independently at individual frequencies, it measures linear
association at the level of the jointly band-limited covariance matrix.

## Definition

Let `S(f)` be a finite Hermitian cross-spectral-density matrix in the
`(frequency, channel, channel)` format used by `multitaper_spectral_matrix`.
For the requested band `[a,b]`, `a < b`, define

\[
S_{[a,b]}=\int_a^b S(f)\,df.
\]

SciFort integrates the **piecewise-linear interpolant** of each matrix entry.
Thus, a partial segment from `f_i` to `f_{i+1}` with clipped endpoints
`u=max(a,f_i)` and `v=min(b,f_{i+1})` contributes

\[
(v-u)\big[(1-m)S(f_i)+m S(f_{i+1})\big],\quad
m=\frac{(u-f_i)+(v-f_i)}{2(f_{i+1}-f_i)}.
\]

The calculation is exact for piecewise-linear spectra and handles arbitrary
strictly distinct, finite input frequencies in **any order**, including
common two-sided FFT order `[0,positive,...,negative,...]`. Nonuniform grids,
overlapping bands, exact grid-boundary endpoints and off-grid band limits
are supported. Bands must lie entirely inside the observed frequency range:
there is no extrapolation.

By default the result represents *integrated spectral power* (covariance
units if the input `S(f)` is a density in power/Hz). With `average=.true.`,
results are divided by `b-a` and retain the input PSD's units. This scaling
is the same for every entry of a given band: multiple/partial coherence is
unchanged, while conditional/explained power changes units accordingly.

**Important:** For PSDs computed with other conventions, such as a
`spectrum` rather than `density` normalization, integration does **not**
automatically produce physically calibrated band variance. The input
frequency grid and PSD normalization must be consistent.

Once the integrated matrix is formed, the existing Hermitian Cholesky
precision solver computes

\[
R_i^2([a,b])=1-\frac{1}{S_{[a,b],ii}\,[S_{[a,b]}^{-1}]_{ii}}
\]

for multiple coherence and

\[
C_{ij\mid\mathrm{rest}}([a,b])=
-\frac{[S_{[a,b]}^{-1}]_{ij}}
{\sqrt{[S_{[a,b]}^{-1}]_{ii}[S_{[a,b]}^{-1}]_{jj}}}
\]

for complex partial coherency. The magnitude squared of the latter is
partial coherence. These measure dependence between *band-limited signals*;
they are **not** weighted averages of frequency-local coherence and do not
imply causality or statistical significance.

## Public Fortran APIs

- `integrate_spectral_bands(spectrum,bands,average)` returns a
  `multitaper_matrix_result` whose first dimension is the frequency band.
- `integrate_tf_spectral_bands(spectrum,bands,average)` returns a
  `multitaper_spectral_matrix_spectrogram_result` with dimensions
  `(band,time,channel,channel)`.
- `band_multiple_coherence_matrix`, `band_multiple_coherence_tf_matrix`,
  `band_partial_coherence_matrix`, and `band_partial_coherence_tf_matrix`
  combine integration with the corresponding previously implemented solver.
- `multitaper_band_multiple_coherence`, `multitaper_tf_band_multiple_coherence`,
  `multitaper_band_partial_coherence`, and
  `multitaper_tf_band_partial_coherence` are direct signal-level generic
  entry points for real or complex rank-two `(sample,channel)` data.

Each `bands` input has shape `(2,nband)`, and each column holds its lower
and upper frequency limits in Hz. These limits can overlap. Result
`%frequencies(ib)` is the **band midpoint** `(lower+upper)/2`, not an FFT bin
or a frequency-weighted centroid. The caller retains the exact band limits.
All time-frequency outputs preserve the original time grid.

Example:

```fortran
use scifort_signal, only : dp, multitaper_spectral_matrix, &
    integrate_spectral_bands, band_multiple_coherence_matrix, &
    band_partial_coherence_matrix, multitaper_matrix_result, &
    multiple_coherence_result, partial_coherence_result
implicit none
real(dp) :: x(512,3)
real(dp) :: bands(2,2)
type(multitaper_matrix_result) :: s, integrated
type(multiple_coherence_result) :: multiple
type(partial_coherence_result) :: partial

! Fill x with multichannel observations sampled at 100 Hz.
bands(:,1) = [5.0_dp,12.0_dp]
bands(:,2) = [12.0_dp,25.0_dp]
s = multitaper_spectral_matrix(x,3.5_dp,kmax=6,fs=100.0_dp, &
    weighting='adaptive')
integrated = integrate_spectral_bands(s,bands)
multiple = band_multiple_coherence_matrix(s,bands,ridge=0.02_dp)
partial = band_partial_coherence_matrix(s,bands,ridge=0.02_dp, &
    return_coherency=.true.)
! integrated%power:   (2,3,3), band-integrated spectral power
! multiple%coherence: (2,3), squared multiple coherence
! partial%coherency: (2,3,3), complex conditional coherency
```

Equivalent direct calls are supported:

```fortran
multiple = multitaper_band_multiple_coherence(x,3.5_dp,bands,kmax=6, &
    fs=100.0_dp,weighting='adaptive',ridge=0.02_dp)
```

`multitaper_tf_band_multiple_coherence(x,nw,nperseg,bands,...)` and the
partial-coherence counterpart work with time-frequency spectral matrices.
Use `axis=2` for `(channel,sample)` input.

### Regularization and limitations

`ridge=` and `shrinkage=` are applied **after band integration** to the
band-level matrix. Both are explicit; `ridge` is trace-relative diagonal
loading and `shrinkage` blends the cross spectrum toward the diagonal.
`return_coherency=.true.` works with the band partial-coherence API.

The existing per-frequency `cv_shrinkage` values and `effective_tapers`
metadata are intentionally **not propagated** to aggregated spectra.
Frequency-local leave-one-taper-out shrinkage and complex-Wishart
sample-count estimates are not statistically calibrated for a multibin
integral. These band APIs therefore offer only manually specified
regularization; neither `auto_shrinkage` nor `cv_shrinkage` is inferred.

An unsuccessful grid, band, or matrix validation returns a spectral result
with the existing nonzero invalid status and zero-size arrays. Ill-conditioned
or singular integrated matrices follow existing partial-coherence rules,
including `%valid` masks and NaN cells; use manual regularization where
appropriate. For strongly overlapping bands, estimates are correlated.
No confidence limits, statistical significance tests, or unqualified
causal interpretation are offered.

## Reproducible tests

```bash
python tools/run_signal_band_coherence.py
python tools/check_signal_band_coherence_public.py
```

The first runs checked (`-fcheck=all`) and optimized (`-O3 -march=native`)
GNU Fortran builds with native regressions plus the independent randomized
NumPy oracle in `tools/signal_reference/verify_band_coherence.py`. The oracle
uses `numpy.interp`, `numpy.trapezoid` and `numpy.linalg.inv`, not SciFort's
quadrature or Cholesky code. The second compiles the complete public-module
source closure and exercises the new APIs via `use scifort_signal`. No new
Fortran runtime dependencies are required.
