# Automatic block length for canonical spectral moving bootstrap

**SciFort extension**; no corresponding `scipy.signal` API. This is an
**exploratory selection heuristic**, not a statistically guaranteed optimal
block length or a way to calibrate confidence intervals automatically.

## API

`select_canonical_block_length(spectra, max_lag, max_block_length)` accepts a
`multitaper_spectral_matrix_spectrogram_result` with nonoverlapping
window-level complex cross-spectral matrices of shape
`(frequency,window,channel,channel)` and increasing times. It returns a
`canonical_block_length_result` with:

- `block_length`, the selected number of complete windows in each bootstrap block;
- `n_windows`, `max_lag`, `lags_used`, `maximum_block_length`;
- `autocorrelations(:)`, the pilot's lag-1 through lag-`max_lag` normalized
  Hermitian Frobenius autocovariances;
- `first_lag_correlation`, `long_run_variance_ratio`, and
  `weighted_lag_ratio`, diagnostics for selection;
- `status`, `partial_coherence_success` (0) or invalid input (1).

**New:** Omitting `block_length=` from
`canonical_moving_block_bootstrap(...)` or either real/complex
`multitaper_canonical_moving_block_bootstrap(...)` generic invokes this
selector automatically. Explicit positive lengths preserve the existing
resampling and RNG behavior *exactly*. An explicitly supplied length of
zero remains an error. The returned bootstrap result adds
`automatic_block_length` and `selection_lags` metadata.

```fortran
use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
    multitaper_canonical_moving_block_bootstrap

type(canonical_moving_block_bootstrap_result) :: uncertainty
real(dp) :: x(4096,4)

! Supply x; this call forms 32 nonoverlapping 128-sample spectral windows.
uncertainty = multitaper_canonical_moving_block_bootstrap( &
    x,3.5_dp,128,[1,2],[3,4], &
    n_resamples=999,seed=2026,kmax=6, &
    weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp)

if (uncertainty%status == 0) then
    print *, 'selected block length:', uncertainty%block_length
    print *, 'lags retained:', uncertainty%selection_lags
end if
```

For precomputed data, the selector may be called directly to inspect the pilot
before choosing or overriding its length. When `bands=` is supplied to the
bootstrap call, the selector uses the *integrated window spectral matrices*
for the requested bands, rather than the original per-frequency matrices.

## Numerical selection rule

For window-level complex Hermitian spectral matrices `S_t(f)`, compute their
window-wise means and centered matrices `X_t(f)`. The lag-`k` pilot
Frobenius autocovariance is

\[
\widehat\gamma_k = \frac{1}{W-k} \sum_{t=1}^{W-k}
    \sum_{f,i,j} \operatorname{Re}\{\overline{X_t(f)_{ij}} X_{t+k}(f)_{ij}\}.
\]

`k=0` uses denominator `W`, and the pilot correlation is
`rho_k = gamma_k/gamma_0`. The sum includes *whole* cross-spectral matrices;
it never assumes that different frequency bins or channels are independent.
To avoid numerical overflow, a common scale is removed before covariance
products. An exactly constant series of matrices selects block length one.

By default, the pilot evaluates up to
`min(W-2,max(1,floor(sqrt(W))))` lags. It retains the initial consecutive
lags with `abs(rho_k) > 2/sqrt(W)`; the first smaller correlation stops
retention. Negative dependence is retained with its sign. If no lag is
retained, select `L=1`.

For `m` retained lags, apply Bartlett weights `w_k = 1-k/(m+1)` and compute

\[
R=1+2\sum_{k=1}^{m} w_k \rho_k,\qquad
G=2\sum_{k=1}^{m} k\,w_k \rho_k.
\]

When `R` is positive, the heuristic selection is

\[
L = \operatorname{clip}_{[1,L_{max}]}\left(
     \operatorname{round}\left[\left(\frac{3W}{2}
       \left(\frac{G}{R}\right)^2\right)^{1/3}\right]\right).
\]

The default `L_max=max(1,floor(W/3))` ensures at least three nominal drawn
blocks. When the estimated long-run variance is nonpositive or numerically
unstable, the algorithm falls back to `L=1`. `max_lag` and
`max_block_length` may be supplied to the standalone selector to investigate
sensitivity, but the integrated automatic bootstrap uses the documented
defaults. The constant `3/2`, the significance cutoff, and the spectral-matrix
Frobenius aggregation are **heuristic engineering choices**, not a theorem
establishing statistical optimality for canonical coherence.

## Statistic-specific influence alternative

For coherence calculations that may be dominated by unrelated high-power
channels, use `select_canonical_influence_block_length` or pass
`selection_method='influence'` to the moving-block bootstrap. The latter
keeps the original `selection_method='spectral'` default unchanged. See the
full [influence-targeted pilot documentation](signal_canonical_influence_block_length.md),
including its delete-one-window formula, mode/frequency target options,
controlled numerical tests and exploratory VAR(1) coverage comparison.

## Statistical limitations

The selector is based on variability of the **entire spectral matrix**, not
on an influence-function estimate targeted to an individual canonical mode.
Spectral components with high power or variance may dominate, and differing
physical channel units can influence the selected length. A nearly uncorrelated
pilot can choose `L=1` even when particular nonlinear spectral statistics
retain temporal dependence. Short time series make pilot autocorrelations
very noisy. Stationarity, block-length sensitivity, mode separation, and
regularization bias all remain important. No automatic selector can
establish coverage without further assumptions and empirical validation.

## Reproducible numerical and coverage tests

```bash
python tools/run_signal_canonical_auto_block.py
python tools/check_signal_canonical_auto_block_public.py
```

The checked and optimized test runners compare independent NumPy pilot
correlations and selected lengths against native Fortran, exercise real and
complex adaptive multitaper and frequency-band APIs, and rerun the earlier
moving-bootstrap and block-jackknife regressions. A reproducible **simulation
study** draws a bivariate Gaussian AR(1) process with common autoregressive
coefficient `phi=0.94` and innovation correlation `rho=0.6`. Its *population*
squared cross-channel coherence is `rho**2=0.36` at every frequency.
For each replicate, 48 nonoverlapping 12-sample tapered windows yield
window-level spectra at frequency zero. The study compares 95% percentile
intervals from 99 bootstrap resamples under automatic and length-one blocks.
It is **a diagnostic, not a calibration certificate**: even with correct
population coherence, both coverage estimates may be considerably below
95% because of finite spectral-window bias, temporal dependence, and the
choice of percentile interval. The script prints the actual Monte Carlo
coverage, widths, and selected block-length distribution rather than
assuming an improvement.

Mathematical background for block bootstrap and n^(1/3) block-length
scaling: Künsch (1989), *The Jackknife and the Bootstrap for General
Stationary Observations*, *Annals of Statistics* 17(3), 1217–1241.
This implementation and pilot cutoff are original SciFort choices, not
translated code from that reference.
