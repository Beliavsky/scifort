# Moving-block bootstrap for canonical spectral coherence

**SciFort extension**; not a function in `scipy.signal`. This is an exploratory
uncertainty estimator for multichannel spectral analysis. Its *pointwise* 95%
percentile intervals (or the opt-in basic and bias-corrected alternatives) are not guaranteed to have 95% frequentist coverage,
particularly when there are few spectral windows, the blocks are too short,
regularization is data-dependent, or canonical modes nearly coincide. Do not
interpret the intervals as significance tests or simultaneous frequency bands.

## API

Optional trailing `interval_method='percentile'` (the unchanged default),
`'basic'`, `'bc'`, or `'studentized'` chooses the pointwise bootstrap interval. See
[bias-aware interval methods](signal_canonical_interval_methods.md) and
[studentized intervals](signal_canonical_studentized.md) for formulas and failure
conditions. The default percentile intervals are not BCa intervals or calibrated statistical tests; see [BCa](signal_canonical_bca.md) for the optional accelerated method.

The spectrum-level API is
`canonical_moving_block_bootstrap(spectra,group_a,group_b,block_length,...)`.
The `block_length` argument is now **optional**: if omitted, SciFort selects a
heuristic length from pilot spectral-window autocovariances as described in
[automatic block length](signal_canonical_auto_block_length.md). The default
`selection_method='spectral'` remains unchanged; alternatively,
`selection_method='influence'` uses delete-one-window sensitivities of the
requested coherence statistic. See
[influence-targeted block selection](signal_canonical_influence_block_length.md).
Optional `selection_mode` (0 for group coherence, positive for a squared
canonical root) and `selection_frequency` (0 for a mean over bins, positive
for a specific frequency or band) only apply to `'influence'` mode.
`spectra` is a `multitaper_spectral_matrix_spectrogram_result` containing
Hermitian cross-spectral matrices with dimensions
`(frequency,window,channel,channel)`, strictly increasing `times`, and
**nonoverlapping** original spectral windows. The existing spectral-window
estimator's `noverlap=0` behavior is required; a positive time-step alone
cannot prove that caller-supplied windows do not overlap.

The `multitaper_canonical_moving_block_bootstrap` generic accepts real or
complex 2-D time-series data and computes those complete, nonoverlapping
spectral windows itself. Existing DPSS `weighting='equal'`, `'eigen'`, and
`'adaptive'`, `axis`, detrending, `fs`, `nfft` and `kmax` options are
forwarded to the SciFort multitaper estimator. Optional `bands(2,nband)`
can integrate each window's **complex cross-spectral matrix** before
resampling and fitting, or `average=.true.` can use its band average.

```fortran
use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
    multitaper_canonical_moving_block_bootstrap

type(canonical_moving_block_bootstrap_result) :: inference
real(dp) :: x(4096,5)

! Populate x; 4096/128 = 32 disjoint spectral-estimation windows.
inference = multitaper_canonical_moving_block_bootstrap( &
    x,3.5_dp,128,[1,2],[3,4],block_length=4, &
    n_resamples=999,seed=2026,conditioned_on=[5], &
    kmax=6,weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp)

if (inference%status == 0) then
    ! Use inference%valid(f,mode) before interpreting rankwise intervals.
    ! inference%lower95 and %upper95 are bootstrap percentile bounds.
    ! inference%coherence_lower95 and %coherence_upper95 are aggregate bounds.
end if
```

For precomputed spectra, use:

```fortran
inference = canonical_moving_block_bootstrap(window_spectra, &
    [1,2],[3,4],block_length=4,n_resamples=999, &
    conditioned_on=[5],seed=2026)
```

An explicitly supplied `block_length` satisfies `1 <= block_length < n_windows` and must
leave at least two drawn blocks; it **need not divide** the number of windows.
`n_windows >= 4`. The number of draws `n_resamples` defaults to 399, with a
minimum of 20 (for useful 2.5% and 97.5% tails, 999 or more is recommended).
`seed` defaults to 12345 and must be from 1 through 2147483646 inclusive.
The random-number stream is local, fixed and deterministic; `random_seed`
and the global Fortran `random_number` sequence are not touched. The seed
is **not cryptographically secure**.

`circular=.false.` (the default) selects overlapping *possible* block
starts from `1:(W-L+1)`, without wrapping. `circular=.true.` permits any
start from `1:W` and wraps the block to the beginning. The latter may
reduce endpoint imbalance under approximate stationarity but asserts a
periodic boundary for resampling, not for the underlying signal.

## Resampling and statistics

For each frequency (or preintegrated band), let `S_w(f)` be the full
cross-spectral matrix from nonoverlapping window `w`. The full estimate is

\[
\widehat\theta(f)=h\left(\frac{1}{W}\sum_{w=1}^{W} S_w(f)\right),
\]

where `h` computes conditional canonical squared coherence roots, in
**descending strength order**, or their determinant-based aggregate
`1-product(1-root)`. The optional ridge and shrinkage parameters are held
fixed in the full fit and every bootstrap refit.

For replicate `r`, draw `ceil(W/L)` starting indices independently,
**with replacement**, and take the corresponding `L` consecutive entire
spectral windows from each start. Concatenate these blocks in draw order;
truncate only the last block to yield exactly `W` resampled windows.
Form their mean spectral matrix and re-evaluate `h`. Within each sampled
block the original window order is preserved. All frequency bins and
channels share the **same** draw indices, preserving their cross-frequency
and cross-channel dependence in a replicate. Original frequency bins and
individual adaptive DPSS tapers are never resampled independently.

For each scalar mode or aggregate statistic with `B` successful refits:

\[
\widehat{\mathrm{SE}}_{boot}
  = \sqrt{\frac{1}{B-1}\sum_{r=1}^B(\theta_r^*-\bar\theta^*)^2},
\qquad
\widehat{\mathrm{bias}}_{boot}=\bar\theta^*-\widehat\theta.
\]

The pointwise 95% interval uses the empirical 2.5th and 97.5th percentiles
of the resampled **squared coherence** values directly, with linear
interpolation at empirical positions `p*(B-1)`. It is not a normal/Fisher
interval, bootstrap-t interval, bias-corrected accelerated interval,
or an assertion of calibrated significance. Raw bootstrap percentile
intervals are automatically within `[0,1]`, but need **not** contain the
full-sample estimate. The new basic and BC options retain that parameter
range; read their separate validity and calibration caveats in
[interval methods](signal_canonical_interval_methods.md).

`inference%valid(f,mode)` is true only if the full-sample root and **every**
resample root are finite, interior and separated from neighboring sorted
roots by `gap_tolerance` (default `1e-6`). This deliberately disallows
individual-mode inference through ambiguous crossings/ties or boundary
roots. Undefined uncertainty fields are quiet NaNs and validity flags are
false. Aggregate `coherence_valid(f)` can remain true even when individual
roots coincide. A singular spectral bin marks the overall status as
`partial_coherence_singular`, while valid bins may still contain usable
results. Invalid input returns `partial_coherence_invalid_input`.

Output metadata include `automatic_block_length`, `selection_lags`,
`selection_method`, `selection_mode`, `selection_frequency`, and optional
`selection_influence(:)` (allocated only with the influence pilot),
`n_windows`, `block_length`, `blocks_per_resample`,
`n_resamples`, `seed` and `circular`. Full-data `squared_modes` and
`coherence` are separate from bootstrap `stderr`, `bias`, `lower95`,
`upper95`, `valid`, `coherence_stderr`, `coherence_bias`,
`coherence_lower95`, `coherence_upper95` and `coherence_valid`.

## Assumptions and limitations

- Windows within a block may be serially dependent. The moving-block
  bootstrap is most defensible for an approximately stationary, weakly
  dependent sequence of spectral-estimation windows. Its consistency
  generally needs a block length that grows with sample size but remains
  small relative to the number of windows. `block_length=1` recovers the
  ordinary window bootstrap; it is **not** the delete-one-window jackknife.
- The default ordinary moving-block resampling can underrepresent windows
  near the ends. Circular blocks are an alternative, not a universal fix.
- The block length, number of resamples and choice of circular geometry can
  materially change the intervals. Examine their sensitivity in applied
  work, preferably with simulation-based coverage tests for the scientific
  process of interest.
- Spectral windows themselves may be temporally correlated even when
  nonoverlapping; multitaper smoothing and adaptive weighting do not make
  nearby frequencies independent. Resampling whole spectral windows
  preserves each window's frequency dependence but cannot guarantee
  exchangeability or calibrated confidence intervals.
- `ridge` and `shrinkage` are **fixed across refits**. If they were selected
  using the original data, the bootstrap does not propagate that selection
  uncertainty. Re-estimating adaptive weights *within* each original window
  is not part of this bootstrap: its unit is the complete per-window matrix.
- The sorted rank-specific roots can exchange identities when modes cross.
  Validity checks reject nearly tied modes, but they are not a full
  solution to mode-label uncertainty; eigenvector/subspace uncertainty
  remains future work.

Standard background: Künsch (1989), *The Jackknife and the Bootstrap for
General Stationary Observations*, *Annals of Statistics* **17**(3),
1217-1241; and Efron and Tibshirani (1993), *An Introduction to the
Bootstrap*. These references supply mathematical context, not code.

## Reproducible validation

```bash
python tools/run_signal_canonical_moving_bootstrap.py
python tools/check_signal_canonical_moving_bootstrap_public.py
```

The checked and optimized Fortran tests cover original versus circular
blocks, non-divisible lengths, deterministic seeding, global RNG isolation,
fixed regularization, conditional groups, band integration, input
validation, degenerate/singular modes, and real/complex adaptive direct
multitaper calls. The independently written NumPy/SciPy reference uses
NumPy Cholesky/SVD, independently generated block indices and quantiles,
and SciPy interpolation/trapezoidal integration; it does not use the
SciFort Fortran eigensolver or its bootstrap implementation.
