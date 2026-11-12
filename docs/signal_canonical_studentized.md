# Studentized moving-block bootstrap intervals for canonical coherence

**SciFort extension**, not a `scipy.signal` function. This is an opt-in,
exploratory *pointwise* 95% confidence-interval construction for squared
canonical coherence modes and determinant-aggregate group coherence, using
**block-jackknife standard errors inside moving-block bootstrap replicates**.
It does not guarantee 95% coverage under serial dependence, model
misspecification, finite sample bias, shrinkage, or mode crossings.

## Fortran interface

Use `interval_method='studentized'` with the existing spectrum-level
`canonical_moving_block_bootstrap` or real/complex
`multitaper_canonical_moving_block_bootstrap` functions. It is a trailing,
optional argument: prior positional calls and the default
`interval_method='percentile'` are unchanged. Other supported methods are `'basic'`, `'bc'`, and
`'bca'` (see [block-aware BCa](signal_canonical_bca.md)).

```fortran
use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
    multitaper_canonical_moving_block_bootstrap

type(canonical_moving_block_bootstrap_result) :: result
real(dp) :: x(4096,5)

! Populate x with five time series; 4096/128 = 32 complete windows.
result = multitaper_canonical_moving_block_bootstrap( &
    x,3.5_dp,128,[1,2],[3,4],block_length=4, &
    n_resamples=999,seed=2026,conditioned_on=[5], &
    kmax=6,weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp, &
    interval_method='studentized')

if (result%coherence_valid(1)) then
    print *, result%coherence_lower95(1),result%coherence_upper95(1)
    print *, 'Original block-jackknife SE:',result%coherence_studentizing_stderr(1)
end if
```

Additional outputs are allocated **only** when `'studentized'` is selected:

- `studentizing_stderr(f,mode)`: original-sample block-jackknife standard
  errors for individual squared canonical roots.
- `coherence_studentizing_stderr(f)`: original-sample block-jackknife standard
  errors for aggregate group coherence.

Existing `stderr`/`coherence_stderr` **still report bootstrap standard
deviations** of the point estimates, not the studentizing jackknife standard
errors. Bootstrap bias, full-sample estimates, sampling diagnostics,
quantiles, seed behavior, and all older interval methods are preserved.
The two kinds of standard error need not agree.

## Mathematics

Let `W` be the number of complete original nonoverlapping spectral-estimation
windows, `L` the moving block length, and `K=ceil(W/L)` the number of
contiguous deletion blocks. All spectra are full complex Hermitian
cross-spectral matrices, not individual FFT bins or taper vectors. For each
frequency/band, the statistic is `theta = h(mean(S_window))`, with fixed
ridge and shrinkage settings. The original and each bootstrap sample
are partitioned into `K` consecutive blocks. The last can be shorter.

For block `k` of size `n_k`, define the delete-one-block refit
`theta_(-k)` using the other `W-n_k` windows. Use weighted pseudovalues
and their jackknife variance:

\[
z_k=\frac{W\widehat\theta-(W-n_k)\widehat\theta_{(-k)}}{n_k},\qquad
w_k=\frac{n_k}{W},\qquad \bar z=\sum_kw_kz_k,
\]
\[
\widehat {\mathrm{SE}}_{\mathrm{block}}^2 =
\frac{\sum_k w_k^2(z_k-\bar z)^2}{1-\sum_kw_k^2}.
\]

With equal block sizes, this is the ordinary delete-one-group jackknife
variance, `sum((z_k-mean(z))**2)/(K*(K-1))`. With an incomplete final block,
the weighted extension is an **approximation**, not an exact finite-sample
variance identity for nonlinear statistics.

For each moving-block resample `b`, calculate `theta*_b` and re-estimate
`SE*_b` by deleting complete **drawn** blocks (not by deleting one
individual window). The studentized pivot and reversed-tail interval are

\[
T_b^*=\frac{\widehat\theta_b^*-\widehat\theta}
                 {\widehat{\mathrm{SE}}_b^*},\qquad
\left[\widehat\theta-\widehat{\mathrm{SE}}_0 Q_T(0.975),
      \widehat\theta-\widehat{\mathrm{SE}}_0 Q_T(0.025)\right],
\]

where `Q_T` is the linearly interpolated empirical quantile. Bounds are
intersected with the coherence parameter space `[0,1]`; an empty intersection
is invalid. The same local Park–Miller RNG, selected blocks, and original
spectral estimator are used by all interval methods.

## Failure conditions and computational cost

- At least **four deletion blocks** are needed for studentization, even
  though ordinary percentile bootstrap calls allow fewer blocks. The
  automatically selected or user-supplied length may not satisfy this.
- When any needed original or replicate block-jackknife standard error is
  zero, numerically unidentifiable, nonfinite, or its refit is singular, the
  corresponding **interval is invalid**, with `NaN` endpoints and a false
  `valid`/`coherence_valid` flag. This does not erase an available bootstrap
  standard deviation or bias estimate.
- Aggregate group coherence also needs to remain in the interior of `(0,1)`
  beyond `gap_tolerance` in the original and delete-block refits.
- Individual squared canonical modes must remain interior to `(0,1)`
  and separated from neighboring modes by `gap_tolerance` in the full
  estimate, every bootstrap replicate, and its delete-block refits.
  Nearly tied or zero modes lack a uniquely trackable rank-specific pivot.
  Group-coherence intervals do not require mode separation, but still need
  positive studentizing variances and interior aggregate strengths.
- `status=partial_coherence_singular` may mean invalid intervals at some
  frequencies, not necessarily a failure of the full-sample estimator.
  Examine individual validity flags and the two kinds of standard error.
- A studentized call evaluates the canonical estimator for each of `K`
  delete-block refits inside each of `B` bootstrap replicates, plus the
  original `K` deletions: roughly `O(B*K)` extra spectral refits. It also
  temporarily stores one cross-spectrum sum per drawn block. This can be
  expensive for many channels/frequencies/windows.

The method assumes stationary dependence that can be reasonably represented
by moving blocks. It does **not** establish independence of adjacent blocks;
overlapping candidate blocks induce dependence in the empirical pool. Fixed
ridge/shrinkage may distort interval coverage and are not refitted during
resampling. The method does not return simultaneous frequency bands,
studentized subspace *orientation* intervals, or calibrated hypothesis tests.
A pointwise interval may not contain the observed point estimate.

## Reproduce numerical tests

```bash
python tools/run_signal_canonical_studentized.py
python tools/check_signal_canonical_studentized_public.py
```

The test runner builds checked and optimized GNU Fortran, runs new and
previous native regressions, and independently checks studentization,
block resampling and older interval methods using NumPy complex matrix
linear algebra and SciPy spectral/quantile references. Tests include
conditional channel groups, multiple modes, frequency integration,
real/complex adaptive multitaper APIs, circular blocks, incomplete last
blocks, constant data, singular mode directions and too few blocks. The
reference regenerates all bootstrap indices independently using the
published integer RNG formula. It does not call the Fortran estimator for
its expected numeric results.

## Exploratory VAR(1) coverage diagnostic

`tools/signal_reference/study_canonical_studentized_var.py` runs the
pre-existing two-scenario correlated three-channel VAR(1) experiment with
**the same 199 circular-block bootstrap draws and 32 datasets per scenario**
for percentile, basic, BC and studentized intervals. It compares intervals
for the smoothed finite-window estimator against the *population
zero-frequency squared coherence*, so bias caused by finite spectral windows
limits interpretability. The observed fraction of intervals containing the
population target was:

| Model | Percentile | Basic | BC | Studentized | Nominal |
|---|---:|---:|---:|---:|---:|
| Persistent target, weak nuisance | 46.9% | 40.6% | 50.0% | 71.9% | 95% |
| Moderate target, persistent nuisance | 96.9% | 78.1% | 90.6% | 93.8% | 95% |

Studentization increased interval widths (approximately 0.58 versus 0.42
in the persistent-target scenario), and did **not** achieve 95% coverage
there. With only 32 trials, these proportions have large Monte Carlo
uncertainty; the comparison is *diagnostic*, not evidence of a reliable
method ranking or a calibrated confidence procedure. Future studies should
vary block length, window duration, VAR dependence and the estimand.

## Background and provenance

For the bootstrap-t/studentized interval, see B. Efron and R. Tibshirani
(1993), *An Introduction to the Bootstrap*. For dependent-data moving-block
resampling, see H. R. Künsch (1989), *The Jackknife and the Bootstrap for
General Stationary Observations*, *Annals of Statistics*, 17(3), 1217–1241.
These are mathematical references; **no implementation code was copied**.
See `CODE_PROVENANCE.md` for the original SciFort code and test provenance.
