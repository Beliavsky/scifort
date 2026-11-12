# Bias-aware moving-block-bootstrap intervals for canonical coherence

**SciFort extension**, not a `scipy.signal` function. The moving-block
bootstrap now supports alternative *pointwise* interval constructions for
conditional canonical squared coherence roots and the aggregate group
coherence. These constructions do **not** guarantee calibrated frequentist
coverage under serial dependence, finite sample bias, estimator smoothing,
mode crossings, or misspecified blocks.

## Interface and compatibility

`canonical_moving_block_bootstrap` and the real/complex
`multitaper_canonical_moving_block_bootstrap` generic accept one additional
**optional trailing argument**:

```fortran
interval_method='percentile'   ! Existing default; unchanged
interval_method='basic'        ! Reverse/bootstrap-basic confidence limits
interval_method='bc'           ! Bias-corrected percentile limits (not BCa)
interval_method='bca'         ! Bias-corrected/accelerated, delete-block jackknife
interval_method='studentized'  ! Bootstrap-t, block-jackknife variances
```

A method name is lowercase and case-sensitive. Unrecognized names (including `'normal'`) yield
`partial_coherence_invalid_input`.
`result%interval_method` records the selected method. All methods use the
**same bootstrap draws**, fixed ridge/shrinkage regularization, quantile
interpolation, existing bootstrap standard errors and bias estimates. The
studentized method additionally computes a delete-block jackknife standard
error for the original estimate and every bootstrap replicate; these are
separately reported in `studentizing_stderr` and
`coherence_studentizing_stderr` when requested. The interval method can
change confidence bounds and, if an interval cannot be constructed, the
corresponding validity flag and status.
Existing calls with no `interval_method` still produce the original exact
percentile results.

For example:

```fortran
use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
    multitaper_canonical_moving_block_bootstrap

type(canonical_moving_block_bootstrap_result) :: inference
real(dp) :: x(4096,5)

! Populate x with 4096 samples in five channels.
inference = multitaper_canonical_moving_block_bootstrap( &
    x,3.5_dp,128,[1,2],[3,4],block_length=4, &
    n_resamples=999,seed=2026,conditioned_on=[5], &
    kmax=6,weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp, &
    interval_method='bc')

if (inference%coherence_valid(1)) then
    print *, inference%coherence_lower95(1), inference%coherence_upper95(1)
end if
```

All pre-existing options remain supported, including optional circular
resampling, automatic or influence-targeted block-length selection,
frequency-band integration, and complex signal input. `n_resamples >= 20` is
required as before. At least 999 replicates are recommended for estimating
tail quantiles, but increasing replicate count does not remove systematic bias.

## Equations

For an original squared mode or group coherence estimate
`theta = h(mean(S_w))` and bootstrap refits `theta_star(1:B)`, let `Q(p)`
denote the linearly interpolated empirical bootstrap quantile at probability
`p` (`p*(B-1)` positions). All these statistics take values in `[0,1]`.
For a nominal two-sided 95% interval:

- **Percentile (default):** `[Q(.025), Q(.975)]`.
- **Basic/reverse:** `[2*theta - Q(.975), 2*theta - Q(.025)]`,
  intersected with `[0,1]`. Truncation may produce an endpoint interval at
  zero or one; if the intersection is empty the bounds become NaN and the
  validity flag is false.
- **Bias-corrected (BC):** calculate
  `p0 = (number(theta_star < theta) + 0.5*number(theta_star == theta))/B`,
  `z0 = Phi_inverse(p0)`, then use
  `[Q(Phi(2*z0-z_0.975)), Q(Phi(2*z0+z_0.975))]`,
  where `Phi` is the standard normal CDF and `z_0.975 = 1.959963984540054`.
  To avoid changing the rank of a resample that reproduces the full estimate
  up to floating-point roundoff, the comparison treats differences within
  `64*epsilon(theta)*max(1,abs(theta))` as ties. The `Phi_inverse` calculation
  uses bounded scalar bisection with `erfc`; no external Fortran dependency
  is required.

BC is distinct from **BCa**: the BC acceleration parameter is zero. The
separate [block-aware BCa method](signal_canonical_bca.md) estimates it from
original delete-block jackknife influences. For bootstrap-t construction with a
block-jackknife denominator, see
[studentized intervals](signal_canonical_studentized.md). The code does not
compute simultaneous or hypothesis-test intervals. If all resamples are strictly
above or below the full estimate (so `p0` equals zero or one), BC is **not
identified**; it returns NaN bounds with a false validity flag, rather than
inventing an infinite normal quantile. Standard error and bias fields remain
available if all resampling fits succeeded. `partial_coherence_singular`
status can therefore mean either a singular fit or an unidentified aggregate
BC interval; use the per-frequency validity flags to interpret results.

Even when the block-bootstrap spectral matrices are valid, rank-specific
intervals require every bootstrap canonical mode to be interior and separated
by `gap_tolerance`; otherwise the existing `valid(f,mode)` remains false.
Aggregate group-coherence intervals may remain valid across repeated roots.
Intervals are not clipped to contain the point estimate; neither percentile
nor BC bootstrap intervals guarantee that property.

## Validation and empirical coverage

`python tools/run_signal_canonical_interval_methods.py` compiles and runs
checked and optimized GNU Fortran implementations, the existing moving-block
bootstrap and jackknife regressions, the previous independent reference, and
`tools/signal_reference/verify_canonical_bootstrap_intervals.py`.
The new independent reference generates its own Park-Miller resampling indices,
computes conditional modes with NumPy/SciPy and uses SciPy normal quantiles.
The checked build tests 120 interval-method configurations (40 base
scenarios times three methods), while the optimized build checks 36
configurations (12 base scenarios times three methods) to keep runtime practical.
The test covers real-valued squared modes from **complex Hermitian spectral
matrices**, conditioned channel groups, ordinary/circular blocks, optional
frequency bands, and all three interval constructions. Repeated resamples and
near-identical point estimates use an explicit rounding tolerance to avoid
spurious midrank discrepancies.

`python tools/check_signal_canonical_interval_methods_public.py` compiles the
complete umbrella `scifort_signal` closure and checks all current public
interfaces including the real and complex direct APIs.

The new diagnostic `tools/signal_reference/study_canonical_interval_var.py`
reuses the bivariate target inside two three-channel VAR(1) settings. It uses
32 datasets per setting, 48 nonoverlapping 12-sample Fourier windows, 199
bootstrap replicates, circular block length 4 and a fixed random seed.
The nominal 95% interval comparison against **population zero-frequency**
squared coherence was:

| Scenario | Percentile coverage | Basic coverage | BC coverage |
|---|---:|---:|---:|
| Persistent target | 46.9% | 40.6% | 50.0% |
| Moderate target | 96.9% | 78.1% | 90.6% |

This experiment is **exploratory**, not a comparative performance claim.
Coverage has high Monte Carlo error with only 32 datasets, and the windowed
single-frequency spectrum has smoothing/finite-window bias relative to the
population zero-frequency spectral target. The selected block length and
regularization are fixed and may be inappropriate. In particular, BC did
not systematically improve coverage, and none of these examples establishes
95% interval calibration. More rigorous studies require more trials,
longer windows, a fixed estimand, multiple dependence strengths, and
studentized or other adjusted constructions.

## Background and provenance

The bias-corrected bootstrap quantile transformation follows the statistical
formulation of Bradley Efron (1987), *Better Bootstrap Confidence Intervals*,
Journal of the American Statistical Association 82(397), 171-185.
Block resampling follows the mathematical setting described in
H. R. Künsch (1989), *The Jackknife and the Bootstrap for General Stationary
Observations*, Annals of Statistics 17(3), 1217-1241.
These are mathematical references: **no implementation code was copied**.
See `CODE_PROVENANCE.md`; SciPy and NumPy are used only in independent tests.
