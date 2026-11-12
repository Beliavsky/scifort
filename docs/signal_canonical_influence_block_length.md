# Statistic-specific block-length selection for canonical coherence

**SciFort extension** (not a `scipy.signal` API). This is an **exploratory**
option for the canonical-coherence moving-block bootstrap, not a consistent
optimal-block selector or a calibrated uncertainty procedure.

## Motivation and APIs

The existing automatic pilot
[`select_canonical_block_length`](signal_canonical_auto_block_length.md)
measures autocovariance of *all* cross-spectral matrix entries, aggregated by
Frobenius products. It can be dominated by a high-power channel that has no
role in the group coherence actually being bootstrapped.

The new `select_canonical_influence_block_length` estimates temporal dependence
from **leave-one-window influence scores for the requested canonical statistic**.
The estimator is available directly:

```fortran
use scifort_signal, only : dp, canonical_block_length_result, &
    select_canonical_influence_block_length

type(canonical_block_length_result) :: pilot

pilot = select_canonical_influence_block_length(window_spectra, &
    [1,2],[3,4],conditioned_on=[5],ridge=0.02_dp,shrinkage=0.10_dp, &
    mode_index=1,frequency_index=2)
if (pilot%status == 0) then
    print *, pilot%block_length, pilot%lags_used
    print *, pilot%influence  ! one scalar for each complete window
end if
```

Or request it from `canonical_moving_block_bootstrap` (precomputed matrices)
or either real/complex `multitaper_canonical_moving_block_bootstrap` overload:

```fortran
result = multitaper_canonical_moving_block_bootstrap( &
    x,3.5_dp,128,[1,2],[3,4],n_resamples=999,seed=2026, &
    conditioned_on=[5],kmax=6,weighting='adaptive',ridge=0.02_dp, &
    shrinkage=0.10_dp,selection_method='influence',selection_mode=1, &
    selection_frequency=2)
```

**Omitting `block_length`** is required for automatic selection. The default
`selection_method='spectral'` is unchanged. `selection_method='influence'`
activates the new pilot, while explicit `block_length` plus `'influence'` is
rejected to avoid suggesting that selection occurred. An unknown method or
`selection_mode`/`selection_frequency` supplied with `'spectral'` is rejected.
All new optional arguments are at the *end* of the original argument list.

Target options:

- `selection_mode=0` (default) selects determinant-based **group coherence**,
  \(1-\prod_m (1-r_m^2)\). A positive mode index selects the indicated
  **squared canonical coherence** in descending-strength order.
- `selection_frequency=0` (default) uses the arithmetic mean of the selected
  statistic over the frequency bins in that call. A positive index chooses a
  single frequency/bin; if `bands=` is supplied, this is a **band index**.
  The bins enter the target statistic together, not as independent windows.
- `conditioned_on`, fixed `ridge` and fixed diagonal `shrinkage` are identical
  for the full estimate and every deleted-window fit. For `bands`, complex
  spectral matrices are integrated **inside each window before** the pilot.
- The standalone selector also accepts `max_lag` and `max_block_length`,
  governed by the same validation and defaults as the original spectral pilot.

The `canonical_block_length_result` adds `mode_index`, `frequency_index`,
`full_statistic` and, for the influence method, `influence(:)`; all previous
fields (`autocorrelations`, `lags_used`, `first_lag_correlation`,
`long_run_variance_ratio`, etc.) are still available. The bootstrap result
adds `selection_method` (`'manual'`, `'spectral'` or `'influence'`),
`selection_mode`, `selection_frequency`, and optional
`selection_influence(:)` in addition to its earlier selection metadata. The
influence array is **unallocated** for explicit and spectral-pilot cases.

## Mathematical construction

Let \(S_t(f)\) be the full complex Hermitian window cross-spectrum, with
\(W\) complete nonoverlapping windows, and let

\[
  \overline S(f)=\frac1W\sum_{t=1}^W S_t(f),\qquad
  \overline S_{(-t)}(f)=\frac{W\overline S(f)-S_t(f)}{W-1}.
\]

Let \(h\) compute the desired conditional group coherence or isolated
canonical squared root from a mean spectral matrix, including the specified
fixed ridge/shrinkage. If several bins are selected, \(h\) denotes the mean
of those **within-frequency** statistics. The delete-one-window scores are

\[
  I_t=(W-1)\left[h(\overline S)-h(\overline S_{(-t)})\right].
\]

These differ from jackknife pseudo-values by the *same constant* for all
\(t\), so their centered lagged autocovariances are identical. They are
numeric delete-one sensitivity estimates, **not** claimed to be exact
population influence functions. A refit is performed for each window, with
fixed regularization and no repeated taper estimation.

The existing Bartlett-weighted plug-in heuristic is then applied to the
**single real series** \(I_t\). In particular,

\[
  \hat\gamma_k=\frac1{W-k}\sum_{t=1}^{W-k}
    (I_t-\overline I)(I_{t+k}-\overline I),\qquad
  \hat\rho_k=\hat\gamma_k/\hat\gamma_0.
\]

See [automatic block length](signal_canonical_auto_block_length.md) for the
pilot-lag cutoff, long-run variance ratio, length formula and caps.
For compatibility, SciFort reuses *exactly* the existing normalized pilot
by presenting the scalar scores as a one-entry pseudo-spectrum; no new
quadrature or statistical constant is introduced. A constant influence
series selects length one.

**Mode validity:** for a specific mode, the selected root must be finite,
strictly interior to \((10^{-6},1-10^{-6})\), and separated from adjacent
sorted roots by more than \(10^{-6}\) in the full and every delete-one fit.
Invalid or singular selected bins are rejected rather than silently changed
by dropping individual windows. For aggregate group coherence, distinct root
orientation is not required. The estimator can be more expensive than the
matrix pilot: it requires \(W+1\) canonical fits across the selected bins.

## Validation and simulation

```bash
python tools/run_signal_canonical_influence_block.py
python tools/check_signal_canonical_influence_block_public.py
```

Native checked/optimized Fortran tests cover analytic complex two-channel
coherence, an irrelevant slow channel, manual-versus-automatic reproducibility,
real/complex direct multitaper calls, adaptive weights, band integration,
mode/frequency targeting and rejected inputs. Independent NumPy calculations
compare **40** conditioned/unconditioned complex scenarios, including first
and second canonical modes, **1,535** leave-one scores, selected block lengths,
full statistics and pilot correlations. Maximum observed absolute deviation
in that deterministic comparison was approximately `9.7e-13`.

A further reproducible **exploratory** three-channel VAR(1) study generates
stationary Gaussian series with a high-variance independent nuisance channel.
The population zero-frequency squared coherence is computed from the VAR
transfer matrix and innovation covariance. For each scenario, 48 Monte Carlo
trials use 48 nonoverlapping 12-sample tapered windows and 99 circular
bootstrap draws per fit. These are small simulations, not coverage guarantees.

| Scenario | Population coherence | Spectral median length | Influence median length | Empirical 95% coverage: spectral / influence / length 1 |
|---|---:|---:|---:|---:|
| Persistent target, white nuisance | 0.6101 | 1 | 2 | 0.500 / 0.542 / 0.500 |
| Weaker target, persistent nuisance | 0.3798 | 2 | 1 | 0.917 / 0.917 / 0.938 |

The first setting illustrates how the influence pilot can find target
dependence that the matrix pilot misses. It also shows **serious
undercoverage**: statistic-specific selection does *not* correct the bias
of short-window spectral coherence or make naive percentile intervals
valid. In the second setting, the influence pilot avoids some nuisance
dependence without improving observed coverage. Different VAR parameters,
window sizes, frequency bands or bootstrap draws may change all results.

The separate study implementation is
`tools/signal_reference/study_canonical_influence_var.py` with its Fortran
driver `canonical_influence_coverage_driver.f90`. A convenience wrapper runs
it with the main test suite.

## Limitations

- **Stationarity** of the window-level matrix sequence and a reasonable
  short-memory approximation are still required for block bootstrap reasoning.
- A scalar pilot targets the selected aggregate or root; a *single* block
  length may not serve all reported modes and frequencies equally well.
  In particular, choosing `selection_mode=1` does not give simultaneous
  inference for mode 2.
- Sorted canonical roots can cross; rank-specific influence estimation
  excludes near ties but does not track evolving canonical directions.
- The 95% percentile intervals are unstudentized and uncalibrated; spectral
  estimation bias, weak modes, long-range dependence and uncertain pilot
  lengths can cause undercoverage.
- Missing-window, overlapping-window and irregular-time dependence are
  outside the estimator's intended assumptions. `times` must increase, but
  increasing timestamps alone do not certify nonoverlapping windows.
- The pilot retains the same engineering cutoff and \(W^{1/3}\) scaling as
  the original selector; it is **not** proven MSE-optimal for canonical
  coherence, especially under adaptive per-window taper weights.

All Fortran code and tests were written for SciFort, without third-party
source imports. NumPy/SciPy are optional independent-test dependencies,
not Fortran runtime dependencies.
