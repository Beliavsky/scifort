# Block-aware BCa intervals for canonical spectral coherence

**SciFort extension, not a `scipy.signal` API.** The optional
`interval_method='bca'` adds bias-corrected and accelerated (*BCa*)
pointwise bootstrap intervals for squared canonical coherence roots and
aggregate conditional group coherence. They are **exploratory**, not
calibrated confidence intervals or hypothesis tests: classical iid BCa
coverage theory does not establish second-order accuracy when windows are
serially dependent, resampled in overlapping moving blocks, or subjected to
spectral smoothing and fixed regularization.

## Interface

Both `canonical_moving_block_bootstrap` (precomputed nonoverlapping-window
spectral matrices) and the real/complex generic
`multitaper_canonical_moving_block_bootstrap` accept the existing trailing
optional `interval_method` keyword. Supported values are
`'percentile'` (default), `'basic'`, `'bc'`, `'bca'`, and `'studentized'`.
All preceding positional calls, legacy interval options, seeded bootstrap
draws, point estimates, standard deviations, and bias fields are preserved.

```fortran
use scifort_signal, only : dp, canonical_moving_block_bootstrap_result, &
    multitaper_canonical_moving_block_bootstrap

type(canonical_moving_block_bootstrap_result) :: result
real(dp) :: x(4096,5)

! Populate x with five channels of observed samples.
result = multitaper_canonical_moving_block_bootstrap( &
    x,3.5_dp,128,[1,2],[3,4],block_length=4, &
    n_resamples=999,seed=2026,conditioned_on=[5], &
    kmax=6,weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp, &
    interval_method='bca')

if (result%coherence_valid(1)) then
    print *, result%coherence_lower95(1),result%coherence_upper95(1)
    print *, 'BCa acceleration: ',result%coherence_bca_acceleration(1)
end if
```

The following arrays are allocated **only** when `'bca'` is requested:

- `bca_acceleration(f,mode)` holds block-jackknife acceleration for an
  individual squared canonical mode, or quiet NaN when it is not identifiable.
- `coherence_bca_acceleration(f)` contains the corresponding acceleration
  for the determinant-aggregate group coherence.

Intervals remain in `lower95`/`upper95` and
`coherence_lower95`/`coherence_upper95`; invalid endpoints are NaN, with
`valid`/`coherence_valid` false. An invalid aggregate interval can set
`status=partial_coherence_singular` while retaining valid bootstrap moments.
The original `studentizing_stderr` outputs remain specific to the separate
`'studentized'` method.

## Mathematics and numerical convention

Let `W` denote the number of original **complete, nonoverlapping spectral
windows**, `L` the moving-block length, and `K=ceil(W/L)` the number of
original contiguous deletion blocks. Full complex Hermitian spectra are
aggregated and the canonical statistic `theta=h(mean(S_window))` is refitted
with the same fixed ridge/shrinkage choices for every deletion. When
`W mod L != 0`, the last block is shorter.

For original deletion block `k` of size `n_k`, calculate

\[
  w_k=\frac{n_k}{W},\quad
  z_k=\frac{W\widehat\theta-(W-n_k)\widehat\theta_{(-k)}}{n_k},\quad
  \bar z=\sum_k w_k z_k,\quad
  u_k=w_k(z_k-\bar z).
\]

The reported acceleration is

\[
 \widehat a=\frac{\sum_k u_k^3}{6(\sum_k u_k^2)^{3/2}}.
\]

For equal deletion blocks, the common weights cancel and this is the
ordinary delete-one-*block* jackknife skewness formula. With an unequal
last block, weighted pseudovalues and contributions provide a **heuristic
cluster-influence approximation**, not a claim of exact block-BCa theory.
Scale-normalizing the `u_k` before accumulating powers reduces numerical
overflow and underflow risk.

For `B` bootstrap replicates `theta_star(b)`, let
`p0=(N_less+N_tied/2)/B`, `z0=Phi_inverse(p0)` and
`z_alpha=Phi_inverse(alpha)`, where `alpha` is `.025` or `.975`.
Ties within `64*epsilon*max(1,abs(theta))` are treated as numerical ties,
matching the existing `'bc'` method. The BCa endpoint probability is

\[
  \alpha_{\mathrm{BCa}}=\Phi\left(z_0+
    \frac{z_0+z_\alpha}{1-\widehat a(z_0+z_\alpha)}\right).
\]

The interval endpoints are the linearly interpolated empirical bootstrap
quantiles at these two adjusted probabilities. This reduces to the
existing `'bc'` formula for exactly zero acceleration.

## Validity and limitations

- At least **four deletion blocks** are required. The full estimate,
  each original delete-block refit, and each bootstrap replicate must be
  numerically usable. Mode-specific intervals also require isolated roots
  away from zero and one. Aggregate intervals require interior group
  coherence but do not require individual mode separation.
- Constant or numerically indistinguishable leave-block pseudovalues give
  unidentified acceleration. Zero or unit bootstrap bias probability,
  nonfinite adjusted probabilities, an acceleration-pole crossing, or
  reversed/degenerate adjusted tail probabilities invalidate the interval.
  Do **not** silently substitute percentile intervals in these cases.
- A finite acceleration is not a measure of interval reliability. The code
  does not estimate effective independent windows, refit adaptive taper
  weights across windows, calibrate mode-subspace orientations, or supply
  simultaneous confidence bands across frequencies.
- Ridge and diagonal-target shrinkage are fixed throughout resampling;
  tuning uncertainty is not included. Spectral band integration is carried
  out before resampling, not applied to already estimated coherence.
- The BCa acceleration requires `K` additional canonical refits on the
  original data. Unlike studentization, it does **not** refit the jackknife
  for every bootstrap replicate; this is substantially less expensive for
  many replicates.

## Tests and exploratory VAR study

```bash
python tools/run_signal_canonical_bca.py
python tools/check_signal_canonical_bca_public.py
```

The new reference regenerates every moving-block index with independent
Python Park-Miller arithmetic, evaluates complex-valued canonical spectra
using NumPy linear algebra, and computes delete-block accelerations and
BCa quantiles using NumPy/SciPy. It checks real/complex native API
forwarding, frequency integration, unequal last blocks, deterministic draws,
zero jackknife influence, singular/unidentifiable intervals, and prior
percentile/BC/studentized compatibility.

A reproducible *small* VAR(1) diagnostic uses 32 simulated datasets per
setting, 199 replicates, circular block length four, and identical
random draws across methods:

| Scenario | Percentile | BC | BCa | Studentized |
|---|---:|---:|---:|---:|
| Persistent target, weak nuisance | 46.9% | 50.0% | 53.1% | 71.0% |
| Moderate target, persistent nuisance | 96.9% | 90.6% | 90.6% | 93.3% |

These are fractions of **valid** intervals containing the population
zero-frequency squared coherence. The point estimator actually targets a
finite-window, frequency-smoothed quantity; this mismatch, along with
Monte Carlo variability, makes the percentages **diagnostic only**. No
claim of correct coverage or superiority should be drawn. In the studentized
column, 31/32 and 30/32 datasets produced valid intervals, respectively;
the other methods produced 32/32 in both cases. The larger and more
useful next experiment should define the population target of the actual
smoothed estimator and vary dependence, window size, block length and
bootstrap replicates over many more simulation trials.

## Provenance

Independent implementation of published mathematical BCa formulas;
no third-party code copied. References: Bradley Efron (1987),
*Better Bootstrap Confidence Intervals*, **JASA** 82(397), 171–185,
https://doi.org/10.1080/01621459.1987.10478410; H. R. Künsch (1989),
*The Jackknife and the Bootstrap for General Stationary Observations*,
**Annals of Statistics** 17(3), 1217–1241. See `CODE_PROVENANCE.md`.
