# Equal-block jackknife for conditional canonical spectral coherence

**SciFort extension**, not a public `scipy.signal` function. This is an
exploratory dependence-aware uncertainty estimator, **not** a calibrated
confidence procedure or statistical significance test.

`canonical_block_jackknife` applies an equal-group, delete-one-block jackknife
to canonical squared spectral coherence. Each block consists of a consecutive
run of **complete, nonoverlapping spectral-estimation windows**. Compared to
deleting single windows, the blocks can preserve short-range temporal
dependence among neighboring windows within the omitted sampling unit.
They do not remove dependence *between* blocks, and the reported standard
errors and intervals have no universal finite-sample coverage guarantee.

## Usage

```fortran
use scifort_signal, only : dp, canonical_window_jackknife_result, &
    multitaper_canonical_block_jackknife

type(canonical_window_jackknife_result) :: result
real(dp) :: x(4096,5)

! Populate x first. 4096/128 = 32 nonoverlapping spectral windows;
! group four neighboring windows per deletion: G = 32/4 = 8 blocks.
result = multitaper_canonical_block_jackknife(x,3.5_dp,128, &
    [1,2],[3,4],4,conditioned_on=[5], &
    kmax=6,weighting='adaptive',ridge=0.05_dp,shrinkage=0.1_dp)

if (result%status == 0) then
    ! result%squared_modes(f,j): full-window-mean canonical strength
    ! result%stderr(f,j), result%bias(f,j): block jackknife estimates
    ! result%lower95(f,j), result%upper95(f,j): heuristic pointwise limits
    ! result%valid(f,j): false for tied or boundary strengths
end if
```

For already estimated spectra:

```fortran
result = canonical_block_jackknife(window_spectra, [1,2], [3,4], &
    block_length=4, conditioned_on=[5], ridge=0.05_dp, shrinkage=0.1_dp)
```

`canonical_window_jackknife` and `multitaper_canonical_window_jackknife` also
accept an optional trailing `block_length=` keyword; the default is `1`, which
reproduces the original delete-one-window estimator exactly. The named block
entry points require a `block_length` argument to make the choice explicit.
The generic `multitaper_canonical_block_jackknife` accepts both real and
complex 2-D signals and uses `noverlap=0` when computing their per-window
multitaper spectra; all existing DPSS weighting choices (`equal`, `eigen`,
`adaptive`) remain available. The **whole window**, not an individual taper or
frequency bin, is resampled. You can additionally pass `bands(2,nband)` and
`average=.true.` to integrate/average complex cross-spectral matrices *within
each window* before fitting canonical modes.

## Estimator

For `W` ordered windows let `L=block_length` and `G=W/L`. The data must contain
**exactly `G` equal-size blocks**: `L >= 1`, `mod(W,L) == 0`, and `G >= 4`.
Each block is the consecutive index range `((g-1)*L+1):(g*L)` in the input
`spectra%times` order. Unused trailing samples in signal-level input do not
create a partial spectral window.

For complex Hermitian spectral matrices `S_w(f)`, the full estimate is

\[
\hat\theta(f) = h\left(\bar S(f)\right), \qquad
\bar S(f)=\frac{1}{W}\sum_{w=1}^{W}S_w(f),
\]

where `h` returns sorted conditional canonical squared coherence roots after
fixed trace-relative ridge and/or diagonal-target shrinkage. For deletion
group `g`, with the index set `B_g`, compute

\[
\hat\theta_{(-g)}(f)=h\left(
\frac{1}{W-L}\sum_{w\notin B_g}S_w(f)\right).
\]

With `G` refitted values and their mean `\bar\theta_{(-\cdot)}`, the
*equal-group* delete-block standard error and bias estimates are

\[
\widehat{\operatorname{SE}}(\hat\theta)=
\sqrt{\frac{G-1}{G}\sum_{g=1}^{G}
(\hat\theta_{(-g)}-\bar\theta_{(-\cdot)})^2},
\qquad
\widehat{\operatorname{bias}}=(G-1)
(\bar\theta_{(-\cdot)}-\hat\theta).
\]

`G`, **not** `W`, appears in these formulas. The same calculation applies to
the determinant-based aggregate group coherence
`1 - product(1 - squared_modes)`. For individual roots, a separate standard
error is also computed on the Fisher scale
`atanh(sqrt(squared_coherence))`; the approximate pointwise 95% interval is
the inverse transformation of the central estimate plus/minus `1.959964`
Fisher-scale standard errors. It is a **normal-approximation interval** and
should not be described as having verified 95% coverage.

## Result fields and validation

The return type remains `canonical_window_jackknife_result` for compatibility,
with two additional metadata fields:

| Field | Meaning |
| --- | --- |
| `n_windows` | Total number `W` of whole windows used for the point estimate |
| `n_blocks` | Number `G` of deletion groups, not the nominal number of independent observations |
| `block_length` | Number `L` of consecutive windows deleted per refit |
| `squared_modes`, `coherence` | Full-window-mean roots and aggregate group coherence |
| `stderr`, `stderr_fisher`, `coherence_stderr` | Delete-group standard errors |
| `bias` | Equal-group jackknife root bias estimate |
| `lower95`, `upper95` | Heuristic pointwise Fisher-scale interval endpoints |
| `valid`, `coherence_valid` | Mode and per-frequency/band validity masks |
| `status` | `0` success, `1` invalid input, `2` singular/unusable spectral cell |

The implementation refuses incompatible block lengths, fewer than four equal
blocks, malformed time grids, invalid groups or spectral matrices, and
invalid bands. For tied/near-tied modes or strengths near zero or one,
ordinary rank-specific inference is nonregular: the corresponding `valid`
entry is false and numerical uncertainty fields are NaNs. The aggregate
group-coherence SE can still be valid when individual roots are tied.

## Statistical limitations

- The windows must be **in temporal order**, with approximately stationary
  dynamics over the analysis interval. Monotonic timestamps are checked; equal
  duration, spacing, stationarity and mixing are **not** verified by the code.
  `noverlap=0` prevents shared sample points but is not proof of independence.
- Contiguous deletion blocks can account for some **short-memory** dependence,
  but adjacent blocks are not guaranteed to be independent, and the theoretical
  quality depends on block length, record length, and decay of dependence.
  A useful large-sample asymptotic regime has both `L` and `G` increasing,
  with `L/W -> 0`; practical analysis should compare several defensible
  block lengths and report sensitivity, not automatically choose a size.
- The method refits **fixed** ridge and shrinkage coefficients. Data-driven
  tuning or adaptive taper estimation *within each spectral window* is not
  re-selected across omitted blocks. Additional uncertainty from tuning is
  therefore omitted.
- Pointwise intervals do not handle multiple comparisons across frequencies,
  and do not provide confidence regions for the individual coefficient vectors
  or tracked subspaces. Nonregular roots remain excluded.
- Groups must divide the available window count evenly; no unequal-block
  approximation or forced truncation is performed. At least four blocks are
  required by the API; substantially more are recommended for interpreting
  standard errors or normal approximations.

This is an exploratory implementation of standard delete-group resampling,
not a numerical port of a public SciPy function. Standard background on
resampling dependent observations includes H. R. Künsch (1989), *The
Jackknife and the Bootstrap for General Stationary Observations*, Annals of
Statistics **17**(3), 1217–1241; Efron and Tibshirani (1993), *An
Introduction to the Bootstrap*. Neither source was imported as code.

## Reproducible validation

```bash
python tools/run_signal_canonical_block_jackknife.py
python tools/check_signal_canonical_block_jackknife_public.py
```

The native tests independently recompute all delete-group means and standard
errors, check block size 1 against the unchanged old estimator, invalid block
sizes, singular/boundary modes, real and complex multitaper calls, adaptive
weights, and frequency bands. The independent NumPy/SciPy reference uses
complex Cholesky/SVD and separately recomputed delete-block means, with
interpolated and trapezoid-integrated frequency bands. The Fortran code does
not acquire any new external dependencies.
