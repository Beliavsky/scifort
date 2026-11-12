# Group-conditional spectral coherence

The `scifort_signal_group_coherence` module implements **determinant group
coherence**, an original **SciFort extension** (not a `scipy.signal` API).
It measures how strongly two *groups* of channels are linearly associated at
a spectral frequency, optionally after conditioning on a third, disjoint group.
It generalizes ordinary pairwise, partial and multiple coherence.

## Definition

Let `S(f)` be a Hermitian, positive-definite multichannel cross-spectral
matrix. Choose two nonempty, disjoint sets of channel indices `A`, `B`, plus
an optional disjoint conditioning set `C`. For `S_D` the principal block of
`S` indexed by `D`, the group statistic is

\[
G^2_{A,B\mid C}(f)=1-\Lambda(f), \qquad
\Lambda=\frac{\det(S_{A\cup B\cup C})\det(S_C)}
{\det(S_{A\cup C})\det(S_{B\cup C})}.
\]

For an empty `C`, define `det(S_C)=1`. The statistic lies in `[0,1]` for
positive-definite matrices. It is symmetric in `A` and `B`, invariant to the
order of channels within each group and invariant to invertible linear
mixing within each group. Numerically the implementation uses four Cholesky
**log determinants**, never forms a determinant by multiplying pivots and
never inverts a matrix to calculate a determinant.

The output component `%log_information` is `-log(Lambda)`, and
`%coherence = 1-exp(-%log_information)`. Computation of the latter uses a
series expansion near zero to avoid cancellation. The name
`log_information` describes an algebraic log determinant ratio: this is
not a p-value, and its interpretation as Gaussian mutual information
requires additional distributional assumptions (and depends on whether
real or proper complex variables are modeled).

An equivalent description uses the canonical coherences `rho_k` between
`A` and `B`, after linearly residualizing both groups against `C`:

\[
\Lambda=\prod_{k=1}^{\min(|A|,|B|)}(1-|\rho_k|^2).
\]

The complementary routines in
[canonical spectral coherence modes](signal_canonical_coherence.md) explicitly
return these individual squared canonical roots.

Unlike **maximum canonical coherence**, the determinant statistic aggregates
all canonical coherence modes. It is *not* the arithmetic average of
pairwise coherence values. No direction of causality is implied.

Special cases:

- `A=[i]`, `B=[j]`, and `C` contains all remaining channels: the result
  equals existing **squared partial coherence** between `i` and `j`.
- `A=[i]`, `B` contains all remaining channels and `C` is empty: the result
  equals the existing **squared multiple coherence** of channel `i`.
- Both `A` and `B` contain one channel, with empty `C`: ordinary magnitude-
  squared coherence.

## Fortran interfaces

All new functions are exported by the top-level `scifort_signal` module.

| Input | Frequency-resolved | Integrated band |
| --- | --- | --- |
| Multichannel CSD | `group_coherence_matrix` | `band_group_coherence_matrix` |
| Time-frequency CSD | `group_coherence_tf_matrix` | `band_group_coherence_tf_matrix` |
| Real or complex signal | `multitaper_group_coherence` | `multitaper_band_group_coherence` |
| Real or complex signal, moving window | `multitaper_tf_group_coherence` | `multitaper_tf_band_group_coherence` |

Matrix results are `group_coherence_result` or
`group_coherence_tf_result`. They contain `%frequencies`, `%coherence`,
`%log_information`, `%shrinkage_used`, `%valid`, and `%status`. Time-frequency
results also contain `%times`; their cell arrays have shape `(frequency,time)`.
The matrix input types are the established `multitaper_matrix_result` and
`multitaper_spectral_matrix_spectrogram_result`. Index vectors are **one-based**
Fortran channel numbers, and can be reordered, but cannot repeat, overlap, be
out of range or have an empty `A` or `B`.

Example, using the standard three-group case:

```fortran
use scifort_signal, only : dp, group_coherence_result, multitaper_group_coherence
implicit none
real(dp) :: channels(512,5)
type(group_coherence_result) :: result

! Fill channels(:,1:5) with finite observations.
result = multitaper_group_coherence(channels, 3.5_dp, [1,2], [3,4], &
    conditioned_on=[5], kmax=6, weighting='eigen', ridge=0.02_dp)
if (result%status /= 0) error stop 'invalid or singular group spectral estimate'
! result%coherence(f) is group coherence of channels (1,2) and (3,4),
! conditional on channel 5, with f in result%frequencies.
```

For band integration, supply `bands(2,nbands)` in Hz. For example:

```fortran
result = multitaper_band_group_coherence(channels, 3.5_dp, &
    reshape([0.0_dp,0.1_dp,0.1_dp,0.3_dp],[2,2]), &
    [1,2], [3,4], conditioned_on=[5], kmax=6, shrinkage=0.2_dp)
```

The band routines integrate the **cross-spectral matrix first** with the
existing exact piecewise-linear algorithm, then compute group coherence.
The frequency coordinates in the result are band midpoints. As with the
existing band coherence routines, frequency-bin shrinkage and effective
sample counts are **not transferred** to integrated bands; these wrappers
accept explicit `ridge` and `shrinkage` only.

The direct multitaper interfaces accept the same `kmax`, `fs`, `nfft`,
`detrend`, `weighting`, `axis`, and `detrend_mode` options as the existing
multitaper CSD functions. Time-frequency variants additionally accept
`nperseg` and `noverlap`.

## Regularization and limitations

Frequency- and time-frequency-level functions support the existing spectral
regularization conventions: diagonal-target `shrinkage` (zeroing a fraction
of off-diagonal cross spectra) plus a `ridge` proportional to the **mean
original diagonal power across all channels**. `auto_shrinkage` uses the
preexisting Gaussian–Wishart plug-in estimator for fixed taper weights,
and `cv_shrinkage` uses the precomputed leave-one-taper-out coefficients,
including adaptive weights. The matrix functions accept
`effective_tapers` as an explicit override for the automatic estimator;
otherwise they use the spectral result's `%effective_tapers` metadata.
The same option conflicts as existing partial coherence apply.

As an implementation choice, **the complete input spectral matrix must be
positive definite after regularization**, even when only some channels
are selected. This lets SciFort reuse the existing validated precision
solver and ensures consistent validity and shrinkage with partial/multiple
coherence, but it means singular *unused* channels may also invalidate a cell.
If one or more cells are singular, `%status` reports the existing singular
status and `%valid` identifies individually valid cells; invalid numeric
outputs are NaN. Structurally malformed inputs return the invalid-input
status with empty result arrays.

Group coherence is descriptive. Small samples, strongly overlapping windows,
unequal taper weights and spectral smoothing affect its sampling distribution.
No confidence intervals, p-values, causality claims or statistical
independence tests are provided by these routines.

## Numerical verification

`python tools/run_signal_group_coherence.py` compiles and tests the targeted
Fortran modules in both checked and `-O3 -march=native` modes. It runs native
regressions for identity reductions, ordering symmetry, regularization,
invalid inputs, real/complex DPSS signals, adaptive cross-validation,
time-frequency and integrated-band processing. Its independent NumPy test
forms conditional block covariances using **Schur complements**, whitens
both groups with Cholesky and compares singular values (not determinant
ratios). This covers 196 randomized Hermitian-matrix configurations and
1,418 frequency/time cells, including strong coupling, unequal channel
scales, unequal group sizes and complex cross-spectrum phases.

`python tools/check_signal_group_coherence_public.py` recompiles the complete
`scifort_signal` dependency closure and exercises the top-level API together
with earlier public coherence regressions. No additional Fortran runtime
libraries are required by the new source module.
