# Conditional canonical spectral coherence modes

**SciFort extension:** `scifort_signal_canonical_coherence` provides the
individual *squared* canonical coherences between two multichannel spectral
blocks, optionally conditioning on a third, disjoint block. SciPy does not
provide a matching top-level `scipy.signal` routine.

Unlike `group_coherence_matrix`, which returns a single determinant summary,
these APIs expose the sorted spectral mode strengths, from largest to smallest.
By default only strengths are returned; optional `return_vectors=.true.`
provides variance-normalized complex canonical coefficient vectors and
flags for non-identifiable or numerically sensitive individual directions.
Optional `return_covariances=.true.` exposes the conditioned, regularized
A/B autospectral matrices used by this normalization, for covariance-aware
subspace tracking without separately recomputing Schur complements.

## Mathematical definition

Let `S(f)` be the Hermitian cross-spectral matrix and choose nonempty disjoint
one-based channel groups `A` and `B`, plus an optional disjoint group `C`.
First form the conditional blocks through Schur complements:

\[
S_{AA\mid C}=S_{AA}-S_{AC}S_{CC}^{-1}S_{CA}, \qquad
S_{BB\mid C}=S_{BB}-S_{BC}S_{CC}^{-1}S_{CB},
\]

\[
S_{AB\mid C}=S_{AB}-S_{AC}S_{CC}^{-1}S_{CB}.
\]

For `C` empty, use the unconditioned blocks. Factor the positive-definite
auto-spectral blocks as `L_A L_A^H` and `L_B L_B^H`. The whitened cross-spectrum
is

\[
W=L_A^{-1} S_{AB\mid C} L_B^{-H}.
\]

The reported `squared_modes` are the **descending eigenvalues** of either
`W W^H` or `W^H W` (whichever matrix is smaller), equivalently the squares
of the singular values of `W`. There are exactly `min(size(A),size(B))`
modes per frequency/time cell; each is in `[0,1]` for valid positive-definite
spectral matrices. The routines use complex Cholesky solves and an original
small-matrix cyclic Hermitian Jacobi eigensolver, with no BLAS/LAPACK runtime
dependency. A nonconvergent or numerically invalid cell is marked singular.

## Canonical variate coefficients (optional)

For a singular value `sigma_k = sqrt(squared_modes(f,k))`, let `u_k`, `v_k`
be paired unit left/right singular vectors of the whitened matrix `W`, with
`u_k^H W v_k = sigma_k >= 0`. The optional coefficient vectors are

\[
a_k=L_A^{-H}u_k, \qquad b_k=L_B^{-H}v_k.
\]

Writing the columns as `Acoef` and `Bcoef`, they satisfy

\[
Acoef^H S_{AA\mid C} Acoef=I, \qquad
Bcoef^H S_{BB\mid C} Bcoef=I,
\]
\[
Acoef^H S_{AB\mid C} Bcoef
=\operatorname{diag}(\sigma_1,\ldots,\sigma_m).
\]

These are *conditional* canonical variates: when `C` is nonempty, the
coefficients multiply the channels of A and B after their frequency-specific
linear projections on C have been removed. They are not regression weights
on the unadjusted original channels, and coefficients for C itself are not
returned. The same regularization is used for both roots and vectors.

The columns of each group are orthonormal with respect to that group's
conditional spectral covariance, **not** necessarily with respect to the
ordinary Euclidean inner product. A common complex phase rotation of a paired
vector leaves its association unchanged; the routine fixes the phase by
making the largest-magnitude entry of each returned A coefficient vector
real and positive. This is a convention within each spectral cell, not an
across-frequency phase continuity guarantee.

For an exactly zero singular value, the opposite-side direction is not
identified by the data. The code chooses a deterministic orthogonal
complement to complete the canonical basis. For repeated singular values,
any unitary rotation of the corresponding mode subspace is valid. In both
cases, individual directions should not be interpreted uniquely.

The existing determinant aggregate satisfies

\[
G^2_{A,B\mid C}=1-\prod_k(1-\rho_k^2),\qquad
-\log\Lambda=-\sum_k\log(1-\rho_k^2).
\]

The original group-coherence estimator remains the source of aggregate
`%coherence` and `%log_information`: these remain stable in cases where
forming a product from rounded eigenvalues would lose precision. For normal
well-conditioned cells the two representations agree to rounding accuracy.

## Tracking across frequencies (optional SciFort extension)

The default canonical results are independently sorted by mode strength at
*each frequency*. To follow the same spectral direction when mode strengths
cross and to align its complex phase, call `track_canonical_modes` on a result
computed with `return_vectors=.true.`. This is a separate postprocessor and
does not modify the original sorted results. See the
[canonical spectral tracking reference](signal_canonical_tracking.md) for the
similarity metric, Hungarian assignment, reset conditions, and outputs.

## API

All eight functions below are exported by `scifort_signal` and directly by
`scifort_signal_canonical_coherence`.

| Source | Frequency-resolved | Integrated frequency band |
|---|---|---|
| Multichannel cross-spectral matrix | `canonical_coherence_matrix` | `band_canonical_coherence_matrix` |
| Time-frequency cross-spectral matrix | `canonical_coherence_tf_matrix` | `band_canonical_coherence_tf_matrix` |
| Raw real or complex multichannel signals | `multitaper_canonical_coherence` | `multitaper_band_canonical_coherence` |
| Raw real or complex, sliding windows | `multitaper_tf_canonical_coherence` | `multitaper_tf_band_canonical_coherence` |

The result types are `canonical_coherence_result` and
`canonical_coherence_tf_result`:

- `%squared_modes(f,k)` for static or band input; `%squared_modes(f,t,k)` for
  time-frequency input. Modes `k` are sorted **descending**. This is a
  magnitude-squared quantity, not a signed or complex-valued coefficient.
- When `return_vectors=.true.`, `%coefficients_a(f,i,k)` and
  `%coefficients_b(f,j,k)` contain conditional canonical coefficients in
  the **input group order**. In the time-frequency result the shapes are
  `%coefficients_a(f,t,i,k)` and `%coefficients_b(f,t,j,k)`. The last index
  `k` matches `%squared_modes` in descending order.
- When `return_covariances=.true.`, `%covariance_a(f,i,j)` and
  `%covariance_b(f,i,j)` return the exact conditional regularized autospectral
  matrices that define the canonical normalization. For TF results the shapes
  are `(f,t,i,j)`. Invalid cells contain complex NaNs. This option may be
  requested independently of `return_vectors` and is disabled by default.
- Also when requested, `%directions_stable(f,k)` or
  `%directions_stable(f,t,k)` flags an isolated, nonzero root. This is a
  conservative numerical heuristic, **not** a confidence interval, a
  statistical stability estimate, or a guarantee of continuity across bins.
  It is false if `sqrt(squared_modes)` is no larger than
  `64*sqrt(epsilon(1.0_dp))`, or if the difference from an adjacent
  squared coherence is no larger than that tolerance.
  This flag is false for exactly and nearly repeated roots.
- `%frequencies` and, for time-frequency results, `%times`.
- `%coherence`, the original determinant-based scalar group coherence per cell.
- `%log_information`, the corresponding negative log determinant ratio.
- `%shrinkage_used`, regularization coefficients inherited from group coherence.
- `%valid`, a per-frequency or frequency/time Boolean mask.
- `%status`, the existing group/partial coherence status code (zero on success,
  invalid-input or singular status otherwise).

For example:

```fortran
use scifort_signal, only : dp, canonical_coherence_result, &
    multitaper_canonical_coherence
implicit none
real(dp) :: observations(512,5)
type(canonical_coherence_result) :: result

! Fill observations(:,1:5) with finite multichannel signal data.
result = multitaper_canonical_coherence(observations, 3.5_dp, &
    [1,2], [3,4], conditioned_on=[5], kmax=6, &
    weighting='eigen', ridge=0.02_dp, return_vectors=.true.)
! result%coefficients_a(f,:,k) weights the residualized A group.
! result%coefficients_b(f,:,k) weights the residualized B group.
! result%directions_stable(f,k) marks isolated, nonzero roots.
! result%squared_modes(f,1:2) are the two sorted squared coherence modes.
! result%coherence(f) aggregates both modes through Wilks' lambda.
```

Groups must be disjoint, nonempty for `A`/`B`, in range, and must not contain
duplicates. An empty `C` (the default) gives unconditional canonical
coherence. Swapping `A` and `B`, or permuting channels within either block,
preserves the mode strengths up to numerical rounding. The coefficient
vectors transform consistently with permutations, up to global phase and
arbitrary unitary rotations within repeated-root subspaces.

With the default `return_vectors=.false.`, coefficient and stability arrays
remain **unallocated**. On an invalid call they also remain unallocated,
including when vectors were requested. For valid input, cells marked
`%valid=.false.` receive complex NaNs in their coefficient arrays, and
`directions_stable=.false.`. Requesting vectors can increase memory and
CPU use, especially for moving-window analysis; existing default performance
and result shapes remain unchanged.

To use the covariance-metric subspace tracker, request **both**
`return_vectors=.true.` and `return_covariances=.true.`, then call
`track_canonical_subspaces(result,covariance_metric=.true.)`. See
[the mathematical and API reference](signal_canonical_covmetric_tracking.md).

## Regularization and frequency integration

The `ridge`, `shrinkage`, `auto_shrinkage`, `effective_tapers`, and
`cv_shrinkage` options follow the existing group-coherence semantics. Fixed
shrinkage acts on the entire cross-spectral matrix, not only the selected
channels. `auto_shrinkage` uses the existing fixed-weight Gaussian–Wishart
plug-in; `cv_shrinkage` uses previously computed leave-one-taper-out
coefficients, including adaptive taper weighting. Those two automatic options
are mutually exclusive. Complete positive definiteness of the regularized
**full spectral matrix** is required, even if not all channels are selected.
Invalid per-cell modes and aggregates are NaN, while the `%valid` mask shows
which cells remain usable.

Band functions first apply `integrate_spectral_bands` or
`integrate_tf_spectral_bands` to the **complex spectral matrix**, then compute
its canonical modes. Pass `bands(2,nbands)` in Hz; band midpoints are returned
in `%frequencies`. The optional `average` returns mean band spectral density
instead of total power. These band APIs accept `ridge` and fixed `shrinkage`
only; they do not propagate frequency-bin automatic or CV shrinkage estimates
or claim effective independent taper counts across frequency bins.

## Numerical and inferential limitations

The procedure gives descriptive multichannel spectral association, not a
causal relationship or a statistical significance test. No p-values or
confidence intervals are supplied. Repeated or near-tied roots can make
individual mode directions unstable even where the strengths are well-defined. Very
nearly singular covariance matrices may need explicit shrinkage or ridge
loading. The Jacobi method is intended for modest group sizes, not as a
replacement for a high-performance general-purpose eigensolver.

## Reproduction

From the repository root (GNU Fortran 14 or a compatible compiler, NumPy):

```bash
python tools/run_signal_canonical_vectors.py
python tools/check_signal_canonical_vectors_public.py
python tools/run_signal_canonical_coherence.py
python tools/check_signal_canonical_coherence_public.py
```

The vector-specific scripts compile checked and optimized native test suites,
compare coefficients against NumPy Cholesky/SVD (150 independent scenarios,
798 spectral cells per build), and rebuild the complete public interface
(96 source modules, 10 regression programs). They check both left and right
variance normalization, phase conventions, rectangular groups, zero and
repeated modes, regularization, and legacy API behavior.

The existing canonical-coherence runner compiles targeted Fortran modules under runtime checks
and `-O3 -march=native`, runs native regressions and compares individual roots
to independent NumPy SVD of whitened, Schur-conditioned spectra. The reference
includes real and complex Hermitian matrices, uneven group sizes, phase,
unequal channel scales, fixed shrinkage and ridge. Its paired public-interface checker builds the full
`scifort_signal` dependency closure and runs the original nine tests.
