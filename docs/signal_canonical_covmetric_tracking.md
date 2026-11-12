# Covariance-metric tracking of canonical spectral subspaces

**SciFort extension (2026-10-30).** The existing
[`track_canonical_subspaces`](signal_canonical_subspace_tracking.md) can now
use a positive-definite, frequency-dependent *conditional spectral covariance*
inner product, instead of comparing canonical coefficient arrays in arbitrary
Euclidean channel coordinates. This is not a public `scipy.signal` routine.

This is **opt-in**. Legacy Euclidean tracking remains the default. The
canonical estimator can return the exact conditional spectral covariance
blocks used to normalize its coefficient vectors, without a separate
recalculation of conditioning, ridge loading, shrinkage, or band integration.

## Usage

```fortran
use scifort_signal, only : dp, canonical_coherence_result, &
    canonical_subspace_tracks_result, multitaper_canonical_coherence, &
    track_canonical_subspaces, canonical_subspace_success
implicit none
real(dp) :: x(1024,5)
type(canonical_coherence_result) :: canonical
type(canonical_subspace_tracks_result) :: tracks

! Populate x(sample,channel) before calling.
canonical = multitaper_canonical_coherence(x,3.5_dp,[1,2],[3,4], &
    conditioned_on=[5],kmax=6,weighting='eigen',ridge=0.02_dp, &
    return_vectors=.true.,return_covariances=.true.)
tracks = track_canonical_subspaces(canonical,gap_tolerance=0.001_dp, &
    min_similarity=0.6_dp,covariance_metric=.true.)
if (tracks%status == canonical_subspace_success) then
    ! tracks%covariance_metric_used is true.
    ! tracks%track_id(f,k) is a persistent BLOCK identifier.
    ! tracks%basis_a and %basis_b are aligned BLOCK BASES.
end if
```

The same two options work in real/complex, static/time-frequency, and
band-integrated canonical-coherence interfaces. In the time-frequency case,
matching proceeds separately along frequency for each time window.

`return_covariances=.true.` is optional on all canonical estimators, and may
also be requested without `return_vectors`. It allocates
`canonical%covariance_a(f,i,j)` and `canonical%covariance_b(f,i,j)` for static
or band results, and `canonical%covariance_a(f,t,i,j)` and
`canonical%covariance_b(f,t,i,j)` for time-frequency results. The matrices
use the input's group ordering. They are conditional on the specified third
group and contain the **same** diagonal shrinkage and trace-scaled ridge
loading as the corresponding canonical modes. When no conditioning group is
provided, they are the regularized principal autospectral blocks. Invalid
spectral cells have complex NaNs. Without the option, the arrays remain
unallocated, preserving legacy allocation behavior.

`covariance_metric=.true.` on the subspace tracker **requires both coefficient
vectors and covariance matrices** with the expected shapes and positive
definiteness at every valid cell. Missing arrays, incorrect shapes, or a
nonfinite, non-Hermitian or non-positive-definite covariance at a valid cell
return `canonical_subspace_invalid_input`. An invalid source cell can contain
NaN covariances because the tracker skips it and restarts after it. The output
`covariance_metric_used` indicates which comparison rule was selected.

## Mathematical definition

Let `A_p,B_p` and `A_q,B_q` be the paired coefficient matrices for two
candidate canonical mode blocks of the same dimension `d`, at consecutive
frequencies `p` and `q`. Let `C_A(p)` and `C_A(q)` be the Hermitian positive-
definite conditional autospectral covariances for group A, and likewise B.
Use a **common midpoint metric** on *both* compared blocks:

\[
  M_A=\tfrac12(C_A(p)+C_A(q)),\qquad
  M_B=\tfrac12(C_B(p)+C_B(q)).
\]

Factor `M_A=L_A L_A^H`, `M_B=L_B L_B^H`. Convert both endpoints to a
common inner-product coordinate system:

\[
  \widetilde A_r=L_A^H A_r,\qquad
  \widetilde B_r=L_B^H B_r,\qquad r\in\{p,q\}.
\]

The original tracker then uses **one Frobenius scale per whole group/block**
for the stacked paired coefficient matrix, rather than per-column scaling:

\[
  X_r=\frac1{\sqrt 2}
  \begin{bmatrix}
    \widetilde A_r / \sqrt{\|\widetilde A_r\|_F^2/d} \\
    \widetilde B_r / \sqrt{\|\widetilde B_r\|_F^2/d}
  \end{bmatrix}.
\]

If `Q_p` and `Q_q` are Euclidean-orthonormal bases for the columns of `X_p`
and `X_q`, the singular values of `Q_p^H Q_q` are the principal-angle
cosines. The minimum cosine is the block similarity. Eligible block pairs
are matched globally by the existing Hungarian assignment. For a matched
block, with `H=X_p^H X_q=U Sigma V^H`, the common right unitary rotation
`R=V U^H` aligns **both** new A/B coefficient blocks. This does not
change their own conditional spectral-variance normalization.

The common metric matters: whitening each endpoint with its *own* covariance
would conflate a change in the metric with a change in the underlying
subspace. The midpoint metric defines one inner product for the comparison.
It is an arithmetic, not an information-geometric, mean; it is a deliberate
numerical comparison rule and not a statistically optimal estimator.

### Invariance to channel units

For fixed, nonsingular linear changes of coordinates `T_A` and `T_B`,
within the respective channel groups,

\[
 C'_A=T_A C_A T_A^H,\qquad A'_r=T_A^{-H} A_r,
\]

with the analogous expressions for B. Every coefficient overlap Gram
matrix formed using the midpoint metric is unchanged under this simultaneous
transformation. Therefore the principal cosines, global matching scores,
block IDs, and common Procrustes right rotations are invariant up to
numerical roundoff. The aligned coefficient bases themselves transform
contravariantly. A **fixed** coordinate change across the tracked frequency
region is essential to this claim; frequency-varying channel transforms do
not generally preserve the same physical subspace trajectory.

This is an invariance property of the **tracking comparison given consistently
transformed covariance matrices and coefficient vectors**. Re-running spectral
estimation on rescaled raw channels can change the output whenever the
estimator uses coordinate-dependent regularization: in particular, SciFort's
trace-scaled ridge and diagonal-target shrinkage are not fully invariant to
arbitrary channel rescalings or mixing. The new option does not make those
upstream estimation steps invariant.

## Limitations and numerical safeguards

The per-group midpoint covariance must be strictly Hermitian and sufficiently
positive definite for a Cholesky factorization. Matrices with numerically
singular pivots or marked non-Hermitian parts are rejected, rather than
silently replaced with identity or an arbitrary ridge. All original
subspace rules still apply: equal block dimensions, nonzero mode strengths,
finite full-rank coefficients, minimum principal cosine, monotone frequency
coordinates, and restarts on FFT wraps or invalid cells.

**Rotated columns of a near-degenerate block are subspace bases, not necessarily
individual canonical eigenvectors.** Original ordered squared strengths
remain available separately. These scores and IDs describe continuity; they
are not calibrated confidence probabilities, inferential uncertainty bounds,
or evidence that two adjacent frequency estimates are independent.

The covariance option increases memory by approximately
`O(n_frequency * (n_A**2 + n_B**2))` complex numbers (times the number of
time windows for TF results). For large channel groups, leave it disabled
unless covariance-aware comparisons are needed.

## Reproduce the tests

```sh
python tools/run_signal_canonical_covmetric_tracking.py
python tools/check_signal_canonical_covmetric_tracking_public.py
```

The native tests verify missing and malformed covariance rejection,
conditioning under ridge/shrinkage, zero-cost defaults, time-window resets,
and invariance to severe complex channel rescaling. The independent
`verify_canonical_covmetric_tracking.py` compares the Cholesky whitening,
principal cosines, complex unitary Procrustes factors, globally assigned
block identities, and aligned coefficients to NumPy and SciPy on seeded
randomized complex scenarios. No external BLAS/LAPACK is required by the
Fortran implementation; Python numerical packages are test oracles only.
