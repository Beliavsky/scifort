# SciFort Procrustes analysis (2026-11-02)

`scifort_spatial` exports real-valued implementations of
`scipy.spatial.procrustes` and `scipy.linalg.orthogonal_procrustes`.
Both use only native Fortran, and depend on no BLAS, LAPACK, or third-party
implementation code. Input matrices use `(observations, dimensions)`, with
one point per row. The methods accept `real(dp)` (binary64).

## Similarity Procrustes: `procrustes`

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : procrustes
implicit none
real(dp) :: source(4,2), target(4,2), disparity
real(dp), allocatable :: standardized(:,:), aligned(:,:)
integer :: status

source(1,:) = [0.0_dp,0.0_dp]
source(2,:) = [1.0_dp,0.0_dp]
source(3,:) = [1.0_dp,1.0_dp]
source(4,:) = [0.0_dp,1.0_dp]
target(:,1) = 4.0_dp*source(:,2)+7.0_dp
target(:,2) = -4.0_dp*source(:,1)+2.0_dp

call procrustes(source,target,standardized,aligned,disparity,status)
! status == 0; disparity approximately zero
```

The subroutine computes:

1. Independently center each matrix, with a subtraction relative to the first
   observation to reduce cancellation for large translations.
2. Normalize both centered configurations to unit Frobenius norm.
3. Find the best rotation/reflection of the second matrix onto the first.
4. Multiply that rotated matrix by the sum of singular values of the
   normalized cross-product matrix, the optimal nonnegative scale.

`standardized` and `aligned` have the same `(n_points,n_dimensions)` shape as
`source` and `target`. The first matrix has mean zero in each column and
squared Frobenius norm one. The second is centered but not necessarily unit
norm. The scalar `disparity` is exactly the sum of squared elementwise
residuals between these two output matrices, up to floating-point rounding.

The disparity is symmetric if `source` and `target` are exchanged, but the
actual returned matrices depend on input order. The method allows reflections
and arbitrary scaling, as does SciPy. Duplicate points are allowed; each row
contributes independently to the objective.

Both input matrices must have identical shape, at least two observations and
one feature, all finite entries, and at least two distinct observations in
each configuration. For failures, both returned matrices have shape `(0,0)`,
`disparity` is NaN, and the optional `status` reports the error. Status values
are `spatial_distance_success=0`, `spatial_distance_invalid_input=1`, and
`spatial_procrustes_no_convergence=2`.

## Orthogonal Procrustes: `orthogonal_procrustes`

```fortran
use scifort_kinds, only : dp
use scifort_spatial, only : orthogonal_procrustes
implicit none
real(dp) :: a(10,3), b(10,3), scale
real(dp), allocatable :: rotation(:,:)
integer :: status

! a and b are populated by the caller.
call orthogonal_procrustes(a,b,rotation,scale,status)
! When status==0, rotation is orthogonal and minimizes
!     sum((matmul(a,rotation)-b)**2)
```

This routine is exposed through `scifort_spatial` for convenience even though
its SciPy counterpart belongs to `scipy.linalg`. It **does not** translate,
center, or scale the input arrays. The returned `rotation` is square with
shape `(dimensions,dimensions)`, and `scale` is the nuclear norm (sum of
singular values) of `transpose(a) @ b`. The routine permits rank-deficient
input, even a zero matrix; when either input is all zeros it returns the
identity mapping and zero nuclear norm.

Invalid matrices yield an empty rotation and NaN scale, with the same status
codes as above. When the mathematical nuclear norm overflows binary64, its
returned value may be +infinity even if the orthogonal rotation remains valid.

## Numerical approach and limitations

Both procedures use a self-contained **one-sided Jacobi** column
orthogonalization on the cross-product matrix. It produces an orthogonal
polar factor from left and right singular directions. Numerically null
singular directions are completed to an orthonormal basis in deterministic
order. Unlike an eigenvalue decomposition of `transpose(C) @ C`, this does
not square the condition number of the input cross-product matrix. Orthogonal
completion is deliberately not guaranteed to match SciPy's choice for
singular or nearly singular inputs; the best-fit objective and disparity
are the compatibility contract in those cases.

The implementation uses a bounded number of Jacobi sweeps. If the method
cannot converge, it returns status 2 without an approximate result. Very
ill-conditioned high-dimensional problems can be more efficiently solved by
an optimized LAPACK SVD; SciFort intentionally keeps this portable fallback.

One important intentional difference: SciFort normalizes each point set
before centering, and anchors centering at a single observation. This avoids
norm overflow/underflow for very large/small data and cancellation from
common offsets. SciPy's high-level routine may overflow when given raw
magnitudes around `1e180`, although the equivalent normalized inputs work.

## Independent validation

```sh
python tools/spatial_reference/run_native_spatial.py --mode both
python tools/spatial_reference/compare_procrustes_scipy.py --mode both
python tools/spatial_reference/compare_procrustes_scipy.py --mode both --seed 23791 --count 400
```

The numerical oracle is SciPy 1.17.0. The checks include randomly rotated
and reflected clouds, both full-rank and rank-deficient data, repeated rows,
under- and overdetermined configurations, large translations, independent
rescalings, and near-rank-one matrices. Nonunique rotations are checked by
orthogonality and optimal objective rather than entrywise matrix agreement.

**Provenance:** Independent MIT-licensed implementation of the mathematical
orthogonal Procrustes problem (Schonemann, *Psychometrika* 31, 1966,
DOI:10.1007/BF02289451), using SciPy's public API only for compatibility
specifications and numerical comparison; no SciPy code was imported.

See also [Hausdorff and spherical geometry](spatial_geometry.md),
[spatial distances](spatial_distance.md), and [KD-tree queries](spatial_kdtree.md).
