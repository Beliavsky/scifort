# Range, null-space and principal-angle routines

`scifort_linalg_subspaces` implements three `scipy.linalg`-style APIs using
SciFort's standalone economy Jacobi SVD (`scifort_linalg_svd`). All calculations
use binary64 real or complex arithmetic; no BLAS/LAPACK dependency is added.

## Public interfaces

| Generic | Argument order | Output |
| --- | --- | --- |
| `orth` | `a, q, status, rcond, max_sweeps` | Allocatable `q(m,r)` spanning the column space |
| `null_space` | `a, z, status, rcond, max_sweeps` | Allocatable `z(n,n-r)` spanning the numerical kernel |
| `subspace_angles` | `a, b, angles, status, rcond, max_sweeps` | Allocatable `angles(min(ra,rb))`, descending, radians |

Here `a` is `m`-by-`n`, `b` has `m` rows, `r` is the numerical rank of
`a`, and `ra`/`rb` are the ranks of `a`/`b`. The output arrays must be
**allocatable** and are allocated by the routine. All input arrays are
unchanged. `orth`/`null_space` provide real and complex overloads;
`subspace_angles` additionally accepts mixed real/complex input matrices.

All three routines use singular values **strictly greater than**
`rcond * max(s)` to define numerical rank. By default,
`rcond=max(m,n)*epsilon(1.0_dp)` for `orth` and `null_space`. The
`subspace_angles` default rank cutoff for each matrix uses its own
`max(rows,columns)*epsilon`. The optional `rcond` on `subspace_angles` is a
**SciFort extension** (SciPy's public routine does not offer a rank cutoff);
it applies the same explicitly supplied value to both inputs.

`max_sweeps` controls the number of Jacobi SVD sweeps (default 120).

## Algorithms

`orth` returns the retained left singular vectors of `A`. `null_space`
returns the discarded right singular vectors, then completes the **economy**
right-singular-vector basis to a full orthonormal basis when `n > min(m,n)`.
The latter uses canonical-vector pivoting and twice-applied modified
Gram-Schmidt; the diagonal of the residual projector is incrementally updated.
This completion is necessary for wide matrices even when `A` has full row rank.
The basis vector phases and signs need not match SciPy, but the projectors
`Q Q^H` and `Z Z^H` represent the same subspaces.

`subspace_angles` first computes bases for both column spaces. Its cosine
information comes from the singular values of `QA^H QB`. For small angles
it uses singular values of the projection residual to obtain sines, avoiding
catastrophic cancellation in `acos(sigma)` when `sigma` is close to one.
Angles are returned from largest to smallest, between zero and `pi/2`.
If either effective rank is zero, the result is an empty array.

Zero-row and zero-column inputs are supported: `orth(0-by-n)` yields an
empty basis, `null_space(0-by-n)` yields the identity of order `n`, and a
subspace with zero dimension yields an empty list of angles.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_subspaces, only : orth, null_space, subspace_angles
real(dp) :: a(3,4), b(3,1)
real(dp), allocatable :: q(:,:), z(:,:), theta(:)
integer :: status

a = 0.0_dp
a(1,1) = 1.0_dp
b = 0.0_dp
b(2,1) = 1.0_dp
call orth(a,q,status)
call null_space(a,z,status)
call subspace_angles(a,b,theta,status)
! q is 3-by-1, z is 4-by-3, theta(1) = pi/2
```

See `example/linalg_subspaces_example.f90` for a complete checked example.

## Status codes and limitations

- `linalg_status_success=0` on success.
- `linalg_status_invalid_shape=1` for inconsistent row counts, negative or
  nonfinite `rcond`, or `max_sweeps<1`.
- `linalg_status_no_convergence=4` for nonfinite matrix elements, Jacobi
  failure, or failed basis completion.
- On invalid input or failure, basis outputs are allocated with zero columns;
  no partially populated basis is published.

The SVD core is an unblocked one-sided Jacobi method. Its numerical rank is
sensitive to near-roundoff singular values and may differ from SciPy's
LAPACK-backed SVD near the selected cutoff. In particular, the existing
Jacobi SVD thresholds extremely small singular directions during completion.
The null-space completion is quadratic in storage and cubic in work; it is
not designed for large sparse operators. As with the underlying SVD, very
ill-conditioned matrices can yield sensitivity to floating-point arithmetic,
and the `max_sweeps` limit may need increasing.

The implementation is not a drop-in Python API: it uses Fortran allocatable
`intent(out)` arrays, integer statuses, and no batch dimensions.

See [the reproducible validation record](linalg_subspaces_validation_2026-11-07.md).
