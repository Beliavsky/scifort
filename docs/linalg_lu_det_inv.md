# General LU decomposition, determinants, and inverses

The independently implemented `scifort_linalg_lu` adds SciPy-style **rectangular**
partial-pivoted LU factorization, including singular matrices. Its routines
are separate from the existing square, nonsingular `lu_factor`/`lu_solve`
interfaces in `scifort_linalg_dense`.

`scifort_linalg_det_inv` adds the operations corresponding to
`scipy.linalg.det` and the general-matrix path of `scipy.linalg.inv`.
All interfaces are generic over `real(dp)` and `complex(dp)` and accept
caller-allocated output arrays and explicit integer status codes.

| Entry point | Arguments (in order) | Output and behavior |
| --- | --- | --- |
| `lu` | `a, p, l, u, status[, check_finite]` | `A=P L U`, `P(m,m)`, `L(m,k)`, `U(k,n)` |
| `lu_indices` | `a, perm, l, u, status[, check_finite]` | `A=L(perm,:) U`, `perm(m)` (one-based) |
| `lu_permuted` | `a, pl, u, status[, check_finite]` | `A=PL U`, `PL(m,k)` |
| `det` | `a, value, status[, check_finite]` | Real or complex determinant of square `A` |
| `inv` | `a, inverse, status[, check_finite]` | Inverse of square nonsingular `A` |

Here `m=size(a,1)`, `n=size(a,2)`, and `k=min(m,n)`.
The permutation index convention differs from Python's zero-based indices:
Fortran uses **one-based** `perm`, so it holds `A=matmul(l(perm,:),u)`.
For a permutation matrix, `p(i,perm(i))=1`; all its other entries are zero.

These routines can factorize singular, tall, wide, square, real, and complex
matrices. `lu` handles zero core dimensions. `det` of a 0-by-0 matrix returns
one; `inv` of a 0-by-0 matrix succeeds with an empty output. The `lu`
convention uses `P(m,m)` even for zero-column input; SciPy 1.17 instead
returns `P(0,0)` in that special case.

`det` returns zero for singular matrices with an exact zero LU pivot and
returns a nonzero status for invalid shapes or nonfinite/overflowed results.
`inv` returns `linalg_status_singular` if its existing dense LU backend finds a
zero pivot, and leaves the output zero on failure. It solves against the
identity matrix using the existing dense LU routines, **not** by computing
cofactors. Tiny nonzero pivots are not a numerical rank test: for severely
ill-conditioned inputs, compare the inverse with a condition estimate.

The `check_finite` keyword defaults to `.true.`. Disabling up-front checks
is unsafe for nonfinite input; the existing dense LU backend still checks
finiteness when `inv` is called. `lu` uses the LAPACK-style complex pivot
selection `abs(real(z))+abs(aimag(z))` for closer agreement with SciPy.

Unlike Python, Fortran requires separate `lu` and `lu_indices` overloads to
change the rank of the permutation output. `lu_permuted` corresponds to
SciPy's `lu(..., permute_l=True)`. Batched arrays and SciPy `inv` options
such as `assume_a`/`lower` are not yet implemented. Only `real(dp)` and
`complex(dp)` inputs are supported; users can cast integer matrices first.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg_lu, only : lu_indices
use scifort_linalg_det_inv, only : det, inv
real(dp) :: a(2,2), lower_factor(2,2), upper_factor(2,2)
real(dp) :: ainv(2,2), determinant
integer :: perm(2), status

a(1,:) = [1.0_dp, 2.0_dp]
a(2,:) = [3.0_dp, 4.0_dp]
call lu_indices(a, perm, lower_factor, upper_factor, status)
call det(a, determinant, status)  ! determinant = -2
call inv(a, ainv, status)
```

See the full [executable example](../example/linalg_lu_det_inv_example.f90)
and the [validation report](linalg_lu_det_inv_validation_2026-11-09.md).

## Reproducing the comparisons

From the SciFort root (GNU Fortran, SciPy and NumPy required):

```sh
gfortran -std=f2018 -O3 -o /tmp/scifort_lu_driver \
  src/kinds.f90 src/linalg.f90 src/linalg_dense.f90 \
  src/linalg_lu.f90 src/linalg_det_inv.f90 \
  reference/linalg_lu_det_inv_driver.f90
python reference/check_linalg_lu_det_inv.py /tmp/scifort_lu_driver
```

Compile the native test similarly, replacing the reference driver with
`test/test_linalg_lu_det_inv.f90`. Existing `fpm test` should discover the new
native test when FPM and its project dependencies are available.
