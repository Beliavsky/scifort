# Complex sparse LSQR least-squares solver

`scifort_sparse_linalg_complex_lsqr` implements the real-valued SciFort `lsqr`
algorithm for **complex binary64** canonical CSR and CSC matrices, without a
BLAS/LAPACK dependency or a dense normal-equations matrix. Both overloads accept
rectangular matrices, including underdetermined systems and rank deficiency.

## Public interface

```fortran
use scifort_kinds, only : dp
use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
use scifort_sparse_linalg_complex_lsqr, only : lsqr_complex
implicit none
type(csr_complex_matrix) :: a
complex(dp), allocatable :: x(:)
complex(dp) :: rhs(3)
integer :: stop, steps
real(dp) :: residual

a = csr_complex_from_dense(reshape([ &
    cmplx(1.0_dp,1.0_dp,dp), cmplx(2.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(1.0_dp,-1.0_dp,dp), cmplx(3.0_dp,0.0_dp,dp) &
], [3,2]))
rhs = [cmplx(2.0_dp,0.0_dp,dp), cmplx(1.0_dp,0.0_dp,dp), cmplx(3.0_dp,1.0_dp,dp)]
call lsqr_complex(a, rhs, x, stop, itn=steps, r1norm=residual, &
    atol=1.0e-12_dp, btol=1.0e-12_dp, iter_lim=100)
```

The signature is `lsqr_complex(a,b,x,istop,itn,r1norm,r2norm,damp,atol,btol,conlim,iter_lim,x0,anorm,acond,arnorm,xnorm,var,calc_var)`.
Only `a`, `b`, `x` and `istop` are required.

- `x0`: optional complex starting vector and center of Tikhonov damping;
  minimizes `||Ax-b||**2 + damp**2 * ||x-x0||**2`.
- `damp`: nonnegative real scalar (default 0).
- `atol`, `btol`: nonnegative real tolerances (default `1e-6`).
- `conlim`: nonnegative condition limit (default `1e8`, 0 disables it).
- `iter_lim`: iteration cap (default `2*min(m,n)`).
- `r1norm`, `r2norm`: *actual* residual norm and augmented residual norm,
  not Krylov estimates; `arnorm` is the actual augmented normal-equation residual.
- `anorm`, `acond`: iterative estimates, not certified bounds; `var` is a
  real-valued diagonal estimate if `calc_var=.true.`.

`istop=0` denotes an already satisfied initial residual; `1` or `2` denotes
LSQR stopping criteria met; `3` a condition-limit stop; `7` exhausted
iterations; `-1` invalid input; `-2` arithmetic breakdown. `istop=2` can
also denote a stationary least-squares problem (including zero columns).
Negative statuses return an empty `x`.

This solver uses `A` and its *conjugate transpose* `A**H`, not `transpose(A)`.
It is a standard Golub–Kahan/Paige–Saunders LSQR implementation, independently
written for SciFort and adapted to Fortran's one-based complex sparse format.
It does not form `A**H*A`, which is generally less numerically stable.

## Validation

- `test/test_sparse_complex_lsqr.f90` verifies tall and wide systems,
  regularization about nonzero `x0`, the normal-equation residual, zero-sized
  matrices, invalid inputs, CSC/CSR agreement, the variance estimate, and a
  6000-by-120 sparse system assembled without dense storage.
- `tools/sparse_reference/compare_sparse_complex_lsqr_scipy.py` compares
  randomized complex systems against both SciPy's `lsqr` and NumPy's dense
  SVD-based least-squares solver, including explicit rank deficiency,
  regularization, zero rows/columns, and nonzero initial centers.

```sh
python tools/sparse_reference/compare_sparse_complex_lsqr_scipy.py --mode both --cases 200
python tools/sparse_reference/run_native_sparse.py --mode checked
```

The core algorithm does not provide preconditioning, arbitrary
`LinearOperator` callbacks, or sparse mixed precision. No new dependencies
are needed.

See the [validation report](sparse_complex_lsqr_validation_2026-11-06.md).
