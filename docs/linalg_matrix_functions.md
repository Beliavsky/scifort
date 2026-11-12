# Dense matrix square root and matrix sign functions

`scifort_linalg_matrix_functions` adds Fortran analogues of
`scipy.linalg.sqrtm` and `scipy.linalg.signm`. The implementation uses
`complex(dp)` for the Schur reduction and real/complex binary64 overloads,
with no new BLAS/LAPACK or other external library dependency.

## API

| Routine | Input | Output | Description |
|---|---|---|---|
| `sqrtm(a,x,status,residual)` | `real(dp)` | `complex(dp)` | Principal complex square root, including negative real eigenvalues |
| `sqrtm(a,x,status,residual)` | `complex(dp)` | `complex(dp)` | Principal complex square root |
| `sqrtm_real(a,x,status,residual)` | `real(dp)` | `real(dp)` | Real square root only; rejects a non-real principal root |
| `signm(a,x,status,tol,max_iter,iterations,residual)` | `real(dp)` | `real(dp)` | Real matrix sign |
| `signm(a,x,status,tol,max_iter,iterations,residual)` | `complex(dp)` | `complex(dp)` | Complex matrix sign |

All arrays are explicit 2-D rank and caller-allocated. Arguments `residual`,
`tol`, `max_iter`, and `iterations` are optional. `a` is read-only; `x` is
initially set to zero and remains zero on failure. Real `sqrtm` inputs produce a
**complex output** because a real matrix can have no real principal square
root, e.g. a matrix with the scalar eigenvalue `-1`. For code requiring a real
output, call `sqrtm_real` instead.

Example:

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_matrix_functions, only : sqrtm, signm
implicit none
real(dp) :: a(2,2), s(2,2), residual
complex(dp) :: r(2,2)
integer :: status

a = 0.0_dp
a(1,1) = 4.0_dp
a(2,2) = -9.0_dp
call sqrtm(a, r, status, residual)
if (status /= linalg_status_success) error stop "sqrtm"
! r(1,1)=2, r(2,2)=(0,3)

call signm(a, s, status, residual=residual)
if (status /= linalg_status_success) error stop "signm"
! s(1,1)=1, s(2,2)=-1
```

## Numerical methods

`sqrtm` uses a complex Hessenberg reduction and shifted QR Schur iteration,
then solves the triangular square-root recurrence

`F(i,j)=(T(i,j)-sum_{k=i+1}^{j-1}F(i,k)*F(k,j))/(F(i,i)+F(j,j))`.

Reconstruct `X=Z*F*Z^H` and verify the true residual `||X*X-A||_F`.
The diagonal complex scalar square root uses the principal branch. For a
**real** input with an exactly negative real eigenvalue, the positive-imaginary
branch is selected, even if Schur QR introduces insignificant imaginary
roundoff at that eigenvalue.

`signm` uses scaled Newton iteration

`X_next = (mu*X + inv(X)/mu)/2`,

with `mu=sqrt(||inv(X)||_F/||X||_F)`. It checks the iterate change,
`||X*X-I||_F`, and commutation with `A`. No `signm` convergence is claimed
solely from small changes in successive iterates.

## Limitations

* Both functions are limited to `matrix_function_max_order=40`. The QR and
  inverse kernels are unblocked and best suited to small/medium dense problems.
* `sqrtm` can return nonconvergence for matrices with non-existent or poorly
  conditioned square roots, including a nonzero nilpotent Jordan block.
* The matrix sign has no analytic definition on eigenvalues with zero real
  part; apart from the explicitly supported zero matrix, such cases can return
  nonconvergence. The routines do not implement SciPy's batch dimensions.
* Near branch cuts, multiple square roots satisfy `X*X=A`. SciPy 1.17's
  `sqrtm` sometimes chooses the negative-imaginary root for real symmetric
  matrices with negative eigenvalues due to roundoff; SciFort chooses the
  positive-imaginary principal branch in those cases.
* The actual returned status codes are `linalg_status_success` (0),
  `linalg_status_invalid_shape` (1), `linalg_status_no_convergence` (4), or
  `linalg_status_singular` (5, for the underlying inversion). For `signm`, an
  imaginary-axis eigenvalue may be reported as nonconvergence rather than
  singular, depending on iteration.

## Build and tests

With `gfortran` and a Python environment containing SciPy and NumPy:

```bash
python tools/linalg_reference/run_matrix_functions_checks.py --mode both --cases 500
```

The script builds both `-O0 -fcheck=all` and `-O3` versions, runs
`test/test_linalg_matrix_functions.f90`, and compares randomized matrices
against SciPy and a symmetric spectral reference on the square-root branch
cut. See `docs/linalg_matrix_functions_validation_2026-11-06.md` for details.
The executable example is `example/linalg_matrix_functions_example.f90`.

No external code was copied: the Schur QR engine was adapted from the
existing MIT-licensed SciFort `src/sparse_linalg_complex_eigs.f90`.
