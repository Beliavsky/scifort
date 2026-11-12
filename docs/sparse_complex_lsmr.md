# Complex sparse LSMR for rectangular CSR and CSC matrices

The `scifort_sparse_linalg_complex_lsmr` module implements the complex
binary64 Golub--Kahan/LSMR least-squares algorithm corresponding to
`scipy.sparse.linalg.lsmr`.  It uses the Hermitian adjoint `A**H`, **not**
plain transpose, and does not form the normal-equations matrix `A**H*A`.

## Public API

```fortran
use scifort_kinds, only : dp
use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
use scifort_sparse_linalg_complex_lsmr, only : lsmr_complex
implicit none
type(csr_complex_matrix) :: a
complex(dp), allocatable :: x(:)
complex(dp) :: b(3)
integer :: istop, itn
real(dp) :: normr, normar

a = csr_complex_from_dense(reshape([ &
    cmplx(1.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), cmplx(2.0_dp,1.0_dp,dp), &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(3.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp) &
], [3,2]))
b = [cmplx(1.0_dp,1.0_dp,dp), cmplx(2.0_dp,0.0_dp,dp), cmplx(1.0_dp,-1.0_dp,dp)]
call lsmr_complex(a,b,x,istop,itn=itn,normr=normr,normar=normar, &
    atol=1.0e-12_dp,btol=1.0e-12_dp,maxiter=100)
```

The generic is `lsmr_complex(a,b,x,istop,itn,normr,normar,norma,conda,normx,
  damp,atol,btol,conlim,maxiter,x0)` with CSR and CSC overloads. Only
`a`, `b`, `x`, and `istop` are required; `x` is allocated with size `n` on
valid solves and size zero on invalid input or arithmetic failure.

- `damp` is nonnegative, with default `0`; the objective is
  `||A*x - b||_2**2 + damp**2 * ||x - x0||_2**2`.
- `x0` is an optional initial guess and center for regularization. If absent,
  the center is zero.
- `atol`, `btol` default to `1e-6`; `conlim` defaults to `1e8`, with zero
  disabling the condition-limit stop; `maxiter` defaults to `min(m,n)`.
- `normr` is the *true augmented* residual norm, `normar` is the *true*
  normal-equation gradient norm, `normx` is the true norm of the returned
  solution; `norma` and `conda` are iterative estimates.
- `istop=0` means zero initial residual/stationary right-hand side; 1 or 2
  represents a tolerance stop; 3 a condition stop; 7 exhausted iterations;
  negative -1 invalid input, -2 arithmetic breakdown. The SciPy-style 4--6
  machine-precision stopping codes are also possible.

The implementation uses one-based canonical CSR/CSC with `complex(dp)` entries.
It is suitable for sparse rectangular complex systems and handles exact rank
deficiency, damping, over/underdetermined systems and a nonzero initial guess.

The recurrence adapts SciFort's existing BSD-3-Clause real LSMR (which follows
SciPy and Fong--Saunders); see
`LICENSES/scipy-lsmr-BSD-3-Clause.txt`. No new runtime dependencies are added.
The implementation currently lacks user-provided matrix-free callbacks,
preconditioning, and iterative reorthogonalization; exceptionally ill-conditioned
systems can still need higher iterations or tolerances.

## Reproduce focused tests

```sh
python tools/sparse_reference/compare_sparse_complex_lsmr_scipy.py --mode both --cases 360
python tools/sparse_reference/run_native_sparse.py --mode checked
```

The numerical checks use SciPy `lsmr` and independent dense NumPy SVD-based
least-squares solutions, and verify the returned exact residual diagnostics.
See `docs/sparse_complex_lsmr_validation_2026-11-06.md` for execution results.
