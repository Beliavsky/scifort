# Matrix-free complex Arnoldi and reusable shift-invert factors

SciFort extends `scipy.sparse.linalg.eigs` with a **Fortran matrix-free**
interface. The `scifort_sparse_linalg_complex_eigs` module accepts callbacks
that evaluate a complex square linear operator without storing a sparse or
dense matrix. It also accepts a matrix-free inverse operator for shift-invert
searches, and provides a validated cached LU factor for ordinary CSR/CSC calls.
The existing CSR/CSC entry points remain source-compatible.

## Matrix-free operator API

```fortran
call eigs_complex_operator(n,k,matvec,eigenvalues,eigenvectors,info, &
    which,tol,ncv,maxiter,v0,iterations,matvecs,residual_norms, &
    sigma,inverse_matvec,norm_estimate)
```

The same overload is available as `eigs_complex(n,k,matvec,...)`. Both
callbacks must use the explicit interface:

```fortran
subroutine apply(x,y,info)
    complex(dp), intent(in)  :: x(:)  ! length n
    complex(dp), intent(out) :: y(:)  ! length n, supplied by caller
    integer, intent(out) :: info     ! zero on success
end subroutine
```

`matvec` **must** compute the original operator `A*x`. In shift-invert mode,
`sigma` and `inverse_matvec` must be given together; the latter must evaluate
`(A-sigma*I)^(-1)*x` (a matrix-free linear solve). The routine ranks Ritz
values of the transformed inverse operator according to `which` (so `LM`
finds eigenvalues near `sigma`) and returns eigenvalues of **A**, via
`lambda = sigma + 1/mu`. It tests each candidate with the actual original
residual `||A*v-lambda*v||_2`. A callback error or a nonfinite callback result
returns `complex_eigs_breakdown = -2` without claiming convergence.

`norm_estimate=` is an optional finite, nonnegative estimate of the operator
infinity norm, defaulting to `1.0_dp`. It scales the residual acceptance
threshold (`tol * max(1, norm_estimate, abs(lambda))`). A bound that is
excessively large can make this convergence test overly permissive; choose a
reasonable bound. The actual residuals are returned in `residual_norms`.
`matvecs` counts Arnoldi operator applications; it does not count original
operator calls made solely for shifted-residual verification.

The `which` choices are `LM`, `SM`, `LR`, `SR`, `LI`, `SI`. As with the CSR
implementation, `SM` requires `ncv=n`; the solver is an explicitly restarted
Arnoldi/Rayleigh–Ritz method, *not* a fully ARPACK-compatible eigensolver.
`0` denotes convergence, `1` a completed budget with unverified Ritz pairs,
`-1` an invalid input, `-2` breakdown, and `-3` a singular shifted LU when
constructing a cache or using a CSR/CSC shift.

## Example: diagonal matrix-free operator

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_sparse_linalg_complex_eigs, only : eigs_complex_operator
    implicit none
    complex(dp) :: d(5), sigma
    complex(dp), allocatable :: lambda(:), vectors(:,:)
    real(dp), allocatable :: residuals(:)
    integer :: info, i
    d = [(cmplx(real(i,dp),0.1_dp*real(i,dp),dp),i=1,5)]
    sigma = d(3)+cmplx(0.2_dp,0.0_dp,dp)
    call eigs_complex_operator(5,1,apply,lambda,vectors,info, &
        sigma=sigma,inverse_matvec=solve_shift,ncv=5, &
        tol=1.0e-10_dp,residual_norms=residuals)
    if (info /= 0) error stop 'matrix-free eigenvalue search failed'
    print *, lambda(1), residuals(1)
contains
    subroutine apply(x,y,status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = d*x
        status = 0
    end subroutine
    subroutine solve_shift(x,y,status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = x/(d-sigma)
        status = 0
    end subroutine
end program
```

This uses no stored matrix and no LAPACK dependency. For a representative
larger operator, see the 2,500-element matrix-free regression
`test/test_sparse_complex_eigs_operator.f90`.

## Reuse one CSR/CSC shift factorization

```fortran
type(complex_eigs_shift_cache) :: cache
! a is an existing canonical complex CSR or CSC matrix; sigma is complex(dp).
call eigs_complex_prepare_shift(a,sigma,cache,info)
if (info /= 0) error stop 'cannot factor shifted matrix'
call eigs_complex(a,2,lambda,vectors,info,sigma=sigma, &
    shift_cache=cache,ncv=20,tol=1.0e-9_dp)
! Repeat with a different v0 or which, keeping exactly the same a and sigma.
```

The cache stores pivoted sparse LU plus a CSR snapshot of the source matrix.
It validates the matrix shape, indices, **all stored values**, shift, and LU
factor validity each time, rejecting stale caches as `-1`. This avoids
accidentally applying factors to a changed matrix. It also increases memory
usage by one CSR matrix copy; factors are reused without copying them for
each call. Cache preparation can use a CSR or CSC matrix; cached solves can
use either representation of the same canonical matrix. Do not pass a cache
without `sigma`, or with a different `sigma`.

## Verification and limitations

- Native regression: `test/test_sparse_complex_eigs_operator.f90` (all six
  eigenvalue-selection modes, analytic inverse, cache reuse/invalidation,
  invalid callback statuses, nonfinite outputs, 2,500-element matrix-free
  interior-spectrum test).
- Randomized reference: run
  `python tools/sparse_reference/compare_sparse_complex_eigs_operator_numpy.py --mode both --cases 90 --seed 20261108 --scipy-samples 36`.
  It checks matrix-free and cached CSR paths against independent complex
  spectra, true residuals, and direct SciPy `eigs` results, in bounds-checked
  and `-O3` builds. See the [validation report](sparse_complex_eigs_operator_validation_2026-11-06.md).
- Linear operators must be deterministic and have dimension `n`; an inverse
  callback must solve the shifted problem consistently, including any internal
  preconditioner. This solver does not estimate approximation error in an
  inexact inverse solve separately from its original-system residual test.
- LU cache is valid for a **fixed** matrix and shift. It does not update
  factorizations for changed matrices or nearby shifts. For changing matrices,
  rebuild the cache or supply an independent `inverse_matvec` callback.
- No generalized eigensystems, left eigenvectors, or ARPACK implicit restarting.
