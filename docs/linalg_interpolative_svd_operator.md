# Matrix-free randomized interpolative SVD

Module: `scifort_linalg_interpolative_svd_operator`

SciPy analogue: [`scipy.linalg.interpolative.svd` applied to a `LinearOperator`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.interpolative.svd.html).

This API computes a thin rank-`k` factorization of a **real or complex implicit
matrix** without ever requiring the dense `m`-by-`n` array. A caller supplies
`matvec` (`A*x`) and `rmatvec` (`A^T*x` for real inputs; `A^H*x` for complex
inputs) callbacks. Outputs satisfy, approximately,

```
A = U * diag(s) * V^H
```

with `(m,k)` left singular vectors, `(k)` singular values in non-increasing
order, and `(n,k)` **right singular vectors** (not `V^H`). Factors are
orthonormal up to numerical accuracy.

## API

The primary generic `randomized_interpolative_svd_operator` supports both
`real(dp)` and `complex(dp)`. A separate alias
`randomized_interpolative_svd_operator_complex` is available for complex
arguments, for consistency with the existing matrix-free ID interfaces.

### Fixed rank

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_linalg_interpolative_svd_operator, only : randomized_interpolative_svd_operator

real(dp) :: u(m,k), s(k), v(n,k)
integer :: status

call randomized_interpolative_svd_operator(m,n,k,matvec,rmatvec,u,s,v,status, &
    seed=18837_int64,oversampling=8)
```

The caller allocates `u`, `s`, and `v` before calling the routine.
`0 <= k <= min(m,n)` is required. Rank zero is valid; no callbacks are
invoked in that case.

### Adaptive relative tolerance

```fortran
integer :: k, status
real(dp), allocatable :: u(:, :), s(:), v(:, :)
real(dp) :: estimated_error

call randomized_interpolative_svd_operator(m,n,1.0e-5_dp,matvec,rmatvec, &
    k,u,s,v,status,seed=18837_int64,oversampling=8, &
    verification_probes=16,estimated_error=estimated_error)
```

The routine selects `k` and allocates the output factors. Relative tolerance
`eps` must satisfy `0 <= eps < 1`. The optional `estimated_error` is the
sampled relative **Frobenius** reconstruction error of the underlying ID.
It is an estimate, not a certified error bound. In particular, it is **not**
SciPy's relative spectral-norm accuracy contract. The rank and random
factors can differ between SciFort and SciPy.

For complex inputs, declare `u` and `v` as `complex(dp)`, while `s` and
`estimated_error` remain `real(dp)`. The adjoint callback must implement a
Hermitian conjugate transpose, not a plain transpose.

### Callback signature

```fortran
subroutine matvec(x, y, status)
    real(dp), intent(in) :: x(:)    !! n-vector input
    real(dp), intent(out) :: y(:)   !! m-vector output
    integer, intent(out) :: status !! zero for success
end subroutine matvec
```

`rmatvec` takes an `m`-vector and returns an `n`-vector. For complex matrices,
replace both real array types with `complex(dp)` and use the conjugate
transpose in `rmatvec`. Nonzero callback status or nonfinite products cause
the routine to return `linalg_status_no_convergence`. Invalid rank,
tolerance, oversampling, probe count, or output shapes produce
`linalg_status_invalid_shape`. On error, fixed-size outputs are reset to
zero, and adaptive outputs are reset to empty arrays.

## Algorithm and costs

1. The existing matrix-free randomized ID collects a short Gaussian row
   sketch using adjoint products and uses pivoted Householder QR to choose
   the skeleton columns.
2. Forward products applied to `k` coordinate vectors recover the selected
   columns `B`, and adjoint products refine the interpolation matrix `P`.
3. QR on `B` and SVD of a `k`-by-`n` core convert `B*P` to orthonormal factors.
   The skeleton is rescaled before QR and SVD; singular values are restored
   to the original units afterward.
4. In adaptive mode, additional independent forward probes assess the
   Frobenius reconstruction error, and the underlying ID grows its sampled
   rank if necessary.

For fixed `k`, the underlying ID invokes exactly `k` forward products,
plus at most `k + min(k+oversampling,min(m,n))` adjoint products (when
`k<n`). The SVD conversion invokes **no additional callbacks**. Dominant
workspace includes an `m`-by-`k` skeleton, `k`-by-`n` interpolation/SVD core,
and a `(k+oversampling)`-by-`n` sketch. The `m`-by-`n` original matrix
is never allocated by this module. Runtime depends on callback costs,
rank, and the existing Jacobi SVD backend; no claim of speedup has been
established by benchmarks.

## Example, tests and limitations

- `example/linalg_interpolative_svd_operator_example.f90` creates implicit
  real and complex rank-three operators without a dense matrix.
- `test/test_linalg_interpolative_svd_operator.f90` checks real/complex
  reconstruction, orthogonality, huge/tiny scales, zero/empty dimensions,
  rank and tolerance interfaces, and error paths.
- `reference/linalg_interpolative_svd_operator_driver.f90` plus
  `test/compare_linalg_interpolative_svd_operator_scipy.py` provide
  independent NumPy and actual SciPy `LinearOperator` comparisons.
- The random generator is a private deterministic Park-Miller/Box-Muller
  sequence rather than NumPy's `Generator`, so identical seeds do not
  imply identical samples or singular vectors.
- The input operator must be linear, and callback output dimensions must
  match `m` and `n`; there is no automatic dimension discovery.
- SVD factors have arbitrary column signs/phases, so test reconstructed
  products and singular values rather than direct equality of vector entries.

See [validation results](linalg_interpolative_svd_operator_validation_2026-11-10.md).
