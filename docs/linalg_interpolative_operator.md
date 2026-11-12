# Matrix-free randomized interpolative decomposition

Module: `scifort_linalg_interpolative_operator`

SciPy analogue: [`scipy.linalg.interpolative.interp_decomp` with a `LinearOperator`](https://docs.scipy.org/doc/scipy-1.17.0/reference/generated/scipy.linalg.interpolative.interp_decomp.html).

Unlike the dense `scifort_linalg_interpolative_random` module, these routines
need **only operator actions**. They do not require an `m`-by-`n` array.
They accept callbacks for both `A*x` and `A^T*x` (real) or `A^H*x`
(complex). Callbacks must return zero status on success. All outputs use
`real(dp)` or `complex(dp)` and one-based column indices.

## API

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_linalg_interpolative_operator, only : randomized_interp_decomp_operator

integer :: k, status
integer, allocatable :: idx(:)
real(dp), allocatable :: proj(:,:), skeleton(:,:)
real(dp) :: err

! Fixed rank: k is input.
k = 3
call randomized_interp_decomp_operator(m,n,k,matvec,rmatvec, &
    idx,proj,status,seed=1188_int64,oversampling=8, &
    skeleton=skeleton)

! Adaptive: k is output, with random-probe verification.
call randomized_interp_decomp_operator(m,n,1.0e-5_dp,matvec,rmatvec, &
    k,idx,proj,status,seed=1188_int64,oversampling=8, &
    verification_probes=16,skeleton=skeleton,estimated_error=err)
```

For complex matrices, use the same argument lists with
`randomized_interp_decomp_operator_complex`, complex coefficients and
skeleton, and *Hermitian* adjoint callbacks. A callback has this shape:

```fortran
subroutine matvec(x,y,status)
    real(dp), intent(in) :: x(:)    !! n-component input
    real(dp), intent(out) :: y(:)   !! m-component output
    integer, intent(out) :: status !! zero on success
end subroutine matvec
```

The adjoint callback receives an `m`-component input and returns an
`n`-component output. The routine invokes it with the conjugate-transpose
convention in the complex case, **not** a plain transpose.

`idx` is a one-based permutation of `1:n`, `proj` has shape `(k,n-k)`, and
optional `skeleton` contains the exact `m`-by-`k` selected columns of the
implicit original operator. Thus

```fortran
A(:,idx(k+1:n)) ≈ matmul(skeleton,proj)
```

A factorization can be applied to a vector `x` *without reconstructing A*:

```fortran
y_approx = matmul(skeleton, &
    x(idx(1:k)) + matmul(proj,x(idx(k+1:n))))
```

The existing `reconstruct_matrix_from_id(skeleton,idx,proj,output,status)`
can reconstruct a full matrix if desired. Even without `skeleton=`, the
routine internally uses only `k` forward products to recover the selected
columns; the caller can re-fetch those columns using `A*e_j`.

## Algorithm and memory use

1. Generate `s` independent Gaussian row probes. An adjoint callback computes
   each sketch row `omega^H A`, giving a small `s`-by-`n` matrix.
2. Pivoted QR of the sketch selects columns, as in standard randomized ID.
3. Call `A*e_j` for the `k` selected columns, obtaining the **exact skeleton**.
4. Orthogonalize that skeleton with existing Householder QR. Call `A^H*q_j`
   for its `k` orthonormal basis vectors to fit interpolation coefficients
   by a small SVD least-squares system. In particular, there is no need for
   forward products for all `n` columns.
5. For a tolerance-based call, apply the original operator to independent
   Gaussian vectors and compare those products to the ID products. Grow
   the sketch and search higher ranks until the sampled residual is small.

For a fixed rank `k`, the implementation requires exactly `k` forward
products plus at most `k+oversampling+k` adjoint products (when `k<n`),
with sketch rows capped at `min(m,n)`. It stores an `s`-by-`n` sketch and
`m`-by-`k` skeleton/basis, **not** an `m`-by-`n` matrix. For larger
requested ranks a sketch may approach full-matrix dimensions; this is not
an asymptotic speed guarantee. The matrix-free advantage is greatest for
expensive implicit matrices with `k` and `s` much smaller than `m` and `n`.

For tolerance-based calls, the default is 16 independent residual probes;
`verification_probes` may be increased for greater statistical confidence.
The acceptance threshold includes a conservative margin to limit
false acceptance due to sampling noise. **This does not certify the error**:
`estimated_error` is a Monte Carlo estimate of the *relative Frobenius*
residual, not a rigorous bound. For critical accuracy, verify independently
with more probes or an exact residual where feasible.

## Edge cases and compatibility

- `k=0` returns `idx=1:n`, `proj(0,n)`, `skeleton(m,0)` without any callbacks.
- Empty dimensions are supported, returning rank zero.
- The rank must be between zero and `min(m,n)`; adaptive tolerance must be
  finite and `0 <= eps < 1`. Oversampling must be nonnegative and the
  probe count must be positive.
- Malformed arguments return `linalg_status_invalid_shape`; callback errors,
  nonfinite products, or failed factorization return
  `linalg_status_no_convergence`. Failed calls return empty `idx`/`proj`.
- Results from a given seed are reproducible within one compiler/build.
  The local Park–Miller/Box–Muller generator is not NumPy's generator.
- SciPy's random `LinearOperator` implementation targets a spectral-norm
  criterion and its pivots/rank may differ. This implementation uses a
  probabilistic relative Frobenius-norm criterion.
- No sparse storage format is required, since callers can supply sparse,
  FFT-based, or other implicit products through callbacks.
- This is an independent implementation and does not copy SciPy source.

## Build and test

Compile the existing `kinds.f90`, `linalg.f90`, `linalg_svd.f90`,
`linalg_qr.f90`, and `linalg_interpolative.f90`, followed by
`linalg_interpolative_operator.f90`. With FPM, automatic source discovery
picks up the new module.

Tests:

- `test/test_linalg_interpolative_operator.f90` — native regression cases
- `reference/linalg_interpolative_operator_driver.f90` — reference driver
- `test/compare_linalg_interpolative_operator_scipy.py` — compares with
  NumPy SVD reconstruction and SciPy's actual `LinearOperator` API
- `example/linalg_interpolative_operator_example.f90` — runs real and complex
  rank-three *implicit operators* without ever building a dense matrix

See [validation details](linalg_interpolative_operator_validation_2026-11-10.md).
