# RQ factorization (`scipy.linalg.rq`)

SciFort module: `scifort_linalg_rq`. The implementation uses the existing
unblocked real/complex Householder QR routines and requires no LAPACK linkage
or new dependencies.

## Interface

```fortran
use scifort_kinds, only : dp
use scifort_linalg_rq, only : rq, rq_r
integer :: status
real(dp) :: a(3,5), r(3,5), q(5,5), re(3,3), qe(3,5)

call rq(a, r, q, status)                    ! full: R(m,n), Q(n,n)
call rq(a, re, qe, status, mode='economic') ! economy: R(m,k), Q(k,n)
call rq_r(a, r, status)                    ! R-only (no Q formed)
```

`k = min(m,n)`. Each subroutine has real(dp) and complex(dp) variants;
`rq` and `rq_r` also accept one leading batch dimension (Fortran layout
`a(batch,m,n)`). Outputs are caller-allocated. All dummy arguments declare
`INTENT`, and status codes are from `scifort_linalg`:

- `linalg_status_success`: valid factorization.
- `linalg_status_invalid_shape`: unsupported mode, wrong output shape,
  or empty batch dimension.
- `linalg_status_no_convergence`: non-finite input or failure in the
  underlying QR factorization.

`rq` modes are `'full'` (default) and `'economic'`. To use SciPy's
`mode='r'` behavior, call `rq_r` with an `m`-by-`n` output. `rq_r` skips
explicit construction of Q. A single batch dimension is available,
but not arbitrary-rank broadcasting. Zero-sized *matrix dimensions*
are supported: in full mode Q is the `n`-by-`n` identity (possibly empty),
while R is an appropriately shaped zero matrix. A zero-sized *batch*
is rejected, matching SciPy's batch convention.

`check_finite` is an optional logical argument, default `.true.`. When set
`.false.`, the up-front finite check is skipped, but **the underlying
SciFort QR kernel still rejects non-finite matrices**; no promise is
made that nonfinite inputs can be processed.

The reconstructed matrix is `A = R Q`. Rows of Q are orthonormal in both
modes (`Q Q^T = I` for real data and `Q Q^H = I` for complex data).
`R` is upper trapezoidal with its nonzero diagonal right-aligned in full
mode. RQ factors are generally non-unique: signs or complex phases may
be different from LAPACK/SciPy even for the same input.

## Implementation

Let `J_d` be the order-reversal permutation of dimension d. Build the
`n`-by-`m` matrix `B = J_n A^H J_m`, compute `B = Q_b R_b`, and set
`R = J_m R_b^H J_k` and `Q = J_k Q_b^H J_n`, where `k=n` in full mode
and `k=min(m,n)` in economic mode. Reversals are implemented using index
operations rather than explicit permutation matrices. R-only mode
extracts the upper-trapezoidal entries from the packed QR array and does
not unpack Q.

As with other unblocked SciFort QR operations, this is currently a
serial dense factorization, not the blocked LAPACK `*GERQF` implementation.
It is intended for correctness and portability, not peak BLAS-3 speed.

See `example/linalg_rq_example.f90`, `test/test_linalg_rq.f90`, and
`reference/check_linalg_rq.py`. Independent comparisons use SciPy 1.17.0.
