# Dense Householder QR and least-squares solvers (`scipy.linalg`)

`scifort_linalg_qr` implements unblocked Householder QR in real and complex
binary64, optionally with column pivoting. The code uses only intrinsic
Fortran operations, no BLAS/LAPACK runtime dependency. Public routines are
subroutines taking caller-owned output arrays, rather than SciPy's
value-returning functions.

| Routine | Signature | Purpose |
| --- | --- | --- |
| `qr` | `(a, q, r, status, pivots)` | Compute reduced/complete QR; `pivots` is optional |
| `qr_factor` | `(a, af, tau, status, pivots)` | Store compact Householder reflectors and R |
| `qr_unpack` | `(af, tau, q, r, status)` | Construct Q and R from stored factors |
| `qr_solve` | `(af, tau, pivots, b, x, status, rank, rcond)` | Reuse pivoted factors for full-column-rank least squares |
| `lstsq` | `(a, b, x, status, rank, rcond, method)` | Full-rank QR and rank-deficient SVD least squares |

`pivots`, `rank`, and `rcond` are optional for `qr`, `qr_factor`, and `lstsq`
where shown. For `qr_solve`, `pivots` is required. Both `lstsq` and `qr_solve`
support a length-`m` vector RHS or an `(m,nrhs)` matrix RHS, where `A` is
`(m,n)`. `nrhs=0` is supported for matrix RHS. `m` and `n` must be positive.

## QR output and pivoting

Let `k=min(m,n)`. For `qr` and `qr_unpack`, the output dimensions select the
mode:

- **Reduced**: `q(m,k)`, `r(k,n)`.
- **Complete**: `q(m,m)`, `r(m,n)`.

If `pivots` is omitted from `qr`/`qr_factor`, the result satisfies
`A = Q R`. If supplied, column pivoting is enabled and the result satisfies
`A(:,pivots) = Q R`, where `pivots` is a **one-based final column permutation**.
It is *not* LAPACK's array of incremental row interchanges. The `R` returned
by `qr` or `qr_unpack` has exact zeros below its diagonal. `Q^H Q=I`; for
real arrays, `H` means transpose. Factor sign/phase conventions may differ
from SciPy's, so compare reconstruction and orthogonality, not raw entries.

`qr_factor` stores `R` in the upper trapezoid of `af(m,n)` and the tails of
its Householder vectors below the diagonal. It returns the vector `tau(k)`.
For this implementation, each elementary reflector is Hermitian:
`H_i = I - tau(i)*v_i*v_i^H`, with `v_i(1)=1`; the `tau` value is real even in
complex arithmetic. Unpacking leaves the input factorization unchanged.

## Least squares: QR for full rank, SVD for deficient rank

For `m >= n` with independent columns, `lstsq` computes the unique minimizer
of `||A*x - b||_2` using column-pivoted Householder QR. For `m < n` with
independent rows, it computes the unique **minimum-Euclidean-norm** solution
to `A*x=b` by factorizing `A^H` (or `A^T` for real input). Neither path forms
normal equations. `qr_solve` is a reusable factor-based method **only for
`m >= n`** and full column rank; its factorization should be produced by:

```fortran
use scifort_linalg_qr, only : qr_factor, qr_solve
! For m >= n, assume a(m,n), af(m,n), tau(n), pivots(n),
! b(m,nrhs), x(n,nrhs), and status have already been declared.
call qr_factor(a, af, tau, status, pivots)
if (status == 0) call qr_solve(af, tau, pivots, b, x, status)
! Repeat qr_solve with other RHS arrays without refactorizing a.
```

`rank` is an optional output indicating the estimated rank (from the
diagonal of the pivoted `R`), and `rcond` is an optional nonnegative threshold.
The default is `epsilon(1.0_dp)*max(m,n)` for `lstsq`, or
`epsilon(1.0_dp)*m` for tall `qr_solve`. An entry is counted if
`abs(R(i,i)) > rcond*maxval(abs(diagonal(R)))`. This is a **heuristic
rank estimate**, not an SVD-based rank determination.

The optional `method` string selects `'auto'` (default), `'qr'` or `'svd'`.
`'auto'` uses QR when it reports full rank; otherwise it falls back to the
Jacobi SVD in `scifort_linalg_svd` and returns an SVD minimum-norm solution
with status `linalg_status_success`. `'qr'` preserves the old
`linalg_status_singular` and zero-output behavior on deficiency. `'svd'`
always uses the SVD and reports the rank based on its singular values.
`qr_solve` remains QR-only and still rejects deficient rank.

`lstsq` does not return singular values or residual sums of squares. For badly
scaled or nearly rank-deficient matrices, forward accuracy is limited. QR
and SVD rank heuristics can differ near the default threshold. In particular,
`'auto'` does not force SVD when QR *estimates* full rank. Specify
`method='svd'` to obtain a singular-value-based rank determination on every
call. See the [SVD API](linalg_svd.md).

## Status, validation and limitations

All failures zero caller-owned numerical outputs (`pivots=0` if applicable).
Success is `linalg_status_success=0`. Invalid shapes, malformed permutations
or `rcond < 0` return `linalg_status_invalid_shape=1`. Nonfinite inputs or
arithmetic return `linalg_status_no_convergence=4`. A rank-deficient `qr_solve` or `lstsq(method='qr')` returns
`linalg_status_singular=5`; default `lstsq` falls back to SVD. Real and complex arrays
receive identical shape and status checks.

Algorithms are unblocked, recompute trailing column norms for pivot selection,
and use `O(m*n)` compact-factor storage plus output storage. They are suited
to portable moderate-sized matrices, not designed to replace highly tuned
LAPACK `GEQRF`, `GEQP3`, `ORMQR`/`UNMQR`, or `GELSD` for large matrices.

See `example/linalg_qr_example.f90` and
`docs/linalg_qr_validation_2026-11-07.md`.
