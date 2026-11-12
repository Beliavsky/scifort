# Dense SVD (economy and full), singular values and least squares

Module `scifort_linalg_svd` implements real and complex binary64 singular-value
decompositions, Moore-Penrose pseudoinverses, and SVD least-squares solves for
rectangular, including rank-deficient matrices. It uses intrinsic Fortran,
without calling or bundling BLAS/LAPACK. The implementation applies one-sided
Jacobi rotations and reuses concepts from SciFort's polar decomposition.

## Public routines

| Generic | Call | Output and scope |
| --- | --- | --- |
| `svd` | `(a,u,s,vh,status,max_sweeps,full_matrices)` | Economy (default) or full square factors |
| `svdvals` | `(a,s,status,max_sweeps)` | Singular values via a dedicated values-only Jacobi iteration |
| `cond2` | `(a,condition,status,max_sweeps)` | Spectral 2-norm condition number (SciFort extension) |
| `matrix_rank` | `(a,rank,status,tol,max_sweeps)` | SVD-based numerical matrix rank (SciFort extension) |
| `diagsvd` | `(s,m,n,d,status)` | Rectangular diagonal matrix `d(m,n)` |
| `pinv` | `(a,ap,status,rank,rcond)` | Moore-Penrose pseudoinverse, shape `(n,m)` |
| `svd_lstsq` | `(a,b,x,status,rank,rcond)` | Minimum-norm solutions, multi-RHS only |
| `scifort_linalg_qr.lstsq` | `(a,b,x,status,rank,rcond,method)` | Vector or multi-RHS, QR/SVD dispatch |

For an `m`-by-`n` input, `k=min(m,n)` and the default `svd`
requires caller-owned `u(m,k)`, `s(k)` and `vh(k,n)` outputs. With
`full_matrices=.true.`, the outputs instead have shapes `u(m,m)`, `s(k)`
and `vh(n,n)`. Both modes reconstruct the matrix as
`A=U(:,1:k) diag(s) Vh(1:k,:)`; `diagsvd(s,m,n,d,status)` supplies the
full rectangular diagonal matrix for `matmul(matmul(u,d),vh)`. The
optional argument is placed **after** `max_sweeps` so existing positional
calls are unchanged. Unlike SciPy's default, SciFort's default remains
**economy SVD** for backward compatibility.

For `m > n`, a complete `m`-by-`m` left orthonormal basis is constructed;
for `n > m`, the right basis is extended to `n` vectors. Completion uses
pivoted canonical vectors with two passes of modified Gram-Schmidt, and
returns `linalg_status_no_convergence` when a stable completion cannot
be constructed. Empty input matrices are still not supported by `svd` or
`svdvals` (an intentional current limitation).

Singular values are nonnegative and sorted descending. Numerically negligible
directions are completed to orthonormal bases, and their singular values
are set to zero. `svdvals` uses the same one-sided Jacobi rotation criterion as `svd`, but
runs a separate singular-values-only iteration. It avoids constructing and
orthogonalizing both singular-vector factors. For real input, it works in
real arithmetic instead of promoting to complex arithmetic; for complex
input, it omits the right-vector rotation updates. This saves temporary
storage and commonly reduces runtime, without changing the `svd` interface.
`svdvals` retains the previous rank-floor convention and zeros singular values
at or below `8*epsilon(1.0_dp)*s(1)`; values below that limit are not resolved.

`cond2` and `matrix_rank` are **SciFort convenience extensions**: they
correspond mathematically to `numpy.linalg.cond(a, p=2)` and
`numpy.linalg.matrix_rank(a)`, not to distinct public `scipy.linalg` functions.
They support rectangular real or complex matrices. `cond2` returns
`s(1)/s(k)` with `k=min(m,n)` and **positive IEEE infinity** for singular
matrices or when the singular-value ratio exceeds the representable range.
The routine reports a success status for singular input; the infinite
condition is the result, not a failure. `matrix_rank` counts singular values
strictly larger than the optional **absolute** `tol`; the default threshold
is `max(m,n)*epsilon(1.0_dp)*s(1)`. Because `svdvals` clips unresolved
singular values at the Jacobi floor, condition numbers of near-singular
matrices may be infinite or differ from LAPACK/SciPy near machine precision.
Both routines reject empty and nonfinite matrices with status codes and set
outputs to zero on failure.
`diagsvd` supports zero-size output matrices, checks the array dimensions
and requires all diagonal values to be finite. It does not require the
diagonal entries to be sorted or nonnegative, matching SciPy's behavior.

`pinv(a,ap,...)` and `svd_lstsq(a,b,x,...)` choose singular values greater
than `rcond*s(1)`; default `rcond=max(m,n)*epsilon(1.0_dp)`. A rank-deficient
matrix is **not** an error: pseudoinverses and least-squares solutions use the
retained singular subspace, with the solution of minimum Euclidean 2-norm.
`rank` is the optional count of retained singular values. `svd_lstsq`
accepts a two-dimensional RHS, including zero RHS columns, and writes
`x(n,nrhs)`. For a vector RHS, use `scifort_linalg_qr.lstsq`.

`scifort_linalg_qr.lstsq` has an optional `method` argument:

- `'auto'` (default): QR for matrices diagnosed as full rank; SVD otherwise.
- `'qr'`: retain the historical QR-only behavior, including `singular` status.
- `'svd'`: use singular values to determine rank for every input.

`qr_solve`, which reuses QR factors for new right-hand sides, remains
QR-only. It still rejects rank-deficient systems.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg_qr, only : lstsq
use scifort_linalg_svd, only : svd, pinv
real(dp) :: a(3,2), b(3), x(2), ap(2,3)
real(dp) :: u(3,2), s(2), vh(2,2)
integer :: status, rank

a = 0.0_dp
a(1,:) = [1.0_dp,1.0_dp]
b = [2.0_dp,1.0_dp,0.0_dp]
call lstsq(a,b,x,status,rank,method='svd') ! x=[1,1], rank=1
call svd(a,u,s,vh,status)
call pinv(a,ap,status,rank)
```

## Validation and limitations

The algorithm normalizes the input magnitude before Jacobi sweeps to reduce
underflow/overflow in sums of squares. Output scale is restored afterward.
`max_sweeps` defaults to 120 and can be overridden; nonconvergence or
nonfinite arithmetic produces `linalg_status_no_convergence=4`. Invalid
shapes, negative/nonfinite cutoffs or sweep limits smaller than one return
`linalg_status_invalid_shape=1`. Numerical outputs are zeroed on failure.

This is an **unblocked, intrinsic-only implementation** intended for small to
medium dense systems. It is slower than tuned LAPACK `GESDD`/`GESVD` for large
matrices and cannot resolve singular values arbitrarily close to floating-point
roundoff. Near the rank cutoff, SciFort and SciPy can classify rank
differently, and forward solutions can differ more than their residuals.
In `'auto'` mode a QR rank estimate of full rank prevents the SVD fallback,
so use `method='svd'` for a consistent SVD rank decision.
Full square factors require O(m^2+n^2) memory and additional O(max(m,n)^3)
orthogonal-completion work, even when the input matrix is very rectangular.
Choose the default economy factors for ordinary least-squares and pseudoinverse
problems.

See [validation results](linalg_svd_validation_2026-11-07.md) and
`example/linalg_svd_example.f90`.

## Full-SVD example

```fortran
use scifort_kinds, only : dp
use scifort_linalg_svd, only : svd, svdvals, diagsvd
real(dp) :: a(4,2), u(4,4), s(2), ss(2), vh(2,2), d(4,2)
integer :: status

a = reshape([1.0_dp,2.0_dp,0.0_dp,0.0_dp, &
             1.0_dp,0.0_dp,3.0_dp,0.0_dp],shape(a))
call svd(a,u,s,vh,status,full_matrices=.true.)
call diagsvd(s,4,2,d,status)
! A is reconstructed as matmul(matmul(u,d),vh).
call svdvals(a,ss,status)
```

See `example/linalg_svd_full_example.f90` and the
[full-SVD validation record](linalg_svd_full_validation_2026-11-07.md).

## Singular-value-only and spectral-diagnostics validation

See [the independent check record](linalg_svd_condition_validation_2026-11-07.md),
`test/test_linalg_svd_condition.f90`,
`tools/linalg_reference/run_svd_condition_checks.py`, and
`example/linalg_svd_condition_example.f90`. The values-only benchmark is in
`tools/linalg_reference/bench_svdvals.f90`.
