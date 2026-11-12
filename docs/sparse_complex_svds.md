# Complex sparse singular-value decomposition (`svds_complex`)

SciFort's `scifort_sparse_linalg_complex_svds` implements a dependency-free
partial SVD for canonical complex binary64 CSR and CSC matrices. The routine
accepts a rectangular matrix and returns the requested singular triplets:

\[
A v_j \approx \sigma_j u_j, \qquad A^H u_j \approx \sigma_j v_j.
\]

```fortran
use scifort_kinds, only : dp
use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_dense
use scifort_sparse_linalg_complex_svds, only : svds_complex
implicit none
type(csr_complex_matrix) :: a
complex(dp), allocatable :: u(:,:), vh(:,:)
real(dp), allocatable :: s(:), residuals(:)
integer :: info, iterations

a = csr_complex_from_dense(reshape([ &
    cmplx(3.0_dp,1.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp), &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(2.0_dp,-1.0_dp,dp), &
    cmplx(0.0_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp) &
], [3,2]))
call svds_complex(a,1,u,s,vh,info,which='LM',tol=1.0e-10_dp, &
    iterations=iterations,residual_norms=residuals)
if (info /= 0) error stop 'complex SVD did not converge'
```

The generic has CSR and CSC overloads:

```fortran
call svds_complex(a,k,u,s,vh,info,which,tol,ncv,maxiter,v0,iterations,residual_norms)
```

- `k`: requested number of triplets, `1 <= k <= min(m,n)`.
- `u`: allocatable complex `(m,k)` left singular vectors.
- `s`: allocatable real `k` singular values, **ascending**, including for
  `which='LM'`.
- `vh`: allocatable complex `(k,n)`, the **conjugate transpose** of the right
  singular-vector matrix; do not treat `vh` as a plain transpose.
- `which`: `'LM'` (default) for largest or `'SM'` for smallest singular values.
- `tol`: relative true two-sided residual tolerance, default `1e-8`, with
  zero clamped to a machine-precision floor.
- `ncv`: maximum Golub--Kahan trial basis width; default
  `min(min(m,n),max(40,8*k+16))` for `'LM'`. The `'SM'` path currently
  **requires `ncv=min(m,n)`**, set automatically when omitted.
- `maxiter`: cap on bidiagonalization steps, bounded above by `ncv`.
- `v0`: optional finite nonzero complex starting vector of length
  `min(m,n)` on the right side of the *tall* orientation: `A` if `m >= n`,
  or `A^H` if `m < n`.
- `iterations`: number of Golub--Kahan steps.
- `residual_norms`: allocatable real `k` array of exact two-sided residuals:
  `hypot(||A*v - s*u||_2,||A^H*u - s*v||_2)`.

`info=0` means that the requested triplets have converged, `info=1`
means the basis/iteration budget was insufficient (the latest computed Ritz
triplets, if any, remain available), `info=-1` indicates invalid arguments,
and `info=-2` denotes numerical breakdown. Invalid inputs initialize empty
outputs. The all-zero matrix is a special case returning orthonormal
coordinate vectors and zero singular values.

## Numerical method and limits

A complex Golub--Kahan iteration with two full reorthogonalization passes
builds a bounded right trial subspace. The small projected operator `A*V`
is factorized by **complex one-sided Jacobi rotations**, and candidate
triplets are checked using both true residual equations. The method does
not construct the full normal matrix `A^H A` or require BLAS, LAPACK,
PROPACK, or ARPACK. On a lucky Krylov breakdown, the iteration can add an
independent coordinate direction to recover missing spectral directions.

When the full shorter-side trial basis is requested, SciFort finishes it
before reporting convergence. This avoids missing repeated dominant
singular values merely because the initial Krylov chain already has several
converged Ritz vectors. For **genuinely partial searches**, however, tiny
true residuals certify the *returned* singular triplets, not necessarily
that every globally largest singular value has been found: invariant
subspaces absent from a starting vector can still be missed.

This implementation is not a restarted PROPACK/ARPACK replacement. It
currently supports complex binary64 CSR/CSC matrices only: no
matrix-free callbacks, alternate SciPy solver backends, nullspace-complete
full SVD beyond `min(m,n)` triplets, or guaranteed fast smallest-value
search on large matrices. The `'SM'` full-basis requirement can need
quadratic memory in the smaller dimension. The projected thin array
`A*V` is dense, but the input matrix `A` is never densified.

Source and tests: `src/sparse_linalg_complex_svds.f90`,
`test/test_sparse_complex_svds.f90`,
`tools/sparse_reference/compare_complex_svds_scipy.py`, and
`tools/sparse_reference/complex_sparse_svds_driver.f90`.
See [validation notes](sparse_complex_svds_validation_2026-11-06.md).
