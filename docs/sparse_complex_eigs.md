# Complex nonsymmetric sparse eigenvalues (`eigs_complex`)

`scifort_sparse_linalg_complex_eigs` adds a dependency-free partial eigensolver
for **nonsymmetric complex binary64** CSR and CSC matrices. Unlike
`eigsh_complex`, which requires Hermitian input, `eigs_complex` computes
complex eigenvalues and right eigenvectors for general square matrices.

```fortran
use scifort_kinds, only : dp
use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_coo
use scifort_sparse_linalg_complex_eigs, only : eigs_complex
implicit none
type(csr_complex_matrix) :: a
complex(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
real(dp), allocatable :: residual_norms(:)
integer :: info

a = csr_complex_from_coo(5,5,[1,2,3,4,5,1],[1,2,3,4,5,3], &
    [cmplx(1.0_dp,0.0_dp,dp), cmplx(2.0_dp,0.5_dp,dp), &
     cmplx(3.0_dp,-0.2_dp,dp), cmplx(4.0_dp,1.0_dp,dp), &
     cmplx(5.0_dp,-1.0_dp,dp), cmplx(0.3_dp,0.8_dp,dp)])
call eigs_complex(a,2,eigenvalues,eigenvectors,info,which='LM', &
    ncv=5,tol=1.0e-10_dp,residual_norms=residual_norms)
if (info /= 0) error stop 'complex eigs failed to converge'
print *, eigenvalues
print *, residual_norms
```

The generic interface is:

```fortran
call eigs_complex(a,k,eigenvalues,eigenvectors,info,which,tol, &
    ncv,maxiter,v0,iterations,matvecs,residual_norms,sigma,shift_cache)
```

`a` may be a canonical one-based `csr_complex_matrix` or `csc_complex_matrix`.
`k` must satisfy `1 <= k < n` for square `n`-by-`n` matrices. The allocated
outputs `eigenvalues(k)` and `eigenvectors(n,k)` contain complex eigenvalues
and **right** eigenvectors, sorted by decreasing selection priority:

| `which` | Selection priority |
| --- | --- |
| `LM` (default) | Largest modulus |
| `SM` | Smallest modulus |
| `LR` / `SR` | Largest / smallest real part |
| `LI` / `SI` | Largest / smallest imaginary part |

`tol` defaults to `1e-8` and is clamped at a machine-precision floor. `ncv`
controls the maximum trial-subspace width, from `k+1` to `n`, defaulting to
`min(n,max(k+min(k+8,n-k),24))`. Unshifted `SM` (shifted or unshifted) requires `ncv=n`; otherwise the call
returns invalid-input status. `maxiter` is a cap on **explicit restart cycles**,
not on sparse matrix-vector products, defaulting to 30. `v0` is an optional
finite nonzero initial vector of length `n`. `iterations` reports completed
restart cycles. In unshifted mode `matvecs` reports sparse matrix-vector applications;
with `sigma` it counts applications of the inverse operator using triangular solves
(not the additional original-`A` products used to verify true residuals).
`residual_norms(j)` is the true Euclidean norm
`||A*eigenvectors(:,j)-eigenvalues(j)*eigenvectors(:,j)||_2`.

The integer status is `0` for convergence, `1` for a search budget exhausted
without certification, `-1` for invalid inputs, `-2` for a projected
Schur failure or numerical breakdown, and `-3` for a singular shifted matrix. On `info=1`, the most recent Ritz
pairs and their residuals are still available. An invalid call leaves empty
outputs. The zero matrix returns zero eigenvalues and coordinate eigenvectors.

## Implementation and limitations

Two-pass Arnoldi orthogonalization constructs a bounded complex trial basis.
A dense **small projected matrix** is solved by a Householder-QR complex Schur
method already used in SciFort's real-matrix Arnoldi solver. The real matrix
`A` is never converted to a dense matrix. For restarts, leading Ritz directions
and a true residual seed are retained explicitly; this is **not** ARPACK's
implicit QR restart. Optional shift-invert uses a single sparse LU factorization of
`A-sigma*I` reused across the Arnoldi iterations. See
[sparse complex shift-invert](sparse_complex_eigs_shift.md).

Convergence uses the actual Ritz residuals with a scale containing the matrix
infinity norm. This certifies computed eigenpairs to the requested residual
criterion, but does not prove that *all* eigenvalues of larger desired modulus
(or real/imaginary part) have been found in a truncated Krylov subspace,
particularly when invariant subspaces or repeated/clustered eigenvalues are
not represented in the starting vector. Repeated or defective eigenvalues are
not guaranteed to yield all independent eigenvectors. `SM` requires a full
basis because partial smallest-transformed-magnitude searches are not yet reliable.

The memory requirement is O(`n*ncv + ncv**2`) plus the sparse input matrix;
very large `ncv` may be expensive. Matrix-free and reusable-factor interfaces are described in
[sparse complex matrix-free eigs](sparse_complex_eigs_operator.md).
There is no generalized problem `A*x=lambda*M*x`, left eigenvectors, or
ARPACK-compatible nonconvergence exception.

Sources: `src/sparse_linalg_complex_eigs.f90`,
`test/test_sparse_complex_eigs.f90`, and
`tools/sparse_reference/compare_sparse_complex_eigs_scipy.py`.
See [validation](sparse_complex_eigs_validation_2026-11-06.md).
