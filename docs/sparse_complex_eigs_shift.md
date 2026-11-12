# Complex shift-invert eigenvalues (`eigs_complex`, `sigma=`)

SciFort now supports an optional **complex shift** when finding eigenpairs of
nonsymmetric complex binary64 sparse matrices with `eigs_complex`.

```fortran
use scifort_kinds, only : dp
use scifort_sparse_complex, only : csr_complex_matrix, csr_complex_from_coo
use scifort_sparse_linalg_complex_eigs, only : eigs_complex
implicit none
type(csr_complex_matrix) :: a
complex(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
real(dp), allocatable :: residual_norms(:)
integer :: info

a = csr_complex_from_coo(6,6,[1,2,3,4,5,6,2], &
    [1,2,3,4,5,6,5], &
    [cmplx(1.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.1_dp,dp), &
     cmplx(3.0_dp,-0.2_dp,dp),cmplx(4.0_dp,0.3_dp,dp), &
     cmplx(5.0_dp,-0.1_dp,dp),cmplx(6.0_dp,0.2_dp,dp), &
     cmplx(0.2_dp,0.1_dp,dp)])
call eigs_complex(a,2,eigenvalues,eigenvectors,info, &
    sigma=cmplx(3.3_dp,-0.1_dp,dp),which='LM',ncv=6, &
    tol=1.0e-10_dp,residual_norms=residual_norms)
if (info /= 0) error stop 'shifted eigs did not converge'
print *, eigenvalues
print *, residual_norms
```

## Algorithm

A pivoted sparse LU factorization of `A-sigma*I` is constructed **once** and
reused to apply the operator `B=(A-sigma*I)^(-1)` during complex Arnoldi. Each
Ritz value `mu` of `B` is mapped to an eigenvalue of `A` with

`lambda = sigma + 1/mu`.

When `sigma` is specified, `which` refers to the **transformed** eigenvalues
`mu`, as in SciPy for complex matrices. Thus `which='LM'` finds eigenvalues
**near** the shift, not eigenvalues with large original magnitude. The six
modes (`LM`, `SM`, `LR`, `SR`, `LI`, `SI`) retain their usual priority definitions
in the transformed spectrum. `SM` still requires `ncv=n`, because a reliable
partial smallest-magnitude search has not been implemented.

The method returns *original* eigenvalues and normalized right eigenvectors.
Convergence is decided with **true original-system residuals**
`||A*v-lambda*v||_2`, scaled by `max(1, ||A||_inf, |lambda|)`.

The existing CSR/CSC generic interface accepts optional `sigma` and
`shift_cache=` arguments. A cache may be prepared once and reused with the
same matrix and shift, avoiding repeated LU factorization. See
[matrix-free eigs and factor caching](sparse_complex_eigs_operator.md).
`matvecs` counts inverse-operator applications, i.e.
LU forward/back substitutions; it does **not** include the original-system
matrix-vector products used to check residuals. `iterations` counts explicit
Arnoldi restart cycles. A `sigma` is validated for finite real/imaginary parts.

## Status codes and limitations

- `0`: converged, with verified true residuals.
- `1`: restart budget exhausted, with current approximate eigenpairs/residuals.
- `-1`: invalid matrix, options, dimensions, or nonfinite `sigma`.
- `-2`: factorization/solve numerical breakdown or projected eigensystem failure.
- `-3`: singular shifted matrix (e.g. `sigma` coincides with an eigenvalue).

This is a **bounded explicit-restart Arnoldi method**, not ARPACK's implicitly
restarted algorithm. The complex LU currently uses unblocked row pivoting
without fill-reducing ordering. Searches using *different* `sigma` values
require separate factorizations; one `shift_cache` may be reused for the same
matrix and shift. Large systems can be much slower or more memory-intensive
than SciPy/SuperLU+ARPACK. Singular or poorly conditioned shifts should be
moved away from the spectrum. Generalized eigenproblems and shift-changing
LU updates are not yet supported. Matrix-free inverse callbacks and reusable
factor caches are [described here](sparse_complex_eigs_operator.md).

Testing: `test/test_sparse_complex_eigs_shift.f90` and
`tools/sparse_reference/compare_sparse_complex_eigs_shift_scipy.py`. See the
[shift-invert validation notes](sparse_complex_eigs_shift_validation_2026-11-06.md).
