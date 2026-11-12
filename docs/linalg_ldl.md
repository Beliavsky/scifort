# Pivoted indefinite symmetric / Hermitian factorization (`scipy.linalg.ldl`)

`scifort_linalg_ldl` provides a dependency-free dense `ldl` decomposition
for real symmetric, complex Hermitian, and complex symmetric matrices.
The algorithm uses **Bunch-Kaufman symmetric diagonal pivoting** with
1-by-1 and 2-by-2 diagonal blocks. It does not require positive definiteness.

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_ldl, only : ldl
real(dp) :: a(n,n), lu(n,n), d(n,n)
integer :: perm(n), status

call ldl(a,lu,d,perm,status)
if (status == linalg_status_success) then
    ! A = lu * d * transpose(lu)
    ! lu(perm,:) is unit lower triangular
end if
```

For a complex Hermitian matrix the factorization is
`A = lu * d * transpose(conjg(lu))`, and for a complex symmetric matrix
it is `A = lu * d * transpose(lu)`.

```fortran
complex(dp) :: a(n,n), lu(n,n), d(n,n)
integer :: perm(n), status

call ldl(a,lu,d,perm,status)                     ! Hermitian (default)
call ldl(a,lu,d,perm,status,hermitian=.false.)   ! complex symmetric
call ldl(a,lu,d,perm,status,lower=.false.)       ! upper input and factor
```

The output convention matches [`scipy.linalg.ldl`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ldl.html):

- `lu` is in original row order and may **not** be triangular directly.
- `perm` contains **one-based** indices; `lu(perm,:)` is unit lower triangular
  by default or unit upper triangular when `lower=.false.`.
- `d` is block diagonal and has blocks of size at most 2. The factors
  reconstruct **the original matrix** directly; no external permutation
  of the input is required to reconstruct it.
- `lower=.true.` (default) reads **only the input lower triangle**, while
  `lower=.false.` reads **only the upper triangle**. The opposite triangle
  need not be initialized. The selected triangle must be finite.
- For complex input, `hermitian=.true.` is the default. To avoid silently
  discarding imaginary diagonal components, a Hermitian input with a non-real
  diagonal returns `linalg_status_invalid_shape`, unlike SciPy's warning.
- Dimension zero is valid. Symmetric matrices with exact zero diagonal
  pivots can be factored, even when singular; they may have zero `d` entries.

## Status and limitations

| Status | Meaning |
|---|---|
| `linalg_status_success` (`0`) | Valid decomposition produced |
| `linalg_status_invalid_shape` (`1`) | Invalid dimensions or imaginary diagonal for Hermitian input |
| `linalg_status_no_convergence` (`4`) | Nonfinite selected input or nonfinite numerical intermediate |
| `linalg_status_singular` (`5`) | Numerically singular 2-by-2 pivot that cannot be inverted stably |

On failure, the output arrays `lu` and `d` are zero and `perm` is zero;
status is never silently ignored. This module is a **new** independent
unblocked implementation; it does not use external BLAS/LAPACK routines.
Time complexity is O(n^3) and working storage is O(n^2). It has not been
optimized for large dense matrices, nor does it offer batched input.

The block determinant is evaluated in scaled coordinates to reduce overflow;
a near-singular 2-by-2 block is rejected rather than returning unbounded
multipliers. It does not provide LAPACK-grade complete pivoting or a condition
number estimate. As with any floating-point factorization, extreme scaling
and near-singular matrices can amplify roundoff.

See [native tests](../test/test_linalg_ldl.f90),
[1,200-case SciPy checks per compiler configuration](../tools/linalg_reference/run_ldl_checks.py),
[example](../example/linalg_ldl_example.f90) and
[validation report](linalg_ldl_validation_2026-11-09.md).
