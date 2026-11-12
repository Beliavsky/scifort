# Dense triangular, LU, and Cholesky solves (`scipy.linalg`)

`scifort_linalg_dense` adds SciPy-style dense linear solves in pure Fortran.
The public interfaces are *subroutines*, not value-returning functions. Both
real binary64 (`real(dp)`) and complex binary64 (`complex(dp)`) are supported.
Vectors and matrices of right-hand sides have separate generic overloads.

| Generic | Meaning | Signature |
| --- | --- | --- |
| `solve_triangular` | Solve a triangular system | `(a, b, x, status, lower, trans, unit_diagonal)` |
| `lu_factor` | Factor `A` with dense partial pivoting | `(a, lu, ipiv, status)` |
| `lu_solve` | Solve using existing dense LU | `(lu, ipiv, b, x, status, trans)` |
| `cholesky` | Clean Cholesky factor of SPD/Hermitian PD matrix | `(a, f, status, lower)` |
| `cho_factor` | Reusable Cholesky factor | `(a, f, status, lower)` |
| `cho_solve` | Solve using a Cholesky factor | `(f, b, x, status, lower)` |

`lower`, `trans` and `unit_diagonal` are optional, keyword-compatible inputs:

- `lower=.false.` (default) reads the upper triangle; `.true.` reads the lower
  triangle. The other triangle has **no influence** on the solution.
- `trans='N'` (default) solves `A*x=b`; `'T'` solves `A^T*x=b`; `'C'` solves
  `A^H*x=b` (the conjugate transpose). For real arrays, `'C'` equals `'T'`.
- `unit_diagonal=.false.` (default) uses the stored diagonal. With `.true.`,
  triangular solves assume *all diagonal entries are exactly one*, regardless
  of their stored finite values. This option applies only to `solve_triangular`.

A direct triangular solve requires an `n`-by-`n` input matrix, plus either
an `n`-element vector or an `(n,nrhs)` matrix of right-hand sides. The output
`x` must have the same shape as `b`. Caller-owned input arrays are unchanged.
The module checks finiteness of the entire input matrix, including the ignored
triangle and an implicit-unit diagonal; it never silently accepts NaNs.

## Example: LU reuse

```fortran
use scifort_kinds, only : dp
use scifort_linalg_dense, only : lu_factor, lu_solve
real(dp) :: a(3,3), lu(3,3), b(3), x(3)
integer :: ipiv(3), status

! Fill a and b with the intended system.
call lu_factor(a, lu, ipiv, status)
if (status == 0) then
    call lu_solve(lu, ipiv, b, x, status)
    ! Reuse the same factors for another right-hand side, or A^T:
    call lu_solve(lu, ipiv, b, x, status, trans='T')
end if
```

`lu_factor` produces the compact combination of unit-lower `L` (strictly
below the diagonal) and upper `U` (on and above the diagonal), with partial
row pivoting. The pivot array has length `n`: `ipiv(k)` is the **one-based**
row exchanged with row `k` at elimination step `k` (`k <= ipiv(k) <= n`).
This differs from SciPy's **zero-based** pivot indices. Algebraically,
`P*A = L*U`, where `P` is the sequence of those row interchanges.

The module supports repeatedly solving against any number of new right-hand
sides, and solves normal, transpose or conjugate-transpose systems from the
*same* `lu` and `ipiv`. The factorization does not modify `a`. See the
complete `example/linalg_dense_example.f90` for real and complex cases.

## Positive-definite Cholesky factorization and reuse

`cholesky(a, f, status, lower)` and `cho_factor(a, f, status, lower)`
accept real symmetric positive-definite or complex Hermitian positive-definite
matrices. The optional `lower` keyword defaults to `.false.` (upper factor).
If `lower=.true.`, then `A = L L^H`; otherwise `A = U^H U`.
For real arrays, `H` is the transpose. Only the selected input triangle is
used in the arithmetic, matching the underlying LAPACK/SciPy convention.
**All** entries are still checked for finite values. For complex matrices,
the imaginary components of the diagonal are ignored by factorization,
as in LAPACK; diagonal factors are real and positive. Neither routine tests
whether the ignored input triangle is Hermitian-consistent.

The returned factor has **zero in its unused triangle**. SciPy's `cho_factor`
leaves this triangle unspecified, so compare only the requested triangle.
Both `cholesky` and `cho_factor` use the same clean output representation.
`cho_solve(f, b, x, status, lower)` solves the SPD/Hermitian system by
triangular substitutions, accepting an `n`-vector or `(n,nrhs)` array of
right-hand sides. It can reuse a stored factor for multiple subsequent solves.
Callers must pass the same `lower` value used for factorization.

```fortran
use scifort_linalg_dense, only : cho_factor, cho_solve
! a(n,n), f(n,n), b(n,nrhs), x(n,nrhs), status declared by caller.
call cho_factor(a, f, status, lower=.true.)
if (status == 0) call cho_solve(f, b, x, status, lower=.true.)
```

`cho_factor` returns `linalg_status_not_positive_definite` (2) if a
Cholesky pivot is nonpositive (including singular/indefinite matrices).
Nonfinite input or arithmetic returns `linalg_status_no_convergence` (4).
`cho_solve` returns `linalg_status_singular` (5) if a supplied factor has
an exactly zero diagonal; malformed sizes return `linalg_status_invalid_shape`
(1). All routines set their output to zero on failure.

An unblocked Cholesky factorization costs `O(n^3)` arithmetic and stores
`O(n^2)` entries. Each reusable solve costs `O(n^2 nrhs)` operations.
The implementation does not do iterative refinement, rank-revealing
factorization, automatic equilibration, or positive-semidefinite pivoting.
Very ill-conditioned matrices can lose forward accuracy even if their
backward errors remain small.

## Status and numerical limitations

- `linalg_status_success` (0): result available.
- `linalg_status_invalid_shape` (1): invalid matrix/RHS/output shape, pivot
  list or `trans` code. Empty matrix order is currently rejected.
- `linalg_status_no_convergence` (4): nonfinite input or nonfinite computed
  result. (The name is shared with other SciFort numerical routines.)
- `linalg_status_singular` (5): an exactly zero pivot or required diagonal.

Outputs are set to zero on failure, including unsuccessful factorization.
Zero RHS columns are permitted when `n > 0`. `lu_factor` and `lu_solve` are
not iterative methods; they do **not** provide condition-number estimates,
iterative refinement, equilibration, or error bounds. Very ill-conditioned
matrices can have large *forward* errors even with tiny backward residuals.

Triangular solves use `O(n^2 nrhs)` arithmetic and `O(n nrhs)` result storage.
The unblocked dense LU factorization costs `O(n^3)` operations and stores
`O(n^2)` coefficients. No BLAS/LAPACK runtime dependency, fixed maximum
matrix size, SciPy batch dimensions or special SciPy overwrite flags are
required or supported. The real and complex interfaces follow SciPy's
mathematical operations but are not binary/ABI compatible with SciPy.

## Testing

Standalone checked/optimized compilation and independent comparisons with
SciPy 1.17.0:

```sh
python tools/linalg_reference/run_dense_checks.py --mode both --cases 600 --seed 218 --max-order 80
python tools/linalg_reference/run_dense_checks.py --mode both --cases 70 --seed 902 --max-order 60 --stress
```

The harness uses GNU Fortran with either `-O0 -g -fcheck=all -fbacktrace` or
`-O3`, checks residuals as well as solutions, and exercises both real and
complex arithmetic, pivoting, multiple RHSs, unit diagonals, and `N/T/C`
operation modes. The native Fortran test also checks invalid and nonfinite
inputs. See `docs/linalg_dense_validation_2026-11-07.md` for results.

## Dense Cholesky testing

```sh
python tools/linalg_reference/run_dense_cholesky_checks.py --mode both --cases 450 --seed 390 --max-order 70
python tools/linalg_reference/run_dense_cholesky_checks.py --mode both --cases 130 --seed 391 --max-order 85 --stress --tolerance 1e-7
```

The native tests cover `n=9`, both triangles, real and complex factors,
vector/matrix right-hand sides, factor reuse, and invalid/nonpositive inputs.
The reference harness compares factor entries (only the selected triangle),
solutions for two right-hand side arrays and normwise backward residuals.
It also tests extreme scales, deliberately inconsistent unused triangles,
and some ill-conditioned matrices. See the detailed validation record.
