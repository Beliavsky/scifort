# Symmetric / Hermitian indefinite linear solves (`scipy.linalg.solve`)

`scifort_linalg_indefinite_solve` implements the **`assume_a='sym'` and
`assume_a='her'` specializations** of `scipy.linalg.solve`, not its general,
banded, triangular, or positive-definite solver modes. It uses SciFort's
Bunch–Kaufman `ldl` decomposition and explicit triangular/block solves for
real symmetric, complex symmetric, and complex Hermitian systems. Matrices
need not be positive definite, but **must be nonsingular** to solve a system.

The module exports two generic Fortran subroutines, both with real and complex
binary64 overloads for vectors `b(n)` and dense RHS matrices `b(n,nrhs)`:

| Generic | Arguments | Meaning |
| --- | --- | --- |
| `solve_indefinite` | `(a, b, x, status, lower, hermitian)` | Factor `A`, then solve `A x = b` |
| `ldl_solve` | `(lu, d, perm, b, x, status, lower, hermitian)` | Reuse existing `ldl` factors for additional RHS |

For real arrays, omit `hermitian`; both real and complex generics accept
`lower` (default **`.true.`**). The complex generics have an optional
`hermitian` (default **`.true.`**); set `hermitian=.false.` for complex
symmetric `A` and transpose instead of conjugate transpose. The default
`lower` choice differs from `scipy.linalg.solve` (which defaults to upper).
Always pass the **same** `lower` and `hermitian` choices to `ldl` and
`ldl_solve`. All dummy arguments have explicit `INTENT` and `status`
returns `linalg_status_*` values.

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_ldl, only : ldl
use scifort_linalg_indefinite_solve, only : ldl_solve, solve_indefinite
real(dp) :: a(n,n), lu(n,n), d(n,n), b(n,2), x(n,2)
integer :: perm(n), status

! Fill symmetric a and two right-hand sides b.
call solve_indefinite(a,b,x,status)

! For repeated calls, factor once; solve as often as needed.
call ldl(a,lu,d,perm,status,lower=.true.)
if (status == linalg_status_success) then
    call ldl_solve(lu,d,perm,b,x,status,lower=.true.)
end if
```

Complex Hermitian systems:

```fortran
complex(dp) :: a(n,n), b(n), x(n)
integer :: status
call solve_indefinite(a,b,x,status,lower=.false.)
```

Complex symmetric (not Hermitian) systems:

```fortran
call solve_indefinite(a,b,x,status,hermitian=.false.)
```

The LDL decomposition convention is `A = lu*D*lu^H` (Hermitian) or
`A = lu*D*lu^T` (symmetric), where `lu(perm,:)` is unit lower or unit upper
triangular and `perm` is a **one-based permutation**. The solver first
permutes `b`, solves through that unit triangular factor, solves each 1x1 or
2x2 diagonal block of `D`, solves through the factor's adjoint/transpose,
and undoes the permutation. It never treats the 2x2 blocks as unrelated
scalar diagonal entries. Both generic interfaces preserve caller input
arrays; `x` is set to zero if `status` indicates failure.

`solve_indefinite` reads only the selected input triangle, just as `ldl`
does. Its unused triangle may contain unrelated data; it is not checked
for finiteness. Reused factors must be finite and have valid, consistent
unit triangular and block diagonal structure. Nonfinite right-hand sides
are rejected. `n=0` and zero RHS columns are supported; factorization of
singular matrices may succeed, but solving them returns singular status.

## Status and limitations

| Status | Meaning |
| --- | --- |
| `linalg_status_success` (`0`) | A solution was computed |
| `linalg_status_invalid_shape` (`1`) | Inconsistent array sizes, invalid permutation, or malformed factor structure |
| `linalg_status_no_convergence` (`4`) | Nonfinite data or a nonfinite numerical result |
| `linalg_status_singular` (`5`) | Zero 1x1 pivot or (near-)singular 2x2 pivot prevents a solve |

For complex Hermitian input, `ldl` rejects non-real diagonal entries with
`invalid_shape`. The solver uses straightforward binary64 arithmetic without
iterative refinement, equilibration, a condition-number estimate, or a
backward error bound. Extremely ill-conditioned or poorly scaled matrices
may yield inaccurate solutions, even when the residual appears small.
The real overloads use temporary complex work arrays to share the same
substitution implementation; this is correct but not yet performance-tuned
for large real systems. The unblocked factorization is `O(n^3)` with
`O(n^2)` workspace; each solve is `O(n^2 nrhs)`.

See the [complete example](../example/linalg_indefinite_solve_example.f90),
[native tests](../test/test_linalg_indefinite_solve.f90),
[streaming SciPy comparison driver](../tools/linalg_reference/driver_indefinite_solve.f90),
[comparison script](../tools/linalg_reference/run_indefinite_solve_checks.py),
and [validation report](linalg_indefinite_solve_validation_2026-11-09.md).
