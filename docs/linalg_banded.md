# Banded dense linear systems (`scipy.linalg`)

`scifort_linalg_banded` implements the following six routines, including
four `scipy.linalg` analogues and two LAPACK-like LU reuse interfaces,
with real and complex binary64 overloads:

| Fortran generic | Purpose | Inputs |
| --- | --- | --- |
| `solve_banded` | General banded Gaussian elimination with local partial pivoting | Packed general band, vector/matrix RHS |
| `solveh_banded` | Symmetric/Hermitian positive-definite banded solve | Packed upper/lower triangle, vector/matrix RHS |
| `lu_factor_banded` | General banded LU with local partial pivoting (similar to GBTRF) | Packed general band; expanded LU band and pivots |
| `lu_solve_banded` | Solve with reusable LU (similar to GBTRS) | LU factors; vector/matrix RHS; optional N/T/C mode |
| `cholesky_banded` | Packed banded Cholesky factorization | Packed upper/lower triangle |
| `cho_solve_banded` | Solve one or more right-hand sides using the packed factor | Cholesky factor, vector/matrix RHS |

Use the module directly, e.g. `use scifort_linalg_banded, only: solve_banded`.
All input and result arrays are caller-allocated except for internal workspace.
Inputs are not modified. This Fortran interface is modeled after SciPy's
*storage and mathematical behavior*, not its function-call argument order:

```fortran
use scifort_kinds, only : dp
use scifort_linalg_banded, only : solve_banded, cholesky_banded, cho_solve_banded
real(dp) :: ab(3,5), b(5), x(5)
integer :: status

ab = 0.0_dp
ab(2,:)   = 4.0_dp
ab(1,2:5) = -1.0_dp
ab(3,1:4) = -1.0_dp
b = [3.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]
call solve_banded(ab, 1, 1, b, x, status) ! x = [1,1,1,1,1]

! Upper-packed SPD storage has 2 rows (one superdiagonal), not 3.
```

For the last operation, prepare the correctly aligned symmetric upper band:

```fortran
real(dp) :: spd_upper(2,5), cb(2,5), b(5), x(5)
integer :: status
spd_upper = 0.0_dp
spd_upper(2,:)   = 4.0_dp
spd_upper(1,2:5) = -1.0_dp
call cholesky_banded(spd_upper, cb, status, lower=.false.)
if (status == 0) call cho_solve_banded(cb, b, x, status, lower=.false.)
```

See the complete executable example at `example/linalg_banded_example.f90`.

## Packed general storage

`solve_banded(ab, lower, upper, b, x, status)` requires
`shape(ab) == [lower+upper+1, n]`, `shape(b) == [n]` or `[n,nrhs]`, and
`shape(x) == shape(b)`. Indexing is **one-based**:

\[
  \mathrm{ab}(\mathrm{upper}+1+i-j,j)=A(i,j),\quad
  \max(1,j-\mathrm{upper})\le i\le\min(n,j+\mathrm{lower}).
\]

The unused corner cells may contain any *finite* data and are ignored.
The lower and upper counts are nonnegative and less than `n`.
There is no fixed order limit. Use this for *nonsymmetric* and
non-positive-definite banded matrices as well as SPD ones.

## Reusable general LU factorization and solve

For repeated general (not necessarily positive-definite) solves, factor once:

```fortran
use scifort_linalg_banded, only : lu_factor_banded, lu_solve_banded
real(dp) :: ab(3,5), lu(4,5), b(5), x(5)
integer :: ipiv(5), status

! Fill ab with the general tridiagonal matrix in the earlier example.
call lu_factor_banded(ab, 1, 1, lu, ipiv, status)
if (status == 0) then
    call lu_solve_banded(lu, 1, 1, ipiv, b, x, status)
    ! Reuse factors for another RHS, a matrix of RHSs, or a transpose solve.
    call lu_solve_banded(lu, 1, 1, ipiv, b, x, status, trans='T')
end if
```

`lu_factor_banded(ab, lower, upper, lu, ipiv, status)` requires
`shape(ab)=[lower+upper+1,n]`, `shape(lu)=[upper+2*lower+1,n]`,
and `shape(ipiv)=[n]`. Output `ipiv(k)` is the **one-based** pivot row
chosen at elimination step `k`; do not convert it to zero-based indices.

With one-based storage and `offset=upper+lower+1`, the upper entries of `U`
and the lower multipliers are stored together as `lu(offset+i-j,j)`.
The first `lower` rows of `lu` provide room for upper-band fill from row
exchanges; no dense matrix is constructed. This is analogous to the expanded
storage used by LAPACK `GBTRF`, but the implementation is independent and
has a distinct Fortran calling convention. Never interpret `lu` as the original
`ab` or pass it directly into `solve_banded`.

`lu_solve_banded(lu, lower, upper, ipiv, b, x, status, trans)` accepts
vectors or matrices of right-hand sides, with `shape(x)=shape(b)`.
The final optional character is `'N'` (default, solve `A*x=b`), `'T'`
(solve `transpose(A)*x=b`) or `'C'` (solve `transpose(conjg(A))*x=b`).
For real matrices `'C'` and `'T'` are equivalent. Factorization and solve
inputs are not modified; the same `lu` and `ipiv` can be used repeatedly.

The existing `solve_banded` calls these same factorization/solve kernels
internally, and its calling convention has not changed. `lu_factor_banded`
uses `O(n*(upper+2*lower+1))` factor workspace while each repeated solve
requires only `O(n*nrhs)` right-hand-side storage and
`O(n*(upper+lower)*nrhs)` operations. The workspace is caller-managed;
subsequent solves avoid refactorizing the coefficient matrix.

On failure to factor, `lu` and `ipiv` are reset to zero; a partial
factorization must not be reused. Invalid pivot arrays, shapes or transpose
flags yield `linalg_status_invalid_shape`. A valid zero pivot yields
`linalg_status_singular`; nonfinite input or factor values yield
`linalg_status_no_convergence`. A solution is zeroed on failure.

## Packed Hermitian/symmetric storage

Let `bw = size(ab,1)-1` denote the half-bandwidth. For `lower=.true.`,
`ab(1+i-j,j) = A(i,j)` when `j <= i <= min(n,j+bw)`. For the default
`lower=.false.`, `ab(bw+1+i-j,j) = A(i,j)` when
`max(1,j-bw) <= i <= j`. Only one triangle is supplied, and complex
inputs are interpreted as Hermitian. Complex diagonal entries must be
real to roundoff; the lower triangular factor has real positive diagonal.

For repeated right-hand sides, factor once:

```fortran
call cholesky_banded(ab, cb, status, lower=.true.)
if (status == 0) then
    call cho_solve_banded(cb, b1, x1, status, lower=.true.)
    call cho_solve_banded(cb, b2, x2, status, lower=.true.)
end if
```

The one-step `solveh_banded(ab, b, x, status, lower=.true.)` is equivalent
to factoring and solving; it does not mutate `ab`.

## Numerical algorithm and complexity

- General systems: reusable locally partial-pivoted Gaussian elimination, searching
  at most `lower` rows below each pivot. Pivoting can introduce at most
  `lower` extra superdiagonals; the work array keeps only this *extended*
  band. Complexity is `O(n*lower*(upper+lower+1))` operations and
  `O(n*(upper+2*lower+1) + n*nrhs)` memory. Factorization costs
  `O(n*lower*(upper+lower+1))`; each subsequent solve costs
  `O(n*(upper+2*lower+1)*nrhs)`.
- Positive-definite systems: a lower-banded Cholesky factor built by
  subtracting previously computed products, followed by banded triangular
  solves. Complexity is `O(n*bw**2 + n*bw*nrhs)` and workspace
  `O(n*(bw+1) + n*nrhs)`.
- Both paths work for lower/upper bandwidth zero (diagonal systems) and
  `n=1`. No dense `n`-by-`n` matrix is assembled, unless the requested
  bandwidth itself is comparable to `n`.

These are independent implementations of standard elimination and Cholesky
formulas, not translations of SciPy/LAPACK kernels. SciPy is used only to
compute independent numerical test values.

## Status and limitations

`status` uses the constants from `scifort_linalg`:

- `linalg_status_success` (0): completed solution/factorization.
- `linalg_status_invalid_shape` (1): wrong dimensions or invalid bandwidth.
- `linalg_status_not_positive_definite` (2): SPD/Hermitian factorization
  cannot proceed, including complex diagonal entries outside the real
  tolerance, invalid values or zero/negative pivots.
- `linalg_status_no_convergence` (4): nonfinite input/arithmetic encountered
  in a general solve (or a nonfinite factor result).
- `linalg_status_singular` (5): exactly zero pivot for general elimination
  or invalid/zero Cholesky pivot supplied to `cho_solve_banded`.

Returned solutions and factors are zero-initialized on failure. There is no
iterative refinement, equilibration, reciprocal-condition estimate or sparse
fill-reducing reordering. Extremely ill-conditioned systems can have large
forward error even when the backward residual is small; pivoting is confined
to the band (as in standard banded LU). Do not call `solveh_banded` on merely
symmetric indefinite matrices: positive definiteness is essential.

## Validation

The standalone test script builds GNU Fortran in checked and `-O3` modes
and compares against SciPy on real and complex matrices with random
bandwidths, multiple right-hand sides, and forced row exchanges. Native tests
also include an order-50,000 tridiagonal matrix, invalid shapes, singular
systems, nonpositive definite matrices, and Cholesky factor reuse:

```bash
python tools/linalg_reference/run_banded_checks.py --mode both --cases 750 --max-order 70 --seed 109 --stress
python tools/linalg_reference/run_banded_checks.py --mode both --cases 650 --max-order 120 --seed 991 --stress
python tools/linalg_reference/run_banded_lu_checks.py --mode both --cases 600 --max-order 85 --seed 171 --stress
```

See `docs/linalg_banded_validation_2026-11-07.md` and
`docs/linalg_banded_lu_validation_2026-11-07.md` for measured results.
The Python/SciPy comparison tests are not runtime dependencies.
