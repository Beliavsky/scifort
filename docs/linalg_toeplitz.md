# Toeplitz matrices and linear solves (`scipy.linalg.toeplitz` / `solve_toeplitz`)

Module: `scifort_linalg_toeplitz` (MIT, originally implemented for SciFort).

## Fortran API

`toeplitz(c, t, status, r)` constructs a real or complex Toeplitz matrix:

- `c(:)` is the **first column**, of length `n`.
- `r(:)` is the **optional first row**, of length `m`. `r(1)` is ignored:
  the diagonal always comes from `c(1)`.
- Without `r`, the first row is `conjg(c)` (the ordinary identity for
  real arrays); a complex `c(1)` remains unchanged on the diagonal.
- The caller supplies output `t(n,m)`, or `t(n,n)` if `r` is absent, and
  `status`. Rectangular and empty constructors are supported.

`solve_toeplitz(c, b, x, status, r, check_finite, used_dense)` solves the
**square** Toeplitz system without explicitly constructing the full matrix
on the normal path. It supports real and complex `c`, vector `b(n)`, matrix
`b(n,nrhs)`, and mixed real/complex coefficient/RHS combinations (a complex
output is required if either input is complex). The output `x` must have the
same dimensions as `b`. The argument `r`, if supplied, must have length `n`.

- `check_finite` defaults to `.true.`. Nonfinite inputs return
  `linalg_status_no_convergence`. With `.false.`, behavior on nonfinite
  coefficients or RHS is unspecified, just as SciPy warns.
- `used_dense` is an optional logical output, a **SciFort extension**:
  it reports that the fast nonpivoting recurrence had to switch to pivoted
  dense LU. A zero leading diagonal can require this even when the complete
  Toeplitz matrix is nonsingular.
- `status` uses constants in `scifort_linalg`:
  `success=0`, `invalid_shape=1`, `no_convergence=4`, `singular=5`.
- The solution array is set to zero on a reported error. Empty systems and
  empty RHS batches are allowed.

**SciPy correspondence:** call SciPy with `solve_toeplitz((c,r),b)`, or
`solve_toeplitz(c,b)` when `r` is omitted. The Fortran interface uses explicit
`x` and integer `status` output arguments. SciPy 1.17+ batch-leading
dimensions, Python dtype coercion, and high-rank broadcasting are not
implemented; the Fortran core accepts one problem and vector/matrix RHS.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg, only : linalg_status_success
use scifort_linalg_toeplitz, only : solve_toeplitz, toeplitz
real(dp) :: c(4), r(4), b(4), x(4)
integer :: status
c = [1.0_dp,3.0_dp,6.0_dp,10.0_dp]
r = [1.0_dp,-1.0_dp,-2.0_dp,-3.0_dp]
b = [1.0_dp,2.0_dp,2.0_dp,5.0_dp]
call solve_toeplitz(c,b,x,status,r)
if (status /= linalg_status_success) error stop 'solve failed'
```

See the full [executable example](../example/linalg_toeplitz_example.f90).

## Algorithm and cost

The fast solver is an **independently written general Toeplitz bordering
recurrence**, mathematically equivalent in complexity to Levinson-style
methods. It does not require positive definiteness or Hermitian symmetry.

Let the leading `k`-by-`k` Toeplitz matrix be `T_k`, its next right-hand
column be `p_k` and next bottom row be `q_k^T`. Maintain
`a_k = T_k^{-1} p_k` and `d_k = T_k^{-T} q_k` (ordinary transpose, **not**
adjoint). The scalar Schur complement is

`S_k = c_0 - q_k^T a_k`.

When extending a solution from order `k` to `k+1`, the new final entry is

`x_{k+1} = (b_{k+1} - q_k^T x_{1:k}) / S_k`,

and the existing solution entries are updated by

`x_{1:k} <- x_{1:k} - a_k x_{k+1}`.

Toeplitz reversal symmetry gives the inverse first-column vectors from
`a_k` and `d_k`, allowing both bordering states to be advanced with O(k)
work. Thus a single RHS costs **O(n^2)**, with **O(n)** auxiliary storage;
for `nrhs` RHS the cost is O(n^2 (1 + nrhs)) and workspace O(n * nrhs).
Scaling the Toeplitz diagonals and RHS by a common coefficient norm avoids
many overflow/underflow failures from uniformly tiny or enormous matrices.

When the current leading Schur complement is too small in scaled binary64
arithmetic (or the fast solution becomes nonfinite), the implementation
uses SciFort's existing partial-pivoted `lu_factor` and `lu_solve` on the
constructed Toeplitz matrix. This fallback costs O(n^3) work and O(n^2)
memory, but handles invertible systems with zero leading principal pivots.
The threshold is `32*epsilon(1.0_dp)` after coefficient normalization.

## Accuracy limitations

Like SciPy's Levinson-Durbin solver, the fast recurrence can be unstable for
severely ill-conditioned or nonnormal Toeplitz systems. The fallback only
catches small leading pivots and nonfinite fast results; it is **not** a
condition estimate, componentwise backward-error check, or iterative
refinement. A successful status does not guarantee a small forward error.
Use a general pivoted factorization (or condition-aware high-precision method)
when the problem is highly ill-conditioned. The dense fallback is not a
superfast Toeplitz solver.

Public API definitions checked against SciPy documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_toeplitz.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.toeplitz.html.
The algorithm was implemented from elementary Toeplitz block identities;
no SciPy, LAPACK or other third-party implementation source was imported.

## Tests

- [Native Fortran regression](../test/test_linalg_toeplitz.f90)
- [Independent SciPy/NumPy reference tests](../tools/linalg_reference/run_toeplitz_checks.py)
- [Comparison driver](../tools/linalg_reference/driver_toeplitz.f90)
- [Validation results](linalg_toeplitz_validation_2026-11-09.md)
