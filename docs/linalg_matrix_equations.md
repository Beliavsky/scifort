# Dense Sylvester and Lyapunov equations

The module `scifort_linalg_matrix_equations` implements real and complex
binary64 equivalents of three `scipy.linalg` matrix-equation routines:

| Generic interface | Equation | SciPy reference |
|---|---|---|
| `solve_sylvester` | `A X + X B = Q` | `scipy.linalg.solve_sylvester` |
| `solve_continuous_lyapunov` | `A X + X A^H = Q` | `scipy.linalg.solve_continuous_lyapunov` |
| `solve_discrete_lyapunov` | `X - A X A^H = Q` | `scipy.linalg.solve_discrete_lyapunov` |

Here `A^H` denotes the conjugate transpose (ordinary transpose for real data).
The interface uses caller-allocated result arrays and an integer status.

## Algorithms (November 7, 2026 update)

All three generics now accept a final optional character argument `method`:

- `method='auto'` (default): original direct Kronecker elimination for at most
  `matrix_equation_max_unknowns=400` scalar unknowns; otherwise Schur if each
  matrix order is at most `matrix_equation_schur_max_order=40`.
- `method='direct'`: force the previous, pivoted dense Kronecker solver and
  retain its 400-unknown limit.
- `method='schur'`: use SciFort's own complex Schur reduction and triangular
  recurrences regardless of problem size (up to order 40 per coefficient
  matrix), including real input by taking the real part of the computed result.

```fortran
use scifort_kinds, only : dp
use scifort_linalg_matrix_equations, only : solve_sylvester, &
    solve_continuous_lyapunov, solve_discrete_lyapunov

! real(dp) :: a(24,24), b(24,24), q(24,24), x(24,24)
! integer :: status
call solve_sylvester(a, b, q, x, status, method='schur')
call solve_continuous_lyapunov(a, q, x, status) ! automatic Schur for 24x24
```

For `solve_sylvester`, compute complex Schur forms `A=U T U^H` and
`B=V S V^H`. After transforming `C=U^H Q V`, solve the triangular equation
`T Y + Y S=C` by backward row / forward column substitution, then recover
`X=U Y V^H` (the Bartels–Stewart approach).

Both Lyapunov routines use **one** Schur decomposition, `A=U T U^H`, then
compute `C=U^H Q U`. For continuous Lyapunov solve `T Y+Y T^H=C`; for
discrete Lyapunov solve `Y-T Y T^H=C`, both by backward row and column
recurrences. The discrete recurrence includes a double off-diagonal sum.
These approaches avoid the original `(m*n) x (m*n)` dense Kronecker matrix.
An additional normalized backward-residual check guards against reporting a
numerically unreliable Schur solution as successful.

The internal Schur implementation uses the dependency-free routines in
`scifort_linalg_matrix_functions`; **no BLAS/LAPACK or SciPy code is bundled**.
The existing direct method and original call signatures remain available.

## Status and numerical caveats

Status constants are exposed by `scifort_linalg`:

- `linalg_status_success` (0): solution computed.
- `linalg_status_invalid_shape` (1): invalid array shapes, `method`, or
  order too large for the requested algorithm.
- `linalg_status_no_convergence` (4): nonfinite input or output, Schur
  factorization failure, or excessive normalized backward residual.
- `linalg_status_singular` (5): a scalar triangular denominator is zero
  to the method's working precision, or a direct Kronecker pivot fails.

The solution is unique only when the associated Sylvester or Lyapunov
operator is nonsingular. Stability of `A` is sufficient for uniqueness in the
standard stable Lyapunov settings, but is **not required** for every valid
mathematical instance.

The native Schur QR algorithm is intended for **small and moderate matrices**;
its order limit is 40, and its performance is not competitive with optimized
LAPACK for large systems. Singular/near-resonant eigenvalue sums or products,
highly nonnormal matrices, and extreme scaling can defeat the triangular
recurrences or cause loss of forward accuracy. A low backward residual is not
an estimate of conditioning. Optional Schur input balancing, iterative
refinement, batched arrays, and higher-order LAPACK-backed solving remain
future work.

## Validation and examples

- Original test: `test/test_linalg_matrix_equations.f90`
- New method-selection, large-system and error-path tests:
  `test/test_linalg_matrix_equations_schur.f90`
- SciPy reference harness: `tools/linalg_reference/run_matrix_equations_schur_checks.py`
- Original direct reference harness:
  `tools/linalg_reference/compare_matrix_equations_scipy.py`
- Runnable example: `example/linalg_matrix_equations_schur_example.f90`

```bash
python tools/linalg_reference/run_native_linalg.py --mode both
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 72 --max-order 24
python tools/linalg_reference/run_matrix_equations_schur_checks.py --mode both --cases 40 --max-order 40 --large
```

Requires GNU Fortran, NumPy and SciPy for the reference tests; no SciPy or
LAPACK dependency is needed to compile the Fortran matrix-equation module.
See `docs/linalg_matrix_equations_schur_validation_2026-11-07.md` for numerical
results and testing limitations.
