# LDL decomposition validation -- 2026-11-09

Implementation: `src/linalg_ldl.f90`. Independent reference:
`scipy.linalg.ldl` from SciPy 1.17.0. Compiler: GNU Fortran 14.2.0.

## Completed validation

- Native tests `test/test_linalg_ldl.f90` passed with `-fcheck=all`.
  Covers real symmetric, complex Hermitian and complex symmetric inputs;
  lower and upper modes; nonadjacent permutations and 2x2 pivots;
  invalid shapes, nonfinite entries, singular and empty matrices; and
  deterministic randomized matrix families with orders 1 through 15.
- `python tools/linalg_reference/run_ldl_checks.py --cases 1200 --mode checked`
  completed **1,200/1,200 comparisons** with bounds checks and floating-point
  traps enabled. Maximum relative infinity-norm reconstruction residual:
  **1.84e-15**.
- `python tools/linalg_reference/run_ldl_checks.py --cases 1200 --mode release`
  completed **1,200/1,200 comparisons** at `-O3` with the same maximum error.
- Both suites compare each input against `scipy.linalg.ldl`, verify SciPy's
  reference reconstruction, and independently check the SciFort identity,
  one-based permutation, unit triangular structure, and block-diagonal shape.
  Cases include orders 0, 1, 2, 3, 4, 5, 8, 12, 20, and 32; symmetric,
  Hermitian, and complex symmetric matrices; positive-shifted, indefinite,
  zero-diagonal, zero, rank-deficient, and scaled matrices (1e-80 and 1e80).

LDL factorizations are nonunique, so the test **does not require identical
pivots or factors** to SciPy. It independently checks the resulting
reconstruction identities instead.

## Reproduce

```bash
python tools/linalg_reference/run_ldl_checks.py --cases 1200 --mode checked
python tools/linalg_reference/run_ldl_checks.py --cases 1200 --mode release
```

The full repository-wide `fpm test` was not run because FPM was not installed
in this environment. These focused tests do not replace full FPM integration.
No third-party implementation source, runtime dependency, or new license was
added.
