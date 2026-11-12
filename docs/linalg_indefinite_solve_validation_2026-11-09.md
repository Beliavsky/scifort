# LDL-based symmetric/Hermitian linear solves — validation (2026-11-09)

New module: `src/linalg_indefinite_solve.f90`, with the existing
`src/linalg_ldl.f90`. Independent numerical oracles:
`scipy.linalg.solve(..., assume_a='sym'/'her')` and NumPy's general dense
`linalg.solve` on the reconstructed selected triangle. SciPy version 1.17.0;
GNU Fortran 14.2.0.

## Completed checks

- Native GNU Fortran test `test/test_linalg_indefinite_solve.f90` passes
  compiled with `-std=f2018 -O0 -fcheck=all
  -ffpe-trap=invalid,zero,overflow`. It exercises 144 deterministic base matrices of
  orders 1–24, each tested for real symmetric, complex Hermitian, and
  complex symmetric systems in both orientations,
  factor reuse, direct solvers and vector/multiple RHS overloads. Additional
  deterministic tests exercise exact 2x2 pivots, remote row swaps, invalid
  shapes and permutations, malformed/near-singular blocks, singular
  systems, nonfinite RHS and zero-size edge cases.
- Existing `test/test_linalg_ldl.f90` passes after adding the new module.
- `python tools/linalg_reference/run_indefinite_solve_checks.py --cases 900
  --mode checked --seed 20261109` passed **900/900** cases.
- `python tools/linalg_reference/run_indefinite_solve_checks.py --cases 1800
  --mode release --seed 20261110` passed **1800/1800** cases; maximum
  combined relative solution/residual error **1.34e-13**.
- `python tools/linalg_reference/run_indefinite_solve_checks.py --cases 1800
  --mode checked --seed 20261111` passed **1800/1800** cases; maximum
  combined relative solution/residual error **3.04e-13**.

**4,500/4,500 independent SciPy/NumPy comparisons passed** across the three
runs. The SciPy and NumPy reference solvers are never linked to or used in
the Fortran numerical implementation. The reference harness sends each
case through both `solve_indefinite` and separately `ldl` + `ldl_solve`,
checking both against the independent oracles, the linear equation residual,
and each other. The unused input triangle contains deliberately conflicting
values so accidental two-triangle reads would be detectable. Test scales
include 1e-80, 1, and 1e80. Orders include 0, 1, 2, 3, 4, 5, 8, 12, 20,
and 32, with real symmetric, complex Hermitian and complex symmetric cases,
multiple right-hand sides, deliberate 2x2 blocks, and upper/lower selection.

To reproduce the independent comparisons:

```bash
python tools/linalg_reference/run_indefinite_solve_checks.py --cases 900 --mode checked --seed 20261109
python tools/linalg_reference/run_indefinite_solve_checks.py --cases 1800 --mode release --seed 20261110
python tools/linalg_reference/run_indefinite_solve_checks.py --cases 1800 --mode checked --seed 20261111
```

These focused checks do not replace full project-level `fpm test` (FPM is
not installed here). No new runtime dependencies, copied third-party code,
or license terms were introduced.
