# Generalized nonsymmetric eigenproblems — validation (2026-11-08)

New module: `src/linalg_generalized_eigen.f90`. Independent reference:
SciPy `scipy.linalg.eig(A,B,left=True,right=True,homogeneous_eigvals=True)`
with OpenBLAS restricted to one thread. Fortran source compiled with GNU
Fortran and no new BLAS/LAPACK dependencies.

## Completed tests

- Main test: `--cases 400 --max-order 40 --seed 91283`, **400/400 in each**
  of checked (`-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow`) and
  optimized (`-O3`) configurations. Maximum homogeneous/projective eigenvalue
  difference: `2.571e-14` checked, `2.554e-14` optimized; maximum normalized
  right residual `9.547e-15`; left residual `6.863e-15`.
- Additional seed: `--cases 600 --max-order 40 --seed 99181`, **600/600**
  in the checked build. Maximum projective spectral difference `7.747e-14`,
  maximum right residual `1.003e-14`, left residual `8.599e-15`. Its combined
  run timed out before finishing the optimized build, which is **not** counted.
- Independent optimized subset: `--cases 250 --max-order 40 --seed 99181`,
  **250/250**, with maximum projective spectral difference `3.510e-14`.
- Existing ordinary nonsymmetric `eig` native checks and **80/80 SciPy
  comparisons in each build** passed (`--seed 1947 --max-order 18`).
- Existing QZ native checks and **16/16 SciPy comparisons in each build**
  passed (`--method iterative --max-order 8 --seed 31209`).
- New generalized eigenproblem native tests passed with runtime checks and
  overflow/invalid/division-by-zero traps. These include an infinite root,
  very large/small physical eigenvalue ratios, the zero-order case, all four
  input-kind overloads, singular pencils, and invalid argument handling.

The new test harness checks the spectrum via optimal matching of projective
homogeneous coordinates (including `beta=0`), normalized right and left
pencil residuals, normalization of eigenvectors, and agreement between `eig`
and `eigvals`. It covers upper/lower triangular matrix pairs, repeated
finite eigenvalues, singular-but-regular pencils, random dense pairs, and
extreme scale ratios. It deliberately does not compare eigenvector entries
or order, because these are not unique.

A densely populated **lower-triangular** random pencil exposed an avoidable
forward spectral error of about `3e-3` when processed through iterative QZ.
The new direct joint-triangular path retains exact diagonal eigenvalues and
reduces that case's projective error to near roundoff. A pencil with
`A=2*B` exposed an overly strict near-zero residual scale; comparing the
homogeneous-pencil residual against its coefficient scale corrected that
numerical false failure without relaxing the absolute residual tolerance.

**Limits:** QZ is currently capped at order 40. Defective, severely
nonnormal, or irregular pencils can return numerical failure or large
forward eigenvalue differences. Passing the residual checks does not guarantee
forward eigenvalue conditioning. The full project FPM test suite was not run
because FPM was unavailable in the execution environment.

Reproduce the principal checks:

```bash
python tools/linalg_reference/run_generalized_eigen_checks.py \
    --mode both --cases 400 --max-order 40 --seed 91283
```
