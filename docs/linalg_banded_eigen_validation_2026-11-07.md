# 2026-11-07 symmetric and Hermitian band eigenvalue validation

The new functions are tested against SciPy 1.17.0 with Python 3 and
GNU Fortran in both checked and optimized modes.

Main reference test:

```sh
python tools/linalg_reference/run_banded_eigen_checks.py \
  --mode both --cases 600 --max-order 96 --seed 1701
```

The initial completed run reported, for each build:

- 600/600 randomized SciPy comparisons passed.
- Maximum normalized eigenvalue difference: `8.533e-16`.
- Maximum normalized eigen-equation residual: `3.346e-15`.
- Maximum normalized eigenvector orthogonality error: `6.601e-16`.
- Native tests passed, including zero and invalid shapes, real and complex
  matrices, NaN rejection and scale extremes from `1e-150` to `1e150`.

As in SciPy, complex Hermitian eigenvalues are real and sorted ascending.
We do not require matching eigenvectors in repeated eigenspaces because the
basis vectors are not unique. Scaled residuals use each matrix's largest
entry as a scaling factor so the comparisons remain meaningful for
very small and large input matrices.

This is an eigenpair-accuracy check, not a benchmark of memory usage or
performance; the general-band implementation remains dense internally.

Extended validation with a second seed:

```sh
python tools/linalg_reference/run_banded_eigen_checks.py --mode both \
  --cases 320 --max-order 128 --seed 2604
```

- Both compiler modes passed 320/320 additional cases each.
- Maximum normalized eigenvalue difference: `1.050e-15`.
- Maximum normalized eigen-equation residual: `3.710e-15`.
- Maximum normalized orthogonality error: `6.060e-16`.

Including the first 120-case smoke test in both compiler modes, the new
solver passed 2,080 Fortran-to-SciPy comparisons before archive extraction.

Neighboring routines passed in checked and optimized builds: 180 Hessenberg,
180 matrix-balancing, 140 SVD-condition, 60 Schur matrix-equation and
120 Procrustes comparisons, a total of 680 regression comparisons.
The complete `fpm test` suite could not be run because `fpm` is not
installed in the current environment.
