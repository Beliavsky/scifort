# Toeplitz construction and solves — validation (2026-11-09)

New source: `src/linalg_toeplitz.f90`; no new FPM dependencies.
GNU Fortran 14.2.0; SciPy 1.17.0; NumPy 2.3.5.

## Native checked build

```bash
gfortran -std=f2018 -O0 -g -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -c \
  src/kinds.f90 src/linalg.f90 src/linalg_dense.f90 src/linalg_toeplitz.f90
gfortran -std=f2018 -O0 -g -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -o test_toeplitz \
  kinds.o linalg.o linalg_dense.o linalg_toeplitz.o \
  test/test_linalg_toeplitz.f90
./test_toeplitz
```

Result: `test_linalg_toeplitz: PASS`. The native test checks real/complex
square Toeplitz constructors, rectangular constructors, complex-conjugate
default row, explicit row ignoring `r(1)`, non-Hermitian systems, vector and
multiple RHS, all mixed-type RHS overloads, deterministic orders 1–40,
real and complex zero-leading-pivot fallback, exactly singular cases,
invalid shapes, nonfinite inputs, zero-size systems and RHS batches, and
scales from 1e-80 to 1e80.

## Independent SciPy and NumPy comparisons

```bash
python tools/linalg_reference/run_toeplitz_checks.py --cases 1600 --mode checked --seed 20261109
python tools/linalg_reference/run_toeplitz_checks.py --cases 1600 --mode release --seed 20261110
```

- Checked: **1,600/1,600 passed**, 64 dense-LU fallbacks;
  maximum relative solution difference to NumPy dense solve **1.994e-14**;
  maximum normalized residual **1.885e-16**.
- Optimized `-O3 -march=native`: **1,600/1,600 passed**,
  64 dense-LU fallbacks; maximum relative solution difference to NumPy
  **1.417e-14**, maximum normalized residual **1.250e-16**.
- Total: **3,200 successful reference comparisons**. Nonfallback cases
  were also independently checked against SciPy's `solve_toeplitz`.
- Matrices include real and complex non-Hermitian/Hermitian Toeplitz,
  explicit/default first rows, mixed real/complex RHS, orders 1–48,
  1/2/4 right-hand sides, deliberately zero first pivots, coefficient
  scales 1e-80 through 1e80. Cases were checked against NumPy's pivoted
  dense solver; standard cases also against SciPy's Toeplitz solver.

These focused tests do **not** replace `fpm test`; FPM is not installed
in this environment. No new external implementation code was copied.
