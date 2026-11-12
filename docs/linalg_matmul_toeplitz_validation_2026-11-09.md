# Toeplitz multiplication — validation (2026-11-09)

New implementation: `src/linalg_matmul_toeplitz.f90`.
No additional third-party dependency. GNU Fortran 14.2.0, SciPy 1.17.0,
NumPy 2.3.5.

## Native Fortran compilation and tests

Compile dependencies in the project root, placing the objects and module
files in a build directory. The tested source compilation order was:

```
src/kinds.f90
src/linalg.f90
src/linalg_dense.f90
src/linalg_toeplitz.f90
src/signal_windows_extra.f90
src/signal_dpss.f90
src/signal_spectral.f90
src/linalg_matmul_toeplitz.f90
```

Compile with `gfortran -std=f2018 -O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow` and link
`test/test_linalg_matmul_toeplitz.f90` against the compiled objects.
Repeat with `-std=f2018 -O3 -march=native`. Both executable runs reported:

`test_linalg_matmul_toeplitz: PASS`

The native tests cover a documented rectangular Toeplitz matrix,
complex explicit and implicit first rows, all eight mixed-kind/rank
combinations, multiple right-hand sides, incorrect shapes, ignored `r(1)`,
nonfinite checks, empty dimensions, zero-column multipliers, complex
orientation, deterministic test sizes 1–128, a 4,096-element structured product
with a prime-length (8,191) FFT embedding, and coefficient scales from
`1e-80` through `1e80`.

The existing `test/test_linalg_toeplitz.f90` also passed in a checked GNU
Fortran build. The new example compiled and ran successfully.

## Reproducible independent Python / SciPy comparisons

From the project root:

```bash
python tools/linalg_reference/run_matmul_toeplitz_checks.py \
  --mode checked --cases 1800 --seed 20261109
python tools/linalg_reference/run_matmul_toeplitz_checks.py \
  --mode release --cases 1800 --seed 20261110
```

Each invocation recompiles the **actual Fortran implementation and a
Fortran test driver**. Every Fortran output is compared against both
`scipy.linalg.matmul_toeplitz` and the independent dense calculation
`scipy.linalg.toeplitz(c, r) @ x`.

| Build | Cases passed | Maximum relative difference vs. SciPy | Maximum relative difference vs. dense |
| --- | ---: | ---: | ---: |
| Checked `-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow` | 1800/1800 | 2.184e-14 | 2.171e-14 |
| Optimized `-O3 -march=native` | 1800/1800 | 1.917e-14 | 1.930e-14 |
| **Total** | **3600/3600** | **2.184e-14** | **2.171e-14** |

The comparisons include real/complex/mixed-type matrices and multipliers,
rectangular and square dimensions, explicit and conjugate-default rows,
large disagreements between the supplied `r(1)` and `c(1)`, one/two/four
right-hand sides, non-power-of-two FFT sizes including primes, zero
coefficients, and coefficient scales 1e-80 to 1e80. Tolerance: normalized
Frobenius difference <= 3e-11 to each independent reference.

`fpm test` was **not** run: FPM is not installed in this environment.
These focused compilations and numerical checks are not a substitute for
full project CI across compilers.
