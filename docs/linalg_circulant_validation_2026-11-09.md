# Circulant construction and FFT solves — validation (2026-11-09)

New source: `src/linalg_circulant.f90`; no new FPM dependencies.
GNU Fortran 14.2.0; SciPy 1.17.0; NumPy 2.3.5.

## Native Fortran builds

The source dependencies were compiled in the following order:

```bash
gfortran -std=f2018 -O0 -g -fcheck=all \
    -ffpe-trap=invalid,zero,overflow -c \
    src/kinds.f90 src/linalg.f90 src/signal_windows_extra.f90 \
    src/signal_dpss.f90 src/signal_spectral.f90 src/linalg_circulant.f90

gfortran -std=f2018 -O0 -g -fcheck=all \
    -ffpe-trap=invalid,zero,overflow -o test_circulant \
    kinds.o linalg.o signal_windows_extra.o signal_dpss.o \
    signal_spectral.o linalg_circulant.o test/test_linalg_circulant.f90
./test_circulant
```

Result: `test_linalg_circulant: PASS`. The same tests also passed
when compiled using `-O3 -march=native` instead of the debug flags.

Native tests cover: SciPy's documented three-element circulant system,
constructor orientation (including complex coefficients), singular
minimum-norm least squares, vectors and multiple RHS, all four
real/complex coefficient/RHS type combinations, rank-two axis positions,
empty matrices and zero-column RHS, negative tolerance, invalid policy,
invalid axes and output shapes, nonfinite inputs, deterministic system
sizes 1–140 (covering radix-2, mixed-radix and Bluestein FFT paths), and
scales 1e-80, 1, and 1e80.

An independent regression run of the previous Toeplitz test passed:
`test_linalg_toeplitz: PASS` in the checked build.

## Independent SciPy/NumPy comparison

From the project root:

```bash
python tools/linalg_reference/run_circulant_checks.py \
    --cases 1200 --mode checked --seed 20261109
python tools/linalg_reference/run_circulant_checks.py \
    --cases 1200 --mode release --seed 20261110
```

The script compiles the actual Fortran code to an executable and compares
its outputs to `scipy.linalg.solve_circulant`, with `scipy.linalg.circulant`
used to compute independent matrix residuals. The tests include 21 sizes
(from 1 to 128, including prime lengths), 1/2/4 right-hand sides,
real/complex/mixed types, full-rank and singular coefficients, explicit and
default tolerances, `'raise'` and `'lstsq'` policies, and scales between
1e-80 and 1e80.

- Checked build: **1,200/1,200 passed**, 148 expected singular detections;
  maximum relative difference **1.646e-14**; maximum normalized residual
  for the nonsingular cases **6.950e-16**.
- Optimized `-O3 -march=native` build: **1,200/1,200 passed**, 148 expected
  singular detections; maximum relative difference **1.513e-14**;
  maximum normalized nonsingular residual **6.253e-16**.
- Overall: **2,400/2,400 successful SciPy comparisons** and **296 expected
  singular detections**.

For near-singular tests with multiple right-hand sides, the harness uses
SciPy's supported trailing-vector axis (`b.T`, `baxis=-1`, `outaxis=-1`)
so its internal singular-mask broadcasting works correctly in SciPy 1.17.0.
This is mathematically equivalent to the corresponding `b(n,nrhs)` call.

The executable example also passed and printed the reference solutions
`[0.75, -0.25, 0.25]` and `[0.25, 1.25, 2.25, 1.25]`.

The full repository `fpm test` command was **not run** because FPM is
not installed in this execution environment. The focused builds and tests
are not a substitute for a complete FPM CI run.
