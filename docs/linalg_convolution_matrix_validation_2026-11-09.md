# Convolution matrix constructor and FFT products — validation (2026-11-09)

New Fortran source:

- `src/linalg_convolution_matrix.f90`
- `src/linalg_matmul_convolution_matrix.f90`
- `test/test_linalg_convolution_matrix.f90`
- `example/linalg_convolution_matrix_example.f90`
- `tools/linalg_reference/driver_convolution_matrix.f90`
- `tools/linalg_reference/run_convolution_matrix_checks.py`

Environment: GNU Fortran 14.2.0, SciPy 1.17.0, NumPy 2.3.5.

## Native regression testing

Build the following dependencies in order:

```
src/kinds.f90
src/linalg.f90
src/linalg_dense.f90
src/linalg_toeplitz.f90
src/signal_windows_extra.f90
src/signal_dpss.f90
src/signal_spectral.f90
src/linalg_matmul_toeplitz.f90
src/linalg_circulant.f90
src/linalg_hankel.f90
src/linalg_matmul_hankel.f90
src/linalg_convolution_matrix.f90
src/linalg_matmul_convolution_matrix.f90
```

All tests below passed with both
`-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`
and `-std=f2018 -O3 -march=native`:

- `test_linalg_convolution_matrix`
- `test_linalg_matmul_toeplitz`
- `test_linalg_toeplitz`
- `test_linalg_circulant`
- `test_linalg_hankel`

The new test exercises **3,990** `(m,n,mode)` combinations with
`m=1..38`, `n=1..35`, and `full`/`same`/`valid` modes. Each is tested
with real and complex coefficients and both vector and multiple-RHS
multiplications. It also checks mixed real/complex coefficient/RHS
combinations; integer, real, and complex batch constructors; correct
cropping when `m>n`; default `full` mode; invalid lengths, modes,
shapes, and nonfinite inputs; factor scaling from `1e-80` to `1e80`;
zero RHS columns; and a 4096-input FFT case.

The executable example prints SciPy's expected `same`-mode vector:

```
 2.000  6.000  -1.000  -12.500  8.000
```

## Independent SciPy and NumPy comparison

From the repository root:

```bash
python tools/linalg_reference/run_convolution_matrix_checks.py \
    --mode checked --cases 1200 --seed 20261109
python tools/linalg_reference/run_convolution_matrix_checks.py \
    --mode release --cases 1200 --seed 20261110
```

The Python script compiles the **actual Fortran modules** into a stream
processing executable, generates random real or complex coefficient
vectors and multiple right-hand sides, and computes independently:

1. `scipy.linalg.convolution_matrix(a,n,mode)` for every matrix entry;
2. `numpy.convolve(a,x,mode)` for each RHS column;
3. `scipy.linalg.convolution_matrix(a,n,mode) @ x` as a cross-check of
   the Python references.

The Fortran constructor and FFT matrix product are each compared to the
corresponding reference. Inputs cover the full type matrix, modes,
both relative coefficient/input lengths, odd/even lengths, multiple
right-hand sides, sizes 1–85, and coefficient scales `1e-80`–`1e80`.

| Build | Cases | Maximum relative constructor difference | Maximum relative FFT product difference |
|---|---:|---:|---:|
| Checked (`-O0`, runtime checks and traps) | 1,200 / 1,200 | 0 | 2.05718e-14 |
| Optimized (`-O3 -march=native`) | 1,200 / 1,200 | 0 | 2.29739e-14 |
| **Combined** | **2,400 / 2,400** | **0** | **2.29739e-14** |

The cross-check scripts enforce a maximum relative numerical tolerance
of `1e-10` on FFT products. The batch-constructor overloads are tested
in native Fortran tests, not in the independent reference driver.

## Limitations

The full `fpm test` suite was not executed because FPM is not installed
in this environment. The focused tests above are not a substitute for
cross-platform CI. Arbitrary-rank SciPy coefficient broadcasting and
single-precision numeric constructors are not yet implemented. The FFT
products use the existing complex FFT backend and may perform more
work than SciPy's optimized `rfft` paths for real input.

Reference API: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.convolution_matrix.html
