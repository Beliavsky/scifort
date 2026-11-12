# Helmert/Hilbert constructor validation — 2026-11-09

## Environment

- GNU Fortran 14.2.0, checked build (`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`)
- GNU Fortran 14.2.0, optimized build (`-O3 -march=native`)
- Python with SciPy 1.17.0 and NumPy 2.3.5, independent reference comparison
- Fortran Package Manager not available; the full `fpm test` suite was **not** run

## Results

| Check | Checked | Optimized |
|---|---:|---:|
| `test_linalg_special_matrices` native test | PASS | PASS |
| Independent SciPy-reference cases | 478/478 | 478/478 |
| Floating matrix entries checked | 584,407 | 584,407 |
| Exact integer/decimal entries checked | 193,652 | 193,652 |
| Expected invalid/overflow statuses | 131 | 131 |
| Maximum scaled absolute difference from reference | 5.807e-15 | 5.807e-15 |
| Existing `test_linalg_companion` | PASS | PASS |

The Python harness includes all three SciPy constructors, full/compact Helmert
variants, floating Hilbert inverses of sizes up to 128, `int64` exact inverse
matrices through order 14, arbitrary-precision decimal inverses through order
64, negative/zero orders, rejected shapes and checked floating overflow at
order 250. For floating inverse entries, the Python harness compares against
SciPy's `invhilbert(n, exact=True)` integer results cast to `float64`, rather
than the less accurate SciPy `exact=False` floating combinatorial intermediate
calculations. Other matrix constructors are compared directly to SciPy's
ordinary floating outputs.

Native tests include full Helmert orthonormality for orders 1 through 96,
zero-sum contrasts, Hilbert definition and symmetry, small-order inverse
reconstruction, known exact values, int64 overflow rejection, zero dimensions,
wrong-shaped outputs and insufficient decimal character lengths.

### Reproduction

From repository root, with a GNU Fortran toolchain installed:

```sh
mkdir -p build_special && cd build_special
gfortran -std=f2008 -g -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow \
    -c ../src/kinds.f90 ../src/linalg.f90 \
    ../src/linalg_helmert.f90 ../src/linalg_hilbert.f90
gfortran -O0 -g -fcheck=all -o native_special kinds.o linalg.o \
    linalg_helmert.o linalg_hilbert.o ../test/test_linalg_special_matrices.f90
gfortran -O0 -g -fcheck=all -o driver_special kinds.o linalg.o \
    linalg_helmert.o linalg_hilbert.o ../reference/linalg_special_matrices_driver.f90
./native_special
python ../tools/compare_linalg_special_matrices.py ./driver_special
```

`fpm test` remains the recommended complete integration check in an
environment with FPM and the package's external dependencies available.
