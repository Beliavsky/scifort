# Eigenvalue interval selection: validation (2026-11-08)

## Build and reference

- Compiler: GNU Fortran 14.2.0, `-std=f2018`.
- Checked: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`.
- Optimized: `-O3`.
- Independent reference: SciPy 1.17.0 and NumPy (installed locally).
- BLAS-related reference-test thread counts constrained to one for
  reproducible performance on small matrices.

The new `tools/linalg_reference/run_eigen_by_value_checks.py` harness compiles
both a native status/error-path test and a batch eigensolver driver in each
compiler configuration. The independent Python comparisons use
`scipy.linalg.eigh`/`eigvalsh` with `subset_by_value`.

## Completed randomized SciPy comparisons

| Command seed | Cases per build | Builds | Total passed | Max normalized eigenvalue difference |
| --- | ---: | ---: | ---: | ---: |
| 20261108, maximum order 24 | 110 | 2 | 220 | `1.230e-15` |
| 1907, maximum order 48 | 320 | 2 | 640 | `1.183e-13` |
| 20261109, maximum order 64 | 210 | 2 | 420 | `3.384e-15` |
| 20261111, maximum order 40 | 100 | 2 | 200 | `1.547e-15` |
| **Total** | | | **1,480** | **`1.183e-13`** |

The maximum normalized residual across these tests was about `5.47e-14`;
the maximum generalized/ordinary orthogonality or normalization error was
about `2.70e-15`. These errors are measured relative to matrix scales rather
than only the selected eigenvalues; such normalization avoids a misleading
large relative error when an eigenvalue is nearly zero.

Cases exercise real and complex types, upper and lower input triangles,
nonfinite values in the unused triangle, generalized types 1/2/3, repeated
spectra, zero matrices, empty selections, IEEE infinite bounds, full-spectrum
selections, and very large/small coefficient scales. Selection endpoints in
randomized tests are placed in spectral gaps to avoid ambiguity from
different last-bit eigenvalue rounding in LAPACK and the native QL solver.
An additional native test verifies exact half-open endpoints and repeated
values, undersized outputs, conflicting selectors, NaN intervals, and zero-size
outputs. A standalone Fortran example also passed in the checked build.

## Existing regression results

- `run_dense_eigen_checks.py --mode both --cases 190 --max-order 40 --seed 184`:
  380 SciPy comparisons plus native tests passed.
- `run_generalized_eigen_checks.py --mode both --cases 190 --max-order 32 --seed 185`:
  380 SciPy comparisons plus native tests passed.
- `run_banded_eigen_checks.py --mode both --cases 70 --max-order 24 --seed 20261111`:
  140 neighboring banded-eigenvalue comparisons plus native tests passed.

## Numerical and infrastructure limitations

This feature filters an already-computed full eigenvalue spectrum and does not
implement a partial-spectrum algorithm. Eigenvalues on an exact interval
boundary can round differently from those computed by SciPy's LAPACK drivers.
Fortran Package Manager (`fpm`) was unavailable in this execution environment,
so `fpm test` was not run; focused GNU Fortran tests with runtime checks were
run instead. The results reported here do not claim a full-project FPM test.
