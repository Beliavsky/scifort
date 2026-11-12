# Randomized interpolative decomposition — validation, November 10, 2026

## Source

New Fortran module `src/linalg_interpolative_random.f90`, with real and
complex fixed-rank and adaptive randomized interpolative decompositions.
There are no changes to existing numerical source modules or FPM dependencies.

## Build and tests

Compiler: GNU Fortran 14.2.0 on Linux, standard free source-form Fortran 2008.
Independent reference: Python 3.13, NumPy 2.3.5, SciPy 1.17.0.

| Test | Checked (`-O0 -fcheck=all`) | Optimized (`-O3 -march=native`) |
|---|---:|---:|
| New native checks | 37,449 | 37,449 |
| NumPy/SciPy cross-language cases | 416 | 416 |
| Successful SciPy randomized reference calls | 260 | 260 |
| SciPy randomized reference exceptions (not SciFort failures) | 20 | 20 |
| Existing deterministic interpolative ID native checks | 14,944 | 14,944 |
| Existing deterministic ID-to-SVD native checks | 32,770 | 32,770 |
| Existing randomized interpolative SVD native checks | 50,408 | 50,408 |
| Existing QR/SVD native regression suites | Passed | Passed |
| New example | Passed | Passed |

The new native tests cover real/complex, fixed rank zero/middle/full, adaptive
Frobenius tolerances 0, 1e-5, 0.2, empty and all-zero matrices, rank-one
matrices, coefficients at 1e-80 and 1e80, invalid ranks/tolerances,
negative oversampling, NaNs, coefficient shapes, column permutations,
seed repeatability, and reconstruction using the pre-existing ID API.

The cross-language cases include independent random rectangular matrices,
low-rank and noisy-low-rank matrices, scale extremes, and both numerical types.
Each valid SciFort decomposition is reconstructed from **original source
columns**, and checked against independent NumPy computations. Where SciPy's
randomized ID accepts a case, the SciPy reconstruction is also validated;
random pivot choices are not expected to be identical.

Observed statistics in both builds:

- Largest relative reconstruction error across *all* fixed-rank inputs was
  `1.0`, occurring when the requested rank is zero (the zero approximation).
- Largest certified relative Frobenius residual in adaptive tests:
  `8.360e-02`, below the requested tolerance for its test case.
- Largest fixed-rank reconstruction error divided by the **optimal**
  truncated-SVD Frobenius error (for cases with a measurable nonzero tail):
  `1.598104`.
- The example reconstructs rank-three real and complex matrices to relative
  errors of order 1e-16.

In reference tests, SciPy's randomized ID raised 20 exceptions per run,
including an internal `IndexError` on some degenerate matrices. These cases
were not removed from the SciFort checks; all SciFort cases passed.

The complete repository-wide `fpm test` suite was **not** run because the FPM
executable is unavailable in this environment. The above checks compile
only the needed SciFort numerical modules manually, without external BLAS.

## Reproduce manually

Compile source modules in order:

```text
src/kinds.f90
src/linalg.f90
src/linalg_svd.f90
src/linalg_qr.f90
src/linalg_interpolative.f90
src/linalg_interpolative_random.f90
```

Then link `test/test_linalg_interpolative_random.f90` or
`validation/reference_drivers/driver_linalg_interpolative_random.f90` with the module objects.
Run independent comparison:

```bash
python test/compare_linalg_interpolative_random_scipy.py ./random_id_driver
```
