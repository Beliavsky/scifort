# Randomized rank estimation: validation (2026-11-09)

## Implementation

New module `scifort_linalg_estimate_rank` implements a randomized range finder
using repeated modified Gram–Schmidt. Dense matrices and matrix-free
operators are supported for real and complex arithmetic. No new
third-party sources or dependencies were introduced.

## Environment

- GNU Fortran 14.2.0, checked: `-std=f2018 -O0 -g -fcheck=all -fbacktrace`.
- GNU Fortran 14.2.0, optimized: `-std=f2018 -O3 -march=native`.
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0.

## Results

| Test | Checked | Optimized |
|---|---:|---:|
| Native rank checks | 5,771 | 5,771 |
| Independent Python/Fortran cases | 472 | 472 |
| Cases with nonempty SVD rank verification | 448 | 448 |
| Previous interpolative-decomposition native checks | 14,944 | 14,944 |
| Previous spectral-norm native checks | 7,293 | 7,293 |

Native tests cover zero matrices, exact low-rank real/complex data, full rank,
rectangular and empty arrays, data scales `1e-90` and `1e90`, seeded
reproducibility, dense/operator agreement, small singular values below a
specified tolerance, invalid tolerances, nonfinite data, and callback errors.

The independent Python harness builds known-spectrum matrices from discrete
sine bases, runs SciFort's real and complex dense and operator estimators,
checks the result against the known numerical rank, and independently checks
nonempty matrices against NumPy SVD rank. SciPy's own `estimate_rank` is run
with separate generators for both array and `LinearOperator` input to
identify the expected stochastic differences. SciPy's dense backend may
report the exact rank plus roughly 7–8 (or more for special small or zero
matrices); SciFort does not add that padding.

The full `fpm test` suite was not run: FPM is not installed here. The
listed native regressions were compiled directly with `gfortran`.
