# Spectral-norm estimators: validation (2026-11-09)

## Scope

New SciFort `scifort_linalg_spectral_norm` routines implement randomized
power-method estimates of operator spectral norms for dense real/complex
matrices and pairs of matrices, and callback-based operators. The algorithms
are original MIT-licensed code, with no added runtime dependencies.

## Compiler and reference versions

- GNU Fortran 14.2.0; `-std=f2018 -O0 -g -fcheck=all -fbacktrace`.
- GNU Fortran 14.2.0; `-std=f2018 -O3 -march=native`.
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0.

## Native checks

`test/test_linalg_spectral_norm.f90` passed **7,293 checks in each build**:

- Real and complex dense matrices of shapes from `0x0` to `12x13`.
- Corresponding matrix-free operator callback results with the same
  starting vector, using forward and conjugate-adjoint products.
- Norm-of-difference estimates using real and complex pairs.
- Zero, rank-one, rectangular, near-rank-deficient and highly scaled matrices.
- Reproducible seeded draws, negative iterations, incorrect sizes,
  nonfinite inputs and failed callback statuses, plus near-cancellation
  accuracy and unrepresentable-output detection.

## SciPy and exact SVD cross-checks

`test/compare_linalg_spectral_norm_scipy.py` passed **840 cases per build**
(210 each for real norm, real difference, complex norm and complex difference).
Of these, 540 had nonzero exact spectral norms.

| Quantity | Checked | Optimized |
|---|---:|---:|
| Native checks | 7,293 | 7,293 |
| Python reference cases | 840 | 840 |
| Largest relative difference from independently reconstructed power iterates | 5.632e-16 | 5.632e-16 |
| Largest relative difference for matrix-free operator vs Python iteration | 5.195e-16 | 4.971e-16 |
| Largest relative error vs exact NumPy SVD norm | 8.883e-03 | 8.883e-03 |

SciPy's public `estimate_spectral_norm` was also run with independent seeded
random streams on each case. Its maximum relative error against the exact
SVD norm over these tests was **2.726e-02**. Differences between independent
randomized estimators are **not defects**, and neither routine promises a
fixed convergence rate for arbitrary spectral gaps.

The full `fpm test` suite was not run because FPM is not installed in the
execution environment. Tests were compiled and run directly with gfortran.

## Compatibility notes

- SciPy accepts dense arrays and `LinearOperator` values. SciFort exposes
  separate dense-array and explicit forward/adjoint callback interfaces.
- Unlike SciPy's `rng`, SciFort's `seed` uses a self-contained PRNG and does
  not reproduce NumPy bit streams. Passing `start` yields deterministic
  cross-language comparisons.
- Complex callback entry points have `_complex` suffix because Fortran
  generic resolution cannot dispatch by callback-procedure signature alone.
- Calls with zero dimensions return zero in SciFort, rather than relying on
  SciPy's behavior for empty operators.
