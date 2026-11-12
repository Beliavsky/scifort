# Complex sparse LSQR validation (2026-11-06)

## Implemented

- Module: `src/sparse_linalg_complex_lsqr.f90`.
- New native test: `test/test_sparse_complex_lsqr.f90`.
- Streaming SciPy/NumPy reference: `tools/sparse_reference/compare_sparse_complex_lsqr_scipy.py`.
- Standalone native sparse runner registers the new source and regression.

## Independent comparisons

- Seed `20261106`: 150 checked random cases, 240 `-O3` cases.
  Maximum solution difference versus SciPy: `4.208e-11`;
  versus dense NumPy SVD: `4.927e-11`.
- Seed `246802`: 250 checked and 250 `-O3` random cases.
  Maximum difference versus SciPy: `6.505e-11`;
  versus dense NumPy SVD: `4.065e-11`.
- Each reference case checks CSR/CSC conversions, exact residual norms,
  damped normal-equation residuals, and an independently computed dense
  least-squares solution. Includes rank-deficient, underdetermined, empty,
  and regularized systems with nonzero `x0`.
- An additional native test constructs a 6000x120 complex sparse matrix
  from 18,000 COO entries without ever allocating a 6000x120 dense array.
  It checks convergence and solution/residual accuracy in both checked and
  full `-O3` builds.

## Regression suite

- The complete standalone sparse regression suite contains **47 programs**.
- All 47 passed with bounds-checked library compilation and `-O0` drivers.
- All 47 passed with library compiled at `-O3` and drivers at `-O0`.
- The new LSQR module *and its focused native driver* were also compiled and
  run with both compiler units at `-O3`.
- The prior complex triangular and Hermitian eigensolver comparison harness
  passed an additional 100 checked cases (100 triangular, 87 eigen) and
  160 optimized cases (160 triangular, 143 eigen).

## Limitations

- `fpm test` for the entire SciFort repository was not run: FPM is unavailable
  in the validation environment. The standalone sparse runner covers only
  sparse modules and does not exercise SciFort's other subsystems.
- LSQR condition estimates and variance estimates are not strict error bounds;
  no preconditioner or matrix-free callback is currently exposed.
- Optimized full-suite drivers used `-O0`; the focused new test and randomized
  reference driver both used `-O3` throughout.
