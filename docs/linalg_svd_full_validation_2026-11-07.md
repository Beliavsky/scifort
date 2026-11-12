# Full dense SVD and singular-value utility validation (2026-11-07)

The new `full_matrices` option in `scifort_linalg_svd.svd`, real/complex
`svdvals`, and rectangular `diagsvd` were evaluated against SciPy 1.17.0
using GNU Fortran 14.2.0 and Python/NumPy 2.3.5.

## Reproduce

From the project root (with SciPy and `gfortran` available):

```bash
python tools/linalg_reference/run_svd_full_checks.py --mode both --cases 300 --max-order 46 --seed 885
python tools/linalg_reference/run_svd_full_checks.py --mode both --cases 250 --max-order 70 --seed 1157
python tools/linalg_reference/run_svd_checks.py --mode both --cases 160 --max-order 25 --seed 911
python tools/linalg_reference/run_subspaces_checks.py --mode both --cases 160 --max-order 25 --seed 1123
python tools/linalg_reference/run_qr_checks.py --mode both --cases 100 --max-order 25 --seed 987
```

The Python harness compiles its drivers with `-std=f2018` and either
`-O0 -g -fcheck=all -fbacktrace` (checked) or `-O3` (optimized).
Standalone native test source: `test/test_linalg_svd_full.f90`.
Runnable example source: `example/linalg_svd_full_example.f90`.

## Results

- Four full-SVD random suites with seeds 805 (160 cases), 885 (300 cases),
  842 (120 cases), and 1157 (250 cases) passed in each build: **1,660**
  checked-and-optimized comparisons in total. The first 460 per mode
  preceded the introduction of explicit `diagsvd` comparisons; the last
  370 per mode also verify `diagsvd` against SciPy.
- Full-SVD reconstruction normalized by input Frobenius norm was at most
  approximately `4.15e-15` in these suites.
- Square-factor unitary/orthogonal error divided by dimension was at most
  approximately `4.4e-16`.
- Singular-value relative normwise error was at most approximately
  `1.28e-15`; rectangular `diagsvd` differences were exactly zero.
- Relative unused-factor null-space residual was below `6.8e-16`; full
  factor subspace projector discrepancies (when the small singular values
  were well separated from zero) were below `5.4e-15`.
- Native tests (both builds) passed for real/complex 1x1 to 8x8 matrices,
  rank-deficient cases, zero matrices, shape checks, invalid sweep limits,
  nonfinite inputs, and zero-length `diagsvd` outputs.
- Legacy economy SVD/pseudoinverse/least-squares (320 comparisons),
  subspace geometry (320 comparisons), and QR (200 comparisons) passed in
  the checked and optimized builds. Existing native subspace tests passed.

## Numerical and build limits

- For full size, an orthogonal complement is constructed by pivoted
  canonical-vector completion and two passes of modified Gram-Schmidt.
  On matrices with extremely ill-conditioned singular spaces, the new
  complement is not guaranteed to have LAPACK's robustness.
- The underlying Jacobi SVD remains intrinsic-only, unblocked, and slower
  than optimized LAPACK on large matrices. `svdvals` currently computes
  singular vectors internally, so it does not save runtime or workspace.
- Nonempty `m,n >= 1` inputs are required by `svd` and `svdvals`.
- For backward compatibility SciFort's `svd` defaults to economy vectors,
  unlike SciPy's default full matrix factors.
- The complete `fpm test` suite was not executed because FPM is not installed
  in the testing environment. The focused tests compile without FPM.
