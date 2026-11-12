# Complex sparse matrix-function actions — validation (2026-11-06)

New module: `src/sparse_linalg_complex_funm.f90`.

- Native bounds-checked sparse suite: **54/54 tests pass**, including
  `test_sparse_complex_funm` and all previously present standalone sparse tests.
- Focused `-O0 -g -fcheck=all -fbacktrace` native test passes.
- Randomized comparison seed `20261106`: **150 random matrices per build**, six
  CSR/CSC/operator vector/block comparisons per case: **900 comparisons per mode**.
  In this batch, 39 public SciPy Krylov outputs per mode also agreed with
  independent dense references.
- Randomized comparison seed `20261107`: **220 matrices per build**, **1,320
  comparisons per mode**; another 63 qualifying SciPy callable cross-checks
  per mode passed.
- Total randomized checks so far: **4,440 individual Fortran/dense-SciPy path
  comparisons** across both modes, plus **204 matching SciPy Krylov outputs**.
  Maximum absolute error: less than `1.2e-14` in these cases.
- Native large-scale matrix-free test: **1,024 complex variables**, no dense
  input matrix; affine action agrees with exact operator application.
- Tests include complex Hermitian Lanczos, restarted nonsymmetric Arnoldi,
  negative/complex scale factors, dense block RHS, nonfinite callback results,
  exhausted budgets, and zero inputs.

The two compiler configurations are GNU Fortran checked `-O0 -fcheck=all`
and optimized `-O3`. The optimized native full-suite attempt passed its first **10/54** programs,
including the new test, before reaching the execution limit. It is not claimed
as a complete optimized suite. The new native test and documentation example
pass with all their library objects and drivers compiled using `-O3`.
The full project-wide FPM test suite has not been run in this environment.

The dense SciPy matrix functions are the primary reference. SciPy 1.17's
public Krylov implementation raises errors for certain zero-input cases;
those are excluded from the direct SciPy Krylov cross-check count but remain
covered by SciFort's native tests and independent dense comparisons.

This is an original MIT-licensed implementation and test harness. It adapts
SciFort's existing real restarted matrix-function Krylov design, not SciPy
source code. No external numerical runtime dependencies were introduced.
