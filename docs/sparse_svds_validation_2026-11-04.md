# Partial sparse SVD: validation and limitations (2026-11-04)

## Design

`scifort_sparse_linalg.svds` constructs a right trial subspace by
Golub--Kahan bidiagonalization of the taller orientation of the CSR matrix.
It performs two reorthogonalization passes and seeds a new independent Krylov
chain after a lucky breakdown. Projecting `A*V` into the resulting thin trial
basis allows a one-sided Jacobi SVD without explicitly forming the full
normal matrix `A^T A`. Ritz triplets are accepted only after computing
**both** true residuals against A and its transpose.

The current algorithm is **not implicitly restarted Lanczos bidiagonalization**.
Its `ncv` limit is a hard cap, and nonconvergence is returned instead of
silently treating unconverged projected triplets as accurate. For `SM`, the
full shorter dimension is required rather than using an unreliable partial
subspace to estimate the global smallest values.

## Reproducing the tests

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_svds_scipy.py --mode both --cases 180 --partial 70
python tools/sparse_linalg_reference/compare_svds_scipy.py --mode optimized --cases 100 --partial 35 --seed 20261105
```

Reference cases include tall and wide matrices, dense and sparse entries,
rank deficiency, repeated real singular values, stored zeros and duplicated
COO coordinates, and partial 55--210-dimensional sparse tests with clustered
background singular values and isolated dominant values.

The standalone native suite includes a **1,200-by-900** CSR matrix, where
only 35 Krylov directions are retained, plus argument validation and a
specific nonconvergence check. Python verifies singular values against
NumPy or SciPy, reconstructs two-sided residuals independently, and tests
singular-vector orthogonality.

## Deliberate API and numerical differences

- Real binary64 CSR only; no complex inputs or generic `LinearOperator` yet.
- SciFort returns sorted singular values and always returns U and V^T.
- `k=min(m,n)` is supported; SciPy's ARPACK backend restricts k further,
  while its PROPACK backend permits full k.
- The default `ncv` is bounded to avoid dense matrix construction for
  large `LM` tasks; `SM` currently needs `ncv=min(m,n)`, which can use
  quadratic memory. This is a limitation, not a fast large-scale SM method.
- Exact zero singular values receive independent left-null directions when
  the trial space spans the shorter side. Extremely ill-conditioned
  singular values at machine precision may be treated as numerical zero.
- Positive `info` indicates a nonconverged but valid partial result;
  convergence reflects individual residuals, not guaranteed discovery of
  the globally requested spectral subspace in an adversarial initial vector.
- Tests are numerical comparisons only. No upstream SciPy, PROPACK, BLAS or
  LAPACK implementation has been imported. The separate pre-existing LSMR
  executable-stack linker warning is unchanged.
