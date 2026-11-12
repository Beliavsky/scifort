# SciFort sparse Lanczos validation — 2026-11-04

## Implementation

The new `eigsh_lanczos` entry point implements real-symmetric, thick-
restarted Lanczos-style Krylov projection with two-pass reorthogonalization.
Selected Ritz vectors are retained with soft locking. Matrix-vector products
use canonical one-based CSR. The method is dependency-free and retains true
residual diagnostics; no dense n-by-n array is needed.

`eigsh(...,method='lanczos')` and `eigs(...,method='arnoldi')` are opt-in;
default `method='subspace'` remains unchanged. The restarted methods do
support `sigma=` shift-invert after the follow-up update.

## Reproduction

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_lanczos_scipy.py --mode both --cases 72
```

The comparison driver independently verifies the selected eigenvalues
against `numpy.linalg.eigvalsh`, true residuals and orthogonality against
the original matrix, and (for well-separated reduced-subspace cases)
`scipy.sparse.linalg.eigsh`. Additional tests cover a 100-by-100
slowly varying tridiagonal spectrum requiring multiple restart cycles,
zero iteration budgets, nonconvergence, symmetry validation, all four
selection modes, and method dispatch.

## Scope

`LM`, `LA`, `SA` support reduced subspaces, but `SM` needs `ncv=n`.
Normal-mode real shift-invert is supported as of the subsequent
restart update. Hard convergence locking remains future work. Some small subspaces stagnate on clustered
spectra; positive `info` signals this without falsely reporting success.

## Recorded results

- GNU Fortran 14.2, `-O0 -g -fcheck=all -fbacktrace`: all 10 native
  sparse regression programs passed, including `test_sparse_lanczos`.
- GNU Fortran 14.2, `-O3`: all 10 native sparse regression programs passed.
- SciPy 1.17 / NumPy comparison: 72 randomized matrices and two difficult
  100-row spectral-edge matrices per build with seed 20261104.
- An additional 48 randomized matrices and the same two difficult matrices
  per build passed with seed 98765.
- Total: **120 randomized matrices and four difficult cases per build**;
  true residuals, spectral ranking, eigenvector orthogonality, and SciPy
  spectral agreement were checked as appropriate. The two difficult cases
  in each comparison batch required multiple restart cycles.
- The existing `eigsh` comparison script passed an additional 12 optimized
  cases. The existing Arnoldi standalone native tests passed in both modes;
  a separate combined legacy comparison command reached the environment's
  execution limit before completing the Arnoldi reference tests.

The older callback-based LSMR code still triggers a GNU executable-stack
linker warning. The full FPM suite was not run because FPM is unavailable
in this environment. These are pre-existing environmental limitations.
