# Complex sparse pivoted LU validation — 2026-11-06

The new `scifort_sparse_linalg_complex_direct` implements complex sparse
Gaussian elimination with partial row pivoting, fill-in and reusable factors,
without BLAS, LAPACK, SuperLU or dense intermediate matrices.

## Test environment

- GNU Fortran 14.2.0, standard `-std=f2018`
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0
- New code requires no runtime dependencies beyond existing SciFort modules.
  NumPy and SciPy are used solely by optional independent comparison tests.

## Focused regression and reference comparisons

```shell
python tools/sparse_reference/run_complex_direct_checks.py --mode both \
    --count 450 --seeds 2667 4113 6577
```

| Configuration | Focused native test | Randomized comparisons |
|---|---|---:|
| `-O0 -g -fcheck=all -fbacktrace` | PASS | 1,350/1,350 |
| `-O3` | PASS | 1,350/1,350 |

**2,700/2,700 total** comparisons succeeded. Maximum absolute difference
from SciPy's solution was about **4.45e-16** across tested systems. The
reference harness independently checked residuals of the normal, transpose
and Hermitian-adjoint systems using NumPy. SciPy factorization and solutions
are computed with `splu`/`spsolve` from the SciPy installation.

Randomized inputs exercise CSR/CSC conversion, unsorted COO assembly, duplicate
coordinates and stored zeros, multiple RHS, nonsymmetric complex values,
explicit row permutations requiring partial pivoting, and all three `trans`
options. The native regression further tests singular and invalid matrices,
a 0-by-0 system, empty RHS blocks, pivot thresholds, valid factor permutations,
and a 256-by-256 tridiagonal sparse matrix, verifying that the LU remains
sparse and is not constructed via a dense matrix.

## Existing sparse regression suite

```shell
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized --fast-tests
python tools/sparse_reference/run_native_sparse.py --mode optimized --fast-tests \
    --start-at test_sparse_lanczos
```

- **45/45** bounds-checked native regression programs passed in one run.
- **45/45** optimized-library regression programs passed across two runs:
  32 through `test_sparse_arnoldi`, then 13 from `test_sparse_lanczos` to
  the end. The full suite compiled the library objects at `-O3`; its test
  drivers were compiled at `-O0` using `--fast-tests` to avoid the execution
  time limit. A single full run with both library and all test drivers at
  `-O3` did not complete within that limit (10 tests passed beforehand).
- The focused complex-direct test was separately compiled **entirely** at
  `-O3` and passed.
- Full project-wide `fpm test` was not executed because FPM is unavailable
  in this environment. These results cover the standalone sparse suite.

## Limitations and future work

The LU algorithm is unblocked sparse Gaussian elimination with partial row
pivoting, not SuperLU's optimized sparse direct implementation. It omits
fill-reducing column reordering, scaling, iterative refinement and supernodes.
Its `pivot_tol` is an **absolute** threshold, unlike SuperLU's relative
diagonal pivot threshold. Extremely large sparse unsymmetric systems may
exhibit substantial fill-in and high runtime. The algorithm is currently
limited to `complex(dp)` one-based CSR/CSC matrices and dense RHS vectors or
blocks.
