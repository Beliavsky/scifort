# SciFort DOK/LIL validation — November 5, 2026

Environment: GNU Fortran 14.2.0, Python 3.13, NumPy, SciPy 1.17.0.

## New numerical and structure comparisons

```console
python tools/sparse_reference/compare_sparse_dok_lil_scipy.py --trials 300 --mode both --seed 20261105
python tools/sparse_reference/compare_sparse_dok_lil_scipy.py --trials 200 --mode both --seed 84723
```

Results:

| Configuration | Randomized cases | Result |
|---|---:|---|
| Bounds-checked `gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace` | 500 | PASS |
| Optimized `gfortran -std=f2018 -O3` | 500 | PASS |

Each randomized case compares both DOK and LIL arrays to SciPy, including
COO construction with duplicates/explicit zeros, 0-by-N and M-by-0 shapes,
assignment/addition/deletion/lookup sequences, `nnz` after every operation,
full dense content, matrix-vector products, and canonical CSR index/data
arrays. In total, 1,000 randomized cases (two formats per case) passed.

## Standalone native tests

```console
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized
```

- Bounds-checked: **34/34 native sparse regression executables passed**.
- Optimized: **34/34 native sparse regression executables passed**.
- The new `test_sparse_dok_lil` additionally exercises 2,600 deterministic
  updates, hash-table growth and deletion, sorted per-row insertion,
  transpose, cross-format conversions, independent deep copies, stored-zero
  elimination, large sparse dimensions, and invalid/nonfinite inputs.
- A final standalone `test_sparse_dok_lil` rerun passed after adding tests
  for nonfinite matrix-vector products, without changes to library code.

The project-wide `fpm test` command could **not** be executed because the
Fortran Package Manager is not installed in this environment. No existing
external dependencies or external numerical source code were changed.
