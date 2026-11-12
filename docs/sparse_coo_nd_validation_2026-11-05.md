# N-dimensional COO validation — November 5, 2026

Environment: GNU Fortran 14.2.0, NumPy and SciPy 1.17.0.

## Randomized comparison with SciPy

```console
python tools/sparse_reference/compare_sparse_coo_nd_scipy.py --mode both --trials 300 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_scipy.py --mode both --trials 250 --seed 274829
```

| Compiler configuration | Cases | Result |
|---|---:|---|
| Bounds-checked `gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace` | 550 | PASS |
| Optimized `gfortran -std=f2018 -O3` | 550 | PASS |

The 1,100 total cases cover ranks 1–7, zero and singleton dimensions,
unsorted duplicate coordinates, explicit zero values, canonical sorted
coordinate/data arrays, C-order dense flattening, C/F-order reshape, axis
permutation, scalar lookup, from-dense roundtrips, and native 2D COO
interoperability. Comparisons were made with SciPy 1.17.0 `coo_array`
operations, not with a reimplementation of SciFort's Fortran recurrence.

## Native sparse regression

```console
python tools/sparse_reference/run_native_sparse.py --mode both
```

- Bounds-checked sparse suite: **35/35 native tests PASS**.
- Optimized sparse suite: **35/35 native tests PASS**.
- The focused `test/test_sparse_coo_nd.f90` also passed after addition of
  5,000-entry sorting stress tests and rank-7 assumed-rank dense conversion.
- The focused test covers invalid shapes, invalid coordinates, nonfinite
  input/aggregation, transactional duplicate summation, allocation caps,
  `int64` product overflow and large sparse shapes without dense allocation.

This runner compiles the sparse subset without external numerical packages.
A full project-wide `fpm test` was not executed because the FPM executable
is unavailable in this runtime.
