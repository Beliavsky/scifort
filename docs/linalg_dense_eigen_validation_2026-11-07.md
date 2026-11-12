# 2026-11-07 dense symmetric/Hermitian eigenproblem validation

The new `eigh` and `eigvalsh` routines were compared with independently
computed `scipy.linalg.eigvalsh` eigenvalues. Eigenvector columns were
checked by residual and orthogonality rather than matching a nonunique
basis in repeated eigenspaces. Tests compiled with GNU Fortran 14.2.0
and used Python NumPy 2.3.5 / SciPy 1.17.0.

## Completed runs

```sh
python tools/linalg_reference/run_dense_eigen_checks.py --mode both \
  --cases 120 --max-order 64 --seed 420
python tools/linalg_reference/run_dense_eigen_checks.py --mode both \
  --cases 500 --max-order 96 --seed 20261108
python tools/linalg_reference/run_dense_eigen_checks.py --mode optimized \
  --cases 260 --max-order 128 --seed 20261109
python tools/linalg_reference/run_dense_eigen_checks.py --mode both \
  --cases 240 --max-order 64 --seed 4401
```

All **1,980 Fortran-versus-SciPy comparisons passed**:

- 240 smoke comparisons (120 in each compiler configuration).
- 1,000 main comparisons (500 in each compiler configuration).
- 260 additional optimized comparisons at maximum order 128.
- 480 independent-seed comparisons (240 per configuration).

In the main test, the maximum normalized eigenvalue disagreement was
`2.664e-14`; eigenpair residual was `6.866e-16` and normalized
orthogonality error `1.776e-15` at worst. Eigenvalue errors are scaled by
the largest absolute matrix entry, not by the smallest selected eigenvalue,
so that subsets near zero are evaluated with a normwise backward-error
criterion rather than artificially inflated relative forward errors.

Both checked (`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`) and
optimized (`-O3`) native tests passed. Tests cover valid and invalid
subsets, ignored matrix triangle, nonfinite used entries, n=0, large
scale, complex Hermitian input, and eigenvector residuals.

## Neighboring regression tests

The following tests passed with **70 cases in each** compiler mode:

- `run_banded_eigen_checks.py` (140 comparisons).
- `run_hessenberg_checks.py` (140 comparisons).
- `run_balance_checks.py` (140 comparisons).

This adds **420 passing regression comparisons**. The full `fpm test`
suite was not attempted because `fpm` is not installed in this environment.
