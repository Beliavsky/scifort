# 2026-11-08 generalized symmetric/Hermitian eigenproblem validation

The Cholesky-reduction extensions to `eigh` and `eigvalsh` were compiled
with GNU Fortran 14.2.0, using NumPy 2.3.5 / SciPy 1.17.0 as an independent
numerical reference. Both checked and optimized builds were exercised:

- Checked: `-std=f2018 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`
- Optimized: `-std=f2018 -O3`
- Python reference calls used one BLAS thread to avoid oversubscription.

## Completed generalized-eigenvalue reference runs

```bash
python tools/linalg_reference/run_generalized_eigen_checks.py --mode both --cases 144 --max-order 24 --seed 102
python tools/linalg_reference/run_generalized_eigen_checks.py --mode both --cases 500 --max-order 40 --seed 20261108
python tools/linalg_reference/run_generalized_eigen_checks.py --mode both --cases 200 --max-order 72 --seed 98741
python tools/linalg_reference/run_generalized_eigen_checks.py --mode both --cases 360 --max-order 64 --seed 20261109
```

**2,408 comparisons passed** in total. The largest observed normalized
absolute eigenvalue difference, using the maximum eigenvalue magnitude of
the full spectrum as scale, was about `4.79e-13`. The largest normalized
eigenvector equation residual was about `1.66e-15`, and the largest
normalized `B` or `B^-1` orthogonality error was about `1.88e-13`.

The comparisons include types 1/2/3, real and complex matrices, upper and
lower storage, subsets, zero and repeated eigenvalues, diagonal and full
positive-definite mass matrices, and scaling from `1e-100` to `1e100`.
The unused triangle was filled with NaNs in selected cases, and complex
imaginary diagonal entries were deliberately poisoned when unused.

The native test also checks positive-definiteness status for singular or
indefinite `B`, invalid shapes/types/ranges/iteration limits, the zero-size
case, and both real/complex eigenvector normalization equations. All native
tests passed in both compiler modes.

## Completed neighboring regression runs

```bash
python tools/linalg_reference/run_dense_eigen_checks.py --mode both --cases 160 --max-order 64 --seed 20261109
python tools/linalg_reference/run_dense_cholesky_checks.py --mode both --cases 60 --max-order 36 --seed 205
python tools/linalg_reference/run_banded_eigen_checks.py --mode both --cases 24 --max-order 18 --seed 210
```

**488** additional comparisons passed (320 dense standard, 120 dense
Cholesky, 48 banded eigen) in checked and optimized builds. A longer
combined regression command exceeded its execution-time limit while
running the banded-eigen comparison; only completed tests are counted.

## Limits

- A full `fpm test` run could not be performed because `fpm` is not
  installed in this environment.
- The Cholesky transformation is numerically sensitive when `B` is
  ill-conditioned; tightly clustered generalized eigenvalues may be
  less accurate than on well-conditioned mass matrices.
- This is an intrinsic-only `O(n^3)` dense algorithm; `subset_by_index`
  computes the complete spectrum first.
- Mixed-kind pairs, batched inputs, `subset_by_value` and LAPACK drivers
  remain unimplemented.
