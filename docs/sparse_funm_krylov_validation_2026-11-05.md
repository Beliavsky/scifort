# Sparse matrix-function Krylov validation (2026-11-05)

## Coverage

- Native regression: `test/test_sparse_funm_krylov.f90`.
- Cross-language randomized reference:
  `tools/sparse_reference/compare_funm_krylov_scipy.py` and
  `tools/sparse_reference/funm_krylov_driver.f90`.
- Native suite: `tools/sparse_reference/run_native_sparse.py`.
- GNU Fortran 14.2.0, SciPy 1.17.0, real binary64.

The native regression covers exact 2x2 rotations, polynomial functions,
restarted symmetric Lanczos and nonsymmetric Arnoldi approximations against
independent `expm_multiply` references, CSR versus matrix-free equality,
forced nonconvergence, zero RHS, empty/rectangular matrices, zero operators,
zero time, oversized dimensions, invalid tolerances, nonfinite matrix-vector
products and matrix-function callbacks, callback error statuses, and a
2048-component matrix-free second-difference operator without a dense input
matrix.

The random driver tests diagonal, symmetric banded and nonsymmetric sparse
matrices, multiple `t` scales, randomly seeded vectors, `exp`, `sin`,
`x**3-x`, and `I+x+0.1*x**2`. For each case, the Fortran CSR and matrix-free
versions are cross-checked, and their result is compared against a dense
SciPy matrix-function evaluation, not just SciPy's Krylov approximation.
Some SciPy 1.17 cases fail with a restarted-breakdown shape exception;
those are reported separately and do not mask independent dense checks.

## Results

| Test | Bounds-checked | Optimized |
|---|---:|---:|
| Random comparisons (seed 20261105) | 400/400 | 400/400 |
| Independent seed 3381 | 200/200 | 200/200 |
| Full native sparse executables | 32/32 | 32/32 |
| Focused 2048-component matrix-free test | PASS | PASS |

GNU Fortran 14.2.0 used `-O0 -g -fcheck=all -fbacktrace` and
`-O3 -march=native`; SciPy was version 1.17.0. The largest absolute
error in the 400-case sample was `4.326e-9` in each build. Maximum
matrix-free/CSR component discrepancy was `4.441e-16` in the checked
build and `8.882e-16` optimized. The alternate-seed sample maximum
absolute error was `1.315e-9`. SciPy itself raised restarted-breakdown
shape errors in 12/400 and 3/200 sampled cases; these are tracked rather
than silently skipped as reference comparisons. Independent dense SciPy
matrix functions were successfully evaluated in all 600 cases per build.

`fpm` is not installed in the validation environment, so the complete
project-wide FPM test suite was **not** run. The native sparse suite was
compiled and run independently.

## Reproduction

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
# After compiling the helper driver with the sparse core modules:
python tools/sparse_reference/compare_funm_krylov_scipy.py \
    /path/to/funm_krylov_driver --seed 20261105 --cases 400
```

The generic Krylov solver is real-only and uses update-based stopping rather
than a certified error bound for arbitrary matrix functions. Convergence
can depend strongly on `f` and on the spectrum of `A`.
