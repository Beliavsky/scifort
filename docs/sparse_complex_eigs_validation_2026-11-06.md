# Complex nonsymmetric sparse eigs validation — 2026-11-06

Implementation: `src/sparse_linalg_complex_eigs.f90`.

## Numerical comparisons

Tested against independently computed dense `numpy.linalg.eigvals`, and
cross-checked a subset with `scipy.sparse.linalg.eigs` (SciPy 1.17.0,
NumPy 2.3.5). These references are **not** dependencies of the Fortran code.
The comparison suite includes random sparse, full, diagonal, and upper
triangular nonsymmetric complex matrices, repeated COO entries, CSR/CSC
storage, optional nonzero starting vectors, and each of the six ordering
modes (`LM`, `SM`, `LR`, `SR`, `LI`, `SI`). The driver independently evaluates
`||A v - lambda v||` using the complex CSR matrix-vector multiplication.

Reproduction:

```sh
python tools/sparse_reference/compare_sparse_complex_eigs_scipy.py --mode both --cases 240 --seed 20261106
python tools/sparse_reference/compare_sparse_complex_eigs_scipy.py --mode both --cases 160 --seed 748365
```

| Configuration | NumPy eigenvalue comparisons | Direct SciPy `eigs` checks | Maximum eigenvalue discrepancy versus NumPy |
| --- | ---: | ---: | ---: |
| Checked `-O0 -fcheck=all` | 400/400 | 66 | 6.56e-14 |
| Optimized `-O3` | 400/400 | 66 | 3.51e-14 |

The independent residual in every randomized test was less than `1e-6`
(typically far smaller). Matching of dense and sparse selected eigenvalues
uses the requested ordering, on randomly generated, well-separated cases.

## Native regression

- **50/50** standalone sparse regression programs passed with bounds checks,
  including the new `test_sparse_complex_eigs`.
- **50/50** passed with library modules built at `-O3` and regression
  drivers built at `-O0` (`--fast-tests`). The new module and its randomized
  reference driver were **also both compiled at `-O3`** separately.
- Focused tests cover full-basis eigenvalue modes, CSC forwarding, verified
  complex eigenvector residuals, repeated dominant eigenvalues, a 48-by-48
  restarted partial Arnoldi search, an intentionally exhausted search budget,
  and invalid options and dimensions.
- The documentation example compiles and runs with bounds checking.
- The full external-dependency `fpm test` suite was not run: FPM is not
  installed in this environment.

## Scope and limitations

This is a bounded explicit-restart Arnoldi implementation, not ARPACK's
implicitly restarted Arnoldi algorithm. True residual tests validate the
returned Ritz eigenpairs but do **not** certify that an excluded invariant
subspace contains no more highly ranked eigenvalue. Results for very close,
repeated, or defective eigenvalues may require a larger basis or a future
block/Schur-locking method. `SM` requires a full basis; complex matrix-free eigen operators are not yet available.
Shift-invert is documented separately in `sparse_complex_eigs_shift.md`.
