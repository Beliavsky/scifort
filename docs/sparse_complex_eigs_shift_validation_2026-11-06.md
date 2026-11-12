# Complex shift-invert eigs validation — 2026-11-06

Implementation: `src/sparse_linalg_complex_eigs.f90`; the new optional
`sigma` argument composes SciFort's existing complex LU factorization with
its complex Arnoldi/Rayleigh--Ritz code.

## Numerical comparisons

NumPy `numpy.linalg.eigvals` independently provides the full complex spectrum
of randomized square matrices. Expected eigenvalues are ordered by each of
the six SciPy shift-invert transformed-eigenvalue rules,
`mu=1/(lambda-sigma)`. A subset of cases also uses
`scipy.sparse.linalg.eigs(...,sigma=sigma)` as an additional reference. The
Fortran driver computes an independent true residual `||A*v-lambda*v||`.

Reproduction:

```sh
python tools/sparse_reference/compare_sparse_complex_eigs_shift_scipy.py --mode both --cases 110 --seed 20261106
python tools/sparse_reference/compare_sparse_complex_eigs_shift_scipy.py --mode both --cases 75 --seed 739184
```

| Build | Randomized NumPy cases passed | Direct SciPy cross-checks passed | Largest eigenvalue error |
| --- | ---: | ---: | ---: |
| Checked (`-fcheck=all`) | 185/185 | 32/32 | 5.425e-14 |
| Optimized (`-O3`) | 185/185 | 32/32 | 5.413e-14 |

These 370 NumPy cases and 64 additional direct SciPy comparisons cover
CSR/CSC storage, genuinely complex nonsymmetric matrices, upper-triangular
and diagonal cases, repeated COO entries, different target shifts, six
ranking modes, and nonzero initial guesses. The full basis (`ncv=n`)
was used in these randomized comparisons to make the complete ordering
independently checkable. The native test additionally exercises a truncated
Krylov search on a 70-by-70 complex banded matrix with `ncv=16`.

## Native Fortran tests

- The checked standalone sparse regression runner completed **51/51** tests
  successfully with library bounds checks, including the new
  `test_sparse_complex_eigs_shift.f90` and all previously present tests.
- With all eight required library modules and the focused regression drivers
  compiled at `-O3`, the new shifted and existing unshifted complex eigs
  native tests passed. The sample program also compiled and ran.
- Original unshifted `eigs_complex` passed an additional 35 randomized
  checked and 30 optimized NumPy comparisons after the optional `sigma`
  integration (6 and 5 SciPy checks, respectively).
- The **complete optimized standalone sparse regression suite was not
  completed**: compilation of the large pre-existing `sparse_linalg.f90`
  exceeded the execution time limit before tests started. Its status must
  not be inferred from the targeted `-O3` results.
- The project-wide dependency-aware `fpm test` was not run because FPM is
  not installed in this environment.

## Scope and numerical limitations

This is a bounded, explicit-restart Arnoldi method using an unblocked,
row-pivoted sparse LU factorization, not ARPACK's implicitly restarted
method or SuperLU. Close shifts can make `A-sigma*I` ill-conditioned;
exactly singular shifts report `-3`. The selection and stopping rules
verify calculated eigenpairs, but cannot guarantee discovery of all
strongly clustered or defective modes in a truncated subspace.

No external source code was copied, and no new Fortran dependencies were
added. NumPy and SciPy are test-only references.
