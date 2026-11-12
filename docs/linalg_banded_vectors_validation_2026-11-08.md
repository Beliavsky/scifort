# Selected Hermitian band eigenvector validation (2026-11-08)

New functionality: `eig_banded` accepts `subset_by_index`, `subset_by_value`,
`num_selected`, `tol`, and `method='auto'|'inverse'|'dense'`, retaining old
positional arguments. Real and complex input, upper/lower SciPy-packed bands,
repeated diagonal spectra, invalid inputs, and output buffer capacities are
tested by `test/test_linalg_banded_vectors.f90`. The real symmetric path now
uses real banded LU/triangular solves (rather than complex LU) and includes
a native real-versus-complex eigenvector agreement check.

## Independent numerical checks

Run:

```sh
python tools/linalg_reference/run_banded_vectors_checks.py \
    --mode both --cases 280 --max-order 160 --seed 15151
```

GNU Fortran checked (`-O0 -g -fcheck=all` with floating-point traps) and
optimized (`-O3`) tests **each passed 280 randomized SciPy comparisons**.
A second independent seed passed 240 comparisons per build, including
matrices of order 256. Across the two seeds (520 comparisons per build),
the maximum normalized eigenvalue difference was `7.333e-15`, the maximum
eigenvector residual was `5.372e-11`, and the maximum normalized
orthogonality difference was `8.882e-16`. The random families include zero and repeated spectra,
well-separated and clustered eigenvalues, scalings up to `1e70` and down
to `1e-70`, and symmetric/Hermitian input stored in either triangle.
These measurements are from the stated test cases, not global bounds.

The checked and optimized native tests also passed a *6,000-order* real
symmetric tridiagonal band with three selected eigenvectors, using only
band and selected-vector workspace. Both builds passed 240 additional
selected band Sturm eigenvalue comparisons and 70 selected-tridiagonal
vector comparisons; the prior 65-case full-band eigenvalue comparison
suite also passed in each build.

## Known limitations

Sturm counting uses unpivoted Hermitian band LDL; the factorization may
become unreliable at near-singular shifts. Shifted band LU is pivoted but
still cannot guarantee converged inverse iteration for strongly clustered
or defective numerical eigenvalue approximations. For `method='inverse'`,
nonconvergence is reported explicitly; `method='auto'` falls back to dense
only for order <= 96. Full-spectrum eigenvectors still require dense
storage. `fpm test` has not been executed in this environment.

## Further validation of real LU specialization

Using GNU Fortran checked and optimized builds, both modes passed a new
160-case seed (`--max-order 96 --seed 25890`) and a 120-case seed
(`--max-order 64 --seed 20261108`). The independent 280-case and 240-case
SciPy batches listed above were repeated and passed in both modes after
introducing the real-only factorization path. Native tests additionally
compare real and complex factorization paths on identical symmetric data.

Regressions: 110 banded-Sturm and 70 full-band eigenvalue SciPy comparisons
passed per compiler mode; 55 dense-eigenvalue comparisons passed per mode.
Two larger combined regression commands reached the execution limit after
those completed tests, so the remaining suites in those runs are **not**
reported as passing. No performance speedup claim is made.
