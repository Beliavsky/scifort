# 2026-11-05 sparse BiCG, CGS, and Kronecker-sum validation

## Scope and provenance

This update adds the standard BiCG and CGS recurrences independently from
published numerical specifications (Barrett et al., SIAM, 1994), plus a
sorted-CSR Kronecker-sum merge implementing SciPy's mathematical definition.
No SciPy source code or third-party solver implementation was imported.

## Build and comparison

GNU Fortran 14.2.0, Fortran 2018, source-compatible with SciFort's existing
`scifort_kinds`, `scifort_sparse`, `scifort_optimize_lsmr`, and
`scifort_sparse_linalg` modules. Reference computations use the NumPy and
SciPy packages in the test environment (SciPy 1.17.0).

```bash
python tools/sparse_reference/compare_sparse_bicg_kronsum_scipy.py --mode both
python tools/sparse_reference/run_native_sparse.py --mode both
```

The standalone reference driver compares 180 randomly generated sparse
Kronecker sums per build against `scipy.sparse.kronsum`, including square
factors of independently varied orders and many implicit or explicit zeros.
For 120 randomly generated nonsymmetric diagonally dominant real systems per
build, it compares Fortran BiCG/CGS with SciPy BiCG/CGS and verifies true
unpreconditioned residuals, nonzero starting guesses and Jacobi preconditioning.

The native tests additionally exercise exact cancellation with stored zeros,
empty matrices, invalid rectangular factors and integer-dimension overflow,
zero RHS, exact starting guesses, early iteration caps, singular breakdown,
missing Jacobi diagonals, NaN and invalid tolerances, and CGS ILU(0).

## Observed results in the supplied snapshot

With GNU Fortran 14.2.0, Python 3.13, and SciPy 1.17.0:

- `compare_sparse_bicg_kronsum_scipy.py --mode checked`: 180/180
  Kronecker sums and 120/120 solver cases passed.
- `compare_sparse_bicg_kronsum_scipy.py --mode optimized`: 180/180
  Kronecker sums and 120/120 solver cases passed.
- The standalone native regression runner passed all 28 checked tests.
  On its optimized run, 27 of 28 passed before the overall tool time limit;
  the remaining `test_sparse_formats_norm` also passed when compiled and run
  separately at `-O3` with the same module sources. No numerical test failure
  was observed.
- The repository's full `fpm test` could not be run because the `fpm`
  executable was not installed in this environment. The standalone sparse
  tests need no network or external FPM dependencies.

## Known limitations

The Fortran solvers operate on canonical binary64 real CSR matrices and do
not support arbitrary matrix-free operators or complex arithmetic. BiCG uses
Jacobi rather than ILU because a correct general ILU transpose solve is not
currently part of the `ilu0_factor` public API. CGS can converge irregularly
or break down, as with the standard method; always inspect `info` and actual
residuals. `csr_kronsum` returns canonical CSR without selecting formats, and
cannot represent dimensions or nnz exceeding default-integer storage limits.

No changes to FPM dependencies or third-party license texts were required.
