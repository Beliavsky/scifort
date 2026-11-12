# Complex sparse CSR/CSC and Krylov validation (2026-11-05)

Environment: GNU Fortran 14.2.0, NumPy 2.3.5, SciPy 1.17.0.
The new code is original MIT-licensed Fortran; neither SciPy nor NumPy is
required by SciFort at runtime.

## Native Fortran regression

`tools/sparse_reference/run_native_sparse.py` compiles the complete native
sparse collection and runs each executable. Results:

| Build | Regression programs | Result |
|---|---:|---|
| `-O0 -g -fcheck=all -fbacktrace` | 43/43 | Passed |
| `-O3` | 43/43 | Passed |

The focused `test_sparse_complex` passes independently in both modes,
including COO duplicates and explicitly stored zeros, dense/CSR/CSC and
N-D COO conversions, conjugation, Hermitian transpose, CSR/CSC matvec,
CSR matmat, Hermitian CG, nonsymmetric BiCGSTAB and GMRES, initial guesses,
Jacobi preconditioning, finite validation and invalid-parameter cases.
The documentation example also compiles and reports a true residual
norm of approximately `4.44e-16`.

## Randomized external references

`tools/sparse_reference/compare_sparse_complex_scipy.py` independently
builds complex sparse matrices using SciPy COO arrays, verifies products
and compares complex solves against NumPy `linalg.solve` and SciPy
`cg`, `bicgstab`, and `gmres`. Randomized cases contain unsorted COO
indices, split duplicate coordinates, explicit stored zeros, varying
starting guesses, Jacobi on/off, CSR/CSC, Hermitian positive-definite
systems and nonsymmetric complex systems.

| Build | Seeds | Cases | Result |
|---|---|---:|---|
| Checked | 20261105, 20261106, 89127 | 2,400 | Passed |
| Optimized | 20261105, 20261106, 89127 | 2,400 | Passed |

The largest Euclidean solution error in those cases was approximately
`7.58e-10`. Additional focused subsets (160 checked, 240 optimized,
with the last using seed 719234) also passed. The test script and
reference driver are included in the source archive.

## Not covered

Full project-wide `fpm test` was **not** run: FPM is not installed in
this environment. The native sparse runner does not represent all SciFort
modules or all external dependency configurations. Complex ILU, complex
sparse direct factorizations and nonbinary64 complex kinds are not present.
