# Sparse LOBPCG validation (2026-11-05)

## Scope

The original independent `scifort_sparse_linalg.lobpcg` implementation
supports real symmetric standard and generalized eigenproblems, Euclidean
or B-orthogonality, constraint directions, multiple initial vectors,
optional SPD CSR preconditioners, and Jacobi preconditioning. No external
BLAS/LAPACK is required.

The Fortran regression covers: largest/smallest eigenvalues; SPD mass
matrices; constraints; repeated eigenvalues; invalid inputs; zero budget;
singletons; orthonormality; and a 3,600-row CSR diagonal eigenproblem.

The randomized driver compares SciFort's eigenvalues with NumPy's
`eigvalsh` and SciPy's dense `eigh` (also in constrained coordinates), and
samples SciPy's sparse `lobpcg`. Independently compute and verify the
returned eigenvectors' true residuals, B-orthogonality and constraint
compliance. Since both methods can return different bases of repeated
eigenspaces, compare invariant quantities rather than vector entries.

The script admits legitimate nonconvergence: a **positive** Fortran status
must correspond to a residual exceeding the requested tolerance. Never
count that as successful convergence. In difficult cases the algorithm
returns the best residual block observed, which need not equal the dense
reference without convergence.

## Reproduce

```bash
python tools/sparse_reference/run_native_sparse.py --mode both
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
  tools/sparse_linalg_reference/compare_lobpcg_scipy.py --mode both --cases 60
```

For the native large CSR regression, use `test/test_sparse_lobpcg.f90`.

## Observed results

Checked (`-O0 -fcheck=all`) and optimized (`-O3`) GNU Fortran 14 builds:

| Check | Checked | Optimized |
|---|---:|---:|
| Native sparse regression programs | 20/20 passed | 20/20 passed |
| Randomized matrices, two seeds (60 each) | 120 | 120 |
| Converged at `1e-7` residual tolerance | 112 | 112 |
| Honestly reported nonconvergence | 8 | 8 |
| Accurate independent dense eigenvalue comparisons for converged cases | 112 | 112 |
| Additional direct SciPy `lobpcg` comparisons agreeing with dense reference | 18 | 18 |
| Independently verified orthogonality, constraints, residual values | 120 | 120 |

The maximum eigenvalue error among converged randomized examples was under
`1.2e-13` in each mode. Cases that did not achieve the requested tolerance
returned **positive** `info`, not false convergence. The independently
validated checks use both ordinary and generalized symmetric eigenproblems
and unconstrained and constrained solves. The standalone native regression
also passes for a 3,600-row sparse matrix and exactly repeated eigenvalues.

The complete FPM suite was not run because `fpm` is unavailable in the
validation environment. The pre-existing GNU executable-stack linker warning
for the LSMR callback code remains; it is unrelated to the new algorithm.

## Limitations

- `B` and the preconditioner action `M` are required to be SPD; the
  implementation checks symmetry of `B` by default, but does not prove SPD.
- The reported constrained residual removes the Lagrange multiplier term.
- The small dense Rayleigh--Ritz step uses the existing Jacobi kernel;
  large block sizes will be more expensive than SciPy's LAPACK backend.
- Only real CSR operators; no callback or complex generalized operators.
- A positive iteration status means the requested residual tolerance was
  not achieved, and is not masked by a good-looking Ritz eigenvalue.
