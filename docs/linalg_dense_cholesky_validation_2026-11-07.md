# Dense Cholesky / reusable positive-definite solve validation (2026-11-07)

This update extends `scifort_linalg_dense` with original Fortran generics
`cholesky`, `cho_factor`, and `cho_solve` for binary64 real symmetric or
complex Hermitian positive-definite matrices. The factorization reads only
its selected triangle, zeros the unused output triangle and leaves input
arrays intact. Both vector and multiple-RHS solutions are supported. No
additional numerical libraries are needed.

## Toolchain

- GNU Fortran 14.2.0, with `-std=f2018`.
- Checked build: `-O0 -g -fcheck=all -fbacktrace`.
- Optimized build: `-O3`.
- Python 3.13, NumPy and SciPy 1.17.0 as numerical references.
- The full `fpm test` suite was not run because `fpm` was unavailable.

## Commands and results

```sh
python tools/linalg_reference/run_dense_cholesky_checks.py --cases 450 --seed 390 --max-order 70 --mode both
python tools/linalg_reference/run_dense_cholesky_checks.py --cases 130 --seed 391 --max-order 85 --mode both --stress --tolerance 1e-7
```

Each command first compiles/runs native Fortran assertions, then compares the
factor and two independently solved right-hand side arrays with SciPy. Cases
alternate between real/complex and upper/lower triangular factorization.
The ignored input triangle is deliberately overwritten, and some complex
cases have imaginary diagonal components, which LAPACK/SciPy ignore.

| Test set | Checked | `-O3` | Worst solution difference | Worst factor difference | Worst relative reconstruction | Worst normwise backward residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Seed 390 | 450/450 | 450/450 | 1.322e-15 | 1.880e-16 | 2.218e-16 | 1.508e-16 |
| Stress seed 391 | 130/130 | 130/130 | 2.985e-10 | 1.921e-13 | 2.950e-16 | 1.231e-16 |
| **Total** | **580/580** | **580/580** | | | | |

The comparison measures the norm of the solution delta relative to
`max(1, ||X_SciPy||_F)`, the factor delta relative to the SciPy factor norm,
`||F^H F-A||_F / ||A||_F` (upper) or `||FF^H-A||_F / ||A||_F` (lower), and
the normwise backward residual
`||AX-B||_F / (||A||_F*||X||_F + ||B||_F)`.
The stress suite includes positive-definite matrices up to order 192 with
condition numbers approaching 1e7, tiny and large uniform scalings, and
complex data. The larger forward differences on ill-conditioned cases are
consistent with small backward residuals.

Native tests cover invalid/non-square factor shapes, mismatched output shapes,
empty RHS-column arrays, positive semidefinite and indefinite failures,
nonfinite input rejection, exact singular stored factors, and clean factor
reconstruction. The runnable example passed in both checked and optimized
builds.

## Other regression runs

```sh
python tools/linalg_reference/run_dense_checks.py --mode both --cases 100 --seed 221 --max-order 35
python tools/linalg_reference/run_banded_lu_checks.py --mode both --cases 65 --seed 172 --max-order 30
python tools/linalg_reference/run_matrix_functions_checks.py --mode both --cases 35 --seed 218
python tools/linalg_reference/run_logarithm_checks.py --mode both --cases 40 --seed 218 --max-order 8
```

All four commands passed in checked and optimized configurations. The banded
LU test includes its native order-50,000 test. An earlier combined shell
command with a larger dense-LU regression count exceeded its execution time
limit after the Cholesky suites had passed; the above separate dense-LU run
completed successfully.

## Numerical limitations

This is an unblocked Cholesky implementation. It does not offer scaling,
iterative refinement, pivoted semidefinite factorization, condition estimates,
error bounds, BLAS-backed optimization, or batch dimensions. The selected
input triangle defines the Hermitian matrix; no check is performed for
Hermitian consistency in the ignored triangle. As in LAPACK, complex
diagonal imaginary parts are ignored. Finite entries in the entire input
matrix are required. Numerical differences on ill-conditioned matrices are
not a sign of a formal forward-error guarantee.
