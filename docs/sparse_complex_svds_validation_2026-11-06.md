# Complex sparse SVD validation — 2026-11-06

## Independent numerical references

The `compare_complex_svds_scipy.py` harness compiles the new Fortran module
with GNU Fortran and compares all returned complex singular triplets with
dense NumPy SVD, including both true residual equations and orthogonality.
It also cross-checks selected `LM` singular values with SciPy's sparse
`svds` (ARPACK backend). The randomized data cover tall/wide matrices,
CSR/CSC formats, duplicated unsorted COO entries, explicitly stored zeros,
complex phases, nonzero complex starting vectors, rank deficiencies,
repeated singular values, zero matrices, and genuinely partial searches.

Reproduce from the project root:

```bash
python tools/sparse_reference/compare_complex_svds_scipy.py --mode both --cases 140 --partial 16 --seed 20261106
python tools/sparse_reference/compare_complex_svds_scipy.py --mode both --cases 150 --partial 12 --seed 20261107
python tools/sparse_reference/run_native_sparse.py --mode checked
```

Observed randomized results:

| Run | Checked | Fully optimized (`-O3`) | Additional SciPy checks per build |
| --- | ---: | ---: | ---: |
| Seed 20261106, 140 + 16 cases | 156 passed | 156 passed | 24 |
| Seed 20261107, 150 + 12 cases | 162 passed | 162 passed | 22 |
| **Total** | **318 passed** | **318 passed** | **46** |

The independent comparisons verify singular values against dense NumPy
SVD and compute `A@v - s*u` and `A^H@u - s*v` directly in Python. In
the final two-seed rerun, the largest two-sided residual was below
`4.60e-12`, within the configured relative convergence tolerance.

The native regression test includes a **4,000-by-120** complex CSR example
with 240 stored entries and a 28-vector partial basis. It also checks
nonconvergence reporting on a deliberately inadequate one-vector basis,
complex transpose handling, zero singular values, and invalid arguments.
The focused native test passed with both all-checking flags and a fully
optimized `-O3` build.

## Full-suite status

The complete standalone sparse regression runner, including the new SVD
test, passed **49/49 tests** in bounds-checked mode. The attempted fully
optimized full-suite runner exceeded the execution limit *while compiling*
the large pre-existing `src/sparse_linalg.f90` module, before starting any
native tests. This is not an observed numerical regression. The new
complex SVD and its dependencies were separately compiled/tested at `-O3`,
as documented above. The full project-wide `fpm test` was not run because
FPM is unavailable in this environment.

## Differences from SciPy and known limits

A complex SciPy PROPACK smallest-value reference returned zero in a case
where dense NumPy SVD gave a nonzero value. Consequently, dense NumPy SVD
and direct residuals are the independent accuracy checks; selected SciPy
`LM` checks use ARPACK. The new routine intentionally requires a complete
smaller-side trial basis for `which='SM'`. For partial `LM` searches,
converged Ritz vectors do not strictly certify that every globally largest
singular value has been found. See the API guide for other differences.
