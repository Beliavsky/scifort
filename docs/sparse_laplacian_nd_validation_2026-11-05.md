# SciFort sparse Laplacian and structural-query verification (2026-11-05)

## Environment

- GNU Fortran 14.2.0 (`gfortran -std=f2018`)
- Python with NumPy and SciPy 1.17.0 (test tools only)
- Bounds-checked mode: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`
- Optimized mode: `-O3`
- All native tests are compiled directly without external BLAS/LAPACK or FPM.

## Randomized SciPy reference checks

Run:

```sh
python tools/sparse_reference/compare_sparse_laplacian_nd_scipy.py \
    --mode both --laplacian-cases 180 --structure-cases 300
```

| Randomized group | Checked | Optimized |
| --- | ---: | ---: |
| N-dimensional Laplacian cases | 180/180 | 180/180 |
| CSR/CSC sparse structural-query cases | 300/300 | 300/300 |

For each Laplacian comparison, verify the matrix-free matvec, canonical CSR
assembly and sparse/dense equivalence, analytical eigenvalues, eigenpair
residuals, eigenvector orthogonality, and CSR/CSC structural metadata.
Test 1–4 grid axes with dimension sizes from 1 through 5 and all three
boundary conditions. Length-one axis cases use SciPy's matrix-free action
because its 1.17.0 sparse converter disagrees with its own matvec.
For truncated eigenvalue checks, compare to the full SciPy eigenspectrum's
largest `m` values. SciPy 1.17.0's optimized `eigenvalues(m)` occasionally
returns phantom modes for certain multidimensional shapes.

Structural queries are compared to SciPy for random **rectangular** CSR and
CSC matrices of up to 45 rows and columns, including explicit stored zeros,
duplicate input coordinates, upper/lower triangular matrices and empty
storage. The Fortran result for empty storage is `(0,0)` bandwidth; SciPy
1.17.0 raises on empty CSR when asked for `spbandwidth`, so this deliberate
extension is independently checked.

## Native regression tests

Run:

```sh
python tools/sparse_reference/run_native_sparse.py --mode both
```

- **Checked:** 31/31 standalone sparse native programs passed.
- **Optimized:** 31/31 standalone sparse native programs passed.
- The checked suite ran in one pass. The optimized suite was run in two
  segments (first nine programs, then the remaining 22 using
  `--mode optimized --start-at test_sparse_matrix`) to avoid exceeding the
  runner's invocation time limit; the same optimized compiler flags and
  source files were used for each segment.
- The new `test/test_sparse_laplacian_nd.f90` additionally checks invalid
  dimensions, constructor overflow, nonfinite matvec inputs, zero-column
  matmat, periodic doubling for two-point axes, singleton consistency,
  full/truncated eigenvectors, and malformed sparse input rejection.

## Scope and attribution

The implementation uses original MIT-licensed Fortran algorithms derived
from mathematical formulas for regular-grid second differences. No SciPy
implementation source code was imported. The new module and structural
queries use SciFort's existing canonical one-based CSR/CSC data types.

Full `fpm test` was **not** run: FPM is unavailable in this environment.
The standalone sparse suite is therefore not a claim that every SciFort
component was rebuilt or tested.
