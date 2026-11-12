# Partial-spectrum Hermitian band eigenvalues: validation (2026-11-08)

## Algorithm

The new `eigvals_banded` `method='sturm'` route uses scaled Hermitian
banded LDL^H factorizations to count negative pivots in `A - lambda I`.
A Gershgorin interval and bisection isolate selected eigenvalues; no full
matrix or full eigenvector array is formed. Workspace is O(n*kd) and the
approximate time complexity is O(n*kd**2*k*log(width/tolerance)) for `k`
selected eigenvalues. Real matrices are handled through the same complex
banded kernel, without loss of mathematical generality.

The default `method='auto'` selects the banded path only when a subset is
requested. Full-spectrum calls retain the existing dense-backed implementation
for backward compatibility. `method='dense'` also supports selection by
filtering the already-computed full spectrum.

## Reproducible tests

Tests were run with GNU Fortran and locally installed SciPy 1.17.0,
using single-threaded BLAS for Python references.

```sh
python tools/linalg_reference/run_banded_sturm_checks.py --mode both \
    --cases 360 --max-order 128 --seed 2093
python tools/linalg_reference/run_banded_sturm_checks.py --mode both \
    --cases 240 --max-order 96 --seed 1092
```

Both runtime-checked and `-O3` builds passed the two independent seeds,
for **1,200 successful SciPy comparisons** in total. The maximum normalized
eigenvalue difference was `5.531e-14` for the first seed, and `1.073e-14`
for the second. The generators cover upper/lower storage, real/complex input,
full and selected spectra, tightly clustered and repeated eigenvalues,
zero/diagonal matrices, and scale factors `1e-150` and `1e150`.

The native test `test/test_linalg_banded_sturm.f90` is compiled and run by
this comparison harness before each randomized batch. It checks one-based
index and half-open value selection, multiple-eigenvalue cases, infinite
bounds, malformed inputs, insufficient capacity, unused storage, nonfinite
stored entries, explicit dense mode and an order-10000 block-banded matrix.
That large test selects only five eigenvalues, comparing them against a
separately computed order-three block reference.

The optimized native program's measured wall time was 0.08 seconds with
peak process RSS of about 2.7 MB on this particular environment, not a
general performance promise.

## Neighboring regressions

The following independent SciPy comparison suites passed in *both* modes:

| Suite | Comparisons across both modes |
|---|---:|
| Existing banded eigensolvers | 200 |
| Selected tridiagonal eigenvectors | 90 |
| Selected tridiagonal eigenvalues | 120 |
| Dense `eigh` | 90 |
| Generalized dense `eigh` | 80 |
| Dense `subset_by_value` | 90 |
| **Total** | **670** |

The native tests integrated into these regression scripts also passed.

## Limitations

The unpivoted banded LDL^H inertia algorithm is less robust than a pivoted
indefinite factorization near singular shifts. The implementation uses
regularized tiny pivots and nonfinite-growth checks, but inaccurate inertia
counts may still occur in challenging clusters. Endpoint membership for
eigenvalues very close to a selection bound can differ from SciPy.
Bisection is also expensive for large requested fractions of the spectrum,
especially when bandwidth is large. The `eig_banded` eigenvector algorithm
still expands bands to a dense matrix. A band-preserving bulge-chasing
reduction for eigenvectors and a more robust inertia factorization remain
future work.

The full `fpm test` suite was not executed because `fpm` was unavailable.
