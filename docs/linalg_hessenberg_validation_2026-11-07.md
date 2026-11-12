# Hessenberg reduction validation (November 7, 2026)

## Independent reference and build

- Tested using GNU Fortran and the installed SciPy 1.17.0 / NumPy.
- Checked compilation: `gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace`.
- Optimized compilation: `gfortran -std=f2018 -O3`.
- Native tests passed in both modes; an additional compile with
  `-ffpe-trap=invalid,zero,overflow` passed the native test.
- No LAPACK, BLAS, or third-party numerical source was added to the library.

## SciPy comparisons

The `run_hessenberg_checks.py` harness streams numerical matrices into a
standalone Fortran executable, then compares independently in Python using
`scipy.linalg.hessenberg(a,calc_q=True)`. Because the reduction is not
unique, tests check invariants and reconstruction instead of exact `H`/`Q`.

| Run | Count | Maximum relative reconstruction error | Maximum orthogonality error |
| --- | ---: | ---: | ---: |
| Seed 20261107, orders through 48, checked | 280 | 1.658e-15 | 4.010e-16 |
| Seed 20261107, orders through 48, `-O3` | 280 | 1.658e-15 | 4.010e-16 |
| Seed 20261108, orders through 96, checked | 480 | 2.472e-15 | 4.140e-16 |
| Seed 20261108, orders through 96, `-O3` | 480 | 2.472e-15 | 4.140e-16 |
| Seed 20261202, orders through 96, checked | 160 | 3.081e-15 | 3.289e-16 |
| Seed 20261202, orders through 96, `-O3` | 160 | 3.081e-15 | 3.289e-16 |
| Additional order-180 test, optimized | 40 | 3.130e-15 | 4.165e-16 |

The test data include empty matrices, real and complex dense matrices,
triangular and already-Hessenberg matrices, Hermitian or symmetric matrices,
rank-deficient outer products, nilpotent/Jordan-like matrices, repeated
eigenvalues, near-triangular matrices and scales up to `1e300` and down to
`1e-300`. The comparison checks exact zeros below the first subdiagonal.
Relative reconstruction errors normalize by the input magnitude before
computing the Frobenius norm to avoid overflow and underflow.

## Limitations

- No SciPy-style batch inputs or optional matrix balancing.
- `Q` can differ from SciPy by signs/phases; this is not a failure.
- An independent LAPACK comparison of factor *entries* is not meaningful;
  compare invariants and reconstruction.
- The full FPM suite is not available in this environment, so the new
  module and its neighboring regression suites are compiled standalone.
