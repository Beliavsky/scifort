# Matrix-free randomized interpolative SVD — validation (2026-11-10)

Baseline: `scifort-linalg-matrix-free-interpolative-2026-11-10.zip`
(3,075 files). New module:
`src/linalg_interpolative_svd_operator.f90`.

## Environment and commands

- GNU Fortran 14.2.0 (Linux), checked `-std=f2008 -O0 -g -fcheck=all -fbacktrace`.
- GNU Fortran 14.2.0, optimized `-std=f2008 -O3`.
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0.
- FPM unavailable: full `fpm test` suite **not** run.

Compile modules in dependency order:

```sh
mkdir -p /tmp/sf-build
cd /tmp/sf-build
for file in kinds linalg linalg_svd linalg_qr linalg_interpolative \
    linalg_interpolative_operator linalg_interpolative_svd_operator; do
    gfortran -std=f2008 -O0 -g -fcheck=all -J. -I. -c "/path/to/scifort/src/${file}.f90"
done
gfortran -std=f2008 -O0 -g -fcheck=all -I. \
    /path/to/scifort/test/test_linalg_interpolative_svd_operator.f90 \
    *.o -o test_linalg_interpolative_svd_operator
./test_linalg_interpolative_svd_operator
gfortran -std=f2008 -O0 -g -fcheck=all -I. \
    /path/to/scifort/reference/linalg_interpolative_svd_operator_driver.f90 \
    *.o -o driver_linalg_interpolative_svd_operator
python /path/to/scifort/test/compare_linalg_interpolative_svd_operator_scipy.py \
    ./driver_linalg_interpolative_svd_operator
```

Repeat with `-O3` instead of `-O0 -g -fcheck=all` for the optimized build.

## Results

| Verification | Checked | Optimized |
|---|---:|---:|
| New native checks | 41,746 | 41,746 |
| Independent NumPy/SciPy harness cases | 536 | 536 |
| SciPy `LinearOperator` SVD returns validated | 140 | 140 |
| SciPy `LinearOperator` SVD calls raising exceptions | 244 | 244 |
| Maximum factor orthogonality error | 3.55e-15 | 2.89e-15 |
| Maximum fixed-rank Frobenius error / optimal SVD error | 1.70604 | 1.70604 |
| Max singular-value mismatch for reconstructed matrix (scaled) | 3.55e-15 | 3.55e-15 |
| Maximum adaptive relative residual | 0.0951475 | 0.0951475 |
| Existing matrix-free ID native checks | 42,911 passed | 42,911 passed |
| Existing deterministic ID native checks | 14,944 passed | 14,944 passed |
| Existing randomized ID native checks | 37,449 passed | 37,449 passed |
| Existing deterministic ID-SVD native checks | 32,770 passed | 32,770 passed |
| Existing randomized ID-SVD native checks | 50,408 passed | 50,408 passed |
| QR and SVD native regression tests | passed | passed |
| Executable implicit rank-three example | passed | passed |

The harness tests fixed ranks `0`, `1`, `3` (where allowed), and full rank;
adaptive tolerances `0.1` and `1e-5`; real/complex rectangular matrices,
zero and low-rank matrices, zero dimensions, coefficients scaled by
`1e-80` or `1e80`, independently calculated optimal SVD truncation errors,
explicit orthogonality, singular values, callback counts and seeded
repeatability. Fixed-rank error ratios exclude cases where the optimal
rank-k error is zero or near floating-point noise.

**SciPy reference qualification:** 244 of the 384 attempted SciPy
`LinearOperator` calls per build raised `ValueError` inside SciPy 1.17.0,
with examples like `could not broadcast input array from shape (3,) into
shape (5,)`. They are counted as reference failures, not as successful
comparisons or SciFort failures. All 536 SciFort cases independently passed
numerical checks against NumPy. SciPy random outputs are not expected to
have the same signs/phases, rank in tolerance mode, or sketch pivots.

The Frobenius-norm residual estimate for adaptive calls is **probabilistic**;
no spectral-norm accuracy guarantee or performance speedup is claimed.
The new module, tests and examples use existing SciFort routines only and
introduce no additional runtime dependencies or changes to existing
Fortran modules.
