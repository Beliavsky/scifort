# RQ factorization validation — 2026-11-09

New files:

- `src/linalg_rq.f90`: real and complex `rq` / `rq_r` plus one-batch overloads.
- `test/test_linalg_rq.f90`: native tests for all shapes `0..14` in both
  dimensions, five matrix variants, both modes, and both numeric types,
  plus rectangular and zero-size batches and invalid-input tests.
- `reference/linalg_rq_driver.f90`: machine-readable Fortran reference output.
- `reference/check_linalg_rq.py`: seeded Python/SciPy independent comparisons.
- `example/linalg_rq_example.f90`: runnable real/complex example.

Environment: GNU Fortran 14.2.0; NumPy 2.3.5; SciPy 1.17.0.

| Test | `-O0 -g -fcheck=all` | `-O3` |
|---|---:|---:|
| Independent SciPy comparison cases | 1,674 passed | 1,674 passed |
| Native RQ cases | 5,382 passed | 5,382 passed |
| Existing QR native regression | passed | passed |
| Existing QR updates regression | 1,872 passed | 1,872 passed |
| Existing QR multiplication regression | 1,296 passed | 1,296 passed |

Maximum scaled reconstruction residual: `1.1575e-15`.
Maximum Q row-orthonormality error: `1.1102e-15`.
Maximum scaled difference in `abs(R)` against SciPy on full-rank cases:
`4.4418e-15`. The R matrix is compared by absolute entry values because
RQ decompositions differ by legitimate row phase factors. Rank-deficient
cases are validated using reconstruction, orthonormality and structural
zeros instead of assuming a unique R.

The SciPy suite contains 558 cases of each of full, economic and R-only
modes, including real/complex matrices, both empty-dimension shapes,
zero matrices, rank-deficient examples, and scales `1e-50` to `1e50`.

For reproducibility, from the repository root:

```sh
mkdir -p /tmp/rq-build
cd /tmp/rq-build
gfortran -std=f2018 -O0 -g -fcheck=all -c \
    /path/to/scifort/src/kinds.f90 \
    /path/to/scifort/src/linalg.f90 \
    /path/to/scifort/src/linalg_svd.f90 \
    /path/to/scifort/src/linalg_qr.f90 \
    /path/to/scifort/src/linalg_rq.f90
gfortran -std=f2018 -O0 -g -fcheck=all -c \
    /path/to/scifort/test/test_linalg_rq.f90 \
    /path/to/scifort/reference/linalg_rq_driver.f90
gfortran -o test_linalg_rq kinds.o linalg.o linalg_svd.o linalg_qr.o \
    linalg_rq.o test_linalg_rq.o
./test_linalg_rq
gfortran -o linalg_rq_driver kinds.o linalg.o linalg_svd.o linalg_qr.o \
    linalg_rq.o linalg_rq_driver.o
python /path/to/scifort/reference/check_linalg_rq.py ./linalg_rq_driver
```

The same files were compiled and tested separately with `-O3`. All three
previous QR regression executables were compiled and run in both builds.

**Scope:** This is a direct compilation of SciFort's core QR dependencies,
not a build of all ~3,000 project files. FPM was unavailable in this
runtime, so the full `fpm test` suite and Windows Intel Fortran builds have
not been verified.
