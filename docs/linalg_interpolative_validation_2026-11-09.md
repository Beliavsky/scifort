# Deterministic interpolative decomposition validation — 2026-11-09

**Environment:** GNU Fortran 14.2.0; SciPy 1.17.0; NumPy 2.3.5.

| Test | Checked compiler (`-O0 -fcheck=all`) | Optimized compiler (`-O3 -march=native`) |
|---|---:|---:|
| Native real/complex ID and reconstruction checks | 14,944 pass | 14,944 pass |
| Independent SciPy reference ID cases | 2,016 pass | 2,016 pass |
| Of which, fixed-rank / adaptive-rank cases | 1,344 / 672 | 1,344 / 672 |
| ID-to-SVD singular-value comparisons (subset) | 288 pass | 288 pass |
| Previous QR / SVD / full SVD / SVD condition regression suites | 4 pass | 4 pass |
| Largest relative difference of SVD singular values | 1.078e-15 | 8.392e-16 |
| Largest relative interpolation coefficient difference on identical pivots | 1.556e-15 | 1.446e-15 |
| Matched deterministic SciPy pivot permutations | 634 / 1,344 | 632 / 1,344 |

Different pivot orders arise in tied or nearly tied columns; exact pivot
agreement is *not required* to produce the same interpolative approximation.
The maximum ratio of SciFort to SciPy relative reconstruction error in the
reference cases is approximately 1.0 (with a `1e-14` denominator floor).
Low fixed-rank approximations can have intentionally large absolute
truncation error; the tolerance-based test additionally verifies the
chosen relative Frobenius reconstruction error.

## Reproduction from project root

Run commands from a fresh directory to avoid mixing Fortran module files:

```sh
mkdir -p /tmp/scifort-id-check
cd /tmp/scifort-id-check
ROOT=/path/to/scifort-main

gfortran -c -std=f2018 -O0 -g -fcheck=all -fbacktrace \
  -J . -I . "$ROOT/src/kinds.f90" "$ROOT/src/linalg.f90" \
  "$ROOT/src/linalg_svd.f90" "$ROOT/src/linalg_qr.f90" \
  "$ROOT/src/linalg_interpolative.f90"
gfortran -std=f2018 -O0 -g -fcheck=all -J . -I . \
  "$ROOT/test/test_linalg_interpolative.f90" ./*.o -o id_native
./id_native

gfortran -std=f2018 -O0 -g -fcheck=all -J . -I . \
  "$ROOT/validation/reference_drivers/driver_linalg_interpolative.f90" ./*.o -o id_reference
python "$ROOT/test/compare_linalg_interpolative_scipy.py" ./id_reference

gfortran -std=f2018 -O0 -g -fcheck=all -J . -I . \
  "$ROOT/example/linalg_interpolative_example.f90" ./*.o -o id_example
./id_example
```

To reproduce the optimized results, compile each source and driver anew in
a clean directory, substituting `-O3 -march=native` for checked options.

The complete `fpm test` suite was **not** run in this environment because
FPM is unavailable. Only the named existing QR/SVD regression programs
were run. No performance benchmark is claimed.
