# Interpolative SVD — validation record (2026-11-09)

Environment: GNU Fortran 14.2.0; SciPy 1.17.0; NumPy 2.3.5.

| Test | Checked `-O0 -fcheck=all` | Optimized `-O3 -march=native` |
|---|---:|---:|
| Native matrix and API checks | 32,770 passed | 32,770 passed |
| Independent numerical reference cases | 1,944 passed | 1,944 passed |
| Direct deterministic SciPy SVD cases | 814 passed | 814 passed |
| Adaptive tolerance cases (NumPy reconstruction) | 592 passed | 592 passed |
| Zero rank / empty dimension cases | 538 passed | 538 passed |
| Existing interpolative decomposition test | Passed | Passed |
| Existing QR, SVD, full SVD, SVD condition tests | Passed | Passed |

Additional observed numerical results:

- Maximum absolute factor orthogonality error: `2.887e-15` in both builds.
- Maximum scaled singular-value difference versus the deterministic SciPy
  ID-SVD reference: `1.136e-15` (checked), `1.152e-15` (optimized).
- Maximum relative reconstruction error across all reference cases: 1.0,
  expected for the intentionally zero-rank approximations of nonzero matrices.
- Full-rank and tolerance-based reconstructions are explicitly checked against
  the requested errors; fixed-rank approximations are compared to SciPy's
  deterministic ID-SVD for approximation quality.
- The executable example produced relative errors `3.055e-16` (real rank-2)
  and `2.286e-16` (complex adaptive rank-2).

The 1,944 cases in each build are **not all direct SciPy SVD comparisons**:
814 directly compare with `scipy.linalg.interpolative.svd(..., rand=False)`;
the remaining 1,130 validate the adaptive/degenerate interfaces against
NumPy and the mathematical ID-SVD identities. Eigenvector phases/signs and
ID pivot choices are not required to match SciPy.

## Reproduce

Run from a fresh directory, with `ROOT` set to the extracted SciFort root:

```sh
mkdir -p /tmp/scifort-id-svd-check
cd /tmp/scifort-id-svd-check
ROOT=/path/to/scifort-main

gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -c -J . -I . \
    "$ROOT/src/kinds.f90" "$ROOT/src/linalg.f90" \
    "$ROOT/src/linalg_svd.f90" "$ROOT/src/linalg_qr.f90" \
    "$ROOT/src/linalg_interpolative.f90" \
    "$ROOT/src/linalg_interpolative_svd.f90"

gfortran -std=f2018 -O0 -g -fcheck=all -J . -I . \
    "$ROOT/test/test_linalg_interpolative_svd.f90" ./*.o -o id_svd_native
./id_svd_native

gfortran -std=f2018 -O0 -g -fcheck=all -J . -I . \
    "$ROOT/validation/reference_drivers/driver_linalg_interpolative_svd.f90" ./*.o -o id_svd_driver
python "$ROOT/test/compare_linalg_interpolative_svd_scipy.py" ./id_svd_driver

gfortran -std=f2018 -O0 -g -fcheck=all -J . -I . \
    "$ROOT/example/linalg_interpolative_svd_example.f90" ./*.o -o id_svd_example
./id_svd_example
```

Repeat from a **new build directory** with `-O3 -march=native` instead of
`-O0 -g -fcheck=all -fbacktrace` for optimized tests. The existing QR/SVD
regressions can be compiled against the same six modules with their
`test_linalg_interpolative.f90`, `test_linalg_qr.f90`, `test_linalg_svd.f90`,
`test_linalg_svd_full.f90`, and `test_linalg_svd_condition.f90` programs.

The entire `fpm test` suite was not run because FPM was not installed in
this environment. No speedup versus direct SVD is claimed.
