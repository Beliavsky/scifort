# `scipy.linalg.pinvh` independent validation — 2026-11-09

## Environment

- GNU Fortran 14.2.0, checked `-std=f2008 -O0 -g -fcheck=all`;
- GNU Fortran 14.2.0, optimized `-std=f2008 -O3 -march=native`;
- Python 3.13, SciPy 1.17.0, NumPy 2.3.5.

## Results

| Test | Checked | Optimized |
|---|---:|---:|
| Independent SciPy comparison matrices | 2,240 passed | 2,240 passed |
| Native SciFort cases | 335 passed | 335 passed |
| Existing dense-eigen regression | passed | passed |
| Maximum relative pseudoinverse discrepancy vs SciPy | 1.087e-11 | 1.956e-12 |

The 2,240 cases per build include real/complex square matrices of orders
0--13, full-rank and rank-deficient positive/negative/indefinite matrices,
unitary rotations, different coefficient scales, one triangle deliberately
inconsistent with the other, and both `lower` settings. The eigenvalue
threshold is varied using absolute and relative cutoffs. The comparison
requires exactly matching effective rank and a relative pseudoinverse
error below 3e-7 for every case (all results were substantially better).

The largest normalized Moore--Penrose `A A+ A - A` residual in this reference
set is 0.01, which is *expected*: some tests intentionally discard nonzero
eigenvalues by requesting a relative threshold of 0.02. Reconstruction is
not exact under intentional truncation, and this residual should not be
misinterpreted as eigensolver accuracy.

Native tests additionally check nonfinite values, invalid shapes, output
cleanup on failure, complex conjugation, batch behavior, and proper treatment
of negative eigenvalues.

## Reproduction

From the root of the extracted SciFort project:

```sh
mkdir -p /tmp/scifort-pinvh-test
cd /tmp/scifort-pinvh-test
src=/path/to/extracted/scifort-main/src
root=/path/to/extracted/scifort-main

gfortran -std=f2008 -O0 -g -fcheck=all -J. -I. -c \
  "$src/kinds.f90" "$src/linalg.f90" "$src/linalg_dense.f90" \
  "$src/linalg_hessenberg.f90" "$src/linalg_banded.f90" \
  "$src/linalg_banded_tridiagonal.f90" "$src/linalg_banded_eigen.f90" \
  "$src/linalg_dense_eigen.f90" "$src/linalg_pinvh.f90"
gfortran -std=f2008 -O0 -g -fcheck=all -I. \
  "$root/test/test_linalg_pinvh.f90" *.o -o test_pinvh
./test_pinvh
gfortran -std=f2008 -O0 -g -fcheck=all -I. \
  "$root/reference/driver_linalg_pinvh.f90" *.o -o driver_pinvh
python "$root/reference/check_linalg_pinvh.py" ./driver_pinvh
```

Recompile the dependencies and both executables with `-O3 -march=native`
to repeat the optimized run. The full `fpm test` suite was not run because
FPM was not installed in the validation environment.
