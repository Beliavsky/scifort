# Dense SVD, pseudoinverse and rank-deficient least squares — validation (2026-11-07)

## Implementation and environment

- Fortran sources: `src/linalg_svd.f90` and the updated `src/linalg_qr.f90`.
- GNU Fortran 14.2.0, with `-std=f2018 -O0 -g -fcheck=all -fbacktrace`
  and separately `-std=f2018 -O3`.
- Python 3.13, NumPy and SciPy used as independent numerical references.
- Builds use `src/kinds.f90`, `src/linalg.f90`, `src/linalg_svd.f90`,
  and `src/linalg_qr.f90` in that order. No BLAS/LAPACK dependency for
  the Fortran implementation.

## Independent SciPy comparisons

The following commands passed in both checked and optimized configurations:

```sh
python tools/linalg_reference/run_svd_checks.py --mode both --cases 260 --max-order 32
python tools/linalg_reference/run_svd_checks.py --mode both --cases 300 --max-order 45 --seed 920 --stress
python tools/linalg_reference/run_svd_checks.py --mode both --cases 400 --max-order 55 --seed 1111 --stress
```

These **960 distinct test matrices in each mode** constitute **1,920
SciPy comparison cases**. The generators include real/complex matrices,
tall/wide/square shapes, exact and near numerical rank deficiency, repeated
columns, zero matrices, one and multiple RHSs, `nrhs=0`, orthogonal rotated
singular spectra down to `1e-10`, and input scaling by `1e-100`/`1e100`.
**511 cases per mode** were rank-deficient under the SciPy SVD cutoff
(1,022 checked/optimized rank-deficient comparisons).

Main comparisons use SciPy `svdvals`, `pinv`, and `lstsq(lapack_driver='gelsd')`
with the same explicit `rcond`. SciPy's LAPACK least-squares method cannot
handle `nrhs=0` in this environment, so the reference solution for those
cases is empty and its numerical rank is determined from the SciPy singular
values. For `svd`, test `U diag(s) Vh` reconstruction, both factor
orthogonalities, and singular values; for `pinv` compare the pseudoinverse;
for `lstsq` compare minimum-norm coefficients and status.

Worst observed results across the three commands:

| Measure | Maximum |
| --- | ---: |
| Relative SVD reconstruction error | `4.90e-15` |
| Unitarity error, dimension-normalized | `6.66e-16` |
| Relative singular-value difference | `1.33e-15` |
| Relative pseudoinverse difference | `3.28e-6` |
| Relative SVD least-squares difference | `5.55e-6` |
| Relative automatic QR/SVD least-squares difference | `8.07e-6` |

The largest forward differences involve matrices with tiny retained singular
values and do not imply reconstruction errors of the same size. Comparison
tolerances are `2e-7` in ordinary cases and `2e-5` in the deliberately
ill-conditioned spectral families; the automatic QR pathway uses an
additional `3e-5` cap for these comparisons.

## Native tests and previous regressions

`test/test_linalg_svd.f90` passed independently in the checked and optimized
builds. It exercises analytic rank deficiency, multiple right-hand sides,
real/complex economy SVD, Moore-Penrose identities, minimum-norm solutions,
custom rank cutoffs, zero matrices, malformed dimensions, negative cutoffs,
invalid method and sweep options, and nonfinite inputs. The runnable
`example/linalg_svd_example.f90` passed in the checked build.

The previous QR regression runner passed **80 cases per build** for its
original full-rank tests (160 comparisons), alongside the existing
`test/test_linalg_qr.f90` in both configurations. A separate QR run of
160 cases per build also passed. An additional larger QR stress run passed
its checked native test before reaching the command's time limit; it is
**not** counted as a completed comparison run.

Other prior regression checks all passed in both configurations:

```sh
python tools/linalg_reference/run_dense_checks.py --mode both --cases 30 --max-order 30 --seed 218
python tools/linalg_reference/run_dense_cholesky_checks.py --mode both --cases 30 --max-order 30 --seed 390
python tools/linalg_reference/run_banded_lu_checks.py --mode both --cases 25 --max-order 30 --seed 171
```

These include the native 50,000-order banded LU regression and 300 banded LU
SciPy comparisons (150 per build). These modules are unchanged in this release.

## Limitations

The one-sided Jacobi implementation is unblocked and less efficient than
LAPACK SVD for large dense matrices. It computes economy factors only.
Numerical rank depends on a relative cutoff and is ambiguous for singular
values close to roundoff. For `lstsq(method='auto')`, the SVD fallback only
occurs when the QR rank estimate declares deficiency; a QR rank estimate
of full rank may not agree with an SVD estimate near the cutoff.

The full `fpm test` suite was not run because `fpm` is not available in this
runtime. The standalone tests and focused regressions listed above did run.
