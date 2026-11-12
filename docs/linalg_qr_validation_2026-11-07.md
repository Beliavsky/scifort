# Validation of dense Householder QR and least-squares routines (2026-11-07)

## Environment and builds

- GNU Fortran 14.2.0, `-std=f2018 -O0 -g -fcheck=all -fbacktrace`
  and `-std=f2018 -O3`.
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0.
- Standalone builds now use `src/kinds.f90`, `src/linalg.f90`,
  `src/linalg_svd.f90` and `src/linalg_qr.f90` (no numerical library dependency).
- Full `fpm test` could not be run because `fpm` is not installed in this
  container. The new native tests and the focused previous-regression tests
  listed below were compiled and run independently.

## Independent SciPy comparisons

Commands (each in both checked and optimized modes):

```sh
python tools/linalg_reference/run_qr_checks.py --mode both --cases 260 --max-order 60 --seed 431
python tools/linalg_reference/run_qr_checks.py --mode both --cases 100 --max-order 55 --seed 495
python tools/linalg_reference/run_qr_checks.py --mode both --cases 95 --max-order 85 --seed 666 --stress
```

All 455 cases in each compiler mode passed (910 comparisons in total).
The `--stress` set includes matrices up to order 180. The runner constructs
random real and complex tall, wide and square systems; some column-scaled
systems force pivoting. It exercises reduced and complete Q, both pivoted
and unpivoted QR, multiple right-hand sides, full-rank least squares,
minimum-norm underdetermined solutions, and additional matrix families with
singular values ranging over a factor of `1e8`. Input scaling includes
`1e-80` and `1e80` factors. Returned solutions are compared with
`scipy.linalg.lstsq` and QR factors checked against the matrix identities.

Maxima over all runs:

| Quantity | Maximum |
| --- | ---: |
| Relative `||Q*R - A(:,pivots)||_F / ||A||_F` | `8.91e-16` |
| Scaled `||Q^H*Q - I||_F` | `4.45e-16` |
| Relative solution difference from SciPy | `1.28e-8` |
| Scaled least-squares optimality/equation residual | `3.79e-16` |

The largest forward solution difference occurred in deliberately
ill-conditioned input; ordinary cases use a `2e-9` comparison threshold,
while the `1e8` condition-number subset has a `3e-6` allowance. Good QR
reconstruction and small residuals do not ensure small forward errors for
ill-conditioned systems.

`test/test_linalg_qr.f90` additionally tests compact factor unpacking, reuse
of the same QR factors for multiple RHSs, correct result types and shapes,
real/complex exact overdetermined and minimum-norm wide solutions, both Q
shapes, and explicit failure paths. Tests reject duplicate permutation indices, rank-deficient systems via
explicit `method='qr'`, invalid `rcond`, incompatible `tau` sizes,
and nonfinite input. The Python test driver also checks `qr_solve` against
`lstsq` after pivoted QR in full-column-rank cases.

## Regression tests

```sh
python tools/linalg_reference/run_dense_checks.py --mode both --cases 80 --max-order 40
python tools/linalg_reference/run_dense_cholesky_checks.py --mode both --cases 80 --max-order 40
python tools/linalg_reference/run_banded_lu_checks.py --mode both --cases 45 --max-order 40
```

Both builds passed the existing dense triangular/LU and dense Cholesky
native tests and 80 independent comparisons per build for each. The banded
LU suite passed 270 SciPy comparisons per build, plus its native
50,000-by-50,000 banded test. None of those solvers was modified by this
release.

## Known limitations

- Rank-deficient least-squares minimum-norm solutions are not implemented;
  the returned status is `linalg_status_singular` with a rank estimate.
- Numerical rank estimation uses pivoted `R` diagonal magnitudes rather
  than singular values, and may differ from SciPy's SVD-based rank.
- Unblocked QR and repeated trailing-column norms are less efficient than
  tuned block LAPACK algorithms on very large matrices.
- The implementation uses finite-arithmetic guards and reports failure for
  nonfinite/overflowed intermediates, rather than implementing an
  exceptional-scaling QR method.
