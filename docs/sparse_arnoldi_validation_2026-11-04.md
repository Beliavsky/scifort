# Restarted Arnoldi eigenpair validation (2026-11-04)

## Environment

- GNU Fortran 14.2.0, `-std=f2018`
- Checked mode: `-O0 -g -fcheck=all -fbacktrace`
- Optimized mode: `-O3`
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0
- No external Fortran libraries or BLAS/LAPACK dependency needed by the tests.

## Commands

```sh
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_linalg_reference/compare_arnoldi_scipy.py --mode both --cases 120
python tools/sparse_linalg_reference/compare_arnoldi_scipy.py --mode both --cases 120 --seed 20261105
```

## Results

- All nine standalone sparse native regression programs passed in checked and
  optimized builds, including the new `test_sparse_arnoldi`.
- New independent reference suite: **240 randomized separated-spectrum cases
  per mode** across two seeds passed SciPy dense eigenvalue comparisons, with additional
  sparse-ARPACK comparisons on LM/LR/SR cases.
- Eight clustered-spectrum similarity-transform tests per mode passed,
  including three eigenvalues separated by only `1e-5` and narrow `ncv=8`.
- **90 out of 248 cases** required more than one Arnoldi cycle per build;
  the maximum was nine cycles.
- Native tests additionally cover a **1,200-row CSR** eigenvalue problem,
  explicit complex spectra, true residual checks, `locked_count`, restart
  accounting, invalid settings, and a zero matrix.
- Returned complex eigenvalues are matched to independent reference values
  without enforcing an arbitrary ordering for near ties; numerical error
  tolerances are documented in the reference script.

## Known limitations

The implementation is explicitly thick restarted with *soft* convergence
locking, not implicitly restarted Arnoldi (IRAM), and can struggle with
highly nonnormal or defective matrices. It uses SciFort's small projected
complex QR Schur eigensolver, which is not a substitute for LAPACK on
ill-conditioned, tightly clustered, or very large projected problems.
`SM` with a reduced search space and complex shift-invert are not yet
supported by `eigs_arnoldi` (but remain available in the original `eigs`).

The full `fpm test` command could not be run because FPM is unavailable in
this execution environment. The pre-existing executable-stack linker warning
from nested Fortran callbacks also remains.

Additional existing optimized SciPy comparison suites were run after this
change: `compare_eigs_scipy.py --cases 25` and
`compare_minres_scipy.py --cases 25`, both passed.

## Subsequent extension (2026-11-04)

Complex `sigma=` shift-invert is now available in `eigs_arnoldi` and
through `eigs(method='arnoldi')`. The inverse-operator Krylov projection
uses the existing complex GMRES kernel, then checks true residuals against
original CSR A. See `sparse_restarted_shift_validation_2026-11-04.md`.
