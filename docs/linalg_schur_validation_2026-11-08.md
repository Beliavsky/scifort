# Validation: public complex Schur and sorted Schur (2026-11-08)

New public `scifort_linalg_schur.schur` for real and complex matrices (complex
Schur form), including the built-in `lhp`/`rhp`/`iuc`/`ouc` selectors and a
Fortran callback for custom eigenvalue selection.

## Completed numerical tests

GNU Fortran, Fortran 2018, in **checked** (`-O0 -g -fcheck=all
-ffpe-trap=invalid,zero,overflow`) and **optimized** (`-O3`) builds:

| SciPy comparison seed | Cases per build | Sorted cases per build | Result |
| --- | ---: | ---: | --- |
| `20261108` | 480 | 217 | All passed |
| `37271` | 384 | 168 | All passed |
| **Total across two builds** | **1728** | **770** | **All passed** |

The maximum normalized errors across those runs were:

| Metric | Maximum |
| --- | ---: |
| Relative Schur reconstruction | `1.32e-14` |
| Relative spectral matching error vs SciPy | `4.70e-15` |
| Unitarity error, normalized | `1.35e-15` |
| Strict lower-triangle values | `0` |
| Selected invariant-subspace residual | `8.51e-15` |

The reference uses `scipy.linalg.schur(..., output='complex')` with the same
sort regions, plus eigenvalue multiset matching against `scipy.linalg.eigvals`.
The callback uses positive imaginary part. Tests cover matrices of order
0–40, real/complex data, triangular and Hermitian forms, scaled data,
repeated and zero spectra, and nonnormal matrices. Cases with eigenvalues
near a selection boundary are tested without sorting because one bit of
roundoff can legitimately change classification. The native tests explicitly
exercise zero-eigenvalue sorting and invalid-input failure statuses.

Existing regression tests also passed in checked and optimized builds:

- `tools/linalg_reference/run_general_eigen_checks.py`: 120 cases per build,
  plus native `eig`/`eigvals` tests.
- `tools/linalg_reference/run_funm_checks.py`: 45 general, 82 stress,
  30 reordered, 30 diagonalizable repeated-spectrum, and 94 adaptive-Taylor
  cases per build, plus native `funm` tests.

No full FPM suite was run (FPM unavailable); results do not establish
robustness on all nonnormal/defective inputs or SciPy's real Schur output.
