# `qr_multiply` validation — 2026-11-09

Source: `src/linalg_qr_multiply.f90`. Reference driver:
`reference/linalg_qr_multiply_driver.f90`. Independent Python oracle:
`reference/check_linalg_qr_multiply.py`.

Environment: GNU Fortran 14.2.0, Python 3.13, NumPy 2.3.5, SciPy 1.17.0.

| Test | Checked `-O0 -g -fcheck=all` | Optimized `-O3` |
|---|---:|---:|
| Independent randomized SciPy cases | 2,480 passed | 2,480 passed |
| Native QR multiplication cases | 1,296 passed | 1,296 passed |
| Existing `test_linalg_qr` | Passed | Passed |
| Existing `test_linalg_qr_updates` | 1,872 cases passed | 1,872 cases passed |

The maximum **relative** numerical deviations in either build were:

- Q multiplication after diagonal phase alignment: `2.62e-15`.
- R entries after diagonal phase alignment: `1.53e-15`.
- Reconstruction of permuted input using aligned Q and Fortran R: `1.53e-15`.

The native tests cover real and complex data, rectangular and square
matrices, rank-deficient input fixtures, left/right modes, conjugated complex
Q, vector and matrix multiplication, column pivoting, mixed real/complex
arguments, and reuse of compact QR factors. Explicit errors include invalid
output and factor shapes, invalid modes, NaN input/factors, and zero-sized
coefficient matrices. Zero multiplier rows are also tested.

For nondegenerate matrices the reference comparison aligns Q column phases
with the ratio of Fortran and SciPy R diagonals (valid QR factors need not
share signs or phases). The script checks the product, reconstruction, and
permutation. The `real A + complex C` branch is compared against SciPy QR
and NumPy matrix multiplication, rather than SciPy 1.17.0 `qr_multiply`,
which discarded the imaginary component of C during our validation.

**Build scope:** These tests compile the implementation and directly related
modules, not every source file in the repository. `fpm test` was not run:
FPM was not installed in the validation environment. No Windows Intel
Fortran tests were performed.
