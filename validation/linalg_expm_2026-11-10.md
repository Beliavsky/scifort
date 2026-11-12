# Validation: dense matrix exponential, 2026-11-10

## Environment

- GNU Fortran 14.2.0 (checked `-O0 -g -fcheck=all` and optimized
  `-O3 -march=native` builds).
- Python 3.13, SciPy 1.17.0, NumPy 2.3.5 (independent references).
- Standalone builds of existing `kinds`, `linalg`, `linalg_dense` and new
  `linalg_expm` modules. FPM was not installed; the full FPM test suite was **not** run.

## Results

| Check | Checked | Optimized |
|---|---:|---:|
| Native `test_linalg_expm` | 999 pass | 999 pass |
| Independent `scipy.linalg.expm` cases | 1,609 pass | 1,609 pass |
| Matrix entries compared to SciPy | 302,801 | 302,801 |
| Cases using one or more squarings | 571 | 571 |
| Worst normalized Frobenius discrepancy | 9.662e-12 | 9.688e-12 |
| Existing `test_linalg_dense` | pass | pass |
| Existing `test_linalg_expm_condition` | pass | pass |
| Existing `test_linalg_matrix_functions` | pass | pass |
| Executable example | pass | pass |

The worst discrepancy occurred in a 2×2 rotation generator with angle
1,000 radians. The reference suite includes real/complex cases, square
orders from zero through 64, different scales, triangular, skew-Hermitian,
Hermitian, normal and nonnormal matrices, rank-deficient inputs, nilpotent
Jordan blocks, large diagonal entries, and high-frequency rotations.

Native tests additionally check one leading batch dimension, malformed output
shapes, nonfinite inputs, initialized failure outputs, the empty case,
Jordan-block identity, and `exp(A)exp(-A) ≈ I` on randomly constructed inputs.

## How to repeat

From the root of the included source tree:

```shell
python tools/linalg_reference/run_expm_checks.py
```

This script compiles and tests in temporary directories; it does not require
FPM, mutate existing Fortran modules, or access the network. It requires
`gfortran`, NumPy and SciPy installed. In addition, the existing condition
number and matrix function regressions were compiled manually using their
source dependencies.

## Limits

This validates the new `expm` subset, not numerical identity with SciPy on
all inputs. The Fortran implementation uses fixed [13/13] Pade rather than
SciPy's variable degree selection; consequently rounding and runtime differ.
There is no `float32`/`complex64` overload, arbitrary batch rank, or BLAS
backend selection. Nonfinite inputs return a failure status; overflows in
`exp(A)` cannot be represented by the output type.
