# Validation: structure detection and SVD diagonal construction (2026-11-09)

Source archive: latest SciFort matrix-constructor snapshot plus two new
self-contained Fortran modules.

## Toolchain

- GNU Fortran 14.2.0, checked build (`-std=f2008 -fcheck=all -fbacktrace`)
- GNU Fortran 14.2.0, optimized build (`-std=f2008 -O3 -march=native`)
- Python 3.13, NumPy 2.3.5, SciPy 1.17.0

## Results

| Validation | Checked | Optimized |
|---|---:|---:|
| Native `test_linalg_structure_diagsvd` | Pass | Pass |
| Independent SciPy comparison cases | 3,624 / 3,624 | 3,624 / 3,624 |
| Individual output values compared | 57,609 | 57,609 |
| Expected invalid inputs detected | 709 | 709 |
| Block-diagonal and Khatri–Rao regression | Pass | Pass |
| Hadamard, Leslie and DFT regression | Pass | Pass |
| Pascal regression | Pass | Pass |
| Companion/Fiedler regression | Pass | Pass |
| Helmert/Hilbert regression | Pass | Pass |

## Test scope

- Integer32, exact integer64 (including values above `2**53`), real64,
  complex64; rank-2 inputs and single-batch rank-3 inputs.
- Sparse and dense matrices, symmetry and Hermitian cases, zero-sized
  cores, rectangular shapes, large and small scales, NaNs and infinities.
- Exact and tolerant comparisons and directed complex conjugation.
- Rectangular SVD diagonal construction with exact integer outputs,
  zero-sized dimensions, multiple batches, and invalid singular-value
  lengths.
- Intentionally invalid square shapes and zero batch sizes.

Native and independent checks passed under both compiler configurations.
No new external package was linked. The **entire** existing FPM test suite
was not run because FPM is not installed in this execution environment.
