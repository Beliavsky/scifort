# DIA, BSR and sparse-norm validation (2026-11-05)

## Implementation

- `src/sparse.f90`: canonical column-aligned DIA and fixed-block BSR formats,
  validation, CSR conversions, dense conversions and direct matrix-vector actions.
  Row pointers, column indices and all matrix coordinates are **one-based**.
- `src/sparse_linalg.f90`: CSR/CSC `norm` for Frobenius, ±1, ±infinity, and
  approximate 2-norm using existing `svds`; one-based `norm_axis` for column
  and row vector norms.
- No imported code, new libraries or changes to the FPM dependency list.

## Reproducible checks

```sh
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_reference/compare_sparse_formats_norm_scipy.py --mode both --cases 120
```

Use `--seed 74213` for an additional independent randomized sample. A
prebuilt comparison executable can be supplied with `--executable PATH` to
avoid repeatedly compiling the full 5,000+ line sparse solver module.

Environment: GNU Fortran with `-std=f2018`, checked flags `-O0 -g
-fcheck=all -fbacktrace`, optimized flags `-O3`, SciPy 1.17.0 and NumPy.

### Verified results

| Check | Checked | Optimized |
|---|---:|---:|
| Original 22 native sparse regression programs, including new format/norm test | 22 passed | 22 passed |
| Randomized comparisons, seed 20261105 | 220 passed | 240 passed |
| Randomized comparisons, seed 74213 | 180 passed | 180 passed |
| Randomized comparisons, seed 2027 (after final input-guard edits) | 120 passed | 120 passed |
| Independently verified spectral 2-norms | 104 | 108 |
| Largest direct DIA/BSR regression | 18,000 rows | 18,000 rows |

The randomized suite compares five matrix norms, column 2-norms and row
1-norms, CSR-to-DIA/BSR-to-dense values, and direct DIA/BSR vector products,
including stored zeros and repeated COO coordinates. Spectral norms are
checked on one case in five to control the runtime of Golub--Kahan SVD.
DIA numerical differences from SciPy/NumPy were no larger than ~1.8e-15
in these comparisons. Both CSC and CSR norm overloads are included in the
native tests, along with arbitrary DIA offset order and invalid arguments.

## Behavior and limitations

- DIA stores all columns of every occupied diagonal, including padding.
  CSR-to-DIA may require much more memory for widely spaced diagonals.
  Back-conversion expands in-bound implicit zeros and may change `nnz`.
- BSR blocks must divide both matrix dimensions. BSR-to-CSR expands full
  stored blocks, including zero padding; partial edge blocks are rejected.
- Matrix norm order 2 reuses the iterative `svds`, and may report
  nonconvergence for difficult matrices. Other unsupported matrix orders
  return an invalid-input status, rather than an invented norm.
- The full FPM suite was not run, because FPM was unavailable in this
  environment. An existing executable-stack warning associated with the
  LSMR callback code remains; no change in that component was attempted.
