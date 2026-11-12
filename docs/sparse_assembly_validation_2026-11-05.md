# Sparse assembly and indexing validation — 2026-11-05

## New interfaces

The new `src/sparse.f90` routines are `csr_hstack`, `csr_vstack`,
`csr_block_diag`, `csr_kron`, `csr_get`, `csr_select` and `csr_slice`.
All use canonical real binary64 one-based CSR storage and reject invalid
arguments via `sparse_invalid_input`.

## Tests completed

GNU Fortran 14.2.0, Python 3.13, NumPy 2.3.5, SciPy 1.17.0:

- `test/test_sparse_assembly.f90` — passed in checked (`-O0 -fcheck=all`)
  and optimized (`-O3`) builds, including 30,000-row concatenation,
  stored-zero preservation, very wide logical matrices (10^9 columns),
  repeated selections, reversed slices, overflow rejection and bad indices.
- Randomized SciPy assembly/indexing comparison driver:
  420 checked cases (seed 13579), 450 checked (seed 93671),
  600 optimized (seed 71881), and 600 checked plus 600 optimized (seed 55091).
  **2,670 new randomized cases passed** in total.
- Existing `test_sparse_matrix` and `test_sparse_csc_diagonals` standalone
  native tests passed in the checked build.
- Earlier `compare_sparse_scipy.py` and `compare_sparse_csc_scipy.py`
  suites each passed 120 cases in **both** checked and optimized builds
  (4,080 numerical matrix comparisons combined).

SciPy 1.17.0 raises `ValueError` for `scipy.sparse.block_diag([])`;
the cross-language driver instead uses a zero-by-zero CSR reference in
that one case. SciFort deliberately returns a valid empty direct sum.
The Fortran slicing API uses inclusive one-based bounds, and reference
indices are translated to Python indexing explicitly.

The full `tools/sparse_reference/run_native_sparse.py --mode checked`
runner exceeded the execution time limit before yielding results. Therefore
**its full 23-program regression suite is not claimed to have passed** for
this update. The full `fpm test` was also not executed (FPM is unavailable).

## Reproduce

```bash
python tools/sparse_reference/compare_sparse_assembly_scipy.py --mode both --cases 600 --seed 55091
python tools/sparse_reference/compare_sparse_scipy.py --mode both --cases 120 --seed 31415
python tools/sparse_reference/compare_sparse_csc_scipy.py --mode both --cases 120 --seed 27182
python tools/sparse_reference/run_native_sparse.py --mode both
```

`run_native_sparse.py` now registers `test_sparse_assembly` as its first
program. The implementation and test code were independently authored for
SciFort; no third-party implementation code was imported.
