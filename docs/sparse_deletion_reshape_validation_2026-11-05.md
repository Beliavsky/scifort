# CSR deletion, reshape and scalable sparse assignment validation

## New capabilities

- `csr_delete_entries`, `csr_delete_rows`, `csr_delete_columns` with
  transactional validation, repeated selectors, and canonical CSR output.
- `csr_reshape`, matching SciPy's row-major (`order="C"`) sparse reshape.
- Sparse RHS `csr_assign` now has O(nnz(target)+nnz(source)+selection length)
  auxiliary storage and *does not* expand a sparse Cartesian selection.

### Reproduce the tests

```sh
python tools/sparse_reference/compare_sparse_deletion_scipy.py --mode both --cases 350 --seed 20261105
python tools/sparse_reference/compare_sparse_deletion_scipy.py --mode both --cases 350 --seed 91327
python tools/sparse_reference/run_native_sparse.py --mode both
```

The native regression additionally checks invalid requests leave input
unchanged, explicit zero preservation, repeated indices, a sparse
50,000 × 50,000 selection and one-billion-column indexing. The randomized
comparison suite validates numerical results and shape against NumPy and
SciPy for five operations, including zero-sized and rectangular shapes.

### Scope and limitations

- Only real-binary64 canonical one-based CSR is supported directly.
- `csr_delete_entries` removes stored entries; it does not shrink dimensions.
- `csr_assign` with scalar or dense source still expands Cartesian updates.
- Dense `reshape` in Fortran has column-major semantics; `csr_reshape` uses
  row-major semantics deliberately for SciPy compatibility.
- New operations do not depend on third-party numerical libraries.

### Verification performed on 2026-11-05

- Two independent seeds with 350 cases per build: **1,400** randomized
  reference cases passed across checked and optimized GNU Fortran builds.
- A third seed with 175 cases per build: **350** further randomized cases
  passed after switching the reshape reference to SciPy CSR `reshape`.
- **Five affected native tests** passed in both checked and optimized modes:
  `test_sparse_deletion_reshape`, `test_sparse_block_assignment`,
  `test_sparse_assembly`, `test_sparse_matrix` and
  `test_sparse_csc_diagonals`.
- Existing `test_sparse_formats_norm`, `test_sparse_linalg`,
  `test_sparse_direct`, and `test_sparse_spilu` passed in checked mode.
- Existing block-assignment and assembly/indexing comparison scripts passed
  **560** cases across the two build modes.
- A broad legacy sparse comparison and optimized sparse-linear-algebra
  integration command exceeded the local execution limit; those checks
  are **not reported as passing**. The full `fpm test` suite was not run.
- The existing `sparse_linalg` executable-stack linker warning persists.
