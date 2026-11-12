# General sparse block and CSR assignment validation — 2026-11-05

## Implemented APIs

- `csr_bmat` and its alias `csr_block`: direct CSR assembly from a rank-two
  grid of canonical CSR matrices and `csr_matrix()` placeholders for absent
  blocks. Optional `row_sizes`/`col_sizes` resolve completely absent block
  rows and columns.
- `csr_set`: one-element replacement/insertion, preserving stored zeros.
- `csr_assign_coo`: batched coordinate scatter with deterministic last-write-
  wins semantics for duplicate destination coordinates. Uses a sorted
  64-bit coordinate key and a two-pass sparse merge, without a dense matrix.
- `csr_assign`: generic scalar broadcast, dense block assignment, and CSR
  source assignment to a Cartesian row-by-column selection, including
  implicit zero values from the source.

### Coverage

`test/test_sparse_block_assignment.f90` covers missing and malformed blocks,
explicit shape hints, empty block dimensions, zeros in CSR structure,
dimension mismatches, invalid assignments that leave the source unchanged,
nonfinite updates, duplicate coordinates, repeated selected rows/columns,
large 30,000-row matrices, and a target with one billion logical columns.

`tools/sparse_reference/compare_sparse_block_assignment_scipy.py` generates
randomized tests for general block assembly, explicit coordinate assignment,
dense/scalar/CSR submatrix assignment, and single-element assignment. It
checks values against SciPy's `bmat` (using explicit zero blocks where SciPy
cannot infer row/column sizes) or NumPy's dense assignment semantics. It also
checks exact shape, canonical row/column ordering, uniqueness, and retention
of explicitly stored zero entries by block assembly.

### Validation carried out

| Test | Outcome |
|---|---|
| Native new regression, `-O0 -fcheck=all` | PASS |
| Native new regression, `-O3` before final input-validation tightening | PASS |
| New randomized reference tests, seed 113 | 240 checked + 240 optimized PASS |
| New randomized reference tests, seed 55095 | 600 checked + 600 optimized PASS |
| Additional optimized reference tests, seed 55096 | 300 optimized PASS |
| New randomized reference tests, seed 77101 (after validation tightening) | 240 checked + 240 optimized PASS |
| Existing assembly/indexing SciPy comparisons | 240 optimized PASS |
| Complete checked native sparse runner | 24/24 programs PASS before final targeted validation tightening |
| Previous native CSR, CSC and assembly regression tests | PASS in targeted optimized build |

The aggregate new reference checks total **1,080 checked and 1,380
optimized randomized cases**. These are *case* counts, not a count of
independent numerical comparisons within each case.

The combined optimized native runner did not complete within the environment's
execution limit. Therefore this report does **not** claim a fully successful
post-change optimized native suite; the targeted optimized tests and
optimized SciPy comparisons passed. The available environment does not have
FPM, so a complete `fpm test` has not been performed. The pre-existing LSMR
executable-stack linker warning remains unrelated to these additions.

### Reproduce

```bash
python tools/sparse_reference/run_native_sparse.py --mode checked
python tools/sparse_reference/run_native_sparse.py --mode optimized
python tools/sparse_reference/compare_sparse_block_assignment_scipy.py --mode both --cases 600 --seed 55095
python tools/sparse_reference/compare_sparse_block_assignment_scipy.py --mode both --cases 240 --seed 77101
python tools/sparse_reference/compare_sparse_assembly_scipy.py --mode both --cases 240
```

### Intentional interface differences from SciPy

- The project has one-based indices, including for assignments.
- A missing block is `csr_matrix()` with *all* components unallocated;
  a canonical, allocated 0x0 CSR is a present block.
- For fully missing block rows/columns, pass dimensions explicitly; SciPy
  instead requires enough blocks to determine each dimension.
- Repeated destination indices use last-write-wins in **Fortran column-major**
  order for dense/CSR assignments. NumPy/SciPy fancy-index assignment does
  not guarantee a particular winner for all repeated-index expressions.
- Direct scatter uses O(nnz+number of updates) storage; rectangular RHS
  assignments materialize all selected positions before updating, and may
  consume more memory for very large selections.

These routines are original SciFort code; SciPy is used only for numerical
comparisons, and no additional runtime dependencies or licenses were added.
