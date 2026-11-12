# Incremental sparse assembly with DOK and LIL

The `scifort_sparse_dok_lil` module provides one-based, real binary64,
2-dimensional sparse matrix formats corresponding to the core construction
and indexing behavior of `scipy.sparse.dok_array` and `lil_array`.

**DOK** uses an open-addressed hash table, offering expected O(1) lookup,
insertion, and removal of individual entries. Deletions use tombstones and
insertion triggers geometric growth or rehashing. **LIL** stores ascending
column coordinates and data in dynamically grown row lists. Lookup within
a row uses binary search; insertion and deletion can take O(nnz in that row)
because entries must be shifted. Both formats are useful for incremental
assembly before converting to CSR or CSC for matrix-vector arithmetic.

The public `dok_matrix` and `lil_matrix` types expose `n_rows` and `n_cols`.
Storage fields are private to preserve structural invariants; do **not**
modify the shape fields after construction. Use `valid_dok` or `valid_lil`
to test the full internal structure and the finite-value constraints.

## API

| Procedure | SciPy behavior / purpose |
|---|---|
| `dok_empty(m,n,info)`, `lil_empty(m,n,info)` | Construct empty, shaped sparse arrays |
| `dok_set(a,i,j,v,info)`, `lil_set(...)` | Replace one coordinate; assigning zero deletes it |
| `dok_get(a,i,j,info)`, `lil_get(...)` | Return one value or implicit zero |
| `dok_add(a,i,j,dv,info)`, `lil_add(...)` | Add increment to one coordinate; zero sum deletes it |
| `dok_delete(a,i,j,info)`, `lil_delete(...)` | Remove a coordinate; absent coordinate is a no-op |
| `dok_clear(a,info)`, `lil_clear(...)` | Remove all stored entries, preserve shape |
| `dok_nnz(a,info)`, `lil_nnz(...)` | Number of stored entries **including explicit zeros** |
| `dok_from_coo`, `lil_from_coo` | Canonicalize duplicates, preserve summed-to-zero coordinates |
| `dok_from_csr`, `lil_from_csr` | Import sorted CSR, preserve stored zeros |
| `dok_from_csc`, `lil_from_csc` | Import CSC, preserve stored zeros |
| `dok_from_dense`, `lil_from_dense` | Import numerical nonzeros only |
| `dok_to_coo`, `lil_to_coo` | Export unique stored coordinates (DOK order is unspecified) |
| `dok_to_csr`, `lil_to_csr` | Export canonical, row-sorted CSR |
| `dok_to_csc`, `lil_to_csc` | Export canonical CSC |
| `dok_to_dense`, `lil_to_dense` | Construct a dense array (only if the dimensions fit memory) |
| `dok_matvec`, `lil_matvec` | Evaluate `a*x` with a vector of length `n_cols` |
| `dok_transpose`, `lil_transpose` | Transpose, retaining explicit zeros |
| `dok_from_lil`, `lil_from_dok` | Convert between assembly formats |
| `dok_eliminate_zeros`, `lil_eliminate_zeros` | Delete stored zero entries |
| `valid_dok(a)`, `valid_lil(a)` | Check one-based coordinates, uniqueness, ordering, and finiteness |

All optional `info` arguments use `sparse_success=0` and
`sparse_invalid_input=1`, provided by `scifort_sparse`. On invalid
construction, the returned matrix is uninitialized and fails `valid_*`.
On invalid element operations, no change is made to the matrix. All stored
values and all valid results are finite; overflowing additions and
matrix-vector results are rejected.

Indices are **one-based**, unlike SciPy's zero-based Python indices.
DOK operations use hash-table order internally; convert to canonical CSR
when deterministic row-major order is needed.

### Explicit zeros

The distinction between *assignment* and *conversion* is important.
As with SciPy's DOK and LIL sparse arrays:

- Setting an element to exactly `0.0_dp` removes it (and does not grow the
  stored-entry count); adding a value that cancels to zero also removes it.
- Converting an existing CSR, CSC or COO matrix **preserves** its explicitly
  stored zero coordinates. COO duplicates are summed first, and their zero
  sums remain stored until the user calls `*_eliminate_zeros` or assigns
  zero to that coordinate.
- `nnz` counts stored entries, including zeros. For example, the numeric
  matrix can be all zero even when `nnz > 0` after COO conversion.

## Example

```fortran
program dok_lil_example
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_to_dense, sparse_success
    use scifort_sparse_dok_lil, only : dok_matrix, lil_matrix, &
        dok_empty, dok_set, dok_add, lil_from_dok, lil_to_csr
    implicit none
    type(dok_matrix) :: d
    type(lil_matrix) :: l
    type(csr_matrix) :: a
    integer :: info

    d = dok_empty(3,3,info)
    if (info /= sparse_success) error stop 'bad shape'
    call dok_set(d,1,1,2.0_dp,info)
    call dok_add(d,1,2,-1.0_dp,info)
    call dok_add(d,2,1,-1.0_dp,info)
    call dok_set(d,2,2,3.0_dp,info)
    if (info /= sparse_success) error stop 'bad entry'
    l = lil_from_dok(d,info)
    if (info /= sparse_success) error stop 'conversion failed'
    a = lil_to_csr(l,info)
    if (info /= sparse_success) error stop 'conversion failed'
    print *, csr_to_dense(a)
end program dok_lil_example
```

## Tests and limitations

`test/test_sparse_dok_lil.f90` exercises empty and rectangular matrices,
structural zeros and duplicate cancellation, sorted CSR/CSC conversions,
transposition, copying, 2,600 deterministic edits, large sparse shapes,
nonfinite inputs, invalid coordinates and arithmetic overflow.

See the [validation log](sparse_dok_lil_validation_2026-11-05.md) for
exact commands and results.

`tools/sparse_reference/compare_sparse_dok_lil_scipy.py` runs randomized
SciPy 1.17 comparisons for repeated assignments, additive updates,
deletions and lookups, checking both stored counts and numerical values at
every operation. CSR output is canonicalized on both sides before comparing
its index/data arrays. Runs can be reproduced with:

```console
python tools/sparse_reference/compare_sparse_dok_lil_scipy.py --mode both --trials 300 --seed 20261105
python tools/sparse_reference/compare_sparse_dok_lil_scipy.py --mode both --trials 200 --seed 84723
python tools/sparse_reference/run_native_sparse.py --mode both
```

Only binary64 real arithmetic and two-dimensional matrices are supported.
There is no Python-style slicing, broadcasting or matrix arithmetic API for
DOK/LIL. Convert to CSR/CSC before using SciFort's sparse solvers.
For broad bulk COO assembly, COO triplets may be more efficient than many
random incremental updates. This module introduces no new dependencies.
