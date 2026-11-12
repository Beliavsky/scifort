# Coordinate sparse matrices (`scipy.sparse.coo_array` subset)

`scifort_sparse_coo` adds two-dimensional real binary64 **COO** storage to
SciFort. Indices are one-based, consistent with the existing `scifort_sparse`
CSR/CSC modules; unlike SciFort's canonical CSR and CSC, COO **allows unsorted
coordinates, repeated (row,column) pairs, and explicitly stored zero values**.
This is useful for assembling finite-element stiffness/mass matrices from
local contributions without making a dense global matrix.

## API

| API | SciPy counterpart | Description |
|---|---|---|
| `coo_matrix` | `scipy.sparse.coo_array` | `n_rows`, `n_cols`, `row(:)`, `col(:)`, `data(:)` |
| `valid_coo(a)` | Structural validation | Require finite values, one-based in-bounds coordinates, and matching buffers |
| `coo_from_triplets(m,n,rows,cols,values,info)` | `coo_array((data,(row,col)),shape=...)` | Copy triplets **without sorting or merging**, including explicit zeros |
| `coo_from_dense(dense,info)` | `coo_array(dense)` | Store only numerically nonzero values |
| `coo_to_dense(a,info)` | `.toarray()` | Sum duplicate coordinates in dense output |
| `coo_to_csr(a,info)` | `.tocsr()` | Canonical CSR; duplicates summed, even when result is zero |
| `coo_to_csc(a,info)` | `.tocsc()` | Canonical CSC; duplicates summed, even when result is zero |
| `coo_from_csr(a,info)` and `coo_from_csc(a,info)` | `.tocoo()` | Retain explicitly stored zeros and shape |
| `coo_transpose(a,info)` | `.transpose()` | Swap coordinates and shape; preserve original triplet order |
| `coo_matvec(a,x,info)` | `a @ x` | Matrix-vector multiplication with duplicate accumulation |
| `coo_get(a,i,j,info)` | Scalar indexing | Return sum over entries at one coordinate, or zero if absent |
| `coo_sum_duplicates(a,info)` | `.sum_duplicates()` | **In-place** canonicalization: sorted row-major coordinates, duplicate sums, zeros retained |
| `coo_eliminate_zeros(a,info)` | `.eliminate_zeros()` | **In-place** removal of explicit zeros without summing duplicates |
| `sparse_find(a,rows,cols,values,info)` | `scipy.sparse.find(a)` | COO, CSR, or CSC overload; canonicalize duplicate entries, then exclude zeros |

Every nontrivial API optionally reports `info`, where `0` is success and
`1` (`sparse_invalid_input`) reports invalid shapes, coordinates, or
nonfinite arithmetic. Failed `coo_sum_duplicates` leaves the input intact.
The data layout has the same explicit-zero and duplicate semantics as SciPy,
except for **one-based** indices.

## Example: sparse finite-element assembly

```fortran
program assemble_coo_example
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_to_dense
    use scifort_sparse_coo, only : coo_matrix, coo_from_triplets, coo_to_csr, sparse_find
    implicit none
    type(coo_matrix) :: coo
    type(csr_matrix) :: csr
    integer, allocatable :: row(:), col(:)
    real(dp), allocatable :: value(:)
    integer :: info

    ! Two overlapping 2-node element matrices [[1,-1],[-1,1]].
    coo = coo_from_triplets(3,3, [1,1,2,2,2,2,3,3], &
        [1,2,1,2,2,3,2,3], &
        [1.0_dp,-1.0_dp,-1.0_dp,1.0_dp, &
         1.0_dp,-1.0_dp,-1.0_dp,1.0_dp], info)
    if (info /= 0) error stop 'invalid element triplets'
    csr = coo_to_csr(coo, info)
    if (info /= 0) error stop 'assembly failed'
    print *, csr_to_dense(csr) ! tridiagonal 3x3 stiffness matrix
    call sparse_find(coo, row, col, value, info)
    if (info /= 0) error stop 'find failed'
    print *, size(value) ! Seven unique nonzero entries
end program assemble_coo_example
```

`coo_to_csr` uses the existing SciFort `csr_from_coo` sorting and summing
implementation, so it does not materialize a dense array. The same applies
to COO-to-CSC, canonicalization, and `sparse_find`. Explicit dense conversion
is the only routine allocating an `n_rows`-by-`n_cols` array.

**Complex matrices, COO-specific slicing/elementwise arithmetic, and
nonfinite stored values** remain outside the two-dimensional COO subset.
For N-dimensional storage, use the separate
[`scifort_sparse_coo_nd` module](sparse_coo_nd.md); mutable DOK and LIL
formats are documented [here](sparse_dok_lil.md). COO matrix-vector multiplication sums in
stored order; floating-point last-bit differences from CSR or SciPy may occur
for unusually ill-conditioned duplicate sums.

## Validation

`python tools/sparse_reference/compare_sparse_coo_scipy.py --mode both --trials 200`

The independent SciPy 1.17 comparison covers randomized rectangular matrices with 0–10 rows and columns, duplicate values/cancellations, explicit zeros, COO to
CSR/CSC and back, dense representation, matrix-vector products, `find`, and
canonicalization. Use `--seed 17` for a second randomized batch.

`python tools/sparse_reference/run_native_sparse.py --mode both` runs the
standalone native regressions including `test_sparse_coo`.

No external numerical dependencies or third-party source code have been added.
