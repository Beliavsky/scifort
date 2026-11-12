# Boolean selection and indexed updates on N-dimensional COO arrays

`scifort_sparse_coo_nd_edit` extends the one-based, real binary64
`scifort_sparse_coo_nd%coo_nd_array` type with boolean selection and mutable
coordinate updates. It depends only on the existing SciFort sparse modules and
the Fortran standard library. The base N-D COO representation is unchanged.

## Public API

| Routine | Behavior |
|---|---|
| `coo_nd_mask_select(a,mask,status,max_nnz)` | Gather entries under a full-rank boolean mask into a **rank-one** sparse array. Keep duplicates and stored zeros. |
| `coo_nd_mask_axis(a,axis,mask,status,max_nnz)` | Select a boolean subset of one axis, preserving the original rank and all other axes. |
| `coo_nd_scatter(a,points,values,status,accumulate,keep_zeros,max_nnz)` | Perform last-write-wins assignment or additive scatter updates at one-based coordinates. |
| `coo_nd_mask_assign(a,mask,value,status,accumulate,keep_zeros,max_nnz)` | Overwrite or add a scalar at every true position in a full boolean mask. |

The full-mask API uses `logical :: mask(:)`, a **C-order flattened** tensor:
the last tensor axis varies fastest. Its length must equal the complete
product of the source dimensions. A mask of shape `(2,3)` corresponds to
positions `(1,1),(1,2),(1,3),(2,1),(2,2),(2,3)` in that order. This is
explicitly different from Fortran's column-major `RESHAPE` order. Convert a
Fortran array to the required order explicitly when constructing the mask.
For large sparse shapes, use pointwise `coo_nd_scatter` instead of a dense mask.

An axis mask has exactly `a%shape(axis)` values, and the output extent for
that axis becomes `count(mask)`. All axis labels and coordinates are **one-based**.

`points` for scatter has shape `(ndim,nupdates)`; `values` has length
`nupdates`. Repeated coordinates are processed **in their supplied order**.
With `accumulate=.false.` (default), the final supplied value wins. With
`accumulate=.true.`, all updates add to the pre-existing value. Negative and
out-of-range coordinates are rejected. No array indexing or broadcasting
syntax is inferred from `points`.

When a successful scatter has at least one update, source duplicates are
coalesced and the resulting COO storage is lexicographically sorted by
coordinates. **Untouched** explicitly stored zeros are retained; touched
positions with a zero final value are removed unless `keep_zeros=.true.`.
An update list of length zero is a no-op that preserves the input's existing
COO storage and duplicates.

All four routines support a configurable `integer(int64)` `max_nnz` (default
1,000,000). For mask assignment, it limits both the number of true mask
positions and the resulting number of stored entries. For other operations,
it limits generated output stored entries. Failed edits are transactional:
invalid shapes, indices, nonfinite values, arithmetic overflow and exceeded
limits return `sparse_invalid_input` and leave the source array unchanged.
Mask selection functions return an unallocated result on error. All successful
operations report `sparse_success` (0).

## Example

```fortran
program sparse_mask_example
    use scifort_kinds, only : dp
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets, coo_nd_get
    use scifort_sparse_coo_nd_edit, only : coo_nd_mask_select, coo_nd_scatter, coo_nd_mask_assign
    implicit none
    type(coo_nd_array) :: a, b
    integer :: info

    a = coo_nd_from_triplets([2,3],reshape([1,2,2,3],[2,2]), &
        [4.0_dp,5.0_dp],info)
    if (info /= 0) error stop 'construct'

    ! C-order mask selects (1,2) and (2,3).
    b = coo_nd_mask_select(a,[.false.,.true.,.false., &
                               .false.,.false.,.true.],info)
    if (info /= 0) error stop 'mask select'
    print *, b%shape              ! 2
    print *, coo_nd_get(b,[1])   ! 4.0

    ! Two updates to the same position; the last assignment wins.
    call coo_nd_scatter(a,reshape([1,2,1,2],[2,2]), &
        [3.0_dp,9.0_dp],info)
    if (info /= 0) error stop 'scatter'
    print *, coo_nd_get(a,[1,2]) ! 9.0

    ! Add 2 at all positions selected by the full mask.
    call coo_nd_mask_assign(a,[.false.,.true.,.false., &
                                .false.,.false.,.true.], &
        2.0_dp,info,accumulate=.true.)
    if (info /= 0) error stop 'mask assign'
    print *, coo_nd_get(a,[1,2]) ! 11.0
end program sparse_mask_example
```

## Complexity and differences from SciPy

`coo_nd_scatter` uses COO coordinate aggregation and a **stable bottom-up
mergesort** to handle `nnz + nupdates` events in approximately
`O((nnz + nupdates) log (nnz + nupdates))` time and proportional temporary
storage. It does not allocate an array proportional to the input shape.
Boolean full-mask routines necessarily require a dense logical mask;
`coo_nd_mask_select` also allocates a prefix map proportional to its length.

SciPy 1.17's N-D COO arrays provide sparse construction and some indexing,
but this mutable scatter interface is a **Fortran-specific extension**, not
a claim that SciPy's COO format supports direct item assignment. Tests use
SciPy's N-D COO construction/dense conversion plus independent NumPy indexing,
assignment and accumulation, including exact stored-entry-count checks.

## Validation

```console
python tools/sparse_reference/compare_sparse_coo_nd_edit_scipy.py --mode both --trials 650 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_edit_scipy.py --mode both --trials 450 --seed 274829 --floating
python tools/sparse_reference/run_native_sparse.py --mode both
```

See [validation results](sparse_coo_nd_edit_validation_2026-11-05.md).

For selection by a **sparse COO mask or coordinate tuples** instead of a
full dense Boolean array, see [sparse predicates and filtering](sparse_coo_nd_filter.md).
