# N-dimensional COO indexing, slicing and pointwise gathering

`scifort_sparse_coo_nd_index` extends the existing real binary64,
positive-rank, one-based `scifort_sparse_coo_nd%coo_nd_array` without altering
its representation. Selection reads only the stored COO entries and never
materializes the full dense input tensor.

The interface is designed to be idiomatic in Fortran, and is **not** a direct
translation of Python's indexing syntax. In particular, all subscripts are
**one-based**, `coo_nd_slice` endpoints are **inclusive**, and when two or
more axes use index arrays, `coo_nd_select` uses the **Cartesian product**
of those selections. NumPy/SciPy simultaneously advanced-indexed arrays can
instead select paired or broadcasted positions. Pointwise coordinate
selection is available through `coo_nd_gather`.

## Public API

| Routine | Purpose |
|---|---|
| `coo_nd_all()` | Select an entire axis, including an empty axis |
| `coo_nd_at(index)` | Select one one-based index, optionally dropping that dimension |
| `coo_nd_slice(first,last,stride)` | Inclusive Fortran-style section; omitted endpoints are based on the sign of the stride |
| `coo_nd_indices(indices)` | Select a one-dimensional array of integer positions, including repeated and unsorted values |
| `coo_nd_select(a,selectors,status,keepdims,max_nnz)` | Apply one selector per axis and return a COO array |
| `coo_nd_gather(a,points,status,max_nnz)` | Select corresponding pointwise indices from a `(ndim,npoints)` coordinate array, returning rank-one COO |

Each selector is of derived type `coo_nd_selector`; combine constructors in a
Fortran array constructor, or declare and populate an array of selectors.
By default, an integer selector **drops** its axis. Use `keepdims=.true.` to
retain it with extent one. When **all** axes are indexed by integers, the
result is a rank-one singleton `[1]`: the underlying COO type does not
represent rank-zero scalars. For single-point scalar reads,
`coo_nd_get(a,indices)` is also available in the base module.

The slice stride defaults to `1`. For a negative stride, omitted `first`
and `last` select the final and first source element, respectively, e.g.
`coo_nd_slice(stride=-1)` reverses an axis. A reversed or empty section
produces an output extent of zero. A stride of zero, an out-of-range
explicit endpoint, or a negative/zero scalar index is rejected rather than
silently clipped to the bounds.

Arrays of indices may be empty, unsorted and repeated. Repetitions replicate
stored entries in the output, including duplicate input entries and stored
zeros. Neither `coo_nd_select` nor `coo_nd_gather` automatically coalesces
duplicates or drops explicitly stored zeros.

## Example

```fortran
program nd_index_example
    use scifort_kinds, only : dp
    use scifort_sparse_coo_nd, only : coo_nd_array, &
        coo_nd_from_triplets, coo_nd_get
    use scifort_sparse_coo_nd_index, only : coo_nd_at, &
        coo_nd_slice, coo_nd_indices, coo_nd_select, coo_nd_gather
    implicit none
    type(coo_nd_array) :: a, b, c
    integer :: status

    a = coo_nd_from_triplets([3,4,5], &
        reshape([3,2,4, 3,2,4, 1,4,3], [3,3]), &
        [7.0_dp,-2.0_dp,9.0_dp],status)
    if (status /= 0) error stop 'construction'

    ! Select A(3,2:4:2,5:2:-1), leaving a (2,4) sparse array.
    b = coo_nd_select(a,[coo_nd_at(3),coo_nd_slice(2,4,2), &
                         coo_nd_slice(5,2,-1)],status)
    if (status /= 0) error stop 'selection'
    print *, b%shape                   ! 2, 4
    print *, coo_nd_get(b,[1,2])      ! 5.0, after summing duplicate entries

    ! Orthogonal indexing, including a repeated first-axis index.
    b = coo_nd_select(a,[coo_nd_indices([3,1,3]), &
                         coo_nd_at(2),coo_nd_at(4)],status)
    if (status /= 0) error stop 'index array'
    print *, b%shape                   ! 3
    print *, coo_nd_get(b,[3])        ! 5.0

    ! Pairwise advanced indexing with two coordinate tuples.
    c = coo_nd_gather(a,reshape([3,2,4, 1,4,3],[3,2]),status)
    if (status /= 0) error stop 'gather'
    print *, c%shape                   ! 2
    print *, coo_nd_get(c,[2])        ! 9.0
end program nd_index_example
```

## Error handling and limits

On success the optional integer `status` is `sparse_success` (0). On invalid
input, incompatible selector rank, out-of-bounds coordinate, or exceeded
output cap, it is `sparse_invalid_input` (1), and the result remains
unallocated. Inputs are never modified. The default `max_nnz` is **one
million** generated COO entries and can be raised with an `integer(int64)`
argument. The limit includes replicated stored zeros and duplicated entries,
so it prevents accidental blow-up caused by repeated index-array selections.

The selection implementation computes the generated entry count before
allocating the result and does not allocate storage proportional to the
product of the source shape. It uses direct coordinate mapping for slices
and linear scans of index arrays for repeated selections. `coo_nd_gather`
currently performs a direct search of stored entries for every requested
point (O(nnz*npoints)); it is intended for moderate point lists, not
billions of random lookups. Scalar reads are already supported by
`coo_nd_get`; future work could accelerate large gathers with a COO index.

## Testing

```console
python tools/sparse_reference/compare_sparse_coo_nd_index_scipy.py --mode both --trials 650 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_index_scipy.py --mode both --trials 450 --seed 274829 --floating
python tools/sparse_reference/run_native_sparse.py --mode both
```

The randomized driver creates SciPy COO inputs with duplicates and stored
zeros and uses NumPy to compute the independent dense reference. It also
cross-checks SciPy's own COO basic indexing for compatible selectors.
See [validation details](sparse_coo_nd_index_validation_2026-11-05.md).

Boolean mask selection and in-place indexed assignment are available in [the N-D COO edit module](sparse_coo_nd_edit.md).
