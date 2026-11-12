# Sparse N-dimensional Boolean broadcasting and reductions

`scifort_sparse_coo_nd_bool_reduce` operates on complement-aware
`coo_nd_bool_mask` objects. The implementation is a Fortran extension:
these are not functions with identical names in `scipy.sparse`. NumPy
broadcasting and reduction results are the independent reference.

A mask has a Boolean `default_value` for its unstored coordinates and a
lexicographically sorted, unique exception list. A mask with implicit `true`
is therefore as sparse as one with implicit `false`. All coordinates and
axis numbers are **one-based**. Matrix/tensor storage is never expanded to
its full shape.

## Public procedures

| Procedure | Result | Semantics |
|---|---|---|
| `coo_nd_bool_broadcast_to(a,shape,status,max_exceptions)` | Boolean mask | Broadcast to `shape`, right-aligning dimensions |
| `coo_nd_bool_and_broadcast(a,b,status,max_exceptions)` | Boolean mask | Broadcast-compatible logical AND |
| `coo_nd_bool_or_broadcast(a,b,status,max_exceptions)` | Boolean mask | Broadcast-compatible logical OR |
| `coo_nd_bool_xor_broadcast(a,b,status,max_exceptions)` | Boolean mask | Broadcast-compatible logical XOR |
| `coo_nd_bool_difference_broadcast(a,b,status,max_exceptions)` | Boolean mask | `a AND NOT b` after broadcasting |
| `coo_nd_bool_count_axis(a,axis,status,keepdims,max_exceptions)` | Integer sparse count array | Exact number of `true` values along one axis |
| `coo_nd_bool_any_axis(a,axis,status,keepdims,max_exceptions)` | Boolean mask | Axis-wise `ANY` |
| `coo_nd_bool_all_axis(a,axis,status,keepdims,max_exceptions)` | Boolean mask | Axis-wise `ALL` |
| `coo_nd_bool_count_get(counts,point,status)` | `integer(int64)` | Retrieve exact count from the count array |
| `coo_nd_bool_any(a,status)` | `logical` | `ANY` over all elements, including implicit ones |
| `coo_nd_bool_all(a,status)` | `logical` | `ALL` over all elements, including implicit ones |
| `valid_coo_nd_bool_count_array(a)` | `logical` | Verify a sparse count array's invariants |

The existing `coo_nd_bool_and`, `_or`, `_xor`, and `_difference` operations
in `scifort_sparse_coo_nd_bool` still require equal shapes. The new
`*_broadcast` variants make broadcasting explicit instead of changing the
behavior of existing code.

### Boolean broadcasting

Dimensions align from the **right**, as in NumPy: source extents must match
the output or equal one. Leading dimensions are implicitly of extent one.
An extent-one source may broadcast to an empty dimension (extent zero). The
background Boolean value is preserved, and only the exception coordinates
are replicated. Binary operations first broadcast to a common compatible
shape, then merge sorted sparse exception sets.

The default maximum number of exceptions is 1,000,000. This cap applies to
each expanded mask and the binary result, **not just** to the final result:
for instance, an expansion larger than the limit can fail even when the
subsequent AND would be uniformly false. Set `max_exceptions` explicitly to
raise the limit when memory permits. No dense Boolean array is allocated.

### Reductions and exact counts

A sum of boolean values along an axis is stored as:

```fortran
type :: coo_nd_bool_count_array
    integer, allocatable :: shape(:)
    integer, allocatable :: coords(:,:)
    integer(int64), allocatable :: values(:)
    integer(int64) :: default_count = 0_int64
end type
```

At every unstored output coordinate, the count is `default_count`; at an
exception coordinate it is the corresponding `values(k)`. For a false
background, this baseline is zero. For a true background, it is the length
of the reduced axis. The implementation projects and aggregates only the
stored exception coordinates, then subtracts from that baseline if needed.
Counts are exact `integer(int64)` integers, not floating-point estimates.

`any_axis` tests whether each count is greater than zero; `all_axis` tests
whether it equals the reduced axis length. For zero-length reduced axes,
`ANY` is false, `ALL` is true, and every count is zero.

The default `keepdims` is false. For rank one, a reduced scalar is represented
by a **length-one array** (shape `[1]`) because SciFort's COO types require
positive rank. Pass `keepdims=.true.` to retain the reduced axis as a
singleton. `coo_nd_bool_any` and `coo_nd_bool_all` return scalar logicals
for full reductions. Both full reductions work even when the *total tensor
size* exceeds `huge(0_int64)`.

The reduction cap applies to the distinct projected exception coordinates
used in the intermediate integer count array. Failure returns an unallocated
result and `sparse_invalid_input`. Inputs are not modified.

## Example

```fortran
program example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_sparse, only : sparse_success
    use scifort_sparse_coo_nd_bool, only : coo_nd_bool_mask, coo_nd_bool_from_points
    use scifort_sparse_coo_nd_bool_reduce
    implicit none
    type(coo_nd_bool_mask) :: a, wide, any_row
    type(coo_nd_bool_count_array) :: counts
    integer :: status

    a = coo_nd_bool_from_points([2,1],reshape([1,1],[2,1]),status)
    if (status /= sparse_success) error stop 'make mask'
    wide = coo_nd_bool_broadcast_to(a,[2,4],status)
    if (status /= sparse_success) error stop 'broadcast'
    counts = coo_nd_bool_count_axis(wide,2,status)
    if (status /= sparse_success) error stop 'count along columns'
    print *, coo_nd_bool_count_get(counts,[1],status) ! 4
    print *, coo_nd_bool_count_get(counts,[2],status) ! 0
    any_row = coo_nd_bool_any_axis(wide,2,status)
    if (status /= sparse_success) error stop 'ANY'
    print *, coo_nd_bool_any(wide,status) ! T
end program
```

## Test and implementation notes

- No extra numerical libraries are required. Broadcasting expands only the
  exception lists; `coo_nd_sum_duplicates` provides stable sorting and sparse
  projection aggregation. The algorithms are original MIT-licensed work.
- The count-array `default_count` is nonnegative and stored coordinates have
  different values from the baseline. All output exception sets are sorted.
- When a rank-one axis is reduced, use `[1]` as the result coordinate.
- Invalid axes, incompatible shapes and exceeded caps return
  `sparse_invalid_input`; failure values are unallocated, not partial.

Run:

```console
python tools/sparse_reference/compare_sparse_coo_nd_bool_reduce_numpy.py --mode both --trials 500 --seed 716
python tools/sparse_reference/run_native_sparse.py --mode both
```

The comparison driver checks every resulting element of small random masks
against NumPy for broadcasting, AND/OR/XOR/difference, axis-wise `any`, `all`,
integer counts, keepdims, rank-one scalar reduction, and full reductions.
Native tests cover enormous extents and stress-sized exception sets without
materializing dense masks.
