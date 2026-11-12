# Sparse predicates and coordinate-based filtering of N-D COO arrays

`scifort_sparse_coo_nd_filter` extends `scifort_sparse_coo_nd%coo_nd_array`
with **sparse Boolean masks** and coordinate-set selection, without constructing
an N-dimensional dense mask. The source type remains the existing one-based,
real binary64 COO type, including duplicate coordinates and explicit stored zeros.

## Routines

| Routine | Operation |
|---|---|
| `coo_nd_compare_scalar(a,op,threshold,status,max_nnz)` | Compare **summed** values at each coordinate and return a canonical sparse mask. Supported `op` strings: `>`, `>=`, `<`, `<=`, `==`, `/=`, `abs>`. |
| `coo_nd_count_nonzero(a,status)` | Count the distinct nonzero coordinates after summing duplicates. Return `integer(int64)`, or `-1` on invalid input/overflow. |
| `coo_nd_filter_sparse(a,mask,status,complement,max_nnz)` | Preserve source entries whose coordinates have a **nonzero** value in the summed COO mask. `complement=.true.` selects the opposite set. |
| `coo_nd_filter_points(a,points,status,complement,max_nnz)` | Filter using one-based integer coordinates of shape `(ndim,npoints)`; duplicate query points are treated as set membership. |
| `coo_nd_filter_value(a,op,threshold,status,max_nnz)` | Apply a scalar predicate to aggregated source values, then retain the original entries at matching coordinates. |

The **Boolean representation** is a real COO `coo_nd_array`: `data=1.0_dp`
means *true*, and absent coordinates mean *false*. A user-supplied mask is
interpreted by first summing any duplicate entries: a coordinate is true
iff its resulting value is nonzero. This deliberately handles cancellation
(`+3,-3` means false) and differs from a logical OR of stored entries.

Because a sparse mask has an **implicit false background**, a comparison is
rejected with `sparse_invalid_input` if it would be true for an implicit zero:
for example `a == 0`, `a /= 1`, and `a >= 0`. This avoids the misleading
result of silently ignoring potentially billions of implicit elements.
For `a /= 0`, `a > 0`, or `abs(a) > 1`, the background remains false and
computing a sparse mask is well-defined. All thresholds and COO data must be
finite. The `abs>` operator means `abs(a) > threshold`.

The comparison result has **unique, lexicographically ordered** coordinates
and only `1.0_dp` values. Filtering preserves source entry order, duplicated
entries and explicitly stored zeros whenever their coordinate is selected;
`coo_nd_filter_value` tests the **coalesced** coordinate value but preserves
all original entries at passing coordinates (including a zero entry paired
with a nonzero one).

An optional `integer(int64) max_nnz` bounds stored output entries and defaults
to 1,000,000. On validation or allocation-cap failure, the result is an
unallocated COO object and `status` is `sparse_invalid_input`. No input is
changed. The query coordinate list and sparse masks must have the source rank
and shape; indices are **one-based**.

## Example

```fortran
program sparse_predicate_example
    use scifort_kinds, only : dp
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets
    use scifort_sparse_coo_nd_filter, only : coo_nd_compare_scalar, &
        coo_nd_filter_sparse, coo_nd_filter_points
    implicit none
    type(coo_nd_array) :: a, mask, b
    integer :: info

    a = coo_nd_from_triplets([2000000000,2000000000], &
        reshape([4,10, 100,2, 4,10],[2,3]), &
        [5.0_dp,-3.0_dp,2.0_dp],info)
    if (info /= 0) error stop 'construct'

    mask = coo_nd_compare_scalar(a,'>',0.0_dp,info)
    if (info /= 0) error stop 'comparison'
    b = coo_nd_filter_sparse(a,mask,info)
    if (info /= 0) error stop 'mask filtering'
    print *, size(b%data)  ! 2 stored entries at (4,10)

    b = coo_nd_filter_points(a,reshape([100,2],[2,1]),info)
    if (info /= 0) error stop 'coordinate filtering'
    print *, size(b%data)  ! 1 stored entry at (100,2)
end program sparse_predicate_example
```

## Complexity, compatibility and validation

`coo_nd_compare_scalar` and `coo_nd_count_nonzero` aggregate COO duplicates
using stable mergesort in `O(nnz log nnz)` time and `O(nnz)` auxiliary space.
Sparse filtering sorts mask coordinates once and uses binary searches to
select each stored source entry: approximately `O(m log m + n log m)` time,
`O(m + nout)` auxiliary space for `m` mask entries and `n` source entries.
A coordinate list is converted to a sparse membership mask. Nothing is
allocated proportional to the product of tensor dimensions.

These are **Fortran-specific sparse selection APIs**, not identical to the
current SciPy `coo_array` item-access or advanced-indexing syntax. The tests
construct SciPy COO reference inputs and compare numerical results, coordinate
ordering, explicit-zero preservation, and exact output stored-entry counts to
independent NumPy/dictionary calculations.

```console
python tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py --mode both --trials 300 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_filter_scipy.py --mode both --trials 450 --seed 274829 --floating
python tools/sparse_reference/run_native_sparse.py --mode checked
```

See [the validation record](sparse_coo_nd_filter_validation_2026-11-05.md).
