# N-dimensional COO reductions, broadcasting and sparse tensor contractions

`scifort_sparse_coo_nd_ops` extends the one-based, real binary64
`scifort_sparse_coo_nd%coo_nd_array` type. It performs reductions, elementwise
broadcasting operations and `tensordot` without forming dense intermediate
tensors. It does not require new libraries or a separate build system.

These are SciFort operations with NumPy-like semantics. SciPy 1.17's
`coo_array.sum(axis=...)` returns dense data, and `coo_array.multiply` does
not yet provide general N-D broadcasting. Numerical tests therefore compare
SciFort's sparse results with NumPy operations on independently constructed
SciPy COO arrays, rather than claiming one-to-one API parity for operations
that SciPy does not expose.

## Procedures

| Procedure | Behavior |
|---|---|
| `coo_nd_sum_all(a, status)` | Sum all stored contributions, returning a real scalar. Duplicate coordinates and explicit zeros contribute as expected. |
| `coo_nd_sum_axis(a, axis, keepdims, status)` | Reduce a single **one-based** axis into canonical COO. |
| `coo_nd_sum_axes(a, axes, keepdims, status)` | Reduce one or more one-based axes; `axes=[]` coalesces duplicates without changing shape. |
| `coo_nd_broadcast_to(a, new_shape, status, max_nnz)` | Expand aligned trailing singleton axes, preserving entry order, duplicates and stored zeros. |
| `coo_nd_add(a, b, status, max_nnz)` | Elementwise COO addition after broadcasting. Canonicalizes duplicates, retaining zeros after cancellation. |
| `coo_nd_multiply(a, b, status, max_nnz)` | Elementwise COO product with broadcasting. Canonicalizes operands and intersects nonzero coordinates. Excludes result zeros. |
| `coo_nd_tensordot(a, b, axes_a, axes_b, status, max_nnz)` | Match and sum over chosen pairs of one-based axes, returning canonical COO and excluding result zeros. |

An optional integer `status` is set to `sparse_success` (0), or to
`sparse_invalid_input` (1) for invalid shapes/axes, arithmetic overflow or
an exceeded storage limit. Failed array-valued operations return an
unallocated result (check `valid_coo_nd`). `sum_all` returns `0.0_dp` on
failure. Input arrays remain unchanged.

The source type must have rank at least 1, as in
`scifort_sparse_coo_nd`. A reduction of every axis requires
`keepdims=.true.` (a rank-preserving sparse result), or the scalar-valued
`coo_nd_sum_all` procedure. A complete tensor contraction is represented as
an N-D COO array of shape `[1]`; use `coo_nd_get(result, [1])` to recover
the scalar. This maintains compatibility with the positive-rank storage type.

`coo_nd_broadcast_to` aligns **trailing dimensions**, adding virtual leading
axes of length one. Its returned COO retains duplicate coordinates and
explicit zeros; binary operations coalesce before returning. A size-one axis
can broadcast to length zero, yielding no entries. Broadcasting to a shape
with fewer axes than the input is invalid.

`max_nnz` is an optional `integer(int64)` storage cap. Its default is
**1,000,000** generated entries, preventing accidental allocation when
expanding singleton dimensions or enumerating tensor-product contributions.
For addition it bounds the total concatenated contributions; for
multiplication it bounds each broadcasted operand; for contraction it bounds
matching nonzero entry pairs *before* output aggregation. Raise the cap
explicitly for large tasks. Intermediate dense arrays are never allocated.

## Example

```fortran
program sparse_tensor_example
    use scifort_kinds, only : dp
    use scifort_sparse_coo_nd, only : coo_nd_array, &
        coo_nd_from_triplets, coo_nd_get
    use scifort_sparse_coo_nd_ops, only : coo_nd_sum_axes, &
        coo_nd_broadcast_to, coo_nd_tensordot
    implicit none
    type(coo_nd_array) :: a, b, reduced, tiled, contracted
    integer :: info

    ! A is shape (2, 3, 2); only A(1, 2, 1) = 4 is stored.
    a = coo_nd_from_triplets([2,3,2],reshape([1,2,1],[3,1]),[4.0_dp],info)
    if (info /= 0) error stop 'A construction'
    reduced = coo_nd_sum_axes(a,[3],status=info)
    if (info /= 0) error stop 'reduction'
    print *, coo_nd_get(reduced,[1,2])   ! 4.0

    tiled = coo_nd_broadcast_to(reduced,[5,2,3],info)
    if (info /= 0) error stop 'broadcast'
    print *, coo_nd_get(tiled,[5,1,2])  ! 4.0

    ! B is shape (2, 1); B(1, 1) = 3.
    b = coo_nd_from_triplets([2,1],reshape([1,1],[2,1]),[3.0_dp],info)
    contracted = coo_nd_tensordot(a,b,[3],[1],info)
    if (info /= 0) error stop 'tensor contraction'
    print *, coo_nd_get(contracted,[1,2,1])  ! 12.0
end program sparse_tensor_example
```

## Computational approach

Reductions reindex existing entries and merge duplicate coordinates with the
existing O(nnz log nnz) stable mergesort. Broadcasting enumerates the minimum
necessary stored-entry replicas (not the full dense shape). Addition
concatenates then canonicalizes; multiplication canonicalizes and merges
sorted nonzero coordinate streams. Tensor contraction groups nonzero entries
by the contracted coordinates, using a stable key sort and merge join. It
enumerates only cross-products within matching key groups, then coalesces
contributions to equal output coordinates. The complexity is dominated by
sorting and the number of matched nonzero pairs, not by the dense array size.

## Validation

```console
python tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py --mode both --trials 900 --seed 20261105
python tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py --mode both --trials 600 --seed 274829
python tools/sparse_reference/compare_sparse_coo_nd_ops_scipy.py --mode both --trials 750 --seed 11111 --floating
python tools/sparse_reference/run_native_sparse.py --mode both
```

Additional validation details: [N-D COO operations test log](sparse_coo_nd_ops_validation_2026-11-05.md).
Limitations: only finite real(dp) data, rank >= 1, and one-based integer
coordinate arrays. There is no general sparse `einsum` parser,
automatic scalar broadcasting, or complex data in this real-valued module. Complex COO tensor operations are now provided by [`scifort_sparse_coo_nd_complex_ops`](sparse_coo_nd_complex.md). COO indexing and
slicing are supported in [`scifort_sparse_coo_nd_index`](sparse_coo_nd_index.md). Dense conversion
is used only in the optional independent test harness.
