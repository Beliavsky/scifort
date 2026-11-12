# Complex N-dimensional COO arrays

`scifort_sparse_coo_nd_complex` adds `coo_nd_complex_array` with
`integer, allocatable :: shape(:), coords(:,:)` and
`complex(dp), allocatable :: data(:)`. Coordinates are **one-based** with
shape `(rank, nnz)`. Duplicates and explicitly stored zeros are supported.
A separate type preserves the existing real-valued COO API.

The storage module provides:

| Procedure | Purpose |
| --- | --- |
| `coo_nd_complex_from_triplets` | Build an array from shape, coordinates and complex values |
| `valid_coo_nd_complex`, `coo_nd_complex_is_canonical` | Validate storage and canonical C-axis order |
| `coo_nd_complex_sum_duplicates`, `coo_nd_complex_eliminate_zeros` | Canonicalize / discard explicit zeros |
| `coo_nd_complex_get` | Read one value (summing duplicate entries) |
| `coo_nd_complex_transpose`, `coo_nd_complex_reshape` | Reorder axes and reshape (C/F orders) |
| `coo_nd_complex_conjugate`, `coo_nd_complex_hermitian` | Conjugate all values / conjugate-transpose rank-two COO |
| `coo_nd_complex_from_dense`, `coo_nd_complex_to_dense` | Assumed-rank dense conversion (ranks 1–7) |
| `coo_nd_complex_from_dense_flat`, `coo_nd_complex_to_dense_flat` | Flat dense array conversion using C/F ordering |
| `coo_nd_complex_from_real`, `coo_nd_complex_to_real` | Convert from real COO or back if each imaginary component is zero |

`scifort_sparse_coo_nd_complex_ops` provides:

| Procedure | Purpose |
| --- | --- |
| `coo_nd_complex_sum_all` | Complex scalar sum of entries |
| `coo_nd_complex_sum_axis`, `coo_nd_complex_sum_axes` | Canonical sparse axis reductions; `keepdims` optional |
| `coo_nd_complex_broadcast_to` | NumPy-like trailing-axis singleton broadcasting |
| `coo_nd_complex_add`, `coo_nd_complex_multiply` | Elementwise arithmetic with broadcasting |
| `coo_nd_complex_tensordot` | Sparse contraction across selected one-based axes |
| `coo_nd_complex_vdot` | Sum of `conjg(a) * b` across corresponding tensor entries |

**Contraction semantics:** `tensordot` does not conjugate its first operand,
matching NumPy; `vdot` does. A contraction to a scalar is represented by
shape `[1]` to preserve positive-rank COO storage. All operations reject
nonfinite complex components (real or imaginary), detect nonfinite results,
and return `sparse_invalid_input` for invalid dimensions, axes or output caps.
As with the real COO implementation, `max_nnz` is an optional `integer(int64)`
cap, defaulting to **1,000,000** generated candidates. No dense data are
created in sparse arithmetic or contraction.

## Example

```fortran
program demo_complex_sparse
    use scifort_kinds, only : dp
    use scifort_sparse_coo_nd_complex
    use scifort_sparse_coo_nd_complex_ops
    implicit none
    type(coo_nd_complex_array) :: a, b, c
    complex(dp) :: inner
    integer :: info

    a = coo_nd_complex_from_triplets([2],reshape([1,2],[1,2]), &
        [cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp)],info)
    if (info /= 0) error stop 'construct'
    b = coo_nd_complex_conjugate(a,info)
    if (info /= 0) error stop 'conjugate'
    inner = coo_nd_complex_vdot(a,a,info)
    if (info /= 0) error stop 'vdot'
    print *, inner  ! (7.0, 0.0)
    c = coo_nd_complex_tensordot(a,b,[1],[1],info)
    if (info /= 0) error stop 'tensordot'
    print *, coo_nd_complex_get(c,[1])  ! (7.0, 0.0)
end program demo_complex_sparse
```

## Validation

```console
python tools/sparse_reference/compare_sparse_coo_nd_complex_scipy.py --mode both --trials 390
python tools/sparse_reference/compare_sparse_coo_nd_complex_scipy.py --mode both --trials 390 --floating
python tools/sparse_reference/run_native_sparse.py --mode both
```

The randomized test compares complex128 dense results from Fortran against
independently built SciPy `coo_array` instances and NumPy tensor operations.
The full regression runner compiles native sparse tests without FPM.
See the [validation report](sparse_coo_nd_complex_validation_2026-11-05.md).

**Limitations:** Only `complex(dp)` (complex128), and positive-rank arrays.
Python indexing, broadcasting of complex scalars and complex CSR/CSC
interoperation are not yet implemented. This module does not change the
real-valued COO types, or provide direct `einsum` notation.
