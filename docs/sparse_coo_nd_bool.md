# Complement-aware sparse Boolean masks

`scifort_sparse_coo_nd_bool` adds logical mask algebra for sparse N-dimensional
COO arrays. **This is a Fortran extension**, not a direct translation of a
named SciPy Boolean-mask class. Its numerical semantics are checked against
independent NumPy/SciPy dense reference arrays.

A `coo_nd_bool_mask` holds:

```fortran
type :: coo_nd_bool_mask
    integer, allocatable :: shape(:)       ! Nonnegative extents, rank >= 1
    integer, allocatable :: coords(:,:)    ! Unique, sorted, one-based exceptions
    logical :: default_value = .false.    ! Value outside the exception set
end type coo_nd_bool_mask
```

Each coordinate in `coords(:,k)` has Boolean value `.not. default_value`;
all other elements have `default_value`. Unlike conventional sparse COO,
this representation can store the logical complement of a sparse set without
enumerating the implicit background. It allocates storage proportional to the
number of exceptions, not the product of the dimensions.

## Public routines

| Routine | Purpose |
|---|---|
| `coo_nd_bool_from_points(shape,coords,status,default_value,max_exceptions)` | Construct a mask from exception coordinates; deduplicate and sort |
| `coo_nd_bool_from_coo(a,status,complement,max_exceptions)` | From real COO: sum duplicates, treat nonzero as true; optionally negate |
| `coo_nd_bool_compare_scalar(a,op,threshold,status,max_exceptions)` | Compare real COO with a scalar, including implicit zero elements |
| `coo_nd_bool_not(a,status)` | Logical NOT in constant additional space |
| `coo_nd_bool_and(a,b,status,max_exceptions)` | Logical intersection |
| `coo_nd_bool_or(a,b,status,max_exceptions)` | Logical union |
| `coo_nd_bool_xor(a,b,status,max_exceptions)` | Exclusive OR |
| `coo_nd_bool_difference(a,b,status,max_exceptions)` | `a AND NOT b` |
| `coo_nd_bool_get(a,point,status)` | Boolean lookup at a one-based coordinate tuple |
| `coo_nd_bool_count_true(a,status)` | Exact count of true entries as `integer(int64)` when representable |
| `coo_nd_bool_filter(a,mask,status,max_nnz)` | Retain real COO entries selected by a compatible Boolean mask |
| `coo_nd_bool_to_coo(mask,status,max_nnz)` | Convert an implicit-false mask to real COO with values `1.0_dp` |
| `valid_coo_nd_bool_mask(mask)` | Validate shape, coordinate bounds, sort order and uniqueness |

Binary mask operations require **identical shapes**; broadcasting is not
performed. Each output uses the result of the Boolean operation on the input
backgrounds as its default. A sorted set merge produces only the coordinates
where the output differs from that default, with complexity `O(n_a+n_b)`
(after validation). There is no allocation proportional to a dense tensor.

`coo_nd_bool_compare_scalar` accepts `>`, `>=`, `<`, `<=`, `==`, `/=`, and
`abs>`, with the last meaning `abs(a) > threshold`. All real data and thresholds
must be finite. First, duplicate COO entries are summed. At unstored
coordinates the numerical value is zero: e.g. `a == 0` yields a default-true
mask with an exception for every coordinate whose summed value is nonzero.
This avoids the restriction of the older
[`coo_nd_compare_scalar`](sparse_coo_nd_filter.md), which can only produce
implicit-false sparse masks.

`coo_nd_bool_from_coo` interprets a real mask's *summed* nonzero values as
true. Thus `2,-2` at one coordinate cancel to false, and an explicitly
stored zero is false. `coo_nd_bool_filter` preserves the source's entry
order, duplicate entries and explicitly stored zeros at selected coordinates.
It does not sum the source entries.

## Example

```fortran
program complement_mask_example
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_sparse_coo_nd, only : coo_nd_array, coo_nd_from_triplets
    use scifort_sparse_coo_nd_bool, only : coo_nd_bool_mask, &
        coo_nd_bool_compare_scalar, coo_nd_bool_filter, coo_nd_bool_count_true
    implicit none
    type(coo_nd_array) :: a, b
    type(coo_nd_bool_mask) :: zero_mask
    integer(int64) :: n
    integer :: status

    a = coo_nd_from_triplets([1000000,1000000], &
        reshape([11,12, 2,9],[2,2]),[7.0_dp,0.0_dp],status)
    if (status /= 0) error stop 'COO construction'
    zero_mask = coo_nd_bool_compare_scalar(a,'==',0.0_dp,status)
    if (status /= 0) error stop 'zero predicate'
    n = coo_nd_bool_count_true(zero_mask,status)
    if (status /= 0) error stop 'count true'
    print *, n  ! 999999999999: only (11,12) is nonzero
    b = coo_nd_bool_filter(a,zero_mask,status)
    if (status /= 0) error stop 'filter'
    print *, size(b%data)  ! 1: explicit zero at (2,9) is retained
end program complement_mask_example
```

## Failure rules and limitations

- All coordinates are **one-based** and use default Fortran `integer`;
  dimensions themselves can be as large as `huge(0)`.
- The default `max_exceptions` and `max_nnz` limits are 1,000,000. Negative
  limits, malformed inputs, nonfinite arithmetic, or exceeded limits return
  `sparse_invalid_input`. Failure results are unallocated; inputs are never
  modified. The caps may be increased explicitly, subject to memory.
- Converting a default-true mask to a conventional COO array is rejected if
  all shape extents are positive; it would require enumerating the complement.
  A zero-element array can be materialized as an empty COO regardless of its
  Boolean default.
- `coo_nd_bool_count_true` detects signed 64-bit overflow for default-true
  masks. For default-false masks, only the stored exception count is needed,
  so even a shape whose total element count overflows 64-bit can be counted.
- Scalar comparisons and real COO conversions aggregate duplicates using
  the existing stable COO sort. A duplicate sum overflowing to nonfinite is
  rejected instead of silently generating an invalid Boolean mask.
- This module provides scalar comparisons, equal-shaped Boolean algebra and
  sparse COO filtering. **Logical broadcasting and Boolean reductions** are
  now available in the separate
  [`scifort_sparse_coo_nd_bool_reduce` module](sparse_coo_nd_bool_reduce.md).
  Full NumPy-style arbitrary Boolean-index assignment is not provided.

## Validation

Run the focused standalone Fortran test in both compiler configurations via
`python tools/sparse_reference/run_native_sparse.py --mode both`, and use

```console
python tools/sparse_reference/compare_sparse_coo_nd_bool_scipy.py --mode both --trials 900 --seed 20261105
```

The Python script compiles a streaming Fortran comparison program and
independently builds full Boolean reference grids from SciPy real COO arrays
and NumPy logical operations, checking exact coordinates, default values,
counts, COO stored-entry order and explicitly stored zeros. See the
[validation report](sparse_coo_nd_bool_validation_2026-11-05.md).
