# Complement-aware N-D Boolean broadcasting and reductions: validation

Date: 2026-11-05

## Scope

New module: `src/sparse_coo_nd_bool_reduce.f90`.

Native regression: `test/test_sparse_coo_nd_bool_reduce.f90`.

NumPy reference: `tools/sparse_reference/compare_sparse_coo_nd_bool_reduce_numpy.py`
using `tools/sparse_reference/sparse_coo_nd_bool_reduce_reference_driver.f90`.

All code is original MIT-licensed Fortran/Python; NumPy is an independent
test reference, not part of the Fortran implementation.

## Checks performed

- GNU Fortran 14.2.0; NumPy 2.3.5.
- Bounds-checked `-O0 -g -fcheck=all -fbacktrace`: 41/41 native sparse
  regression executables passed.
- Optimized `-O3`: 41/41 native sparse regression executables passed.
- Two independent random seeds (`716`, `9276`), 500 cases/seed/configuration:
  1,000 NumPy comparisons per configuration; 2,000 total, all passed.
- Reference modes include: unary broadcasting; broadcasted Boolean
  AND/OR/XOR/difference; axis-wise `ANY`, `ALL`, and exact `integer(int64)`
  counts; `keepdims`; full Boolean reductions; ranks one through five;
  repeated exception coordinates, implicit-true/false backgrounds, and
  zero-sized dimensions.
- Native tests additionally exercise: 12,000 exception coordinates,
  a tensor with shape `[2000000000,2000000000,2000000000]`, true-background
  counts from a two-billion-length reduced axis, full logical reductions
  with total element count too large for int64, incompatible dimensions,
  output caps, and invalid axes.

## Reproduction

```console
python tools/sparse_reference/compare_sparse_coo_nd_bool_reduce_numpy.py --mode both --trials 500 --seed 716
python tools/sparse_reference/compare_sparse_coo_nd_bool_reduce_numpy.py --mode both --trials 500 --seed 9276
python tools/sparse_reference/run_native_sparse.py --mode both
```

The full project-wide `fpm test` was not run because `fpm` is unavailable
in this environment. The standalone native sparse regression suite covers
all existing sparse regression executables and the new test in both modes.

## Documented limitations

- Arrays have positive rank and use one-based default-integer coordinates.
  A rank-one axis reduction returns shape `[1]` rather than rank zero.
- Broadcasting operates on sparse exception coordinates; an expansion cap
  (default 1,000,000) prevents accidental huge intermediate allocations.
- These are Fortran-specific sparse Boolean APIs, not direct equivalents
  of public SciPy sparse Boolean mask routines.
- The new integer count arrays support exact single-axis counts only;
  multiaxis reduction and more general sparse integer arithmetic remain
  possible future extensions.
