# Block-diagonal and Khatri–Rao matrix operations (2026-11-09)

Two new MIT-licensed modules extend SciFort's `scipy.linalg` coverage:

| SciPy operation | SciFort module | Fortran interface |
|---|---|---|
| `scipy.linalg.block_diag(*arrs)` | `scifort_linalg_block_diag` | `as_block` and `block_diag` |
| `scipy.linalg.khatri_rao(a, b)` | `scifort_linalg_khatri_rao` | `khatri_rao` |

Both use caller-allocated outputs and the existing integer status conventions:
`linalg_status_success` (0) means success, and
`linalg_status_invalid_shape` (1) means an invalid dimension, unsupported
output conversion, or incompatible batch dimensions. **Outputs are zero-filled
on failure.** No existing modules or third-party libraries are changed.

## `block_diag`: heterogenous blocks with `as_block`

Since Fortran does not support Python-style variadic arguments with mixed
shapes and types, the API accepts an array of `type(linalg_block)` objects.
Construct each block with the overloaded pure function `as_block`:

```fortran
use scifort_kinds, only : dp
use scifort_linalg_block_diag, only : as_block, block_diag
real(dp) :: a(2,2), b(1,3), result(3,5)
integer :: status
! initialize a, b before calling
call block_diag([as_block(a), as_block(b)], result, status)
```

A scalar is normalized to a 1×1 block, a rank-1 vector is a 1×N row,
a rank-2 array is used unchanged, and a rank-3 array represents one
leading batch dimension `(batch, rows, cols)`. Blocks in an array can
combine integers, real numbers, and complex numbers (including scalar,
vector, and matrix blocks). To preserve precision, integer blocks are
stored as `integer(int64)` until copied into the selected output type.

Supported output kinds are `integer(int32)`, `integer(int64)`, `real(dp)`,
and `complex(dp)`, each with rank 2 or 3. The **caller chooses the output
kind**, equivalent to selecting the promoted NumPy dtype. An output of
insufficient type (e.g., complex input into a real output, or real input
into an integer output) is rejected. Narrowing `int64` to `int32` is
range checked. Like SciPy, an empty rank-1 vector is still a 1×0 block,
while a rank-2 matrix of shape `(0,0)` contributes zero rows and columns.
An empty list of blocks gives a 1×0 output matrix.

When an input block is rank 3, the output must also be rank 3. Different
batch dimensions broadcast if they are equal or one is 1; one batch
dimension may be zero, following NumPy's broadcasting rules. Unbatched
blocks are reused in each batch. Unlike SciPy's arbitrary leading batch
shapes, this API currently supports only **one leading batch dimension**.
An explicit singleton-batch input cannot be passed to a rank-2 output.

## `khatri_rao`: column-wise Kronecker product

Given `a(n,k)` and `b(m,k)`, the result has shape `(n*m,k)` and entries

\[
  \mathrm{out}((i-1)m+j,\ell) = a(i,\ell)b(j,\ell)
  \quad\text{for}\quad 1\le i\le n,\;1\le j\le m.
\]

```fortran
use scifort_kinds, only : dp
use scifort_linalg_khatri_rao, only : khatri_rao
real(dp) :: a(2,3), b(3,3), result(6,3)
integer :: status
! initialize a, b before calling
call khatri_rao(a, b, result, status)
```

`khatri_rao` implements all 16 input-kind combinations among
`integer(int32)`, `integer(int64)`, `real(dp)`, and `complex(dp)`,
with output promotion to the widest required kind (and an explicitly typed
output). Each input may have rank 2 or 3, yielding four rank combinations
per kind pair. Batched outputs have shape `(batch,n*m,k)` and support
broadcasting a singleton batch. In contrast to `block_diag`, a zero-sized
**batch** is rejected, matching SciPy. Zero row or zero column counts in
unbatched inputs are valid. Input column counts must agree, and product
row counts are checked against default-integer indexing bounds. Integer
multiplication itself is performed with the selected output integer kind;
large operands may overflow, as with ordinary Fortran integer arithmetic.

The algorithms require O(sum of block sizes) operations for `block_diag`
plus zeroing its output, and O(nmk) operations for `khatri_rao`.
They do not build temporary Kronecker-product matrices.

## Example and verification

- Runnable example: [`example/linalg_block_diag_khatri_rao_example.f90`](../example/linalg_block_diag_khatri_rao_example.f90)
- Native tests: [`test/test_linalg_block_diag_khatri_rao.f90`](../test/test_linalg_block_diag_khatri_rao.f90)
- Seeded cross-checks: [`tools/compare_linalg_block_diag_khatri_rao.py`](../tools/compare_linalg_block_diag_khatri_rao.py)
- Reference driver: [`reference/linalg_block_diag_khatri_rao_driver.f90`](../reference/linalg_block_diag_khatri_rao_driver.f90)
- Results: [`docs/linalg_block_diag_khatri_rao_validation_2026-11-09.md`](linalg_block_diag_khatri_rao_validation_2026-11-09.md)

Reference API documentation:
[`block_diag`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.block_diag.html)
and [`khatri_rao`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.khatri_rao.html).
