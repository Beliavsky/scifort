## N-dimensional coordinate sparse arrays

`scifort_sparse_coo_nd` supports arbitrary-rank real COO arrays with one-based
coordinate vectors, C/F-order reshaping, transpose, duplicate canonicalization
and dense/2D-COO conversion. See [N-dimensional COO API](sparse_coo_nd.md).

## Incrementally editable DOK and LIL matrices

`scifort_sparse_dok_lil` provides DOK hash-table and LIL sorted-row
storage as counterparts to `scipy.sparse.dok_array` and `lil_array`.
They permit efficient single-coordinate updates, preserve stored zeros when
converting from compressed/coordinate formats, and support CSR/CSC output.
See [complete DOK/LIL API and tests](sparse_dok_lil.md).

# Generic rectangular sparse matrices (`scipy.sparse` subset)

`scifort_sparse` supplies real-binary64, **one-based, canonical CSR and CSC**
matrix types and the elementary operations needed to build sparse numerical
algorithms without dense storage. It is independent of the earlier
`scifort_sparse_csgraph` graph module. No new runtime dependencies or
third-party implementation code are required.

## Kronecker sum: `csr_kronsum`

`csr_kronsum(a, b, status)` takes canonical **square** real CSR matrices
`a` (m-by-m) and `b` (n-by-n) and constructs their SciPy-style Kronecker sum

`kron(I_n, a) + kron(b, I_m)`

as a canonical `(m*n)`-by-`(m*n)` CSR matrix. Indices and row pointers remain
one-based. No dense `m*n`-by-`m*n` intermediate or temporary Kronecker matrices
are constructed. Instead, it merges the sorted row streams of the two terms,
coalescing overlapping diagonal locations and preserving explicitly stored
zeros, including zeros resulting from exact cancellation.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense, csr_kronsum, sparse_success
implicit none
type(csr_matrix) :: ax, ay, lap2d
integer :: status
ax = csr_from_dense(reshape([2.0_dp,-1.0_dp,-1.0_dp,2.0_dp],[2,2]))
ay = csr_from_dense(reshape([3.0_dp,-1.0_dp,-1.0_dp,3.0_dp],[2,2]))
lap2d = csr_kronsum(ax, ay, status)
if (status /= sparse_success) error stop 'invalid Kronecker sum'
```

Invalid inputs include nonsquare operands, invalid or nonfinite CSR data,
unrepresentable result dimensions or stored-entry counts, and nonfinite sums
of overlapping diagonal entries. An empty 0-by-0 factor produces an empty
result. Complexity is O(n*nnz(a) + m*nnz(b) + m*n), with row-pointer storage
O(m*n) and output/work storage bounded by the sum of the input terms' stored
entries. This function always returns CSR; unlike SciPy, it does not select
an output format or handle complex matrices.

Validation: `test/test_sparse_kronsum.f90` and
`tools/sparse_reference/compare_sparse_bicg_kronsum_scipy.py`.

## Sparse random matrices and `diags_array`-style construction (2026-11-05)

The new **`scifort_sparse_constructors`** module is layered on top of
`scifort_sparse` and SciFort's existing `scifort_random` generator, without
changing existing CSR or CSC APIs. It provides real-binary64 canonical CSR
matrices with one-based indexing:

| API | Description |
|---|---|
| `csr_diags_array(diagonals,offsets,n_rows,n_cols,status)` | Build matrices from unpadded, ragged diagonals, with optional square-shape inference and scalar broadcasting when an explicit shape is supplied |
| `csr_random(n_rows,n_cols,density,state,status,distribution,data_values)` | Sample exactly `round(density*n_rows*n_cols)` distinct coordinates, filling values with uniform or normal variates, or caller-supplied finite data |
| `csr_rand(n_rows,n_cols,density,state,status)` | Uniform-data shorthand for `csr_random` |
| `csr_random_array(...)` | Alias of `csr_random`, emphasizing sparse-array usage |
| `sparse_diagonal` | Wrapper for a one-dimensional allocatable diagonal, usable with an array of derived-type constructors |

### Example

```fortran
use scifort_kinds, only : dp
use scifort_random, only : rng_state, rng_seed
use scifort_sparse, only : csr_matrix, sparse_success
use scifort_sparse_constructors, only : sparse_diagonal, csr_diags_array, csr_rand
implicit none
type(rng_state) :: rng
type(csr_matrix) :: a, r
integer :: status

! One-based CSR banded matrix. Scalars broadcast with an explicit shape.
a = csr_diags_array([sparse_diagonal([1.0_dp]), &
                     sparse_diagonal([-2.0_dp]), &
                     sparse_diagonal([1.0_dp])], [-1,0,1], 100,100,status)
if (status /= sparse_success) error stop 'invalid diagonal input'
call rng_seed(rng,12345)
r = csr_rand(1000,2000,0.001_dp,rng,status)  ! 2000 distinct stored cells.
if (status /= sparse_success) error stop 'invalid random matrix input'
```

**Semantics:** Diagonal offsets are column minus row, as in SciPy. Without a
shape, the result is square with side length equal to the first diagonal's
length plus the absolute value of its offset. With a shape, singleton input
diagonals broadcast to the required diagonal length; longer diagonal vectors
are truncated, and insufficient non-scalar diagonal vectors are rejected.
Repeated offsets and non-finite values are rejected. Explicit zero-valued
diagonal entries remain stored.

The random constructor chooses coordinates **without replacement**, using
Floyd sampling and an open-addressed hash set. Temporary memory is
O(nnz + number of rows), not O(number of rows × number of columns). The
number of stored coordinates is calculated with nearest-integer, ties-to-even
rounding, as in SciPy. Random generation consumes an explicit `rng_state`,
uses SciFort's MRG32k3a PRNG, and commits state only after success. `density`
0 and 1 are supported, as are empty matrices and large logical column counts.
The uniform distribution uses `(0,1)` and the normal distribution uses
Box–Muller. `data_values` is an array of explicit finite values in sampled
coordinate order; it is **not** a Python callable/data sampler.

**Differences from Python:** The routines always produce CSR (not a selectable
sparse format), and require the explicit SciFort `rng_state`. An identical
numeric seed does not produce NumPy/SciPy's exact random coordinates or values.
Scalar broadcasting requires both `n_rows` and `n_cols`. No complex/integer
sparse dtype or custom callback-based value sampler is provided.

**Verification:** `test/test_sparse_constructors.f90` checks structural
integrity, reproducibility, explicit zeros, invalid-input RNG transactionality,
large sparse shapes, uniform/normal distribution sanity, scalar broadcasting,
shape inference, and boundary cases. The randomized comparison driver
`tools/sparse_reference/compare_sparse_constructors_scipy.py` checks
independently generated diagonals against SciPy's `diags_array`, and validates
random-matrix count, uniqueness, and value properties. Run:

```bash
python tools/sparse_reference/compare_sparse_constructors_scipy.py --mode both --cases 250
python tools/sparse_reference/run_native_sparse.py --mode checked
```

## Structural deletion, sparse-to-sparse assignment, and CSR reshape (2026-11-05)

The following routines extend canonical CSR indexing without dense work arrays:

| API | Operation |
|---|---|
| `csr_delete_entries(a,rows,columns,status)` | Remove *stored* elements at paired one-based coordinates; missing coordinates and repeated deletion requests are permitted |
| `csr_delete_rows(a,rows,status)` | Remove selected complete rows and renumber the remaining rows; repeated row indices are accepted |
| `csr_delete_columns(a,columns,status)` | Remove columns, renumber retained columns, and preserve their relative order without allocating an array of length `a%n_cols` |
| `csr_reshape(a,n_rows,n_cols,status)` | Return a canonical CSR view of the same row-major flattened elements in a different shape (new allocation; the source is unchanged) |
| `csr_assign(a,rows,columns,sparse_values,status)` | Efficient sparse-to-sparse assignment, overwriting implicit source zeros by removing existing target elements |

These APIs validate the entire request before changing an input matrix.
They preserve explicitly stored zeros in retained entries, and preserve the
number and value of explicit source zeros in sparse assignment unless their
selected coordinates are overwritten by a later repeated index. Deleting
entries does **not** change matrix shape, whereas deleting whole rows/columns
does. `csr_reshape` uses **C-style row-major** flattening (to match
`scipy.sparse.csr_matrix.reshape(order='C')`), not Fortran column-major
`reshape`. Reshape supports empty dimensions when element counts agree.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_coo, csr_delete_columns, &
    csr_assign, csr_reshape, sparse_success
implicit none
type(csr_matrix) :: a, rhs, b
integer :: info
a = csr_from_coo(2,4,[1,2],[1,4],[2.0_dp,3.0_dp])
rhs = csr_from_coo(2,2,[2],[2],[5.0_dp])
call csr_assign(a,[1,2],[1,4],rhs,info)  ! Removes old (1,1), sets (2,4)=5
if (info /= sparse_success) error stop 'sparse assignment failed'
call csr_delete_columns(a,[2,3],info)   ! A now has shape (2,2)
b = csr_reshape(a,1,4,info)            ! Row-major reshape to (1,4)
```

Standalone verification is available in
`test/test_sparse_deletion_reshape.f90` and
`tools/sparse_reference/compare_sparse_deletion_scipy.py`.

## General block assembly and CSR assignment (2026-11-05)

These new APIs extend the existing assembly and indexing operations:

| API | Operation |
|---|---|
| `csr_bmat(blocks,status,row_sizes,col_sizes)` | Assemble a rank-two array of CSR blocks, with `csr_matrix()` representing an absent (zero) block |
| `csr_block(...)` | Alias for `csr_bmat` (all supported arguments are identical) |
| `csr_set(a,row,column,value,status)` | Insert or replace a single one-based entry in-place |
| `csr_assign_coo(a,rows,columns,values,status)` | Apply coordinate replacements in order; duplicate coordinates take the **last** value, unlike the summation convention of `csr_from_coo` |
| `csr_assign(a,rows,columns,values,status)` | Cartesian selected-row/column assignment from a **scalar**, a real `values(:,:)`, or a canonical `csr_matrix` |

`csr_bmat` infers each block-row height and block-column width from present
matrices, checking consistency. Missing blocks must be represented by a default
`csr_matrix()` with *unallocated* components, not by a valid zero-sized CSR
matrix (which is treated as present). If a complete block row/column is absent,
pass explicit `row_sizes`/`col_sizes` so its dimensions are unambiguous. Explicit
sizes, when provided, must agree with every present block. Dimensions and total
stored entries are checked for default-integer overflow. The result is canonical
CSR, built directly in row order and preserving stored zero entries, with no
dense intermediate.

For assignment, all indices are one-based; `csr_assign` uses the Cartesian
product of `rows` and `columns` (not paired/fancy indexing). The dense source
must be shape `(size(rows), size(columns))`; the sparse source must have those
same dimensions. Each selected cell is replaced, including implicit zeros from
sparse sources. Repeated target indices are resolved in **Fortran column-major**
source traversal order (last update wins); SciPy does not promise the same tie
order for advanced indexing. Explicit source zeros remain stored; **implicit
sparse-source zeros remove any preexisting target entry without allocating a
stored zero**. Invalid values, shapes, indices or malformed inputs produce
`sparse_invalid_input` without changing `a`. Coordinate scatter uses O(nnz +
updates) memory, **independent of logical column count**. Sparse-to-sparse
Cartesian assignment now uses O(nnz(target) + nnz(source) + number of selected
rows and columns) temporary storage rather than expanding the Cartesian
product. Dense and scalar Cartesian assignment still materialize
`size(rows) * size(columns)` updates, so their very large selections should
be chunked.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_bmat, csr_from_dense, &
    csr_assign, csr_set, csr_to_dense, sparse_success
implicit none
type(csr_matrix) :: a, b, result, grid(2,2)
integer :: info
a = csr_from_dense(reshape([1.0_dp,2.0_dp],[2,1]))
b = csr_from_dense(reshape([3.0_dp,4.0_dp],[1,2]))
grid = csr_matrix()               ! Unallocated CSR = absent block
grid(1,1) = a
grid(2,2) = b
result = csr_bmat(grid,info)
if (info /= sparse_success) error stop 'block assembly failed'
call csr_set(result,1,3,5.0_dp,info)
call csr_assign(result,[2,3],[1,2],0.0_dp,info)
if (info /= sparse_success) error stop 'sparse assignment failed'
```

The new routines use no external solver or matrix dependencies. Native and
randomized SciPy/NumPy comparisons are at:
`test/test_sparse_block_assignment.f90` and
`tools/sparse_reference/compare_sparse_block_assignment_scipy.py`.

## Sparse matrix assembly and indexing (2026-11-05)

All operations below accept canonical `csr_matrix` inputs and return a new
canonical CSR matrix. The constructors operate on sparse arrays directly,
without allocating dense intermediate arrays. They preserve stored zeros;
invalid shapes, coordinates, nonfinite Kronecker products, and dimension/nnz
overflow yield `sparse_invalid_input` and an unallocated CSR result.

| API | Operation |
|---|---|
| `csr_hstack(blocks,status)` | Concatenate nonempty `type(csr_matrix)` arrays along columns; each block has the same number of rows |
| `csr_vstack(blocks,status)` | Concatenate nonempty arrays along rows; each block has the same number of columns |
| `csr_block_diag(blocks,status)` | Put rectangular CSR blocks along the diagonal; an empty array produces a valid `(0,0)` result |
| `csr_kron(a,b,status)` | Sparse Kronecker product (including explicit zero-valued structural products) |
| `csr_get(a,row,column,status)` | Return one element, zero if unstored; validates one-based coordinates |
| `csr_select(a,rows,columns,status)` | Cartesian selection of arbitrary row/column index arrays, including reordering and repetition |
| `csr_slice(a,row_first,row_last,col_first,col_last,status,row_step,col_step)` | Inclusive one-based rectangular slices, with optional nonzero signed strides |

Fortran bounds in `csr_slice` are **inclusive** (unlike Python's exclusive
slice stop). Positive strides follow `first:last:stride`; negative strides
reverse the traversal. A positive-step empty slice may use `(n+1,n)` or
`(1,0)`; an empty dimension itself can also be sliced. For general index
arrays `csr_select` is more flexible. A repeated index creates a distinct
output row or column. The selector does not allocate an auxiliary array
proportional to the *logical* input column count; this matters for huge,
mostly empty sparse matrices.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_dense, csr_hstack, csr_kron, &
    csr_slice, csr_select, csr_to_dense, sparse_success
implicit none
type(csr_matrix) :: a, b, joined, product, picked
integer :: info
a = csr_from_dense(reshape([1.0_dp,2.0_dp,0.0_dp,3.0_dp],[2,2]))
b = csr_from_dense(reshape([4.0_dp,5.0_dp],[2,1]))
joined = csr_hstack([a,b],info)         ! Shape (2,3)
if (info /= sparse_success) error stop 'hstack failed'
picked = csr_slice(joined,2,1,3,1,info,row_step=-1,col_step=-1)
product = csr_kron(a,b,info)           ! Shape (4,2)
picked = csr_select(joined,[2,1,2],[3,1,3],info)
```

Only real binary64 CSR matrices are supported directly. Convert to/from CSC
with `csc_to_csr` and `csr_to_csc` if that storage format is desired.
`csr_hstack([])` and `csr_vstack([])` return an invalid-input status;
`csr_block_diag([])` instead returns a valid `(0,0)` CSR matrix. This last
case follows the intended empty direct-sum identity even though SciPy 1.17
raises an exception for its own `block_diag([])` implementation.

## Diagonal (DIA) and block CSR (BSR) formats

`dia_matrix` and `bsr_matrix` are additional binary64 representations,
independent of the existing canonical CSR and CSC types.

| API | Operation |
|---|---|
| `csr_to_dia(a,status)`, `dia_to_csr(a,status)` | CSR/DIA conversion without a dense matrix |
| `dia_from_data(n_rows,n_cols,offsets,data,status)` | Construct from SciPy-style column-aligned diagonal data |
| `dia_from_dense`, `dia_to_dense`, `dia_matvec`, `valid_dia` | DIA construction, conversion, multiplication, validation |
| `csr_to_bsr(a,block_rows,block_cols,status)`, `bsr_to_csr(a,status)` | CSR/BSR conversion with rectangular dense blocks |
| `bsr_from_dense`, `bsr_to_dense`, `bsr_matvec`, `valid_bsr` | BSR construction, conversion, multiplication, validation |

DIA offsets are *column minus row*, sorted strictly increasing, and the
Fortran layout `data(k,j)` stores element `(j-offsets(k),j)` at **column** `j`.
Invalid/off-matrix padding is ignored; repeated diagonal offsets are rejected.
The DIA storage size is `number_of_diagonals * n_cols`, which can be much
larger than the original CSR when the occupied diagonals are widely spaced.
Conversions expand every in-bound diagonal slot, including zero values, into
CSR storage; they preserve the matrix values, not necessarily the count of
explicitly stored zeros.

BSR uses `block_rows` × `block_cols` dense blocks in `data(:,:,k)`,
with one-based block-column indices and one-based block-row pointers.
Block columns are strictly increasing within each row. The input shape must
be divisible by the corresponding block dimensions (partial edge blocks are
not yet supported). Converting BSR to CSR expands all entries in stored
blocks, including their zero padding. BSR block products use dense block
matvecs directly rather than reconstructing CSR.

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, dia_matrix, bsr_matrix, &
    csr_from_dense, csr_to_dia, csr_to_bsr, dia_matvec, bsr_matvec
implicit none
type(csr_matrix) :: a
type(dia_matrix) :: d
type(bsr_matrix) :: b
real(dp), allocatable :: y(:)
a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2]))
d = csr_to_dia(a)
b = csr_to_bsr(a,2,2)
y = dia_matvec(d,[2.0_dp,3.0_dp])
y = bsr_matvec(b,[2.0_dp,3.0_dp])
```

## Example

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csr_from_coo, csr_transpose, &
    csr_matmul, csr_matvec, csr_to_dense, sparse_success
implicit none

type(csr_matrix) :: a, a_transpose, gram
real(dp), allocatable :: y(:), gram_dense(:,:)
integer :: status

! Rectangular 3-by-2 matrix; repeated (1,2) coordinates are summed.
a = csr_from_coo(3,2,[1,1,3],[2,2,1], &
    [2.0_dp,3.0_dp,4.0_dp],status)
if (status /= sparse_success) error stop 'invalid coordinates'
y = csr_matvec(a,[1.0_dp,2.0_dp],status)       ! [10,0,4]
a_transpose = csr_transpose(a,status)
gram = csr_matmul(a_transpose,a,status)
gram_dense = csr_to_dense(gram,status)          ! diagonal([16,25])
```

## Representation and invariants

The public `type(csr_matrix)` stores `n_rows`, `n_cols`, `indptr(:)`,
`indices(:)` and `data(:)`. CSR offsets and coordinates are one-based:

- `size(indptr) = n_rows+1`, `indptr(1) = 1`, and
  `indptr(n_rows+1) = size(data)+1`.
- A row `i` occupies positions `indptr(i):indptr(i+1)-1`; its columns
  `indices` are **strictly increasing** and lie in `1:n_cols`.
- Duplicate coordinates are combined by addition during COO construction.
  An explicitly stored zero is permitted and remains stored until removed
  by `csr_eliminate_zeros`.
- All stored values must be finite binary64 numbers. NaN, infinities,
  bad dimensions, noncanonical storage, and nonfinite computed results
  are rejected with `sparse_invalid_input` (1); success is `sparse_success` (0).
- Zero dimensions, including `(0,n)` and `(m,0)`, are supported.
- A failed constructor or arithmetic function returns a default,
  unallocated `csr_matrix`. Invalid vector/matrix conversions return
  empty arrays. This allows callers to check the status before using output.

The public `type(csc_matrix)` has the same components, but its pointer array
has length `n_cols+1`; `indptr(j):indptr(j+1)-1` stores one column whose
**row** indices are strictly increasing. `valid_csc` checks analogous
canonical-storage invariants. `csr_to_csc` and `csc_to_csr` run in
O(nnz + n_rows + n_cols) time without dense intermediates and preserve
explicit stored zeros and rectangular dimensions.

This type models numeric **matrices**, not adjacency graphs: missing entries
mean zeros. `scifort_sparse_csgraph` remains the API for weighted graph
edges, where a stored zero can denote a real edge and parallel arcs may
remain separate. Do not interchange these types without deliberately
choosing the desired semantics.

## Functions

| API | Purpose |
|---|---|
| `csr_from_coo(n_rows,n_cols,rows,columns,values,status)` | Construct from unsorted triplets; sum duplicates, sort columns |
| `csr_from_dense(dense,status)` | Store the exact nonzero entries of a finite rank-2 array |
| `csr_to_coo(a,rows,columns,values,status)` | Extract triplets, including explicitly stored zeros |
| `csr_to_dense(a,status)` | Allocate an ordinary rank-2 Fortran array |
| `valid_csr(a)` | Check the full structural and finiteness contract |
| `csr_eye(n_rows,n_cols,offset,status)` | Unit-valued diagonal, optionally rectangular or shifted |
| `csr_transpose(a,status)` | Canonical sparse transpose |
| `csr_matvec(a,x,status)` | Sparse matrix-vector product |
| `csr_matmul(a,b,status)` | Sparse matrix-matrix product |
| `csr_add(a,b,status,alpha,beta)` | `alpha*a + beta*b`; scales default to 1 |
| `csr_elementwise_multiply(a,b,status)` | Entrywise product, not matrix multiplication |
| `csr_scale(a,factor,status)` | Scale all stored values; preserve stored zeros |
| `csr_eliminate_zeros(a,status)` | In-place removal of entries exactly equal to zero |
| `csc_from_coo(n_rows,n_cols,rows,columns,values,status)` | Construct canonical CSC from unsorted triplets, summing duplicate coordinates |
| `csc_from_dense(dense,status)` | Construct CSC, omitting numerical zeros |
| `csr_to_csc(a,status)` and `csc_to_csr(a,status)` | Sparse CSR/CSC conversion, preserving stored zeros |
| `csc_to_dense(a,status)` | Convert canonical CSC to dense values |
| `csc_matvec(a,x,status)` | Matrix-vector product from CSC storage |
| `csc_transpose(a,status)` | Transpose a CSC matrix without a dense intermediate |
| `valid_csc(a)` | Check full CSC structure and finiteness |
| `csr_diagonal(a,offset,status)` | Extract one diagonal as a dense vector, filling unstored entries with zero |
| `csr_diags(n_rows,n_cols,offsets,diag_ptr,values,status)` | Assemble multiple packed diagonals with varying lengths |
| `csr_setdiag(a,values,offset,status)` | Replace an entire diagonal in place; preserve stored zeros |
| `csr_triu(a,offset,status)` and `csr_tril(a,offset,status)` | Extract upper and lower triangular CSR matrices at the specified diagonal |

These function names are idiomatic native Fortran entry points rather
than Python class methods; `indices`, `rows`, and `columns` use **one-based**
indexing (SciPy uses zero-based indices). Integer statuses are returned
instead of Python exceptions. All operations are restricted to `real(dp)`;
complex numbers, multiple numeric precision types, and broadcasting
are not yet supported directly; sparse solvers and factorizations are supplied
separately by `scifort_sparse_linalg`.

Addition, Hadamard multiplication and matrix multiplication omit exact
zero-valued outputs automatically. Constructors and transpose preserve
explicit stored zeros; `csr_scale` also preserves them.

## CSC and diagonal operations example

```fortran
use scifort_kinds, only : dp
use scifort_sparse, only : csr_matrix, csc_matrix, csr_diags, &
    csr_diagonal, csr_setdiag, csr_to_csc, csc_matvec
implicit none
type(csr_matrix) :: a
type(csc_matrix) :: c
real(dp), allocatable :: d(:), y(:)

! Construct 3x4 with a main diagonal (3 entries) and an upper +2
! diagonal (2 entries). Packed slices are [ptr(1):ptr(2)-1], etc.
a = csr_diags(3,4,[0,2],[1,4,6], &
    [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp])
d = csr_diagonal(a)             ! [1,2,3]
call csr_setdiag(a,[9.0_dp,8.0_dp,7.0_dp])
c = csr_to_csc(a)
y = csc_matvec(c,[1.0_dp,1.0_dp,1.0_dp,1.0_dp])
```

`csr_diags` uses a **packed ragged-diagonal** API rather than forcing
unequal-length diagonals into a rectangular array. The one-based `diag_ptr`
has length `size(offsets)+1`, starts at 1, ends at `size(values)+1`, and
its adjacent slices specify each diagonal's values in **increasing row**
order. Each diagonal must contain exactly the number of positions inside
`(n_rows,n_cols)`; negative offsets denote lower diagonals. Duplicate
offsets, incomplete diagonals and nonfinite values are rejected. Stored
zeros remain explicit after construction or `csr_setdiag`.

`csr_setdiag` accepts the **complete** target diagonal, unlike SciPy's more
flexible scalar/broadcasting and partial update variants. Invalid input
leaves the original matrix unchanged. `csr_triu` and `csr_tril` preserve
explicit stored zeros for the entries that survive the triangular mask.

These are native Fortran APIs, not identical Python call signatures.
Potentially huge dense conversions and dimension extremes remain bounded
by ordinary Fortran integer and address-space limits.

## Storage and complexity

Let `z` be the number of stored nonzeros and `m`, `k`, `n` the matrix sizes.

| Operation | Complexity / storage |
|---|---|
| COO to CSR/CSC | O(z log z) worst-case sorting, O(z+m+n) storage |
| CSR to CSC and back | O(z+m+n) work, O(z+m+n) storage |
| CSR/CSC diagonal extraction and triangular masks | O(z+m) work (scan), O(max(m,n)) diagonal vector |
| Dense to CSR | O(mn) work, O(z+m) CSR storage |
| Transpose | O(z+m+n) work, O(z+n) output |
| Matvec | O(z+m) work, O(m) output |
| Addition | O(z_A+z_B+m), at most z_A+z_B output entries |
| Elementwise product | O(z_A+z_B+m), at most min(z_A,z_B) entries |
| Matrix multiplication | O(sum over A(i,k) of nnz(B(k,:)) + row-wise sorting), O(n+m+nnz(C)) auxiliary/output storage |

`csr_matmul` uses a dense **one-row** accumulator of length `n` (not a
`m` by `n` matrix). Its two-pass row assembly produces canonical, sorted
column indices with no dense matrix intermediate. It is suitable for large
sparse problems when `n` is moderate, but an extremely wide matrix can
still make the length-`n` workspace expensive.

## Testing

Native regression, including a 12,000-by-12,000 sparse matrix product:

```text
python tools/sparse_reference/run_native_sparse.py --mode both
```

Randomized, independent comparisons against SciPy/NumPy (SciPy 1.17 was
used during development):

```text
python tools/sparse_reference/compare_sparse_scipy.py --mode both --cases 430
```

Checked mode uses GNU Fortran `-O0 -g -fcheck=all -fbacktrace` and optimized
mode `-O3`. These standalone tests require GNU Fortran and (for the reference
comparison) Python, NumPy, and SciPy; none is a runtime dependency of the
Fortran module.

CSC and diagonal regression (including a 25,000-row sparse conversion),
checked and optimized builds:

```text
python tools/sparse_reference/run_native_sparse.py --mode both
python tools/sparse_reference/compare_sparse_csc_scipy.py --mode both --cases 350
```

The latter checks CSC conversion, round trips, transpose, packed diagonals,
triangular masks, diagonal extraction and replacement against SciPy, including
random rectangular and empty shapes, duplicate COO entries, and stored zeros.

See [sparse code reuse and attribution](sparse_reuse_repositories.md) for an
upstream-license survey relevant to future `scipy.sparse.linalg` work.
