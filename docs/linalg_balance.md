# Matrix balancing (`scipy.linalg.matrix_balance`)

The independent `scifort_linalg_balance` module implements a real/complex
binary64 similarity balancing transform, following the classical principles
of Parlett--Reinsch balancing. It computes

    B = inv(T) A T    (equivalently A T = T B)

with `T = P D`, where `P` is a permutation matrix and `D` is a positive
diagonal matrix whose entries are exact integral powers of two. The ordering
can isolate uncoupled rows/columns. Subsequent binary-power scaling tries
to decrease the combined off-diagonal row/column 1-norm within the remaining
coupled block. Eigenvalues, determinant, and trace are preserved by the
similarity in exact arithmetic, although floating-point eigensolvers may
produce different rounded eigenvalues.

## Public Fortran interface

```fortran
use scifort_kinds, only : dp
use scifort_linalg_balance, only : matrix_balance
real(dp) :: a(4,4), b(4,4), t(4,4), factors(4)
integer :: status, indices(4)
! Assign finite square input a.
call matrix_balance(a, b, status, t, &
    permute=.true., scale=.true., &
    scale_factors=factors, permutation=indices)
```

`matrix_balance` is a generic over `real(dp)` and `complex(dp)`. Inputs and
matrix outputs have the same type. The `b` and `t` arrays must have shape
`(n,n)`; `scale_factors` and `permutation` must have length `n`. Input and
output arrays are supplied by the caller and are not reallocated. Both
`permute` and `scale` default to `.true.`. Each optional output may be omitted.
The empty `(0,0)` matrix is accepted.

| Argument | Intent | Meaning |
| --- | --- | --- |
| `a` | `in` | Finite, square real or complex input `(n,n)` |
| `b` | `out` | Balanced matrix `(n,n)`, with `B = inv(T) A T` |
| `status` | `out` | `linalg_status_*` code |
| `t` | optional `out` | Full similarity matrix `(n,n)` |
| `permute` | optional `in` | Isolate uncoupled rows/columns; default true |
| `scale` | optional `in` | Apply power-of-two scaling; default true |
| `scale_factors` | optional `out` | Positive scale entries `D(j,j)` |
| `permutation` | optional `out` | **One-based** original index of column `j` after reordering |

The two separate outputs reconstruct the full similarity matrix as
`T(permutation(j),j) = scale_factors(j)`, with all other entries zero.
`B(i,j) = A(permutation(i),permutation(j)) * scale_factors(j) / scale_factors(i)`.
SciPy's `separate=True` permutation is **zero-based**, so add one when
comparing its indices with SciFort's. When both flags are false, `B=A`,
`T=I`, and the scale/permutation outputs are identity. The input is never
overwritten. SciPy's `overwrite_a`, `check_finite=False`, and batch inputs
are not supported.

## Status, algorithm, and limits

- `linalg_status_success` (0): a finite result and similarity transform
  were produced, including the empty array case.
- `linalg_status_invalid_shape` (1): an input or output shape is wrong.
- `linalg_status_no_convergence` (4): nonfinite input, nonconvergent
  balancing iterations, or numerically unrepresentable output.

On invalid shape or nonfinite input, outputs are initialized to zero.
The core first isolates zero-coupled eigenvalues through symmetric row/
column swaps. It then normalizes the working matrix and performs up to 100
sweeps of off-diagonal balancing, accepting a scale only when the predicted
row/column norm decreases by at least 5%. Scale changes are restricted to
integer powers of two. Numeric guards skip scaling moves that could make
scale factors or matrix entries too large or too small. This can yield a
valid similarity transform without achieving perfect row/column balance.

The real overload reuses the complex-valued core and therefore allocates
complex working storage. This is an independent, unblocked `O(n^2)` storage
and approximately `O(n^2)` per sweep algorithm, not an exact reimplementation
of LAPACK's `xGEBAL`. Both the permutation order and the accepted binary
scale factors can differ from SciPy. Particularly ill-conditioned or
extremely heterogeneously scaled matrices may suffer underflow in the
working normalization; balancing is heuristic and is **not** guaranteed
to improve eigenvalue conditioning. It may be harmful on some matrices.

## Testing

The native tests cover real/complex similarity identities, optional output
modes, invalid shapes and nonfinite inputs, empty matrices, and extreme
uniform scales. The separate SciPy comparison script checks independent
`scipy.linalg.matrix_balance` results for a wide variety of matrices,
validates the full similarity identity and the power-of-two scales, and
reports relative off-diagonal 1-norm quality without asserting identical
factor arrays.

```bash
python tools/linalg_reference/run_balance_checks.py \
    --mode both --cases 650 --max-order 96 --seed 20261108
```

See `example/linalg_balance_example.f90` and
`docs/linalg_balance_validation_2026-11-07.md`.

Reference: B. N. Parlett and C. Reinsch, "Balancing a Matrix for
Calculation of Eigenvalues and Eigenvectors", *Numerische Mathematik* **13**
(1969), 293--304, DOI: 10.1007/BF02165404. The public SciPy contract is
used for cross-library consistency and as an independent numerical reference.
No LAPACK, SciPy, or other third-party implementation source was incorporated.
