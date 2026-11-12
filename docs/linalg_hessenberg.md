# Upper Hessenberg reduction (`scipy.linalg.hessenberg`)

`scifort_linalg_hessenberg` adds an independent, MIT-licensed, real/complex
binary64 Hessenberg reduction based on Householder similarity transformations.
For a finite square matrix `a(n,n)`, it computes an upper Hessenberg matrix
`h(n,n)` and, optionally, an orthogonal or unitary matrix `q(n,n)` satisfying

    A = Q H Q^H

where `Q^H` means conjugate transpose (ordinary transpose for real matrices).
The elements `h(i,j)` for `i > j+1` are set to exact zero. Unlike QR,
Hessenberg reduction is a similarity transformation and preserves eigenvalues.
It is useful as a preliminary step before nonsymmetric eigensolvers and Schur
factorization.

## Fortran interface

```fortran
use scifort_kinds, only : dp
use scifort_linalg_hessenberg, only : hessenberg
real(dp) :: a(8,8), h(8,8), q(8,8)
integer :: status
! Set a.
call hessenberg(a, h, status, q)
! Optional q: when no explicit basis is needed, call hessenberg(a,h,status).
```

The generic accepts both real(dp) and complex(dp) matrices; input and output
kinds must match. The dummy arguments are, in positional order:

| Argument | Intent | Meaning |
| --- | --- | --- |
| `a` | `in` | Finite square input; `(0,0)` is valid |
| `h` | `out` | Caller-allocated Hessenberg result, same shape and type as `a` |
| `status` | `out` | A `linalg_status_*` code |
| `q` | optional `out` | Caller-allocated orthogonal/unitary similarity matrix, same shape and type |

This follows the mathematical contract of `scipy.linalg.hessenberg` but uses
SciFort's caller-allocated arrays and status argument rather than SciPy's
`calc_q` flag and return values. The input is never overwritten; there is
no `overwrite_a`, `check_finite=False`, or batched-array option.

### Status and edge cases

- `linalg_status_success` (0): calculation succeeded. For zero-size inputs,
  zero-size results are returned. For an all-zero matrix, `h=0`, `q=I`.
- `linalg_status_invalid_shape` (1): non-square `a` or incorrectly shaped
  `h`/`q`. Caller-owned outputs are initialized to zero on entry.
- `linalg_status_no_convergence` (4): nonfinite input or unrepresentable
  intermediate/output values. Unlike an iterative Schur solver, this routine
  has a fixed number of Householder steps; its numerical-failure status
  is not a QR iteration limit.

The implementation normalizes the input by its largest matrix component to
avoid overflowing Householder norms. Within each column, it scales the
Householder vector before forming the reflector. Complex reflectors use a
unit phase to avoid cancellation in their first component. Lower entries
are explicitly zeroed after each similarity step. The original scaling is
restored in the returned `h`; exceptionally large results may overflow even
when all input entries are finite, in which case a numerical-failure status
is returned. Values below binary64 representability may underflow.

This is an unblocked `O(n^3)` algorithm with `O(n^2)` working memory
(including the output); it is not expected to match optimized LAPACK on
large matrices. Eigenvectors, Schur factors, reordering and matrix balancing
are **not** computed by this routine. The decomposition is not unique: do
not compare `H` or `Q` element-by-element with SciPy's implementation.

## Validation

A standalone native test checks real and complex reconstruction, orthogonality,
structural zeros, optional-output behavior, empty arrays and error status.
Independent SciPy comparisons check `Q H Q^H`, `Q^H Q=I`, exact structural
zeros, trace and Frobenius-norm invariants, and reconstructability relative to
SciPy's `hessenberg(...,calc_q=True)`. The reference is used for tests only;
no SciPy or LAPACK source is included in SciFort.

```sh
python tools/linalg_reference/run_hessenberg_checks.py \
    --mode both --cases 480 --max-order 96 --seed 20261108
```

A small stand-alone example is `example/linalg_hessenberg_example.f90`.
See `docs/linalg_hessenberg_validation_2026-11-07.md` for test results and
additional test commands.
