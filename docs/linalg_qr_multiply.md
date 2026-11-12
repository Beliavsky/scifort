# Householder QR multiplication (`scipy.linalg.qr_multiply`)

`scifort_linalg_qr_multiply` adds `qr_multiply` to compute the economic QR
factorization of a matrix and multiply its Q factor by one or more vectors,
**without explicitly constructing Q**. It also adds `qr_apply`, a SciFort
extension for reusing compact QR factors when multiplying multiple inputs.

Both functions use the existing self-contained `scifort_linalg_qr` Householder
implementation. There is no new LAPACK, BLAS, or FFT dependency.

## API

```fortran
use scifort_kinds, only : dp
use scifort_linalg_qr, only : qr_factor
use scifort_linalg_qr_multiply, only : qr_multiply, qr_apply

! Given A(m,n), let k = min(m,n).
! Left multiplication: C(k,p) -> result(m,p), R(k,n).
call qr_multiply(a, c, result, r, status, mode='left')

! Right multiplication: C(p,m) -> result(p,k), R(k,n).
call qr_multiply(a, c, result, r, status, mode='right')

! Omit mode for the SciPy-compatible 'right' default.
! Passing pivots(n) activates column-pivoted QR.
call qr_multiply(a, c, result, r, status, mode='left', pivots=pivots)

! Reuse precomputed Householder factors (SciFort extension):
call qr_factor(a, af, tau, status)
if (status == 0) call qr_apply(af, tau, c, result, status, mode='left')
```

Overloads are provided for `real(dp)` and `complex(dp)` coefficient matrices
and multipliers of rank one (vectors) or rank two (matrices). Mixed
real/complex arguments are supported by `qr_multiply`; the resulting
`product` and `r` must both be complex. Reusable `qr_apply` accepts real/real
or complex/complex factors and multipliers; convert mixed types explicitly.
Outputs are caller-allocated. All dummy arguments use Fortran semantics:
input arrays are unmodified; there are no `overwrite_a`/`overwrite_c`
options. The API does not yet accept leading batch dimensions.

## Mathematical conventions and array shapes

The economic Q has shape `(m,k)`, and R has shape `(k,n)`, for positive
`m` and `n`. Omitting `pivots` yields `A = Q R`; requesting `pivots(n)`
enables column pivoting, so `A(:,pivots) = Q R`. Permutations use **one-based
Fortran indices**, rather than SciPy's zero-based indices.

| Mode | Multiplier matrix | Output matrix | Vector input and output |
|---|---|---|---|
| `'left'` | `C(k,p)` | `Q C` with shape `(m,p)` | `c(k)` -> `result(m)` |
| `'right'` (default) | `C(p,m)` | `C Q` with shape `(p,k)` | `c(m)` -> `result(k)` |

`conjugate=.true.` applies the **elementwise complex conjugate** of Q (not
its conjugate transpose); it has no effect for a real Q. This matches the
meaning of SciPy's `conjugate` argument. The default is `.false.`.

`qr_apply(af,tau,c,product,status,mode,conjugate)` accepts exactly the
compact Householder factors written by `scifort_linalg_qr::qr_factor`,
including those produced with pivoting. It does not take a pivot vector,
since the reflectors already define the Q factor. These routines support
any number of multiplier columns/rows, including **zero right-hand sides**.

## Algorithm and limitations

For `A(m,n)` and `k=min(m,n)`, `qr_factor` computes the compact
Householder reflectors `H(1),...,H(k)`, with `Q = H(1)...H(k)` restricted
to the first `k` columns. A left multiplication starts with `[C; 0]` and
applies the reflectors in reverse order. A right multiplication applies
the reflectors to C in forward order, then keeps its first `k` columns.
Both operations use the original compact factors and avoid allocating Q.

The multiplication requires `O(m*k*p)` arithmetic for `p` multiplier
vectors, plus a temporary of shape `(m,p)` or `(p,m)`. Fresh `qr_multiply`
also incurs the usual `O(m*n*k)` QR factorization work and stores `af(m,n)`.
`qr_apply` avoids that factorization overhead when factors are reused.

The Householder sign/phase conventions may differ from LAPACK/SciPy while
producing equally valid decompositions. Compare reconstructions, Q
orthogonality, and products after aligning column phases, not raw Q or R
entries. Exact rank-deficient inputs may have less uniquely determined Q
columns. Invalid matrix sizes, output sizes, modes, or permutation sizes
produce `linalg_status_invalid_shape`; nonfinite factors or multipliers
produce `linalg_status_no_convergence`. All outputs are initialized to zero
when a call fails.

**SciPy 1.17.0 compatibility note:** Its `qr_multiply` implementation can
silently discard the imaginary component of a complex multiplier when the
matrix A is real. SciFort instead promotes this combination to complex and
preserves the mathematically correct complex product; the independent
comparison script uses SciPy `qr` and NumPy multiplication as its oracle for
this particular combination. This is an intentional correctness difference,
not a claim of byte-for-byte agreement with that SciPy release.

## Reproducible tests

From the project root, with GNU Fortran and SciPy installed:

```sh
mkdir -p build/qr_multiply_check
cd build/qr_multiply_check
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -c \
  ../../src/kinds.f90 ../../src/linalg.f90 ../../src/linalg_svd.f90 \
  ../../src/linalg_qr.f90 ../../src/linalg_qr_multiply.f90
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -o test_qr_multiply \
  *.o ../../test/test_linalg_qr_multiply.f90
./test_qr_multiply
gfortran -std=f2018 -O0 -g -fcheck=all -fbacktrace -o qr_multiply_driver \
  *.o ../../reference/linalg_qr_multiply_driver.f90
python ../../reference/check_linalg_qr_multiply.py ./qr_multiply_driver
```

The independent reference compares both left and right products against
SciPy QR with appropriate column phase alignment, checks R reconstruction,
and compares direct SciPy `qr_multiply` whenever its result is well-defined.
See `docs/linalg_qr_multiply_validation_2026-11-09.md` for results.
