# Helmert and Hilbert matrices (2026-11-09)

SciFort implements three matrix constructors from `scipy.linalg` using two new
MIT-licensed Fortran modules. The routines are stand-alone apart from SciFort's
existing kind and status declarations, and have **no additional dependencies**.

| SciPy function | Fortran module and call |
|---|---|
| `scipy.linalg.helmert(n, full=False)` | `scifort_linalg_helmert`: `call helmert(n, h, status, full=...)` |
| `scipy.linalg.hilbert(n)` | `scifort_linalg_hilbert`: `call hilbert(n, h, status)` |
| `scipy.linalg.invhilbert(n, exact=False)` | `scifort_linalg_hilbert`: `call invhilbert(n, h, status)` |
| `scipy.linalg.invhilbert(n, exact=True)` | `invhilbert(n, h_int64, status)` for `n<=14`; `invhilbert_exact_decimal(n, h_string, status)` at any feasible order |

All matrix outputs are **caller-allocated** and status is an `integer`. The
routines initialize the output before validating shapes. A successful call
returns `linalg_status_success` (0), a bad shape/order or unsupported exact
`int64` order returns `linalg_status_invalid_shape` (1), and overflow in
floating-point `invhilbert` returns `linalg_status_no_convergence` (4), leaving
the output zero-filled. Exact decimal output that has insufficient character
length similarly returns status 1 and is blank-filled.

## Helmert

`helmert` produces an orthonormal matrix of Helmert contrasts. For an order
`n>=1`, `full=.true.` requires `h(n,n)`; the default `full=.false.` requires
`h(n-1,n)` and omits the constant first row. The first full row equals
`1/sqrt(n)` in all columns. For one-based `k=2..n`, row `k` has its first
`k-1` entries equal to `1/sqrt(k*(k-1))`, its `k`-th entry equal to
`-(k-1)/sqrt(k*(k-1))`, and all remaining entries zero.

`n=1` is valid: compact output has shape `(0,1)`, while the full output is
`[1]`. We reject `n<=0` with a status rather than reproducing NumPy's
exception for SciPy's zero-order Helmert constructor.

## Hilbert

`hilbert(n, h, status)` accepts `n>=0`, with output shape `(n,n)`. Each
one-based entry is `h(i,j)=1/(i+j-1)` as `real(dp)`; `n=0` yields `(0,0)`.

`invhilbert` **does not** numerically invert a nearly singular Hilbert
matrix. It uses the closed-form combinatorial identity

\[
  (H_n^{-1})_{i,j}=(-1)^{i+j}(i+j-1)
  {n+i-1\choose n-j}{n+j-1\choose n-i}
  {i+j-2\choose i-1}^{2}.
\]

This generates the inverse's enormous integer coefficients without solving
an ill-conditioned dense linear system. The floating-point version multiplies
binomial factors with overflow checks, returning `real(dp)` values. For high
`n`, simply storing the inverse in double precision may overflow; status 4
is then returned, without partial output.

### Exact integer options

- For `0<=n<=14`, supply an `integer(int64)` output matrix to the **same**
  `invhilbert` generic interface. SciPy's exact integer entries fit signed
  64 bits in this range. For `n=15` or higher, this overload returns status
  1, avoiding silent overflow.
- At any practical size, call `invhilbert_exact_decimal(n, text, status)`
  with a `character(len=...) :: text(n,n)` array. The implementation constructs
  **exact** integers through prime factorization and base-10^9 multiplication.
  Each output string is the integer's decimal representation, including a
  minus sign if applicable. It does **not** require compiler-specific
  128-bit integers or a third-party big-integer package. Ensure sufficient
  fixed character length for every entry; otherwise the operation fails
  atomically. Limited by memory, integer loop-index range and runtime.

The string overload is a *SciFort extension*: Python's `exact=True`
produces NumPy `int64` for `n<=14`, and arbitrary-precision Python objects
for `n>14`; Fortran has no corresponding standard integer object dtype.

## Fortran example

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_kinds, only : dp
use scifort_linalg_helmert, only : helmert
use scifort_linalg_hilbert, only : hilbert, invhilbert, invhilbert_exact_decimal
integer :: status
real(dp) :: contrasts(3,4), h4(4,4), h4_inv(4,4)
integer(int64) :: h4_exact(4,4)
character(80) :: h16_decimal(16,16)

call helmert(4, contrasts, status)     ! 3x4; choose full=.true. for 4x4
call hilbert(4, h4, status)
call invhilbert(4, h4_inv, status)     ! floating real(dp)
call invhilbert(4, h4_exact, status)   ! exact integer(int64)
call invhilbert_exact_decimal(16, h16_decimal, status)
! trim(h16_decimal(8,8)) == '42475099528537378560'
```

See the [complete runnable example](../example/linalg_special_matrices_example.f90)
and [numerical validation report](linalg_special_matrices_validation_2026-11-09.md).

## Compatibility details

- Fortran arrays use one-based indices. Python returns newly allocated
  arrays; Fortran takes a preallocated output matrix and an explicit status.
- All constructors are real-valued (`real(dp)`), just like SciPy's
  floating-point Helmert and Hilbert matrices. The integer exact path
  explicitly uses `iso_fortran_env:int64`.
- Exact arbitrary-precision entries use **strings**, not directly arithmetical
  Fortran integers; callers can pass the strings to an arbitrary-precision
  package of their choice when needed.
- Input validation rejects negative `n` uniformly; SciPy's `hilbert(-1)`
  currently returns an empty array, but negative orders are not meaningful.
- Matrix dimensions scale as `O(n²)`, while constructing every exact decimal
  inverse element uses prime factorizations and is intended for moderate `n`.
  Avoid directly multiplying ill-conditioned Hilbert matrices at large `n` to
  validate their inverses: the residual is dominated by double-precision
  round-off even when the constructed inverse coefficients are accurate.

SciPy public APIs:
[helmert](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.helmert.html),
[hilbert](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.hilbert.html), and
[invhilbert](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.invhilbert.html).
