# Pascal and inverse Pascal matrices (2026-11-09)

The `scifort_linalg_pascal` module implements `scipy.linalg.pascal` and
`scipy.linalg.invpascal` for symmetric, lower-triangular and upper-triangular
matrices. It builds exact integer matrices without performing unstable dense
inversions. No external dependency is added.

## API

| Fortran call | Meaning |
|---|---|
| `call pascal(n, p_real, status, kind='symmetric')` | Pascal coefficients as `real(dp)` |
| `call pascal(n, p_int64, status, kind='lower')` | Exact Pascal coefficients as signed `integer(int64)` |
| `call invpascal(n, p_real, status, kind='upper')` | Inverse Pascal coefficients as `real(dp)` |
| `call invpascal(n, p_int64, status, kind='symmetric')` | Exact inverse coefficients as signed `integer(int64)` |
| `call pascal_exact_decimal(n, p_text, status, kind='symmetric')` | Exact Pascal coefficients as decimal strings |
| `call invpascal_exact_decimal(n, p_text, status, kind='symmetric')` | Exact inverse coefficients as decimal strings |

The arguments are:

- `n`: nonnegative order; `n=0` produces an empty `(0,0)` matrix.
- `p_real`, `p_int64`, or `p_text`: caller-allocated rank-two output of shape
  `(n,n)`, storing real64, signed int64, or fixed-length character strings.
- `status`: integer output, `linalg_status_success` (0) on success,
  `linalg_status_invalid_shape` (1) for bad order/kind/dimensions or too-short
  decimal strings, `linalg_status_no_convergence` (4) if the requested
  numeric representation overflows. An unsuccessful output is cleared.
- `kind`: optional, exactly `symmetric` (default), `lower`, or `upper`.

### Mathematical definition

With one-based indices and binomial coefficients `C(a,b)`, the three Pascal
matrices are

- symmetric: `P(i,j) = C(i+j-2,i-1)`;
- lower: `L(i,j) = C(i-1,j-1)` for `i >= j`, otherwise zero;
- upper: `U = transpose(L)`.

The symmetric Pascal matrix factors as `P = L * transpose(L)`.
The lower/upper inverses are signed Pascal matrices with entries
`L_inv(i,j) = (-1)^(i-j) C(i-1,j-1)` for `i>=j`. Therefore,

```
P_inv(i,j) = (-1)^(i+j) sum(k=max(i,j),n) C(k-1,i-1)*C(k-1,j-1)
```

This formula is evaluated using **positive products and sums** and a final
sign, avoiding the severe cancellation encountered with generic inversion.
All numeric constructors reject overflow instead of silently wrapping or
producing infinities. For the symmetric inverse, the simple summation is
`O(n^3)`; the other constructors are `O(n^2)`.

### Representations and limits

`integer(int64)` can represent every symmetric Pascal entry through order
34, and every symmetric inverse Pascal entry through order 35. For lower or
upper Pascal matrices **and their inverses**, `int64` works through order 67.
Larger orders return status 4 when the first coefficient overflows.

For arbitrary orders where machine integers are insufficient, supply a
`character(len=...) :: p_text(n,n)` array to the `*_exact_decimal` routines.
They use an internal base-10^9 unsigned big-integer representation, adding
signs for inverse matrices. Every returned decimal string is exact. The
implementation checks character length, and clears the entire output if any
entry does not fit. These string-based entry points are **SciFort extensions**;
SciPy uses Python arbitrary-precision integer objects instead.

The floating overloads use recurrence and finite-range checks. Because of
floating-point rounding, entries may differ slightly from SciPy's
`exact=False` binomial implementation, but represent the same matrices.

## Example

```fortran
use, intrinsic :: iso_fortran_env, only : int64
use scifort_linalg_pascal, only : pascal, invpascal, pascal_exact_decimal
integer :: status
integer(int64) :: p(5,5), p_inv(5,5)
character(100) :: big(50,50)

call pascal(5,p,status)
call invpascal(5,p_inv,status)
call pascal_exact_decimal(50,big,status)
! trim(big(50,50)) == '25477612258980856902730428600'
```

See the [runnable example](../example/linalg_pascal_example.f90) and
[validation instructions](linalg_pascal_validation_2026-11-09.md).

## SciPy compatibility

SciPy documents the functions at
[pascal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.pascal.html)
and [invpascal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.invpascal.html).
SciFort uses one-based indices, caller-allocated arrays and an explicit error
status. The default `real(dp)` or `integer(int64)` output type is selected by
the caller, instead of SciPy's `exact=` Boolean; arbitrary-precision integer
results use decimal strings, not Python object arrays. Fortran's signed
`int64` has a smaller positive range than SciPy's unsigned `uint64` output
for exact Pascal matrices, and SciFort will reject values that do not fit.
