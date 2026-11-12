# Polynomial companion and symmetric Fiedler matrices

SciFort adds three `scipy.linalg` matrix constructors using original,
MIT-licensed Fortran source. No LAPACK, FFT, Python, or other library is
needed by these two modules.

| Fortran module | Public interface | SciPy analog |
|---|---|---|
| `scifort_linalg_companion` | `companion(a,c,status)` | `scipy.linalg.companion` |
| `scifort_linalg_companion` | `fiedler_companion(a,c,status)` | `scipy.linalg.fiedler_companion` |
| `scifort_linalg_fiedler` | `fiedler(a,f,status)` | `scipy.linalg.fiedler` |

All three are **subroutines**, with caller-allocated output and status codes:
`linalg_status_success` (0) or `linalg_status_invalid_shape` (1). Every
output is initialized to zero before shape/leading-coefficient validation;
failure leaves a zero-filled output. The routines are `pure` and support
one leading batch dimension.

## Example

```fortran
use scifort_kinds, only : dp
use scifort_linalg_companion, only : companion, fiedler_companion
use scifort_linalg_fiedler, only : fiedler
real(dp) :: poly(5), c(4,4), fc(4,4), d(3,3)
integer :: status

poly = [1.0_dp,-16.0_dp,86.0_dp,-176.0_dp,105.0_dp]
call companion(poly,c,status)
call fiedler_companion(poly,fc,status)
! Both c and fc have eigenvalues 1, 3, 5, and 7.
call fiedler([1.0_dp,4.0_dp,12.0_dp],d,status)
! d(1,:) = [0,3,11].
```

See the [runnable example](../example/linalg_companion_example.f90).

## Shapes and types

For input `a(:)` of length `n`, `companion` takes a matrix `c(n-1,n-1)`;
`n>=2` and `a(1)/=0` are required. The first row is `-a(2:n)/a(1)`;
entries immediately below the main diagonal are 1, and all other entries
are zero. Coefficients are in **descending polynomial power order**.

`fiedler_companion` has the same shape for `n>=2`, but is pentadiagonal.
For `n=0` or `n=1`, it accepts an empty `c(0,0)` and reports success.
For degrees one and two, it coincides with the ordinary companion matrix.
The coefficients are scaled by `a(1)`; for very large coefficient ratios,
both companion constructions can lose floating-point accuracy. Monic
coefficients are preferable. Nonzero leading coefficients are mandatory
for positive-degree polynomials.

For `fiedler`, the output `f(n,n)` satisfies
`f(i,j)=abs(a(i)-a(j))`, with `f(i,i)=0`; n=0 is allowed. Complex input
returns real, nonnegative distances.

There are overloads for default `integer`, `real(dp)` and `complex(dp)`.
Integer inputs produce a `real(dp)` companion and Fiedler companion,
but an integer ordinary Fiedler distance matrix. Real inputs produce real
matrices and complex polynomial inputs produce complex companion matrices.
For complex Fiedler input, **the output is real(dp)**.

The batch overloads accept `a(batch,n)` and
`c(batch,n-1,n-1)` / `f(batch,n,n)`; they do not yet support an arbitrary
number of leading dimensions. If any batch polynomial has a zero leading
coefficient, the whole call returns invalid-shape status and zero output.
For `batch=0`, the regular dimension checks still apply.

The constructor does not attempt to compute roots. Use an eigensolver (or
polynomial root solver) on the returned companion matrix when roots are
needed. Both companion matrices have the same eigenvalues but generally
**different entries and eigenvectors**.

## Differences from SciPy

- Python returns arrays directly and raises `ValueError`; SciFort uses a
  preallocated output and explicit status.
- SciPy supports arbitrary batch ranks; SciFort supports one leading batch
  dimension only.
- `fiedler_companion` with `n<2` is represented by a `(0,0)` Fortran array,
  while SciPy returns an empty one-dimensional ndarray for the unbatched case.
- SciPy accepts additional NumPy dtypes; SciFort supports default integer,
  `real(dp)`, and `complex(dp)`. Integer distance arithmetic must remain
  representable in the default Fortran integer type.
- No preemptive NaN/Inf rejection is performed; inputs with nonfinite
  coefficients follow ordinary IEEE arithmetic and may yield nonfinite
  output (as with SciPy's constructors). No normalization is applied beyond
  division by the leading polynomial coefficient.

See the [validation report](linalg_companion_validation_2026-11-09.md).
The mathematical API was checked with the public SciPy documentation:
[companion](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.companion.html),
[fiedler_companion](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.fiedler_companion.html), and
[fiedler](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.fiedler.html).
