# Dense matrix trigonometric and hyperbolic functions

The `scifort_linalg_trigonometric` module provides six SciPy-style matrix
functions, each accepting both real and complex binary64 input matrices:

| Procedure | Math | SciPy correspondence |
| --- | --- | --- |
| `sinm` | sin(A) | `scipy.linalg.sinm` |
| `cosm` | cos(A) | `scipy.linalg.cosm` |
| `tanm` | sin(A) cos(A)^(-1) | `scipy.linalg.tanm` |
| `sinhm` | sinh(A) | `scipy.linalg.sinhm` |
| `coshm` | cosh(A) | `scipy.linalg.coshm` |
| `tanhm` | sinh(A) cosh(A)^(-1) | `scipy.linalg.tanhm` |

These are **matrix functions**, not elementwise array functions. As with the
underlying matrix exponential, a matrix need not be diagonalizable. The
module is self-contained beyond existing MIT-licensed SciFort modules and
adds no required dependencies.

## Fortran API

```fortran
use scifort_kinds, only : dp
use scifort_linalg_trigonometric, only : sinm, cosm, tanm, sinhm, coshm, tanhm
real(dp) :: a(3,3), x(3,3)
complex(dp) :: b(3,3), y(3,3)
integer :: status

call sinm(a, x, status)
call tanm(b, y, status)
call coshm(a, x, status)
```

Each routine has the signature `call function_name(a, x, status)`. The
input/output element type must be the same (`real(dp)` or `complex(dp)`).
Both arrays must be square and of matching shape, with order **1 to 40**.
The caller owns and allocates output arrays. Status values come from
`scifort_linalg`:

- `linalg_status_success` (0): computation succeeded;
- `linalg_status_invalid_shape` (1): invalid or mismatched dimensions;
- `linalg_status_no_convergence` (4): nonfinite data, excessive magnitude,
  numerical overflow, or convergence/reconstruction failure;
- `linalg_status_singular` (5): matrix cosine/cosh is singular or numerically
  too close to singular for a reliable tangent solve.

On any failure, the function output is zero. The routines do not call
`stop` or `error stop` and do not allocate hidden persistent state.

## Numerical implementation

An independent matrix Taylor-series implementation computes sine and cosine
**together**, or their hyperbolic counterparts. The 1-norm of A is reduced
to at most 0.5 using a power-of-two scaling, then both series are evaluated
using powers of A^2. Repeated double-angle identities recover the original
argument:

- `sin(2A) = 2 sin(A) cos(A)`;
- `cos(2A) = cos(A)^2 - sin(A)^2`;
- `sinh(2A) = 2 sinh(A) cosh(A)`;
- `cosh(2A) = cosh(A)^2 + sinh(A)^2`.

This avoids the relative loss of accuracy in subtracting `exp(A)-exp(-A)`
(or `exp(iA)-exp(-iA)`) when A is tiny. Tangent functions use SciFort's
existing pivoted small-dense inverse and verify the resulting right solve.
Fortran complex arithmetic handles analytic continuation automatically.

## Limitations

- This is a small dense implementation without LAPACK acceleration or
  Schur/Parlett evaluation. It is intended for matrices of order <= 40;
  behavior for highly nonnormal or very large-norm matrices is not guaranteed.
- The scaling iteration limit is 60 and the initial Taylor limit is 40 terms.
  Inputs large enough to overflow intermediate results return a failure
  status rather than IEEE infinities as a successful result.
- Matrices whose cosine/cosh is at or near a pole of tangent/tanh return a
  singular status rather than an arbitrarily huge, poorly conditioned value.
- Matrix function evaluation can be ill-conditioned even when all operations
  are finite. For critical applications, compare against a well-tested
  Schur/Padé implementation.

## Tests

Requires GNU Fortran, NumPy and SciPy for the randomized reference tests:

```sh
python tools/linalg_reference/run_trigonometric_checks.py --mode both --cases 1200
```

The script builds and runs the standalone native regression suite and compares
randomized real/complex matrices against SciPy 1.17.0. All six operations are
included. It uses temporary compilation directories and makes no changes to
the project's FPM dependencies. For an executable demonstration, build
`example/linalg_trigonometric_example.f90` with the listed module sources.

See `docs/linalg_trigonometric_validation_2026-11-07.md` for the recorded
validation commands and results.
