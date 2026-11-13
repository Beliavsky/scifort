# Lebedev quadrature on the unit sphere

`lebedev_rule` supplies Cartesian points and surface-integration weights for
all 32 rules supported by SciPy. It is exported by `scifort_integrate` and
the implementation module `scifort_integrate_lebedev`.

```fortran
call lebedev_rule(n, x, w, status)
```

`n` is the algebraic order, not the number of points. The pure subroutine
allocates binary64 `x(3,m)` and `w(m)`, where `m` is the rule's node count.
Each column of `x` is a unit-sphere point. The weights sum to `4*pi`, up to
roundoff, so `sum(w*f(x))` approximates the surface integral directly.
For a surface average, divide the weighted sum by `4*pi`.

## Example

This complete program integrates `exp(x_coordinate)`. Its exact surface
integral is `4*pi*sinh(1)`.

```fortran
program lebedev_example
    use scifort_integrate, only : dp, lebedev_rule, integrate_status_success
    implicit none
    real(dp), allocatable :: x(:,:), w(:)
    real(dp) :: integral, reference
    integer :: status

    call lebedev_rule(17, x, w, status)
    if (status /= integrate_status_success) error stop 1
    integral = sum(w*exp(x(1,:)))
    reference = 4.0_dp*acos(-1.0_dp)*sinh(1.0_dp)
    print '(a,es24.16)', 'Surface integral: ', integral
    print '(a,es12.4)', 'Absolute error:   ', abs(integral-reference)
end program lebedev_example
```

Complex integrands use the same real points and weights with a complex
weighted sum. Vector-valued integrands can be summed component by component.
For a sphere with center `c(3)` and positive radius `r`, evaluate the integrand
at `c+r*x(:,j)` and multiply the weighted sum by `r**2`.

## Available rules

| Algebraic order `n` | Nodes `m` |
| --- | --- |
| 3 | 6 |
| 5 | 14 |
| 7 | 26 |
| 9 | 38 |
| 11 | 50 |
| 13 | 74 |
| 15 | 86 |
| 17 | 110 |
| 19 | 146 |
| 21 | 170 |
| 23 | 194 |
| 25 | 230 |
| 27 | 266 |
| 29 | 302 |
| 31 | 350 |
| 35 | 434 |
| 41 | 590 |
| 47 | 770 |
| 53 | 974 |
| 59 | 1202 |
| 65 | 1454 |
| 71 | 1730 |
| 77 | 2030 |
| 83 | 2354 |
| 89 | 2702 |
| 95 | 3074 |
| 101 | 3470 |
| 107 | 3890 |
| 113 | 4334 |
| 119 | 4802 |
| 125 | 5294 |
| 131 | 5810 |

Only these orders are accepted; for example, 33 is unsupported. No rounding
to a neighboring rule occurs. The same output variables can be reused on
successive calls; `intent(out)` replaces their previous allocations.

## Status and numerical contract

`status` is optional and uses existing integration constants:

| Status | Meaning |
| --- | --- |
| `integrate_status_success` (0) | Rule successfully allocated and generated. |
| `integrate_status_invalid_argument` (2) | Unsupported integer order. |
| `integrate_status_numerical_failure` (3) | Output allocation failed. |

Both output arrays are unallocated on failure, including when previously
allocated arrays were supplied. With no status argument, check allocation
before accessing the outputs. The library performs no printing or input and
does not stop execution. Immutable tables require no runtime data files.

The rules integrate Cartesian polynomials of total degree at most `n` exactly
in exact arithmetic. Binary64 coefficient and arithmetic errors remain.
Octahedral symmetry includes reflections and permutations of coordinates.
The implementation preserves the published coefficients and their original
point ordering; the current ordering agrees with SciPy 1.18.1.

**Orders 13, 25, and 27 contain negative weights**, as in the source rules.
All other supported rules have positive weights. Negative weights are
intentional and may cause cancellation. The routine does not clamp or
renormalize weights. Accuracy for general, nonsmooth, or sharply localized
functions is not guaranteed by the polynomial order; compare orders or use
problem-specific error analysis. No adaptive integrator or error estimate is
provided by this node/weight routine.

## Sources and validation

The tables and orbit expansion adapt SciPy v1.18.1's BSD-2-Clause
[`_lebedev.py`](https://github.com/scipy/scipy/blob/v1.18.1/scipy/integrate/_lebedev.py),
Copyright (c) 2010, Robert Parrish. Parrish's implementation credits
Christoph van Wuellen's Fortran translation of Dmitri Laikov's original C
routines. The requested scientific reference is:

V. I. Lebedev and D. N. Laikov, "A quadrature formula for the sphere of the
131st algebraic order of accuracy", *Doklady Mathematics* 59(3), 477-481 (1999).

See [provenance](../CODE_PROVENANCE.md), the
[full upstream license](../LICENSES/scipy-lebedev-BSD-2-Clause.txt), and the
[validation record](../validation/lebedev_2026-11-13.md). Every rule is checked
against SciPy and independent analytic moments, including symmetry and signed
weights. Native Fortran coverage is the current scope; C ABI and wrappers
remain deferred.
