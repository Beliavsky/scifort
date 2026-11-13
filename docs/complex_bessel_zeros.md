# Complex Bessel Y zeros

`y0_zeros`, `y1_zeros`, and `y1p_zeros` in `scifort_special` and
`scifort_bessel_zeros` accept either real or complex allocatable output arrays.
Their output types select the overload; no `complex=True` switch is needed.
The existing real routines and their argument ordering are unchanged.

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_special, only : y0_zeros
    implicit none
    complex(dp), allocatable :: roots(:), derivatives(:)
    integer :: status
    call y0_zeros(4,roots,derivatives,status)
    if (status == 0) print *, roots(1), derivatives(1)
    ! Root approximately -2.403016631993776 + 0.539882313042664i
    ! Y0' approximately -0.100747689288557 + 0.881967710099448i
end program example
```

All three overloads are pure subroutines with the same argument pattern:

```text
call y0_zeros(nt, zeros, derivatives[, status])
call y1_zeros(nt, zeros, derivatives[, status])
call y1p_zeros(nt, zeros, values[, status])
```

Both output arrays must be `complex(dp), allocatable`. The complex overloads
require `1 <= nt <= 1000`. Invalid counts return status
`bessel_zeros_status_invalid_argument` (1) and empty arrays. Numerical failure
returns `bessel_zeros_status_numerical_failure` (2) and empty arrays, never a
partial table. Success is `bessel_zeros_status_success` (0). Status is optional.

## Root sequence and associated values

The routines return the principal-branch upper-half-plane sequence corresponding
to SciPy's complex mode. Imaginary parts are positive; real parts decrease
strictly. Lower-half-plane conjugates are not included.

| Routine | Roots | Associated values |
| --- | --- | --- |
| `y0_zeros` | `Y_0(z)=0`, starting near `-2.4030+0.5399i` | Actual derivative `Y_0'(z)=-Y_1(z)` |
| `y1_zeros` | `Y_1(z)=0`, starting near `-0.5027+0.7862i` | `Y_1'(z)=Y_0(z)-Y_1(z)/z` |
| `y1p_zeros` | `Y_1'(z)=0`, starting near `+0.5768+0.9040i` | `Y_1(z)` |

The first complex `Y_1'` root has **positive** real part. Subsequent roots
have negative real part. This follows SciPy's observed root sequence despite
the SciPy documentation's description of negative real parts for all roots.

SciFort returns the actual `Y_0'` derivative, preserving its real-routine
convention. SciPy 1.15.3's `y0_zeros` companion instead equals `+Y_1`, the
opposite sign from its documented derivative. Negate that SciPy companion
when comparing it with SciFort. The other two companion conventions agree.

## Method and limitations

The implementation uses damped Newton refinement with existing complex Bessel
kernels. Derivatives come from the Bessel recurrences and differential equation.
Initial string phases and the imaginary asymptote `log(3)/2` follow DLMF;
the exceptional order-one roots use coarse geometric initial guesses.
Root-specific boxes prevent crossing onto another string, and acceptance
requires both a small Newton correction and a small residual relative to
the local function amplitude. Refinement is bounded to 60 iterations per root,
with at most 24 backtracking trials per iteration. Storage is linear in `nt`.

The 1000-root ceiling bounds work and the intended validation domain. It does
not apply to the existing real overloads. These routines are floating-point
root finders rather than interval-certified complete zero enumerators; their
accuracy also depends on the existing complex Bessel kernels. Large tables can
be substantially slower than small ones because complex function evaluations
use angular continuation. No shared kernel or build setting was changed.

See the [validation record](../validation/complex_bessel_zeros_2026-11-13.md)
for focused tests, root residuals, independent multiprecision checks, and
observed SciPy limitations. No C ABI overloads were added.

## Mathematical sources

The implementation is original MIT-licensed SciFort code based on
[DLMF recurrences](https://dlmf.nist.gov/10.6), the
[Bessel equation](https://dlmf.nist.gov/10.2.E1),
[analytic continuation](https://dlmf.nist.gov/10.11.E6),
[large-argument phases](https://dlmf.nist.gov/10.17), and
[complex zero geometry](https://dlmf.nist.gov/10.21#ix).
See [provenance](../CODE_PROVENANCE.md) for implementation and reference details.
