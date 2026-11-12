# Complex-valued local Wiener filtering

The `wiener` generic in `scifort_signal` and `scifort_signal_denoise`
now accepts `complex(dp)` input arrays of ranks 1, 2, and 3. This extends
the existing real and default-integer overloads without changing their
signatures or result types.

```fortran
use scifort_signal, only : dp, wiener
complex(dp) :: samples(128), image(12, 24)
complex(dp), allocatable :: filtered(:), smooth(:,:)
integer :: status

filtered = wiener(samples, kernel_size=5, status=status)
smooth = wiener(image, kernel_shape=[3, 5], &
    noise=cmplx(0.02_dp, 0.01_dp, dp), status=status)
```

`kernel_size=` applies one size to all axes; `kernel_shape=` accepts one
positive length per axis. The default is three per axis. Even window sizes
and windows larger than the data are permitted. Filtering uses zero padding,
including at the boundaries, and returns the input array's shape. It accepts
noncontiguous array sections.

## Noise options and SciPy compatibility

`noise=` can be omitted (estimated from the mean local variance), a
`real(dp)` scalar or a `complex(dp)` scalar. A complex explicit noise value
allows exact modeling of SciPy's complex-valued noise option. The result
always has type `complex(dp)`.

**Important:** SciPy 1.17 computes the complex *pseudo-variance*

\[
V = E[X^2] - (E[X])^2
\]

instead of the conventional magnitude-squared variance
`E[abs(X)**2] - abs(E[X])**2`. Consequently both its locally estimated
variance and its default estimated noise may be complex. NumPy's
comparison `V < noise` uses the real component first, then the imaginary
component to break ties. SciFort reproduces this convention so that the
result agrees numerically with the public `scipy.signal.wiener` API.

When the local pseudo-variance compares below the noise estimate, the
output equals the local mean; otherwise the gain is
`1 - noise / V`. An exact zero pseudo-variance with zero noise yields a
quiet NaN, consistent with the indeterminate gain in SciPy. This legacy
behavior should not be assumed equivalent to a physically motivated
complex Wiener estimator based on signal power.

Nonfinite inputs, nonfinite explicit noise values, zero or negative window
sizes, conflicting window arguments or an incorrectly sized
`kernel_shape=` result in `denoise_status_invalid_input` (1) and a
shape-preserving complex array filled with quiet NaNs. The default
success status is `denoise_status_success` (0). In contrast, existing
real/integer implementations retain their previous failure conventions.

## Validation

Run from the project root with GNU Fortran, Python, NumPy and SciPy installed:

```bash
python tools/run_signal_wiener_complex.py
```

The runner compiles with runtime checks and `-O3`, executes the native
regressions including strided arrays and input validation, compares 27
complex configurations/1,161 output samples independently against SciPy
1.17, and reruns the 1,050 existing real/integer median/Wiener reference
cases. The greatest observed complex absolute error is `2.066e-15` in
both build modes. The full FPM suite was not run during this iteration.
