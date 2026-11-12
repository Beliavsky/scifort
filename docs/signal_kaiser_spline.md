# Gaussian spline approximation and Kaiser FIR design

SciFort exposes four independent MIT-licensed implementations of the
corresponding `scipy.signal` numerical algorithms through `scifort_signal`:
`gauss_spline`, `kaiser_beta`, `kaiser_atten`, and `kaiserord`. No new external
Fortran dependencies are required. All real arguments and results use `dp`.

```fortran
use scifort_signal, only : dp, gauss_spline, kaiser_beta, kaiser_atten, &
    kaiserord, kaiser_order_result, firwin, design_status_success
real(dp) :: knots(3), beta, atten
real(dp), allocatable :: taps(:)
type(kaiser_order_result) :: design
integer :: status

knots = gauss_spline([-1.0_dp, 0.0_dp, 1.0_dp], 3)
beta = kaiser_beta(65.0_dp)  ! 6.20426
atten = kaiser_atten(211, 0.0375_dp)  ! approximately 64.481 dB

design = kaiserord(65.0_dp, 24.0_dp/500.0_dp)
if (design%status == design_status_success) then
    ! A symmetric Kaiser FIR window with SciPy-compatible shape.
    taps = firwin(design%numtaps, [0.35_dp], &
        window='kaiser', window_parameter=design%beta, status=status)
end if
```

The `kaiserord` result has `numtaps`, `beta` and `status` components. This is
SciFort's derived-type equivalent of SciPy's `(numtaps, beta)` tuple. Widths
are normalized to Nyquist frequency (one corresponds to pi radians/sample).
The Kaiser formulas are empirical design approximations, not guaranteed
frequency-response bounds. `firwin` still returns an allocatable array and
sets an optional status code; adding `window_parameter` at the end preserves
existing positional calls. The optional parameter is rejected for window
families without a shape argument.

`gauss_spline(x, n)` is elemental: a scalar or an array of any supported Fortran
rank can be supplied as `x`, and the spline order `n` is nonnegative. It
calculates the normal approximation with variance `(n+1)/12`. For a negative
order, the Fortran function returns a quiet NaN rather than raising a Python
exception. Large-magnitude knots are evaluated without intermediate overflow.

`kaiser_beta` and `kaiser_atten` are also elemental. For `kaiserord`, SciFort
returns `design_status_invalid_input` when attenuation magnitude is below 8 dB,
width is nonpositive, either argument is nonfinite, or the inferred tap count
overflows the default Fortran integer. This is deliberately stricter than
SciPy for negative or zero widths.

## Reproducible verification

Run the native regression:

```sh
python tools/signal_reference/run_signal_suite.py --jobs 6
python tools/signal_reference/run_signal_suite.py --optimized --jobs 6
```

The SciPy comparison script compiles a small independent Fortran driver:

```sh
python tools/signal_reference/compare_signal_kaiser_spline.py
```

With GNU Fortran 14.2.0 and SciPy 1.17.0 the new functions and Kaiser
`firwin` path passed 2,388/2,388 numerical comparisons, maximum absolute
error `7.77e-16`. The full checked and staged optimized signal suites passed
102/102 tests each. The new standalone routines also compiled at `-O3` with
IEEE invalid, zero-divide, and overflow traps enabled. Project-wide FPM
integration remains unverified here.

Mathematical sources: Oppenheim and Schafer, *Discrete-Time Signal
Processing* (Kaiser empirical formulas); Bouma et al. (2007), *Fast and
Accurate Gaussian Derivatives Based on B-Splines* (Gaussian approximation).
Public SciPy documentation was consulted for API conventions, and SciPy 1.17
was used for independent numeric reference comparisons. No SciPy
implementation code was copied or linked.
