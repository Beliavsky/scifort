# Cubic and quadratic spline interpolation

SciFort now provides four `scipy.signal` B-spline interfaces through
`scifort_signal`:

- `cspline1d(samples, lamb=0, status=...)` computes cubic spline coefficients.
- `qspline1d(samples, lamb=0, status=...)` computes quadratic spline coefficients.
- `cspline1d_eval(coeff, newx, dx=1, x0=0, status=...)` evaluates cubic splines.
- `qspline1d_eval(coeff, newx, dx=1, x0=0, status=...)` evaluates quadratic splines.

The coefficient functions accept rank-one `real(dp)` or `complex(dp)` samples.
Evaluation accepts rank-one real locations and returns either a rank-one array
or a scalar, matching the coefficient kind. `dx` and `x0` describe the knot
spacing and origin. Both evaluation functions handle out-of-range positions
through mirror symmetry, including positions many knot intervals beyond either
endpoint. Invalid calls return an empty allocatable array and the nonzero
`spline_status_invalid_input` code; scalar evaluation returns a quiet NaN.

```fortran
use scifort_signal, only : dp, cspline1d, cspline1d_eval, &
    qspline1d, qspline1d_eval, spline_status_success
real(dp), allocatable :: cubic(:), quadratic(:), out(:)
real(dp) :: sample(6), times(9)
integer :: status, i
sample = [1.0_dp, 0.25_dp, 2.0_dp, -1.0_dp, 0.75_dp, 3.0_dp]
times = [(real(i - 2, dp), i=1,9)]
cubic = cspline1d(sample, status=status)
if (status == spline_status_success) then
    out = cspline1d_eval(cubic, times, status=status)
end if
quadratic = qspline1d(sample)
```

For nonzero cubic smoothing, pass `lamb` to `cspline1d`. SciPy's smoothed
cubic recurrence is defined for positive `lamb` but becomes unstable for some
small values; SciFort rejects cases with a non-decaying recursive pole
(`abs(rho) >= 1`), as well as `lamb <= 1/144` and inputs with fewer than
three samples. This deliberate safeguard differs from SciPy 1.17's unstable
output for some small smoothing coefficients. A smooth cubic calculation with
`lamb=0.2_dp` is supported. As in SciPy, `qspline1d` does **not** implement
nonzero smoothing and reports an invalid-input status when requested.

For single-sample coefficient arrays, the coefficient routines preserve SciPy's
existing special-case normalization. Evaluation of a singleton coefficient
array is defined for any finite coordinate in SciFort instead of allowing
SciPy's recursive mirror-boundary handling to recurse indefinitely. SciFort
also supports complex quadratic evaluation by directly evaluating the complex
coefficients; SciPy 1.17 currently fails for this case owing to a real-valued
scratch array.

## Verification

```sh
python tools/signal_reference/compare_signal_spline_scipy.py
python tools/signal_reference/run_signal_suite.py --jobs 5
```

A deterministic set of 480 cases compares spline coefficients and evaluation
against SciPy 1.17.0, including real and complex inputs, cubic smoothing,
nonunit spacing, arbitrary origins, and mirror boundaries. For complex
quadratic splines, SciPy is evaluated separately on the real and imaginary
parts. Maximum coefficient difference was `1.351e-15`; maximum evaluation
difference was `2.798e-14`.

The added `test/test_signal_spline.f90` is also included in the standalone
signal regression suite. GNU Fortran's runtime-checked suite passed 103/103
tests; the focused spline test also passed at `-O3` with IEEE traps enabled.
An attempted complete optimized rebuild exceeded the execution limit, and
project-wide `fpm test` could not be run here because FPM was unavailable.

The recursion and mirror-boundary algorithm are adapted from SciPy 1.17's
BSD-3-Clause `scipy.signal._spline_filters` implementation. Its upstream
copyright and license terms are recorded in `THIRD_PARTY_LICENSES.md` and
`CODE_PROVENANCE.md`; no external dependencies are introduced.
