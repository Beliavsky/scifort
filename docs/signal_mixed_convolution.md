# Mixed real/complex one-dimensional convolution and correlation

`scifort_signal` now exposes both operand orders (`real(dp), complex(dp)`
and `complex(dp), real(dp)`) for the array-returning rank-one APIs:

| Fortran function | Output | Supported options |
|---|---|---|
| `convolve` | `complex(dp)[:]` | `mode`, `method`, `status` |
| `correlate` | `complex(dp)[:]` | `mode`, `method`, `status` |
| `fftconvolve` | `complex(dp)[:]` | `mode`, `status`, `axes=[1]` |
| `oaconvolve` | `complex(dp)[:]` | `mode`, `status` |
| `choose_conv_method` | `'direct'` or `'fft'` | `mode`, `measure`, `times`, `status` |

Inputs need not be manually converted to `complex(dp)`. **For correlation,
the second operand is conjugated** (even when the first is real). As with
SciPy's 1-D correlation, `mode='same'` uses the first signal's length and
the rank-one centering convention, which differs by one index from the
rank-two `correlate2d` convention for some even-length filters.

```fortran
use scifort_signal, only : dp, convolve, correlate, fftconvolve, oaconvolve
real(dp) :: x(4) = [1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp]
complex(dp) :: h(2) = [cmplx(1.0_dp,2.0_dp,dp), &
                        cmplx(0.0_dp,-1.0_dp,dp)]
complex(dp), allocatable :: y(:), r(:)

y = convolve(x, h, mode='same', method='auto')
r = correlate(h, x, mode='full', method='fft')
y = fftconvolve(x, h, mode='full', axes=[1])
y = oaconvolve(x, h, mode='valid')
```

The existing real-real and complex-complex APIs are unchanged, including
positional argument order. The new overloads reuse the established complex
kernels, so supported modes (`full`, `same`, `valid`), algorithm selection,
FFT slicing, shape rules, and status/error behavior are identical after
promoting the real argument. Strided rank-one sections are accepted. A bad
mode, bad axis, or empty signal returns an empty result and a nonzero `status`.

## Numerical verification

Run `python tools/run_signal_mixed_convolution.py`. With GNU Fortran, NumPy,
and SciPy 1.17 this runs native regression tests and 960 deterministic
independent SciPy comparisons in a checked build and a targeted optimized
build. The reference driver checks 182,304 complex values across both
operand orders, all three modes, direct/FFT/auto convolution and correlation,
and FFT/overlap-add convolution. The maximum error was `5.177e-13` in the
checked build and `5.154e-13` in the targeted `-O3` build.

The optimizations apply to the two modified convolution modules and test
programs; large unchanged FFT dependencies use `-O0` in the optimized build.
For a broader source check, `python tools/check_signal_mixed_public.py`
compiles the 91-source public `scifort_signal` dependency closure in
Fortran module order and runs the existing `test_signal_convolve` and
`test_signal_fftconvolve` programs, plus a public-interface compilation and
runtime smoke test of all five new generics. These also passed with GNU
Fortran 14.2. The public checker uses `build/signal_mixed_public` unless
`--build-dir` specifies another destination.

The focused runner and public checker do not replace the full `fpm test`
suite, which has not been run here.
