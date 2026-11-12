# Arbitrary numerator/denominator filters in `decimate`

SciFort's `decimate` and `decimate_into` now accept digital transfer-function
coefficients via **trailing** `b=` and `a=` keywords. This extends the
previous Hamming FIR, Chebyshev-I IIR, caller-provided FIR `taps=`, and
real second-order-section `sos=` choices. No existing positional arguments
have changed.

The coefficient vectors are in increasing lag order:

\[
H(z) = \frac{b_0 + b_1 z^{-1} + \cdots}{a_0 + a_1 z^{-1} + \cdots},
\quad a_0 \neq 0.
\]

Each vector must be nonempty, finite, and have a nonzero leading denominator.
Both `b` and `a` must be supplied together; they cannot be combined with
`taps=`, `sos=`, `ftype=`, or `n=`. SciFort reports
`resample_status_invalid_input` if these requirements fail.

| Input | Coefficient types | Output |
| --- | --- | --- |
| `real(dp)` | both `real(dp)` | `real(dp)` |
| `complex(dp)` | both `complex(dp)` | `complex(dp)` |
| `real(dp)` | both `complex(dp)` | `complex(dp)` |

For complex data with real-valued transfer functions, convert the coefficient
arrays using `cmplx(b,0.0_dp,dp)` and `cmplx(a,0.0_dp,dp)`. The pre-existing
`taps=` and the original `sos=` overloads remain real-valued for real or
complex input. Additional complex-SOS overloads are documented in
[`signal_decimate_complex_sos.md`](signal_decimate_complex_sos.md).
The Fortran generic requires the two coefficient arrays to have the same kind
and real/complex type, rather than silently mixing types.

```fortran
use scifort_signal, only : dp, decimate, decimate_into
implicit none
real(dp) :: x(120), b(3), a(3), channels(2,120)
real(dp) :: output(2,40)
real(dp), allocatable :: y(:)
complex(dp) :: cb(3), ca(3), complex_output(2,40)
integer :: i, status

x = [(sin(0.1_dp*real(i,dp)),i=1,120)]
b = [0.22_dp,0.09_dp,0.06_dp]
a = [1.0_dp,-0.56_dp,0.18_dp]
y = decimate(x,3,b=b,a=a,status=status)
! Explicit causal IIR, then downsampling:
y = decimate(x,3,b=b,a=a,zero_phase=.false.,status=status)

channels(1,:) = x
channels(2,:) = 2.0_dp*x
call decimate_into(channels,output,3,b=b,a=a,axis=-1,status=status)

! Complex coefficients may yield complex results from real input.
cb = cmplx(b,[0.06_dp,0.04_dp,-0.02_dp],dp)
ca = cmplx(a,[0.00_dp,0.12_dp,-0.04_dp],dp)
call decimate_into(channels,complex_output,3,b=cb,a=ca, &
                   zero_phase=.false.,axis=-1,status=status)
```

When `size(a)==1`, `decimate` treats the coefficients as a FIR filter,
normalizing taps by `a(1)` before applying the existing `upfirdn` or
polyphase delay-compensated algorithm. Longer denominators select a direct
transfer-function IIR filter and reuse SciFort's `lfilter` for causal mode
or `filtfilt` for zero-phase mode. In all cases output length is
`ceiling(size(x,axis)/q)` and `zero_phase` defaults to `.true.`.
For IIR filters, an adequately long input is required for the default
forward/backward padding; as with SciPy, the caller is responsible for
filter stability, alias rejection and sensible cutoff selection.

The interface is available for rank-one, rank-two and rank-three
array-returning `decimate` calls, and for rank-one through rank-15
assumed-rank `decimate_into` calls. Real-to-complex overloads use `b=` and
`a=` as required keywords; the existing same-type overloads retain optional
`b=` and `a=` with error handling for missing counterparts. Axes are
one-based or negative (default is the final axis). `decimate_into` evaluates
all lanes into a staged buffer and does not alter the caller's output on
failed validation, incorrect output shape or failed filtering.

These explicit `b/a` options provide a Fortran equivalent of the numeric
transfer function of SciPy's `signal.dlti` supplied to `signal.decimate`;
there is no Python `dlti` object or implicit conversion from poles/zeros.
SciFort uses direct `filtfilt` for all `b/a` IIR systems. SciPy may instead
convert real-pole DLTI systems to SOS internally, so floating-point
roundoff can differ. Long or ill-conditioned high-order filters are better
represented using `sos=` when the coefficients are real.

## Validation

From the repository root, with GNU Fortran, NumPy and SciPy available:

```sh
python tools/run_signal_decimate_transfer_function.py
```

The checked and targeted `-O3` builds each compare 576 independent SciPy
1.17 test cases (16,320 output samples) spanning real and complex inputs,
real and complex IIR coefficients, real/complex FIRs, nonunit `a0`, both
phase modes, multiple factors and signal lengths. Maximum absolute error is
`8.951e-16` for the seeded comparison. Native regressions exercise real and
complex rank-1 to rank-3 paths, rank-15 assumed-rank output, strided
array sections, mismatched filter arguments, nonfinite coefficients and
transactional output preservation. The runner also reruns the previous 336
custom FIR/SOS SciPy cases and the 53-case and 184-case multirate suites.

Optimized compilation covers both modified Fortran modules and the test
drivers; unchanged dependencies are compiled with runtime checks. The full
`fpm test` suite is not verified in the validation environment because FPM
is unavailable. The Fortran implementation adds no third-party dependency.
