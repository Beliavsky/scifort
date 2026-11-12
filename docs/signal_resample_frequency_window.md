# Fourier resampling with caller-provided spectral windows

**New:** Named-window shape keywords (`window_parameter`,
`window_parameter2`, `window_coefficients`) are documented separately in
[Parameterized Fourier resampling windows](signal_resample_parameter_windows.md).


SciFort implements the frequency-domain weighting in `scipy.signal.resample`
with **FFT-order array windows or callable windows** (zero frequency first, then positive and
negative frequencies). The public `scifort_signal` module exports both
array-returning `resample` (real/complex, ranks 1–3) and caller-output
`resample_into` (real/complex, ranks 1–15).

The optional arguments, appended after existing parameters for backward
compatibility, are:

- `frequency_window(:)` — real-valued FFT-bin weights for **real or complex**
  input samples.
- `complex_window(:)` — complex-valued FFT-bin weights for **complex**
  input samples.
- `window_function` — Fortran callback returning **real** FFT-bin weights.
- `complex_window_function` — Fortran callback returning **complex** FFT-bin
  weights; only available for complex input signals.

Pass at most one window form: `window='hann'`, `frequency_window=...`,
`complex_window=...`, `window_function=...`, or
`complex_window_function=...`. Window arrays must have exactly the input sample count
along the selected axis. Invalid lengths, nonfinite values, or conflicting
window arguments produce `fourier_resample_status_invalid_input`. The
caller-allocated output of `resample_into` is unchanged after errors.

## Examples

```fortran
use scifort_kinds, only : dp
use scifort_signal, only : resample, resample_into, &
    fourier_resample_status_success
real(dp) :: x(7), weights(7), data(2,7), output(2,11)
real(dp), allocatable :: y(:)
integer :: status

! Supply FFT bin weights in FFT order (not frequency-sorted order).
x = [1.0_dp,2.0_dp,-1.0_dp,4.0_dp,0.0_dp,1.0_dp,3.0_dp]
weights = [0.6_dp,1.1_dp,0.1_dp,1.3_dp,0.8_dp,0.4_dp,1.7_dp]
y = resample(x,11,frequency_window=weights,status=status)
if (status /= fourier_resample_status_success) error stop 'resample failed'
data(1,:) = x
data(2,:) = 2.0_dp*x
call resample_into(data,output,11,axis=2,frequency_window=weights,status=status)
```

A complex signal can use either real weights via `frequency_window=` or
complex weights via `complex_window=`, including when using precomputed
complex Fourier coefficients with `domain='freq'`.

## Callable window example

Like SciPy's `signal.resample(..., window=callback)`, SciFort evaluates a
user-defined function at all **FFT-bin frequencies** and applies its result
before discarding or interpolating Fourier coefficients. Callback frequencies
are in **cycles per input sample** (`fftfreq(n)`), not Hertz; time-grid `t`
changes the output times, not the callback frequencies. For even `n`, the
Nyquist bin is supplied as **-0.5**. For odd `n`, frequencies end with
`-1/n`. Frequencies are in native FFT order, not sorted order.

```fortran
module my_windows
    use scifort_kinds, only : dp
    implicit none
contains
    function frequency_response(f) result(w)
        real(dp), intent(in) :: f(:)
        real(dp) :: w(size(f))
        w = exp(-25.0_dp*f**2)
    end function frequency_response
end module my_windows

program demo
    use scifort_kinds, only : dp
    use scifort_signal, only : resample, resample_into
    use my_windows, only : frequency_response
    implicit none
    real(dp) :: x(8), multichannel(8,3), output(12,3)
    real(dp), allocatable :: y(:)
    integer :: status

    x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,3.0_dp,2.0_dp,1.0_dp,0.0_dp]
    y = resample(x,12,window_function=frequency_response,status=status)
    multichannel(:,1) = x
    multichannel(:,2) = 2.0_dp*x
    multichannel(:,3) = 3.0_dp*x
    call resample_into(multichannel,output,12,axis=1, &
        window_function=frequency_response,status=status)
end program demo
```

A complex-valued callback uses the same argument list, with a
`complex(dp) :: w(size(f))` result; pass it as
`complex_window_function=my_complex_response` for complex input.
Both array-returning `resample` (ranks 1–3) and caller-output
`resample_into` (ranks 1–15) evaluate the callback **once per invocation**,
not once per channel. The weight array is reused for every lane, including
noncontiguous array sections. A NaN or infinity in a returned window is an
invalid input; caller-output arrays remain unchanged. Callback exceptions
(such as an explicit Fortran `error stop`) cannot be caught by SciFort.

## Real-input convention

SciPy's real-input resampler uses a one-sided Fourier transform, so it
averages each pair of FFT-order weights indexed `k` and `-k mod n` before
multiplying the spectrum. SciFort does the same for both named and custom
real windows. This can differ from applying the original asymmetric weights
to a complex cast of the data; that complex interpretation remains available
through the complex-input overload.

One-based positive axis indices or negative axis indices are accepted on
rank-two/three array-returning overloads and rank-independent `resample_into`.
Existing time-grid arguments (`t`, `output_times`) and `domain` continue to
work. Note `complex_window=` and `domain='freq'` are restricted to the
complex-input interfaces.

## Testing and limitations

`test/test_signal_resample_frequency_window.f90` exercises the public API,
including noncontiguous outputs, rank-15 arrays and errors without partial
writes. Run `python test/compare_signal_resample_frequency_window_scipy.py`
for 900 SciPy 1.17 comparison cases (32,664 output values) and `SCIFORT_OPT=-O3`
for the optimized kernel. The observed maximum absolute difference is
`5.76e-14` with GNU Fortran 14.2.

`test/test_signal_resample_callable_window.f90` checks odd/even frequencies,
noncontiguous outputs, rank-15 inputs, one-shot evaluation and error recovery.
`python test/compare_signal_resample_callable_window_scipy.py` tests 450
independent SciPy 1.17 callable-window cases (12,722 values) across real and
complex signals, complex/real callbacks, multiple resampling ratios, and
`domain='freq'`. Maximum absolute difference in GNU Fortran 14.2 checked
and optimized tests was `5.69e-14`.

SciFort uses radix-2 FFTs for power-of-two lengths, a short direct DFT
for non-power-of-two lengths up to 16, and Bluestein FFTs for larger
non-power-of-two lengths. An explicit reusable `fft_plan` is available for
repeated transforms, but `resample` currently uses the unplanned FFT path.
No extra Fortran dependencies are needed.
