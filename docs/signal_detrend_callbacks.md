# Per-segment spectral detrending callbacks

SciFort's legacy spectral estimators support a trailing optional
`detrend_callback=` procedure argument. The callback runs once for each
**unwindowed** segment, before taper multiplication, Fourier transformation
and spectrum averaging. It provides the equivalent of SciPy's callable
`detrend` argument with an idiomatic Fortran in-place procedure interface.

The option is available in the real and complex vector implementations of
`periodogram`, `welch`, `csd`, `coherence`, `stft`, and `spectrogram`, as well
as the rank-independent `periodogram_into`, `welch_into`, `csd_into`,
`coherence_into`, `stft_into`, and `spectrogram_into` interfaces. The former
support vectors; the latter use the existing selected-axis and output-shape
rules (spectral outputs through rank 15, STFT/spectrogram input through rank 14).

## Callback signatures

The public `scifort_signal` module exposes
`spectral_real_detrend_callback` and `spectral_complex_detrend_callback`:

```fortran
abstract interface
    subroutine spectral_real_detrend_callback(segment)
        import dp
        real(dp), intent(inout) :: segment(:)
    end subroutine
    subroutine spectral_complex_detrend_callback(segment)
        import dp
        complex(dp), intent(inout) :: segment(:)
    end subroutine
end interface
```

Provide an ordinary procedure with the appropriate signature. The caller's
original signal is not modified: the library passes a private copy of each
segment. The real callback is selected for real-only estimators, and the
complex callback for complex-input estimators, including mixed real/complex
CSD and coherence. The callback must return finite values. A callback
producing any NaN or infinity causes the estimator to return
`spectral_status_invalid_input`. For `*_into`, the caller's preallocated
output is unchanged if any segment fails. Callbacks may be called several
times for a single signal (for example, coherence computes several spectra);
avoid relying on callback side effects or global state.

## Example: remove a quadratic trend

```fortran
program spectral_quadratic_example
    use scifort_signal, only : dp, spectral_result, welch, spectral_status_success
    implicit none
    real(dp) :: x(128), t
    integer :: i
    type(spectral_result) :: spectrum

    do i = 1, size(x)
        t = real(i, dp)
        x(i) = sin(0.2_dp*t) + 0.002_dp*t*t + 0.05_dp*t
    end do
    spectrum = welch(x, window='hann', nperseg=32, noverlap=16, &
        detrend_callback=remove_quadratic)
    if (spectrum%status /= spectral_status_success) error stop 'Welch failed'
    print *, spectrum%frequencies(1:5), spectrum%power(1:5)
contains
    subroutine remove_quadratic(segment)
        real(dp), intent(inout) :: segment(:) !! segment to detrend
        real(dp) :: t(size(segment)), q(size(segment)), av, slope, curvature
        integer :: k, n
        n = size(segment)
        do k = 1, n
            t(k) = (real(k,dp) - 0.5_dp*real(n+1,dp))/real(n,dp)
        end do
        q = t*t
        q = q - sum(q)/real(n,dp)
        av = sum(segment)/real(n,dp)
        slope = sum(segment*t)/sum(t*t)
        curvature = sum(segment*q)/sum(q*q)
        segment = segment - av - slope*t - curvature*q
    end subroutine remove_quadratic
end program spectral_quadratic_example
```

For a real multi-channel array with a channel-specific time axis, pass the
same procedure to `welch_into` and specify `axis=1` (one-based). SciFort
applies it independently to each channel; there is no need to reshape or
copy the input manually. The `*_into` procedures still require correctly
sized caller-owned output arrays.

## Validation and defaults

`detrend_callback=` and the existing character-valued `detrend=` option
are mutually exclusive. Specifying both returns invalid-input status;
neither takes precedence. When no callback is provided, existing constant,
linear and no-detrend behavior remains unchanged. The default detrend value
is still chosen independently by each estimator: periodogram/Welch/CSD/
coherence use a constant, STFT uses none, and spectrogram uses a constant.
With a callback, the estimator skips that built-in detrending step entirely.
No new windowing/scaling conventions are introduced.

This is an **API adaptation** rather than a direct signature translation:
SciPy accepts a Python callable through `detrend=`, which returns a transformed
array. SciFort's new `detrend_callback=` has an in-place `intent(inout)`
procedure and is deliberately separate from the string-based `detrend=`.

Run the new numerical checks with:

```sh
python tools/run_signal_detrend_callbacks.py
```

The standalone runner compiles seven core spectral source modules in checked
and `-O3` modes, runs a native test covering both output ranks, mixed-type
cross spectra, noncontiguous arrays, invalid specifications, nonfinite callback
outputs and preservation of caller-owned buffers, and checks 48 independent
SciPy 1.17 reference cases comprising 2,098 spectral values and both coordinate
grids. Maximum absolute differences in this environment were `3.122e-13`
(checked) and `7.264e-13` (optimized). The previously existing 440-case
spectrogram-mode and 192-case two-sided-spectrum SciPy suites also passed in
checked mode. The top-level `scifort_signal` module and its 91-source dependency
closure passed a compiler syntax check. The full FPM test suite was not run.

Provenance: original SciFort Fortran implementation using the existing FFT,
window, channel-walking and detrending kernels. Algorithmic reference:
SciPy's documented public callable detrend support,
https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html .
SciPy and NumPy are only used in the external comparison script.
