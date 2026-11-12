# Mixed real/complex cross-spectral analysis (2026-10-28)

SciFort now accepts one real and one complex input in the **rank-one**
`csd`, `coherence`, and `short_time_fft%cross_spectrogram` interfaces. The
existing same-type overloads and rank-independent `csd_into` /
`coherence_into` interfaces are unchanged; the latter already support mixed
channels. All mixed spectra use **two-sided FFT bins**, including the negative
Nyquist bin for an even `nfft`. A requested `return_onesided=.true.` on mixed
`csd` is ignored just as in the complex/complex overload.

```fortran
use scifort_signal, only : dp, csd, coherence, cross_spectral_result, &
    spectral_result, short_time_fft, short_fft_success
real(dp) :: x(256)
complex(dp) :: y(256)
type(cross_spectral_result) :: cross
type(spectral_result) :: consistency
type(short_time_fft) :: analysis
complex(dp), allocatable :: cross_tf(:,:)
integer :: status

! Supply the samples to x and y before these calls.
cross = csd(x, y, nperseg=64, noverlap=32, nfft=128,
            window='hann', average='median')
consistency = coherence(y, x, nperseg=64, noverlap=32)
call analysis%initialize_window('hann', 32, 8, fs=100.0_dp,
                                mfft=64, fft_mode='twosided')
if (analysis%status /= short_fft_success) error stop 'plan initialization'
call analysis%cross_spectrogram(x, y, cross_tf, status, padding='zeros')
if (status /= short_fft_success) error stop 'cross spectrogram'
```

`csd` preserves the standard cross-spectrum order: `conjg(FFT(x))*FFT(y)`;
`short_time_fft%cross_spectrogram` uses its own existing convention:
`STFT(x)*conjg(STFT(y))`. Interchanging inputs therefore conjugates the
result. Both mixed-order overloads are available and accept noncontiguous
one-dimensional array sections.

`csd` supports the same `fs`, `window`, `nperseg`, `noverlap`, `nfft`,
`detrend`, `scaling`, `average`, `window_parameter`, `window_parameter2`,
`window_coefficients`, `window_values`, and `detrend_callback` options as
the complex/complex implementation. The callback receives a **complex**
rank-one segment even for an originally real channel; it cannot be combined
with a string `detrend` argument. Mixed `coherence` shares the matching
options (except the `average` and `scaling` options, which are not in the
legacy SciPy `coherence` signature). The `short_time_fft` method accepts
`p0`, `p1`, and padding (`zeros`, `edge`, `even`, `odd`). It requires a
**two-sided or centered** FFT plan because one of the inputs is complex;
onesided/onesided2X plans reject mixed cross-spectrograms with
`short_fft_invalid_input`, returning an empty output. Unequal input
lengths are also rejected by the short-time cross-spectrum method.

## SciPy validation note

SciPy 1.17's *legacy* `signal.csd` and `signal.coherence` currently have
order-dependent errors when given mixed real/complex input. The reference
driver therefore promotes the real operand explicitly to `complex128` for
both legacy estimators. This tests the intended numerical formula against
SciPy without depending on its mixed-dtype dispatcher. The newer
`ShortTimeFFT.spectrogram(x, y)` accepts mixed input directly and is
compared without casting. See `tools/compare_signal_mixed_cross.py`.

Validation uses **384 independent SciPy 1.17 comparisons** with **31,008**
complex or real coefficients, covering both input orders, periodic Hann and
Gaussian windows, mean/median averaging, none/linear detrending, spectrum
scaling, frequency grids, boundary padding and ShortTimeFFT phase shifts.
Maximum observed error is below `5.2e-13` in checked and targeted optimized
GNU Fortran builds. Native regression also covers a complex detrend
callback, input striding and invalid requests. Run the standalone suite:

```bash
python tools/run_signal_mixed_cross.py
```

The runner checks the modified spectral and short-time FFT modules at `-O3`;
the unchanged window/DPSS/FFT dependencies are compiled at `-O0` for that
optimized variant. It also runs the existing native `test_signal_short_time_fft`,
`test_signal_short_time_fft_helpers`, and `test_signal_spectrogram_modes` tests.
The full FPM suite is separate and is not implied by these results.
