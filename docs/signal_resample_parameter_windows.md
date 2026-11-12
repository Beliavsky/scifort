# Parameterized Fourier-domain resampling windows

`scifort_signal::resample` now supports the same named-window shape parameters
as SciFort's existing spectral estimators and `resample_poly`, using additional
optional Fortran keywords:

- `window_parameter`: first real shape parameter, such as Gaussian standard
  deviation, Kaiser beta, Tukey alpha or DPSS time-bandwidth product.
- `window_parameter2`: second real parameter for generalized Gaussian
  (`power`, `width`) and Taylor (`nbar`, `sll`) families.
- `window_coefficients(:)`: cosine-series coefficients for `general_cosine`.

These keywords work with real and complex vector `resample` calls, rank-2/3
axis-selecting `resample` functions, and rank-1 through rank-15
caller-allocated `resample_into`. Existing positional arguments and the
traditional default (no spectral window) remain unchanged.

```fortran
use scifort_kinds, only : dp
use scifort_signal, only : resample, resample_into
real(dp) :: values(16), channels(16,2), enlarged(25,2)
real(dp), allocatable :: result(:)
integer :: status

result = resample(values, 25, window='gaussian', &
                  window_parameter=3.5_dp, status=status)
call resample_into(channels, enlarged, 25, axis=1, &
                   window='kaiser', window_parameter=4.5_dp, status=status)
```

The named taper uses `get_window(..., fftbins=.true.)` to construct the
**periodic window in the frequency domain**. As in SciPy, its samples are
FFT-shifted to zero-first FFT-bin order before multiplying the input FFT.
This is *not* time-domain FIR windowing. For real input, the positive- and
negative-frequency taper values are averaged. For complex inputs, all FFT
bins retain their own weights, including with `domain='freq'`.

Named-window shape parameters require `window='family'` and cannot be
combined with `frequency_window`, `complex_window`, `window_function` or
`complex_window_function`. A parameter supplied for a family that does not
accept it, or a second parameter for a one-parameter family, is rejected.
Invalid requests return `fourier_resample_status_invalid_input`; for the
`resample_into` subroutine the caller-provided output is unchanged.

For a general-cosine window, use e.g.

```fortran
result = resample(values, 25, window='general_cosine', &
                  window_coefficients=[0.6_dp, 0.3_dp, 0.1_dp])
```

Verification is reproducible with:

```bash
python tools/run_signal_resample_param_windows.py
```

This runs checked and targeted optimized GNU Fortran builds, 252 cases against
SciPy 1.17 (3,696 values; maximum absolute discrepancy `1.523e-14`) and
existing 53-case/184-case multidimensional/complex multirate comparisons.
SciPy 1.17 rejects `general_cosine` parameter tuples in `resample`'s array-API
path, so for that family the test independently evaluates SciPy's
`windows.general_cosine(..., sym=False)` and passes its FFT-shifted values
as a frequency-bin array. All other families use SciPy's `resample` with
native window tuples. NumPy and SciPy are test-only dependencies; no external
Fortran runtime dependency was added. The FPM-wide suite was not run here.
