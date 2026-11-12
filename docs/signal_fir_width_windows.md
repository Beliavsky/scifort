# FIR filter-design transition widths and parameterized windows

SciFort extends two existing `scipy.signal`-compatible FIR interfaces. Both
retain their original positional argument order and use `dp` real coefficients.

## `firwin`: transition-width Kaiser designs

`firwin(numtaps, cutoff, ..., window_parameter, width)` now accepts an
optional, trailing real `width` keyword. `width` is an approximate transition
width **in the same frequency units as `fs`**, not normalized to the sampling
frequency. When present, it selects a symmetric Kaiser window with

```text
attenuation = kaiser_atten(numtaps, width / (fs / 2))
beta        = kaiser_beta(attenuation)
```

As with SciPy, a supplied transition width **overrides** `window` (including
`window_parameter`) rather than combining with it. A width must be finite and
strictly positive. The default sample frequency remains `fs=2.0_dp`.

```fortran
use scifort_signal, only : dp, firwin, design_status_success
real(dp), allocatable :: taps(:)
integer :: st

taps = firwin(41, [350.0_dp], fs=4000.0_dp, width=120.0_dp, &
              pass_zero='lowpass', status=st)
if (st /= design_status_success) error stop 'firwin failed'
```

For an explicit Kaiser beta rather than an approximate transition width, use
`window='kaiser', window_parameter=beta` as before. `width` does not change
existing band-edge or parity restrictions: an even number of taps cannot have
nonzero response at Nyquist.

## `firwin`: scalar cutoffs, logical modes, and odd-edge multiband designs

The `firwin` generic now accepts real scalar or rank-1 cutoff arguments. For
both forms, `pass_zero` accepts a logical or a character choice, and all
previous positional arguments (including `window_parameter`, `width`, and
`status`) keep their order. In Fortran the overload is selected by type, not
by an optional flag: calls using logical `pass_zero` resolve automatically.

```fortran
use scifort_signal, only : dp, firwin, design_status_success
real(dp), allocatable :: taps(:)
integer :: st

! Single cutoff is a scalar; lowpass is the default.
taps = firwin(31, 0.35_dp, status=st)
if (st /= design_status_success) error stop 'scalar FIR design failed'

! A 3-transition high/low/high multiband response (Nyquist passband).
! This needs odd numtaps.
taps = firwin(31, [0.12_dp, 0.34_dp, 0.67_dp], &
              pass_zero=.false., window='kaiser', &
              window_parameter=4.0_dp, status=st)
if (st /= design_status_success) error stop 'multiband FIR design failed'

! Complementary 3-transition design starts with a DC passband and ends
! with a stopband at Nyquist; an even tap count is allowed.
taps = firwin(30, [0.12_dp, 0.34_dp, 0.67_dp], &
              pass_zero=.true., width=0.035_dp, status=st)
if (st /= design_status_success) error stop 'even-length multiband FIR failed'
```

For `m` cutoff frequencies, the Nyquist passband flag is
`(mod(m,2)==0) .eqv. pass_zero`: the DC passband selection flips at each edge.
With `pass_zero='bandpass'` the DC band is rejected; with
`pass_zero='bandstop'` the DC band is accepted. The string names `lowpass`
and `highpass` still require exactly one cutoff. String `bandpass` and
`bandstop` now accept any cutoff count of at least two, not just pairs.
An even-length filter is valid only if its Nyquist band is a stopband.

The native test and independent SciPy 1.17 reference runner are:

```bash
python tools/run_signal_firwin_multiband.py
```

The checked and optimized builds each verify 443 valid reference designs
(10,370 coefficients, maximum absolute error `4.44e-16`) plus 362 intentionally
invalid cases with correct error status, as well as native symmetry, scaling,
parity, scalar-dispatch, and Kaiser-width regression checks.

## `firwin2`: parameterized symmetric windows

`firwin2` now accepts four trailing optional keywords:

| Keyword | Meaning |
| --- | --- |
| `window_parameter` | First named-window shape, such as Kaiser beta, Gaussian sigma or Tukey alpha |
| `window_parameter2` | Second shape for `general_gaussian` or `taylor` |
| `window_coefficients` | General cosine coefficient array when `window='general_cosine'` |
| `window_values` | Explicit real symmetric window samples, one per tap |

SciFort's existing `spectral_named_window` validation is shared with the
spectral estimators. In `firwin2`, the window is generated with
`fftbins=.false.` to match SciPy's **symmetric** FIR-design convention.
Nonfinite values and incompatible shape parameters are rejected by returning
`fir_advanced_status_invalid_input` and an empty coefficient array.

```fortran
use scifort_signal, only : dp, firwin2, fir_advanced_status_success
real(dp), allocatable :: taps(:)
integer :: st

taps = firwin2(29, [0.0_dp, 0.25_dp, 0.7_dp, 1.0_dp], &
               [1.0_dp, 0.8_dp, 0.2_dp, 0.0_dp], &
               window='gaussian', window_parameter=5.0_dp, status=st)
if (st /= fir_advanced_status_success) error stop 'firwin2 failed'
```

`window_values` is a **SciFort extension**, not an argument offered by SciPy's
`firwin2`. It may not be combined with a named window or any parameter. Its
length must equal `numtaps`; the values are multiplied directly into the
unwindowed FIR coefficients. `window='none'` remains the existing unwindowed
option and cannot be given shape parameters.

## Validation

Run from the project root:

```bash
python tools/run_signal_fir_width_windows.py
```

The runner compiles with GNU Fortran 2018, first with `-fcheck=all` and then
with `-O3`, executes native invalid-option regressions, and compares 320
filter designs (6,736 coefficients) to SciPy 1.17 with absolute errors below
`4.45e-16`. The comparison includes lowpass, highpass, bandpass and bandstop
Kaiser-width designs; FIR types I–IV for `firwin2`; finite symmetric windows,
two-parameter windows and explicit window arrays. SciPy 1.17 cannot pass
`general_cosine` through `firwin2` due to its `xp`/`device` forwarding: the
reference checks this family against SciPy's unwindowed `firwin2` result
multiplied by a separately generated symmetric `general_cosine` window.

Existing verification scripts were also rerun:

```bash
python tools/signal_reference/verify_design.py
python tools/signal_reference/verify_fir_advanced.py
python tools/run_signal_twosided.py --checked-only
```

SciPy/NumPy are reference-only Python dependencies. The Fortran library adds
no external dependencies, no third-party implementation code, and no new
exported ABI symbols.
