# Additional SciPy-compatible spectral window families

`scifort_signal.get_window` now supports 25 named window families, including
coefficient-vector `general_cosine` and DPSS single-taper windows. The elementary windows are independent; the
three advanced constructions are adapted from SciPy 1.17 (BSD-3-Clause), with
license attribution. There are no new Fortran dependencies. The returned one-dimensional array is `real(dp)`.

```fortran
use scifort_kinds, only : dp
use scifort_signal, only : get_window, spectral_status_success
real(dp), allocatable :: w(:)
integer :: status

! Periodic by default (appropriate for spectral estimation).
w = get_window('blackmanharris', 128, status=status)
! Symmetric window (typically for FIR filter design).
w = get_window('kaiser', 129, fftbins=.false., parameter=8.0_dp)
w = get_window('general_gaussian', 129, parameter=2.0_dp, parameter2=16.0_dp)
w = get_window('tukey', 128, parameter=0.25_dp)
! Dolph-Chebyshev attenuation in dB.
w = get_window('chebwin', 128, parameter=80.0_dp)
! Taylor's nbar and sidelobe attenuation in dB.
w = get_window('taylor', 128, parameter=5.0_dp, parameter2=70.0_dp)
! KBD is symmetric-only and requires an even number of samples.
w = get_window('kbd', 128, fftbins=.false., parameter=8.0_dp)
! Arbitrary origin-centered cosine-series coefficients.
w = get_window('general_cosine', 128, coefficients=[0.42_dp, 0.5_dp, 0.08_dp])
! DPSS time-half-bandwidth product NW: single, peak-normalized taper.
w = get_window('dpss', 128, parameter=3.5_dp)
```

Window names (aliases supported) are `boxcar`, `hann`, `hamming`, `blackman`,
`bartlett`, `flattop`, `tukey`, `triang`, `parzen`, `bohman`,
`blackmanharris`, `nuttall`, `barthann`, `cosine`, `lanczos`,
`exponential`, `kaiser`, `gaussian`, `general_gaussian`, and
`general_hamming`, `chebwin`, `taylor`, `kaiser_bessel_derived`,
`general_cosine`, and `dpss` (alias `slepian`).
The first seven were previously present.

| Name | `parameter` | `parameter2` |
| --- | --- | --- |
| `tukey` | alpha (default 0.5) | Not used |
| `exponential` | decay `tau` (default 1) | Not used |
| `kaiser` | beta (required) | Not used |
| `gaussian` | standard deviation (required, positive) | Not used |
| `general_gaussian` | exponent `p` (required, positive) | width `sig` (required, positive) |
| `general_hamming` | coefficient alpha (required) | Not used |
| `chebwin` | attenuation in dB (required) | Not used |
| `taylor` | integer `nbar` (default 4) | sidelobe suppression dB (default 30) |
| `kaiser_bessel_derived` / `kbd` | Kaiser beta (required) | Not used |
| `general_cosine` | `coefficients(:)` (required, not a scalar `parameter`) | Not used |
| `dpss` / `slepian` | Time-half-bandwidth product `NW` (required) | Not used |
| Other names | Not used | Not used |

`fftbins=.true.` (default) selects periodic windows; `fftbins=.false.`
selects symmetric windows. Following SciPy, the periodic family is computed
by extending the corresponding symmetric family by one element and dropping
the final sample, except for constant/one-element windows. The KBD window is symmetric-only
and requires an even positive length (or zero). With `get_window`, you must
set `fftbins=.false.` for KBD, unlike SciPy's other periodic-window defaults.
One behavioral
correction is that the *default* Tukey alpha is now 0.5 (previously 0.25);
call with `parameter=0.25_dp` to reproduce the old choice.

An invalid name, missing required parameter, invalid width, or negative
length returns a zero-length array and sets `status` to
`spectral_status_invalid_input`. SciFort continues to accept zero-length windows for the other families
(returning an empty array and success); SciPy 1.17's `get_window` rejects
`Nx=0` before dispatch. The new DPSS dispatch follows SciPy and rejects
zero-length windows. The direct `general_cosine` accepts an empty coefficient
vector, returning all zeros when `n > 1`. Like SciPy, windows of length 1
return `[1]` even if the coefficient sequence is empty.

## Direct advanced-window API

The following procedures are reexported by `scifort_signal` and
`scifort_signal_spectral`:

```fortran
use scifort_signal, only : chebwin, taylor, kaiser_bessel_derived, general_cosine
w = chebwin(129, 95.0_dp)    ! symmetric by default
w = chebwin(128, 95.0_dp, sym=.false.)
w = taylor(129, nbar=6, sll=80.0_dp, norm=.false.)
w = kaiser_bessel_derived(128, 12.0_dp)
w = general_cosine(129, [0.35875_dp, 0.48829_dp, 0.14128_dp, 0.01168_dp])
```

Direct functions default to symmetric windows (`sym=.true.`), as in SciPy's
`scipy.signal.windows` namespace. `get_window` defaults to periodic windows
(`fftbins=.true.`), also as in SciPy. Use direct `taylor` for `norm=.false.`;
`get_window('taylor', ...)` always normalizes the peak. `get_window` accepts
aliases `cheb`, `taylorwin`, `kbd`, and `kaiser bessel derived`.

All direct advanced-window functions return an empty array and an invalid-input
status for invalid shapes or arguments (if the caller requests `status`), rather than
raising exceptions. Chebyshev synthesis uses a log-scaled direct DFT, which
is quadratic in length; it prioritizes a self-contained, accurate reference
implementation rather than FFT efficiency. Taylor limits `nbar` to 1000 and
sidelobe suppression to 5000 dB. These numerical/parameter limits are
SciFort-specific protections.

## Verification

`python tools/signal_reference/compare_windows_extra.py` compiles the
included Fortran driver with GNU Fortran and checks 1,092 cases against
SciPy 1.17. The reference covers real window values, parameter variants,
symmetric and periodic windows, and lengths 1–127; largest observed absolute
difference was 1.43e-14. `test/test_signal_windows_extra.f90` and
`test/test_signal_windows_public.f90` add native and public-interface tests.

The advanced-window reference script `python test/compare_signal_windows_advanced.py`
compiles a Fortran driver and compares 1,718 deterministic/seeded cases with
SciPy 1.17, covering direct and dispatcher calls, odd/even lengths, symmetric
and periodic modes, Taylor peak/DC normalization, KBD constraints, and edge
lengths. Maximum absolute error was 4.82e-13. Native regressions are
`test_signal_windows_advanced.f90` and
`test_signal_windows_advanced_public.f90`.

The additional test `python test/compare_signal_cosine_dpss.py` validates
1,092 origin-centered general cosine and DPSS cases against SciPy 1.17
(maximum absolute error 2.64e-12). Targeted Fortran regressions are
`test_signal_windows_cosine_dpss.f90` and
`test_signal_windows_cosine_dpss_public.f90`.

`get_window('dpss', n, parameter=NW)` returns a single taper, like SciPy's
`get_window(('dpss', NW), n)`; use the separate `dpss` function for several
tapers and concentration ratios. `get_window('general_cosine', ..., coefficients=a)`
uses an array argument rather than SciPy's Python tuple. Existing DPSS and
cosine-series numerical kernels are reused; only the public API is extended.
The build order for standalone GNU Fortran invocations is `kinds`,
`signal_dpss`, `signal_windows_extra`, then `signal_spectral`.
For all families except DPSS, zero-length windows continue to follow the
older SciFort API rather than SciPy's stricter generic dispatcher.

## Kaiser filter design helpers

`kaiser_beta(attenuation)`, `kaiser_atten(numtaps, width)`, and
`kaiserord(ripple, width)` are also available through `scifort_signal`.
Use `firwin(numtaps, cutoff, window='kaiser', window_parameter=beta)` to
apply the resulting shape parameter. See
[Gaussian splines and Kaiser FIR design](signal_kaiser_spline.md) for complete
examples and differences in invalid-input behavior.
