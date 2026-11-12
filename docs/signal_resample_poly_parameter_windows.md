# Parameterized polyphase resampling windows

SciFort's `resample_poly` now accepts shape parameters for its FIR prototype
window, following the named-window conventions of SciPy's
[`scipy.signal.resample_poly`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample_poly.html).
The new keywords are appended to existing argument lists, preserving the
existing positional calls:

| New optional keyword | Purpose | Examples |
| --- | --- | --- |
| `window_parameter` | First shape parameter | Gaussian `std`, Kaiser `beta`, Tukey `alpha`, generalized Gaussian `power`, Taylor `nbar` |
| `window_parameter2` | Second shape parameter | Generalized Gaussian `width`, Taylor sidelobe level `sll` |
| `window_coefficients` | Real cosine-series coefficients | `general_cosine` family |

These arguments are supported by all real and complex vector `resample_poly`
overloads, the existing axis-selected rank-2/3 array-returning overloads,
and the caller-allocated, assumed-rank `resample_poly_into` overloads through
rank 15. Complex FIR tap overloads retain their existing behavior.

```fortran
use scifort_kinds, only : dp
use scifort_signal_resample, only : resample_poly
implicit none
real(dp) :: x(64)
real(dp), allocatable :: y(:)
integer :: st

x = 1.0_dp
y = resample_poly(x, 3, 2, window='gaussian', &
                  window_parameter=4.0_dp, status=st)
! Equivalent window shape to SciPy's window=('gaussian', 4.0).
```

For two parameters, pass `window='general_gaussian'`,
`window_parameter=1.3_dp`, and `window_parameter2=4.4_dp`. General cosine
uses, for example, `window='general_cosine'` and
`window_coefficients=[0.42_dp, 0.50_dp, 0.08_dp]`.

The implementation generates a symmetric prototype window, multiplies it by
the ideal low-pass sinc, normalizes the taps to unit DC gain, and uses the
existing polyphase filtering and endpoint padding implementation. The default
unparameterized Kaiser window and legacy `beta=` behavior are unchanged.
The new parameters must not accompany explicit `taps` or `ctaps`; `beta=`
may not be combined with them. Parameters not applicable to the selected
window family are rejected rather than ignored. On invalid input,
rank-one and array-returning calls return an empty output and set
`resample_status_invalid_input`; `resample_poly_into` leaves the supplied
output unchanged.

The same interfaces support all previously implemented boundary modes,
including `constant`, `mean`, `median`, `line`, and `reflect`, without
changing their defaults.

## Verification

Run from the project root (requires `gfortran`, Python, NumPy, and SciPy):

```sh
python tools/run_signal_poly_windows.py
```

The standalone runner tests the native Fortran regression, and 256
independently generated real/complex SciPy 1.17 comparisons covering 6,976
output samples (maximum absolute error `2.987e-15`). It runs a checked build
and a **targeted** `-O3` build of the two modified resampling modules;
unchanged FFT/spectral dependencies use checked compilation in both builds.
Existing independent resampling regressions also passed: 53 multidimensional,
184 complex-filter, and 738 rank-4-to-15 SciPy cases. The full FPM test suite
was not run.

**SciPy reference caveat:** SciPy 1.17 currently rejects
`window=('general_cosine', coefficients)` when provided directly to
`resample_poly` through its internal `firwin` window dispatcher. For this
family, the comparison harness independently constructs the symmetric
window with SciPy `windows.general_cosine`, multiplies SciPy's unwindowed
`sinc` FIR prototype, normalizes the coefficients, then passes these explicit
coefficients to `resample_poly`. No SciPy code is required at runtime by
SciFort.
