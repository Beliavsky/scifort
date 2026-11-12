# Continuous-order band-stop filter objective

SciFort's `scifort_signal` now exports `band_stop_obj`, a Fortran port of the
standalone numerical function `scipy.signal.band_stop_obj`. The function
supports the original Butterworth (`'butter'`), Chebyshev (`'cheby'`), and
elliptic (`'ellip'`) analog band-stop prototype families. It computes the
**continuous** (non-rounded) filter order as one candidate passband edge
moves; this is an objective calculation, **not** a filter design.

```fortran
use scifort_signal, only : dp, band_stop_obj, band_stop_success
real(dp) :: order
integer :: status
order = band_stop_obj(0.5_dp, 0, [0.5_dp,4.0_dp], &
                      [1.5_dp,2.5_dp], 1.0_dp, 40.0_dp, 'ellip', status)
if (status /= band_stop_success) error stop 'invalid design'
```

Inputs are angular frequencies in rad/s. `ind` is deliberately **0-based**
(`0` or `1`) for parity with SciPy; the other parameters are a candidate
passband edge `wp`, original two passband edges `passb`, fixed stopband
edges `stopb`, positive `gpass` and `gstop` in dB, and family name. The
function returns a quiet NaN and optional
`band_stop_invalid_input` status when inputs or results are invalid.
The Fortran API requires strictly positive ordered passband edges and
finite inputs. A Butterworth objective can be negative for a nonphysical
trial with `nat < 1`, matching SciPy; Chebyshev and elliptic objectives
reject such trials rather than returning complex results.

Elliptic integrals use SciFort's existing independent arithmetic-geometric
mean implementation. The complementary integrals are evaluated without
subtracting nearly equal floating-point numbers. Consequently the Fortran
answer can be *more accurate* than SciPy's unmodified `band_stop_obj`
formula at extreme stopband attenuations; tests also use the independent
`scipy.special.ellipkm1` kernel as a stable high-attenuation oracle.

## Reproducible checks

```bash
python tools/run_signal_band_stop_obj.py
```

This checks both checked `-fcheck=all` and optimized `-O3` builds,
1,102 independently generated numerical reference cases, including
922 direct `scipy.signal.band_stop_obj` comparisons, and invalid-input
checks. Requires GNU Fortran, Python, NumPy and SciPy; no FPM or system
BLAS/LAPACK library is needed for this focused test.
