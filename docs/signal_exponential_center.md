# Exponential and Poisson windows with an arbitrary center

The existing `get_window('exponential', n, parameter=tau)` and
`get_window('poisson', n, parameter=tau)` remain backward-compatible:
`tau` is the decay constant and the window is centered at the usual
symmetric position.

`get_window` now also accepts `parameter2=center` to choose an **arbitrary**
window center while retaining SciFort's first-parameter convention:

```fortran
use scifort_signal, only : dp, get_window
real(dp), allocatable :: w(:)
w = get_window('exponential', 8, fftbins=.true., &
               parameter=2.25_dp, parameter2=0.0_dp)
```

Here `w(k+1) = exp(-abs(k-center)/tau)` for `k=0,...,n-1`.
As in `scipy.signal.windows.exponential`, specifying an explicit center
is only valid for **periodic** windows (`fftbins=.true.`); it is rejected
for `fftbins=.false.`. Unlike SciPy's positional tuple
`('exponential', center, tau)`, SciFort deliberately preserves its
existing keyword parameter order: **`parameter=tau, parameter2=center`**.

The same behavior is supported by `scifort_signal_windows_extra.window_samples`
and `scifort_signal_spectral.spectral_named_window`, so other spectral
estimators can use the new window.

```bash
python tools/run_signal_exponential_center.py
```

This runs native assertions and compares 400 randomly generated periodic
windows, including zero/one-point lengths and negative/noninteger centers,
with SciPy in checked and optimized GNU Fortran configurations.
