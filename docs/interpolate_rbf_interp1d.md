# Scattered-data radial basis interpolation and legacy interp1d

This increment adds two **independently written**, MIT-licensed Fortran modules:

- `src/interpolate_rbf.f90` for SciPy-style `RBFInterpolator`.
- `src/interpolate_legacy.f90` for SciPy's legacy `interp1d` interface.

Both modules are exported by `scifort_interpolate`, in addition to being usable
directly. Neither adds a required FPM dependency or calls BLAS/LAPACK.
SciPy 1.17 was used as a **black-box numerical test reference**; its source
was not copied or translated.

## Radial basis interpolation

`rbf_interpolator(y, d, ...)` accepts `y(npoints,ndim)` and real or complex
`d(npoints)` or `d(npoints,nchannel)`. It returns either
`rbf_interpolator_type` or `complex_rbf_interpolator_type`, according to
`d`'s type. Each object supports `%valid()`, `%dimensions()`, `%channels()`,
`%value(q[,channel])`, and `%evaluate(q)` where `q(ndim)` yields one value per
channel, or `q(nquery,ndim)` yields `(nquery,nchannel)`.

Kernel names: `linear`, `thin_plate_spline`, `cubic`, `quintic`,
`multiquadric`, `inverse_multiquadric`, `inverse_quadratic`, and `gaussian`.
The module additionally implements `matern1_2`, `matern3_2`, `matern5_2`,
which are not supported by the *installed SciPy 1.17 reference package*;
those kernels have separate analytic native tests.

Optional construction arguments:

- `kernel`: default `thin_plate_spline`.
- `epsilon`: positive scale factor, required for kernels other than
  `linear`, `thin_plate_spline`, `cubic`, `quintic`.
- `degree`: polynomial augmentation degree `-1` (none) through 5; the
  default is the usual minimum (`linear`/`multiquadric`: 0,
  `thin_plate_spline`/`cubic`: 1, `quintic`: 2, others: 0).
- `smoothing`: common nonnegative regularization, default zero; OR
  `smoothing_values(npoints)` for heterogeneous regularization.
- `neighbors`: optional count of nearest points to refit for each query.
  Omitting this option fits a single global RBF. Local fits construct a
  reusable balanced KD-tree, pruning distant site groups rather than
  scanning all observations for each query. A local-subproblem coefficient
  cache is not implemented.
- `status`: `rbf_status_success=0`, `rbf_status_invalid_input=1`,
  or `rbf_status_singular=2`.

The polynomial-augmented dense system is solved with partial-pivoting Gaussian
elimination using no external libraries. The global path takes `O(n^3)`
work and `O(n^2)` storage. Local evaluation uses KD-tree neighbor selection
(typically sublinear in n in low dimensions, but linear worst case) plus a
roughly `O(k^3)` dense solve for `k` neighbors. This is intended for moderate dataset sizes.
Polynomial terms use centered/scaled coordinates to improve conditioning.
Near-singular kernel systems can still lose precision; callers should
consider smoothing, epsilon, or local fitting as appropriate.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: rbf_interpolator_type, rbf_interpolator

type(rbf_interpolator_type) :: rbf
real(dp) :: sites(4,2), observations(4), prediction

sites(1,:) = [0.0_dp, 0.0_dp]
sites(2,:) = [1.0_dp, 0.0_dp]
sites(3,:) = [0.0_dp, 1.0_dp]
sites(4,:) = [1.0_dp, 1.0_dp]
observations = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
rbf = rbf_interpolator(sites, observations, kernel='thin_plate_spline')
prediction = rbf%value([0.3_dp, 0.7_dp])
```

## interp1d

`interp1d(x, y, ...)` accepts real abscissas `x(npoints)` and real or
complex ordinates `y(npoints)` or `y(npoints,nchannel)`. Multichannel data
therefore use **axis 0** (not the Python `interp1d` default last axis).
The return type is `interp1d_type` or `complex_interp1d_type`.

Supported `kind` strings:

`linear`, `nearest`, `nearest-up`, `previous`, `next`, `zero`,
`slinear`, `quadratic`, `cubic`.

The `slinear`/`quadratic`/`cubic` modes use SciFort's existing native
B-spline collocation. The other modes use direct interval searches; at
nearest midpoints `nearest` rounds toward the lower knot and `nearest-up`
toward the upper knot.

Options include `bounds_error` (true by default without extrapolation),
`fill_value` (default NaN), independent `fill_below`/`fill_above`,
`extrapolate`, `assume_sorted` (false by default), and `status`.
Duplicate knots are rejected. `copy=False`, arbitrary Python `axis`,
array-valued fill broadcasting, and numeric `kind` selectors are not yet
implemented. The API offers `%valid()`, `%channels()`,
`%value(x[,channel])`, and `%evaluate(x)` for scalar or vector queries.

```fortran
use scifort_kinds, only: dp
use scifort_interpolate, only: interp1d_type, interp1d

type(interp1d_type) :: f
real(dp) :: prediction

f = interp1d([0.0_dp, 1.0_dp, 2.0_dp], [0.0_dp, 1.0_dp, 4.0_dp], &
             kind='quadratic', extrapolate=.true.)
prediction = f%value(1.5_dp)
```

## Verification

Both comparison drivers compile checked (`-fcheck=all`) and optimized
(`-O3`) GNU Fortran executables and compare against SciPy/NumPy. Two
independent seeds were run in both compiler configurations:

- `compare_rbf_scipy.py`: **220 cases per seed/build**, 3,850 complex
  output entries per run; all eight SciPy-1.17 kernels, up to 4 coordinates,
  global/local fitting, smoothing, real/complex channels.
- `compare_interp1d_scipy.py`: **360 cases per seed/build**, 34,862 complex
  output entries per run; all nine modes, unsorted knots, bounds and
  extrapolation, real/complex channels, exact nodes and midpoint ties.
- `test/test_interpolate_rbf.f90` verifies polynomial reproduction, local
  and global fits, analytic Matérn values, error handling, and complex data.
- `test/test_interpolate_legacy.f90` verifies knots, tie conventions, fills,
  extrapolation, bounds errors, and real/complex public APIs.

Reproduce with:

```bash
python tools/interpolate_reference/compare_rbf_scipy.py
python tools/interpolate_reference/compare_rbf_scipy.py --seed 8675309
python tools/interpolate_reference/compare_interp1d_scipy.py
python tools/interpolate_reference/compare_interp1d_scipy.py --seed 8675309
python tools/interpolate_reference/run_all_native_interpolation.py
```

The full FPM test suite was **not** run locally because FPM is unavailable.
