# Native periodic scalar spline validation (2026-11-13)

The native `splrep(per=.true.)` path supports interpolation and automatic-knot
smoothing for degrees 1-5. The last ordinate and weight are ignored. The
Fortran signature adds optional `per` after the existing `status` argument.
Periodic fixed knots and continuation remain unsupported; the TCK evaluator
retains its existing extrapolation and clipped-integration conventions.

## Sources and method

The adapted kernels and retained BSD-3-Clause license are recorded in
[provenance](../CODE_PROVENANCE.md). The reference is an independently
compiled, unchanged perazz/FITPACK core at commit
`1907c7bad9f5c4557970e906bdd4909eb52b2d4c`. No downloaded source or compiled
artifact is committed. `splrep_periodic_upstream_driver.f90` calls its
`percur`, `splev`, `splder`, and `splint` directly.

`check_splrep_periodic.py` uses deterministic random seed 20261114: 180 fits
with nonuniform sites, degrees 1-5, varying weights, deliberately unequal
closing ordinates, interpolation, small positive targets, and constant limits.
One additional cubic interpolation fixture uses 200 samples of
`cos(x)+sin(2*x)` on `[0,2*pi]`, reproducing the geometry and `.01` relative
plus `.0001` absolute analytic tolerance of upstream `test_periodic_fit`.
Only that fresh-fit geometry is reproduced; the complete upstream test suite
and its extrapolating periodic object interface were not run.

Each fit checks the residual bound, periodic endpoint derivatives, and
evaluation of native tuples with independent SciPy `BSpline`. Positive
adaptive fits also check the cyclic constrained normal equations using an
independently constructed basis and kth-derivative jumps. The checker
estimates a nonnegative penalty multiplier and verifies scaled stationarity
to `2e-8`. When the reference has the same knots, the squared weighted
prediction difference is bounded by the absolute residual difference plus
`1e-10*max(1,sum((w*y)**2))` roundoff allowance.

## Results

GNU Fortran (GCC) 17.0.0 20260510 (experimental), Windows; Python SciPy 1.15.3.
Checked `-O0` and optimized `-O3` builds both use `-g -fbacktrace -fcheck=all
-ffpe-trap=invalid,zero,overflow -Werror=implicit-interface`.

- Native periodic regression: 115 checks passed.
- Full `fpm test`: passed, exit code 0. A checked smoke test through the
  public `scifort_interpolate` module also passed.
- Existing scalar smoothing (108 checks), parametric smoothing (125 checks),
  parametric interpolation, TCK calculus, and fixed-knot spline regressions pass.
- Original pinned core: 181 fits, 120 adaptive; 45,738 coefficient, value,
  derivative, and integral comparisons. Maximum normalized difference
  `9.91144e-11`, below the `3e-8` comparison limit.
- Maximum relative adaptive residual error: `0.000852808`, below `0.0011`.
- Maximum independent scaled stationarity error: `7.27735e-14`.
- SciPy diagnostic comparison: 181 fits, 45,449 comparable numbers; 13 fits
  have different knot sets and 78 compared arrays exceed a `3e-8` normalized
  difference. Maximum normalized difference is `60.9693`. These diagnostics
  are recorded, not asserted to be close elementwise comparisons.

The SciPy comparisons assert fitting contracts, seam continuity, independent
representation evaluation, normal equations, and the same-knot prediction
bound. Tight elementwise agreement is asserted against the pinned original
core. Exact SciPy knots, coefficients, and higher derivatives are not promised:
adaptive knot selection and finite-tolerance smoothing can follow different
paths, particularly near interpolation. Residual agreement provides no
general higher-derivative error bound. The numerical contract and these
compatibility limitations are documented in the [TCK guide](../docs/interpolate_fitpack_tck_backend.md).

Native tests also cover zero data, invalid active inputs, ignored nonfinite
closing entries, negative/NaN targets, minimal sample counts, seam continuity,
tuple padding, calculus identities, default/nonperiodic behavior, large common
coordinate offsets, large/small scales, subnormal weights, and explicit
numerical failure for unrepresentable normalized targets.

## Reproduction

Compile the native driver with the current core and compatibility modules and
the remaining SciFort library. Compile the upstream driver against the exact
original core separately, without linking the adapted native core. Then run:

```text
python validation/check_splrep_periodic.py build/periodic-debug/driver.exe
python validation/check_splrep_periodic.py build/periodic-debug/driver.exe --upstream-executable build/periodic-debug/upstream.exe
python validation/check_splrep_periodic.py build/periodic-release/driver.exe
python validation/check_splrep_periodic.py build/periodic-release/driver.exe --upstream-executable build/periodic-debug/upstream.exe
fpm test
```

The standalone native programs were compiled with the flags above, linking
fresh `src/interpolate_curfit_core.f90`, `src/interpolate_fitpack_compat.f90`,
and `src/interpolate.f90` objects with the existing archive for unchanged
support modules. The full FPM run checks the ordinary project integration.
