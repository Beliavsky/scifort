# Rectangular-grid spherical interpolation

`rect_sphere_bivariate_spline` is a native binary64 counterpart to SciPy's
`RectSphereBivariateSpline` for interpolation (`s=0`). It uses FITPACK's
periodic bicubic grid solver without an optional backend or new dependency.
Adaptive smoothing remains unsupported.

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate, only : rect_sphere_bivariate_spline, &
        sphere_bivariate_spline_type
    implicit none
    type(sphere_bivariate_spline_type) :: fit
    real(dp) :: u(6), v(8), r(6,8), t
    integer :: i, j, status
    do i = 1, size(u)
        u(i) = scifort_pi*real(i,dp)/7.0_dp
        t = u(i)/scifort_pi
        do j = 1, size(v)
            r(i,j) = 1.0_dp+2.0_dp*(3.0_dp*t*t-2.0_dp*t*t*t)
        end do
    end do
    do j = 1, size(v)
        v(j) = 2.0_dp*scifort_pi*real(j-1,dp)/real(size(v),dp)
    end do
    fit = rect_sphere_bivariate_spline(u,v,r, &
        pole_continuity=[.true.,.true.], pole_values=[1.0_dp,3.0_dp], &
        pole_flat=[.true.,.true.],status=status)
    if (fit%valid()) print *, fit%value(scifort_pi/2.0_dp,0.7_dp), status
end program example
```

The example represents a cubic latitude function exactly; its output is 2
to rounding accuracy. The constructor and type are exported by both
`scifort_interpolate` and `scifort_interpolate_sphere`.

## Inputs and pole options

```text
fit = rect_sphere_bivariate_spline(u,v,r,s,pole_continuity,pole_values, &
                                  pole_exact,pole_flat,status)
```

- `u` contains at least four strictly increasing colatitudes in `(0,pi)`.
- `v` contains at least four strictly increasing longitudes, with
  `-pi <= v(1) < pi` and `v(last) < v(1)+2*pi`. Omit the duplicate seam
  endpoint. All coordinates are finite and measured in radians.
- `r` is finite, with shape `(size(u),size(v))`.
- Optional `s` defaults to zero. Positive values return
  `sphere_spline_unsupported_smoothing` (5); negative or nonfinite values
  return `sphere_spline_invalid_input` (1).
- The optional pole arrays each have length two, ordered north then south.
  Unlike SciPy's scalar-or-pair options, these are always arrays.
  `pole_continuity` defaults to false and requests first-derivative continuity
  at the corresponding pole. `pole_flat` defaults to false and requests zero
  slope there; a true flat flag requires the corresponding continuity flag.
- `pole_values` supplies finite pole observations. Omission leaves both
  unspecified; a quiet NaN entry leaves just that pole unspecified. Infinity
  is invalid. Missing values and required pole derivatives are estimated by
  FITPACK's preliminary coarse-grid least-squares procedure.
- `pole_exact` defaults to false. With `s=0`, supplied pole values are always
  interpolated exactly regardless of this flag, following FITPACK. The flag
  has no effect at an unspecified pole.

## Results and limits

The result is the shared `sphere_bivariate_spline_type`. Its scalar `%value`,
paired `%values`, Cartesian `%grid`, knot, and coefficient methods follow
the [spherical spline evaluation API](interpolate_sphere.md#evaluation-and-metadata).
Partial derivative orders are independently 0, 1, or 2. Query colatitudes
must be in `[0,pi]`; finite longitudes wrap into the period beginning at
`v(1)`. Automatic longitude wrapping is a **SciFort extension**.

`%residual()` is zero for a successful interpolation. `%rank()` returns -1:
the banded solver does not estimate numerical rank. `%degrees_of_freedom()`
returns `size(u)*size(v)+2+2*count(pole_continuity)`, the constrained tensor
space dimension before supplied/estimated pole values and flatness are fixed.
It is neither an estimated rank nor statistical residual degrees of freedom.
These two diagnostic methods are **SciFort extensions**.

Status 0 means success; 1 means invalid input; 3 means numerical/range failure;
4 means the workspace guard was exceeded; and 5 means unsupported smoothing.
Failures return an invalid object. The LSQ-only rank-warning status 2 is not
used by this constructor.

The distance of each extreme colatitude from its pole, all adjacent coordinate
gaps, and the longitude seam gap must exceed `128*epsilon(dp)*2*pi` radians
(about `1.79e-13`). Smaller separations return status 3. The solver also caps
the aggregate of its scratch workspace, tensor coefficient storage, and input
copy at 2,000,000 entries. This is not a total process-memory limit. Unlike
the scattered-data SVD interface, it has no 256-coefficient cap.
Data and pole values are scaled internally; unrepresentable physical
coefficients are rejected. These guards do not guarantee accuracy for every
ill-conditioned grid. The current minimum is four sites on each axis even
where SciPy accepts fewer latitude sites with particular pole constraints.

One observed compatibility difference concerns unspecified pole estimation.
For a constant grid with a different supplied opposite pole, SciPy 1.15.3's
legacy initialization can leave the missing pole at its initial estimate.
SciFort includes both pole values in the optimization step range, as in
modern FITPACK. The resulting estimate agrees with an independent constrained
least-squares calculation. Both implementations still interpolate the input
grid and the supplied pole. See the [validation record](../validation/rect_sphere_2026-11-13.md).

## Provenance

The grid driver and banded solver are adapted from Federico Perini's
[perazz/FITPACK](https://github.com/perazz/fitpack), based on Paul Dierckx's
FITPACK. BSD-3-Clause notices and the full
[license](../LICENSES/perazz-fitpack-BSD-3-Clause.txt) are retained.
The exact source pin and significant changes are recorded in
[CODE_PROVENANCE.md](../CODE_PROVENANCE.md). No C ABI is added.
