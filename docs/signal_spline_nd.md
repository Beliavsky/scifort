# Rank-independent spline interpolation (ranks 1-15)

The `scifort_signal` module provides **Fortran 2018 assumed-rank** spline
filtering and coordinate interpolation for `real(dp)` and `complex(dp)` arrays.
The implementation is independent of the older rank-1 `cspline1d` and rank-2
`cspline2d`/`qspline2d` APIs; those interfaces and their smoothing semantics
are preserved.

| API | Operation |
|---|---|
| `spline_filter_nd(x, coeff, order, status, mode)` | Cardinal B-spline coefficients along every axis |
| `spline_filter1d(x, coeff, axis, order, status, mode)` | Coefficients along one **one-based** axis |
| `cspline_nd(x, coeff, status)` | All-axis order-3 spline coefficients |
| `qspline_nd(x, coeff, status)` | All-axis order-2 spline coefficients |
| `map_coordinates_nd(coeff, coords, values, order, status, mode, prefilter, cval)` | Interpolate at arbitrary points |
| `map_coordinates_grid_nd(x, coords, y, order, status, mode, prefilter, cval)` | Coordinate grids whose trailing extents are the output shape; defaults to prefiltering |
| `geometric_transform_nd(x, mapping, y, order, status, mode, prefilter, cval)` | Callback mapping from output to input coordinates; output rank may differ |
| `prepare_spline_plan_nd(x, plan, order, status, mode, prefilter, cval)` | Cache coefficients for repeated coordinate sampling (SciFort extension) |
| `evaluate_spline_plan_nd(plan, coords, values, status)` | Sample a cached plan at flat one-based coordinates |
| `evaluate_spline_plan_grid_nd(plan, coords, y, status)` | Sample a cached plan into a shaped destination |
| `evaluate_spline_plan_with_gradient_nd(plan, coords, values, gradient, status)` | Analytic coordinate derivatives and interpolated values in one pass |
| `evaluate_spline_plan_grid_with_gradient_nd(plan, coords, y, gradient, status)` | Real/complex values and shaped derivatives on a coordinate grid (output ranks 1-14) |
| `evaluate_spline_plan_with_hessian_nd(plan, coords, values, gradient, hessian, status)` | Analytic gradient and full symmetric coordinate Hessian; flattened `(source_rank, source_rank, points)` |
| `evaluate_spline_plan_grid_with_hessian_nd(plan, coords, y, gradient, hessian, status)` | Grid-shaped values, gradients, and Hessians (output ranks 1-13) |
| `evaluate_spline_plan_displacement_grid_with_gradient_nd(plan, field, y, gradient, status, matrix, offset, batch_size)` | Grid-shaped displacement warp and per-component sensitivity (output ranks 1-14) |
| `evaluate_spline_plan_displacement_with_hessian_nd(plan, displacements, y, gradient, hessian, status, matrix, offset, batch_size)` | Flat-field displacement warp and first/second displacement sensitivities (output ranks 1-15) |
| `evaluate_spline_plan_displacement_grid_with_hessian_nd(plan, field, y, gradient, hessian, status, matrix, offset, batch_size)` | Grid-shaped displacement Hessians with two leading source axes (output ranks 1-13) |
| `evaluate_spline_plan_affine_with_hessian_nd(plan, matrix, offset, y, dmatrix, doffset, hmatrix, hcross, hoffset, status, batch_size)` | Analytic matrix-matrix, matrix-offset and offset-offset Hessians of a cached affine transform |
| `evaluate_spline_plan_affine_with_gradient_nd(plan, matrix, offset, y, dmatrix, doffset, status, batch_size)` | Batched affine warp and derivatives with respect to its affine parameters |
| `evaluate_spline_plan_tensor_nd(plan, axis_grid, y, status)` | Evaluate a Cartesian product of axis coordinates via separable 1-D contractions |
| `evaluate_spline_plan_affine_nd(plan, matrix, offset, y, status, batch_size)` | Batched nonseparable affine transform of cached coefficients |
| `evaluate_spline_plan_geometric_nd(plan, mapping, y, status, batch_size)` | Point-by-point callback transform of cached coefficients |
| `evaluate_spline_plan_geometric_batch_nd(plan, mapping, y, status, batch_size)` | **SciFort extension:** vectorized coordinate callback, called once per block |
| `evaluate_spline_plan_displacement_with_gradient_nd(plan, displacements, y, gradient, status, matrix, offset, batch_size)` | Warped image and its derivative with respect to each displacement component |

The arguments `x`, `coeff` and the destination of `spline_filter1d` are
assumed-rank. The caller allocates the output with the **same rank and shape**
as the input. The numerical implementation flattens data in Fortran array
order and processes strided lines along each axis. It uses a narrowly banded
spline collocation solver (and an O(n) periodic IIR spline prefilter),
without implementing 15 different algorithms.
Input and output may be noncontiguous array sections.

```fortran
use scifort_signal, only: dp, spline_filter_nd, map_coordinates_nd, &
    spline_status_success
real(dp) :: samples(10,12,8), coeff(10,12,8)
real(dp) :: coordinates(3,2), values(2)
integer :: status

! Initialize samples first, then calculate cubic coefficients.
call spline_filter_nd(samples, coeff, order=3, status=status)
if (status /= spline_status_success) error stop 'spline filter failure'

coordinates(:,1) = [2.0_dp, 3.25_dp, 7.0_dp]
coordinates(:,2) = [11.5_dp, 5.0_dp, 2.75_dp]
call map_coordinates_nd(coeff, coordinates, values, order=3, &
    mode='mirror', status=status)
if (status /= spline_status_success) error stop 'interpolation failure'
```

## Affine-parameter Hessians and second-order displacement sensitivities (2026-10-24)

`evaluate_spline_plan_affine_with_hessian_nd` evaluates a cached affine
warp together with analytic first and second derivatives of every pixel
with respect to its affine matrix and offset. All affine coordinates follow
`source = 1 + offset + matrix * (output - 1)`. Its outputs are:

- `y` has the caller-supplied output shape and may have ranks 1-15.
- `dmatrix(input_rank,output_rank,npixels)` and
  `doffset(input_rank,npixels)` are the existing first derivatives.
- `hmatrix(input_rank,output_rank,input_rank,output_rank,npixels)` gives
  matrix-matrix second derivatives. For example, `hmatrix(a,j,b,k,p)`
  differentiates pixel `p` with respect to `matrix(a,j)` and `matrix(b,k)`.
- `hcross(input_rank,output_rank,input_rank,npixels)` gives
  matrix-offset second derivatives, and
  `hoffset(input_rank,input_rank,npixels)` gives offset-offset derivatives.

The Hessians follow directly from the coordinate Hessian `H` at each output
site: `hoffset(a,b,p)=H(a,b,p)`,
`hcross(a,j,b,p)=(output_j-1)*H(a,b,p)`, and
`hmatrix(a,j,b,k,p)=(output_j-1)*(output_k-1)*H(a,b,p)`.
These are derivatives of the interpolation result with respect to parameters;
they are **not** image-space resampling Jacobians. Computation takes place in
bounded batches without finite differencing. Because the complete Hessian
contains quadratic numbers of affine parameters, the existing
8-million-element guard is scaled by `input_rank**2*(output_rank+1)**2`.

`evaluate_spline_plan_displacement_with_hessian_nd` accepts a flattened
`displacements(input_rank,npixels)` and returns the corresponding warp,
`gradient(input_rank,npixels)`, and
`hessian(input_rank,input_rank,npixels)`. Optional `matrix` and `offset`
provide the affine base map. With no affine base, source and destination
ranks must agree. `evaluate_spline_plan_displacement_grid_with_hessian_nd`
accepts `field(input_rank,shape(y))`, returning gradient and Hessian tensors
with one and two leading source-rank axes respectively. Its maximum output
rank is 13; the flattened interface supports rank 15. Both real and complex
images, noncontiguous output sections, spline orders 0-5, all eight boundary
modes, and configurable batch sizes are supported.

```fortran
real(dp) :: matrix(2,2), offset(2), y(6,5)
real(dp) :: dm(2,2,30), db(2,30), hm(2,2,2,2,30)
real(dp) :: hc(2,2,2,30), hb(2,2,30)
type(spline_interpolation_plan) :: plan
integer :: status
! Prepare plan from a real image and initialize matrix/offset.
call evaluate_spline_plan_affine_with_hessian_nd(plan,matrix,offset, &
    y,dm,db,hm,hc,hb,status)
```

The test `test/test_signal_spline_parameter_hessians.f90` verifies the chain
rule, rank-changing affine maps, rank-15 inputs, strided grid output, and
parameter finite differences across all eight boundary modes and orders
0-5. `tools/signal_reference/compare_spline_parameter_hessians_scipy.py`
compares the result with independently differentiated SciPy interpolation.
At nondifferentiable spline knots the Hessian follows SciFort's symmetric
knot convention, and short-array reflect initialization may differ from
SciPy as previously documented.

## Analytic spline Hessians and grid displacement sensitivities (2026-10-24)

`evaluate_spline_plan_with_hessian_nd` evaluates each spline interpolant, its
coordinate gradient, and its complete **symmetric Hessian** in one traversal
of the spline support. `coords` has shape `(source_rank, npoints)`; the
caller supplies `values(npoints)`, `gradient(source_rank,npoints)`, and
`hessian(source_rank,source_rank,npoints)`. Both `real(dp)` and `complex(dp)`
plans work, with spline orders 0-5, all eight boundary modes, noncontiguous
array sections, and source ranks 1-15.

For orders 2-5, the implementation uses the analytic second derivative of a
centered cardinal B-spline, along with products of derivatives on distinct
axes for the mixed Hessian entries. First-order *same-axis* second
derivatives vanish inside spline pieces, but **mixed derivatives need not
vanish** for tensor-product bilinear or multilinear interpolation. Order-zero
Hessians vanish away from jumps. At nondifferentiable knots and boundary
reflection corners, classical derivatives may be undefined; as in the
existing gradient API, numerical conventions are applied at exact knots.
Hessians differentiate cached spline coefficients with respect to *query
coordinates*, not the original samples or the prefiltering operation.

```fortran
real(dp) :: image(20,18), coords(2,6), values(6), gradient(2,6)
real(dp) :: hessian(2,2,6)
type(spline_interpolation_plan) :: plan
integer :: status
! Initialize image and coords.
call prepare_spline_plan_nd(image,plan,order=3,mode='mirror',status=status)
if (status/=0) error stop 'plan preparation failed'
call evaluate_spline_plan_with_hessian_nd(plan,coords,values,gradient,hessian,status)
if (status/=0) error stop 'Hessian evaluation failed'
```

For a shaped destination `y(n,m)`, the grid API uses
`coords(source_rank,n,m)`, `gradient(source_rank,n,m)`, and
`hessian(source_rank,source_rank,n,m)`. The two additional leading
Hessian dimensions limit grid-shaped outputs to rank 13, even though the
flattened Hessian API supports source ranks up to 15. Both APIs enforce
the existing eight-million-element workspace guard, scaled by the
number of stored Hessian components.

`evaluate_spline_plan_displacement_grid_with_gradient_nd` accepts a dense
field `field(source_rank,n,m)` and returns `y(n,m)` and
`gradient(source_rank,n,m)`. The latter is the derivative of each sampled
value with respect to its own displacement components, **not** the spatial
gradient of the displacement field or a derivative with respect to image
samples. Optional `matrix`, `offset`, and `batch_size` use the conventions
of the flattened displacement interface. The grid-shaped form supports
output ranks 1-14, while the original flattened API remains available for
rank-15 destinations. Both real and complex images are supported.

Validation is provided by `test/test_signal_spline_hessians.f90` and
`tools/signal_reference/compare_spline_hessians_scipy.py`, which compares
values, gradients, and second derivatives to SciPy's independent
`ndimage.map_coordinates` finite differences. Short-array `reflect`
prefilter differences previously documented in this guide remain unchanged.

## Coordinate conventions

`coords` has shape `(rank(coeff), npoints)` and `values` is length `npoints`.
Unlike SciPy's `ndimage.map_coordinates`, **SciFort coordinates are one-based**:
`[1, 1, 1]` identifies the first element of a rank-3 array. Axis arguments
are also one-based. The default mode is `mirror` (whole-sample symmetry),
matching the old SciFort spline evaluators. `mode='constant'` returns the
supplied `cval` for any coordinate outside the data bounds. Additional modes
`reflect`, `grid-mirror`, `nearest`, `wrap`, `grid-wrap`, and `grid-constant`
are now supported. See [Boundary modes and geometric transforms](signal_spline_transforms.md)
for precise semantics, prefilter limits, and affine, shift, zoom and rotation APIs.

The `mode` argument to `spline_filter_nd` and `spline_filter1d` is
optional and defaults to `mirror`. It follows SciPy's *filter* conventions:
`reflect`, `grid-mirror`, and `nearest` use half-sample symmetry;
`mirror`, `constant`, `grid-constant`, and `wrap` use whole-sample symmetry;
`grid-wrap` uses periodic spline filtering in O(n) time per line. Unlike
SciPy's recursive `reflect` initializer on short lines, SciFort solves the
finite collocation equations exactly. See the boundary-mode notes below.

`map_coordinates_nd` takes *prefiltered coefficients* by default
(`prefilter=.false.`). Pass `prefilter=.true.` to filter an ordinary sample
array before interpolation. This differs from SciPy's default `prefilter=True`.
The default interpolation order is 3; valid orders are 0 through 5.
Filtering for orders 0 and 1 is a copy operation (SciPy's separate
`ndimage.spline_filter` API itself accepts only orders 2–5).

The coefficient routines implement ordinary B-spline **interpolation**, not
cubic smoothing (`lamb > 0`). The older `cspline1d` and `cspline2d` interfaces
continue to provide their respective smoothing options.

## Shape-preserving coordinate grids

For a SciPy-style output shape, call `map_coordinates_grid_nd`. Its coordinate
array has leading extent `rank(input)` and trailing extents equal to the
caller-allocated output shape. All coordinates remain **one-based**. This new
interface defaults to `prefilter=.true.` (unlike the older flat-coordinate
interface), and handles real/complex input arrays, noncontiguous sections,
and source ranks 1-15. The output rank is 1-14 because standard Fortran has
a maximum array rank of 15 including the leading coordinate dimension.

```fortran
real(dp) :: input(40,30), coordinates(2,12,9), result(12,9)
integer :: status
! Initialize input and coordinates; coordinates(:,i,j) supplies a 2-D position.
call map_coordinates_grid_nd(input,coordinates,result,order=3, &
    mode='grid-wrap',status=status)
if (status /= spline_status_success) error stop 'grid mapping failed'
```

The new test `test/test_signal_spline_periodic_grid.f90` exercises rank-15
inputs, noncontiguous output slices, invalid shapes and mode values, and
periodic knot reconstruction. The reference script
`tools/signal_reference/compare_spline_periodic_grid_scipy.py` compares
periodic prefilter coefficients and shaped coordinate outputs with SciPy 1.17.

## Analytic spline-coordinate derivatives (2026-10-24)

`evaluate_spline_plan_with_gradient_nd` returns interpolated values and the
analytic first derivative with respect to **each one-based source coordinate**.
For an input rank `r` and `n` sample locations, `coords` and `gradient` have
shape `(r,n)` and `values` has length `n`. Real and complex plans are supported;
`gradient` has the same numeric type as the cached image. The routine uses
B-spline derivative weights and computes values and derivatives together,
without `2*r` additional interpolation passes.

```fortran
use scifort_signal, only : dp, spline_interpolation_plan, &
    prepare_spline_plan_nd, evaluate_spline_plan_with_gradient_nd
real(dp) :: image(40,30), coords(2,3), values(3), gradient(2,3)
type(spline_interpolation_plan) :: plan
integer :: status
! Initialize image and coords before calling these routines.
call prepare_spline_plan_nd(image,plan,order=3,mode='mirror',status=status)
if (status/=0) error stop 'spline preparation failed'
call evaluate_spline_plan_with_gradient_nd(plan,coords,values,gradient,status)
if (status/=0) error stop 'spline gradient failed'
! gradient(1,:) and gradient(2,:) are derivatives along image axes 1 and 2.
```

### Grid-shaped gradients and affine-parameter sensitivities (2026-10-24)

The original `evaluate_spline_plan_with_gradient_nd` takes flat coordinates
of shape `(source_rank, npoints)`. For multidimensional destinations, use
`evaluate_spline_plan_grid_with_gradient_nd` with coordinates and gradient
of shape `[source_rank, shape(y)]`. The routine supports real and complex
plans and noncontiguous sections. Standard Fortran limits `rank(y)` to 14
because the component axis adds one dimension. All source ranks 1-15 remain
supported; a rank-15 destination can use the existing flat-coordinate API.

For repeated affine transformations, `evaluate_spline_plan_affine_with_gradient_nd`
returns the transformed array and **parameter derivatives** in one batched
spline evaluation. The affine mapping for an output coordinate `p`, converted
to one-based Fortran convention, is `q = 1 + offset + matrix * (p - 1)`.
The output derivatives satisfy:

- `doffset(i,k) = dI(q_k) / dq_i`;
- `dmatrix(i,j,k) = (p_j - 1) * doffset(i,k)`.

Here `k` is the output pixel number in **Fortran column-major order**.
The caller allocates `dmatrix(source_rank, output_rank, size(y))` and
`doffset(source_rank, size(y))`, even when `y` is multidimensional.
The algorithm calculates coordinate gradients per batch (1-8192 points)
and applies the chain rule. Both input and output ranks may reach 15,
and the input and output ranks may differ. Output shapes and work limits
are checked; invalid calls return a nonzero status.

```fortran
real(dp) :: image(40,30), matrix(2,2), offset(2), y(12,9)
real(dp) :: coordinates(2,12,9), sampled(12,9), grad_grid(2,12,9)
real(dp) :: dmatrix(2,2,12*9), doffset(2,12*9)
type(spline_interpolation_plan) :: plan
integer :: status
! Initialize image, matrix, offset and coordinates.
call prepare_spline_plan_nd(image,plan,order=3,mode='mirror',status=status)
if (status/=0) error stop 'preparation failed'
call evaluate_spline_plan_grid_with_gradient_nd(plan,coordinates,sampled,grad_grid,status)
if (status/=0) error stop 'grid gradient failed'
call evaluate_spline_plan_affine_with_gradient_nd(plan,matrix,offset,y,dmatrix,doffset,status)
if (status/=0) error stop 'affine sensitivity failed'
```

These routines differentiate the **interpolant**, not the original sampled
image or the plan-preparation step. As for other spline gradients, order-zero
and first-order knots are not differentiable in the classical sense.
The constant-fill mode returns zero derivatives outside the input domain.
For validation, use `test/test_signal_spline_grid_affine_gradients.f90` and
`tools/signal_reference/compare_spline_grid_affine_gradients_scipy.py`.
The latter uses SciPy's `map_coordinates` and independent central differences
of each matrix and offset parameter.

`evaluate_spline_plan_displacement_with_gradient_nd` returns the warped image
plus the derivative at each output pixel with respect to the **displacement
vector at that pixel**. The caller supplies a flat displacement array with
shape `(source_rank, size(output))`, following Fortran column-major order.
The gradient also has shape `(source_rank, size(output))`, avoiding a rank-16
result for rank-15 image arrays. Optional `matrix`, `offset`, and `batch_size`
follow `evaluate_spline_plan_displacement_nd`, with the same one-based source
coordinates and zero-based affine offset convention.

```fortran
real(dp) :: offsets(2,40*30), warped(40,30), sensitivity(2,40*30)
! Fill offsets, e.g. offsets = 0.0_dp for an identity warp.
call evaluate_spline_plan_displacement_with_gradient_nd( &
    plan,offsets,warped,sensitivity,status,batch_size=256)
if (status/=0) error stop 'warp and sensitivity failed'
```

These are derivatives of the *interpolant* with respect to coordinate offsets,
not derivatives with respect to image samples, and not the Jacobian of an
affine transformation or the gradient of the displacement field itself.
Spline order zero has zero derivative between nearest-neighbor jumps. At
knots where a first-order spline is nondifferentiable, the routine returns
the average of the one-sided slopes (the symmetric finite-difference limit).
The mathematical derivative at a corner is not unique; evaluate away from
those knots for derivative-sensitive applications. The derivative of constant-fill
interpolation is zero outside the original data domain. Both APIs preserve
all eight boundary modes, cached prefilter semantics, and source ranks 1-15.

Native tests are in `test_signal_spline_plan_gradients.f90` and
`test_signal_spline_displacement_gradients.f90`; independent central-difference
SciPy checks are in `compare_spline_plan_gradients_scipy.py`. A reproducible
speed comparison against finite differences is in
`benchmarks/bench_ndimage_spline_gradients.f90`.

## Reusable interpolation plans (2026-10-23)

The **SciFort extension** `spline_interpolation_plan` owns precomputed spline
coefficients, the source shape, boundary mode, order, and exterior `cval`. This
is useful for repeated interpolation of the **same source image** at many
different coordinates: `map_coordinates_nd(...,prefilter=.true.)` computes
and packs the source coefficients again on *every call*, whereas a plan does
that work **once**. This optimizes repeated interpolation without changing
the B-spline evaluation or boundary-mode mathematics.

```fortran
use scifort_kinds, only: dp
use scifort_signal_spline_nd, only: spline_interpolation_plan, &
    prepare_spline_plan_nd, evaluate_spline_plan_nd
real(dp) :: source(128,128), points_a(2,100), points_b(2,100)
real(dp) :: results_a(100), results_b(100)
type(spline_interpolation_plan) :: plan
integer :: status
! Initialize source, points_a, and points_b first.
call prepare_spline_plan_nd(source, plan, order=3, mode='reflect', status=status)
if (status /= 0) error stop 'preparation failed'
call evaluate_spline_plan_nd(plan, points_a, results_a, status)
if (status /= 0) error stop 'first query failed'
call evaluate_spline_plan_nd(plan, points_b, results_b, status)
if (status /= 0) error stop 'second query failed'
```

`prepare_spline_plan_nd` defaults to order 3, mode `mirror`, and
`prefilter=.true.`; it expects **ordinary samples**, unless
`prefilter=.false.` is specified to store caller-supplied spline coefficients.
Real and complex source arrays have separate overloads. The cached plan owns
its data and is unaffected by later changes to the source. Repreparing a plan
replaces its contents; a failed preparation leaves it invalid, not stale.

Flat queries use `coords(rank_source,npoints)` and one-based indexing.
`evaluate_spline_plan_grid_nd` also accepts a shaped coordinate tensor with
shape `(rank_source,shape(output))`, subject to the usual rank-15 Fortran
limit for the coordinate array. Real plans require real destinations and
complex plans require complex destinations; mismatches return invalid-input
status. The same resource guards, singleton-axis optimizations, and known
short-array `reflect` differences as `map_coordinates_nd` remain in force.
The plan is read-only during evaluations and may be reused for several calls.

`test/test_signal_spline_plan.f90` exercises repeated queries, original-source
independence, real/complex values, rank-15 and strided arrays, shaped grids,
and failed rebuilds. The new SciPy reference script
`tools/signal_reference/compare_spline_plan_scipy.py` covers 1,440
output values in multiple boundary modes and orders. The workload-specific
benchmark `benchmarks/bench_ndimage_spline_plan.f90` compares one preparation
and multiple cached queries to repeated direct interpolation; it reports
CPU time and verifies numerical equality.

## Fast tensor-product interpolation grids (2026-10-23)

For a Cartesian-product output grid, `evaluate_spline_plan_tensor_nd` avoids
forming the full `(rank, output_shape)` coordinate array and avoids evaluating
the complete multidimensional kernel at each pixel. Instead, it applies a
one-dimensional interpolation transform along each axis. This is especially
useful for large rectilinear resampling grids, where the same input axis
coordinates appear across many output pixels.

```fortran
use scifort_kinds, only : dp
use scifort_signal_spline_nd, only : spline_interpolation_plan, &
    spline_axis_grid, prepare_spline_plan_nd, evaluate_spline_plan_tensor_nd
real(dp) :: source(80,60), resampled(128,96)
type(spline_interpolation_plan) :: plan
type(spline_axis_grid) :: axis_grid(2)
integer :: i, status

! Initialize source first. These are one-based source coordinates.
axis_grid(1)%coordinates = [(1.0_dp+0.6_dp*real(i-1,dp),i=1,128)]
axis_grid(2)%coordinates = [(1.0_dp+0.6_dp*real(i-1,dp),i=1,96)]
call prepare_spline_plan_nd(source,plan,order=3,mode='reflect',status=status)
if (status/=0) error stop 'spline preparation failed'
call evaluate_spline_plan_tensor_nd(plan,axis_grid,resampled,status)
if (status/=0) error stop 'tensor evaluation failed'
```

`axis_grid(d)%coordinates` is an allocatable `real(dp)` vector of one-based
positions for output axis `d`. There must be exactly one vector per source
axis, and `shape(y)` must match their lengths. Thus input and output have the
same rank (up to 15), but their axis lengths may differ. Both real and complex
plans and noncontiguous destinations are supported. Invalid or nonfinite
coordinates, missing axes, or shape mismatches return invalid-input status.
The plan retains its boundary mode, degree, `prefilter` choice, and exterior
value. In particular, `constant` outside samples return the specified fill;
`grid-constant` and the special filtered singleton axes are handled directly.

For large grids the complexity is proportional to the number of values in
successive 1-D contractions times `order+1`, rather than the complete
`(order+1)**rank` product for every pixel. Peak tensor working buffers and
output are capped at 8,000,000 elements; larger queries return
`spline_nd_status_too_expensive`. Arbitrary, nonseparable coordinate mappings
should still use `evaluate_spline_plan_nd` or `evaluate_spline_plan_grid_nd`.

`test/test_signal_spline_tensor.f90` checks real/complex outputs, orders 0-5,
all eight boundary modes, short singleton axes, noncontiguous destinations,
rank-15 inputs and invalid arguments. The independent SciPy reference
`tools/signal_reference/compare_spline_tensor_scipy.py` checks 22,464
real/imaginary values; short-array `reflect` initializers retain the earlier
numerical exception. The benchmark
`benchmarks/bench_ndimage_spline_tensor.f90` compares the two methods on a
regular grid with the **same prefiltered plan** (so prefilter cost is excluded).

## Batched nonseparable transforms with cached coefficients (2026-10-23)

`evaluate_spline_plan_affine_nd` and `evaluate_spline_plan_geometric_nd` let
callers apply many nonseparable transforms of a fixed image without
recomputing the spline prefilter. They reuse **the same prepared plan** as
`evaluate_spline_plan_nd`, including its interpolation order, boundary mode,
exterior fill, and prefilter choice.

```fortran
use scifort_kinds, only : dp
use scifort_ndimage, only : spline_interpolation_plan, &
    prepare_spline_plan_nd, evaluate_spline_plan_affine_nd
real(dp) :: source(100,80), destination(70,60), transform(2,2), offset(2)
type(spline_interpolation_plan) :: plan
integer :: status

! Initialize source, transform, and offset before these calls.
call prepare_spline_plan_nd(source,plan,order=3,mode='reflect',status=status)
if (status/=0) error stop 'could not prepare source'
call evaluate_spline_plan_affine_nd(plan,transform,offset,destination,status)
if (status/=0) error stop 'could not evaluate affine transform'
```

The affine matrix has shape `(rank(source),rank(destination))`, so rectangular
matrices allow the output rank to differ from the input rank. Its elements
map *zero-based* destination indices to *zero-based* source positions:
`source_index = matrix * (destination_index - 1) + offset + 1` in one-based
Fortran coordinates. This convention matches SciPy's affine transform offset.

The geometric callback receives **one-based integer destination coordinates**
and writes **one-based real input coordinates**, using the existing
`ndimage_coordinate_map` interface. Callback calls occur once per output
pixel in Fortran array order, independent of the batch size. As with other
user-provided callbacks, it must not retain references to the temporary
coordinate arrays. The plan and its coefficients are read-only during a call.

Both overloads support `real(dp)` and `complex(dp)` plans, noncontiguous
output sections, and input and output ranks 1-15. The optional `batch_size`
defaults to 256 and can range from 1 through 8192. Coordinate workspace is
bounded by `rank(source)*batch_size`, not the output pixel count. The
allocated flat output remains proportional to the output size, with an
8,000,000-element resource guard; larger outputs return
`spline_nd_status_too_expensive`. Invalid shapes, mismatched output types,
nonfinite mappings, nonfinite matrix/offset arguments, and invalid batch
sizes return status 1. On an error, do not use the destination.

`test/test_signal_spline_plan_transforms.f90` checks all boundary modes,
orders 0-5, real/complex paths, arbitrary and rectangular mappings,
noncontiguous outputs, rank-15 sources, and failure handling. The SciPy
reference `tools/signal_reference/compare_spline_plan_transforms_scipy.py`
checks 19,008 scalar affine and callback values. The workload-specific
`benchmarks/bench_ndimage_spline_plan_transforms.f90` measures the benefit of
reusing a cached plan over repeated prefiltered affine transformations.

## Vectorized coordinate callbacks for cached spline plans (2026-10-24)

`evaluate_spline_plan_geometric_batch_nd` is a SciFort extension that
uses the same spline coefficient plan and interpolation kernel as
`evaluate_spline_plan_geometric_nd`, but maps up to `batch_size` destination
points with **one callback invocation**. This is useful when the mapping is a
vectorizable displacement field, nonlinear warp, or external Fortran model.
It is not an additional SciPy API.

The `ndimage_coordinate_batch_map` callback has the interface:

```fortran
subroutine mapping(output_coords, input_coords)
    use scifort_kinds, only : dp
    integer, intent(in) :: output_coords(:,:) !! (destination rank, batch count)
    real(dp), intent(out) :: input_coords(:,:) !! (source rank, batch count)
end subroutine mapping
```

For example, to map a rank-2 image through a nonlinear warp:

```fortran
use scifort_ndimage, only : spline_interpolation_plan, &
    prepare_spline_plan_nd, evaluate_spline_plan_geometric_batch_nd
! source(:,:) and output(:,:) should be allocated by the caller.
call prepare_spline_plan_nd(source, plan, order=3, mode='reflect', status=code)
call evaluate_spline_plan_geometric_batch_nd(plan, mapping, output, &
    status=code, batch_size=256)
```

Both arrays in the callback use **one-based coordinates**, with one point
per column. Columns follow Fortran array order; the last batch may have
fewer columns than `batch_size`. The callback must populate every element
of `input_coords` and must not retain a reference to temporary arrays.
The routine supports real and complex interpolation plans, ranks 1-15,
noncontiguous outputs, and rank-changing mappings. Existing scalar callbacks
continue to work without changes.

`batch_size` defaults to 256 and must be in `1:8192`. Coordinate storage
scales as `(source rank + destination rank) * batch_size`, independently
of output size; a final flat result buffer is still required. The output
is limited to 8,000,000 elements as in the earlier plan evaluators.
Invalid batch sizes, types, and nonfinite callback coordinates return status
1; the destination is unspecified on error.

The test `test/test_signal_spline_plan_batch_mapping.f90` compares the new
mapping API with pointwise interpolation over eight modes, orders 0-5,
complex inputs, rank changes and rank-15 arrays. The independent SciPy
reference is `tools/signal_reference/compare_spline_plan_batch_mapping_scipy.py`.
`benchmarks/bench_ndimage_spline_batch_mapping.f90` compares scalar and
vectorized callbacks with a shared prepared plan. Speedups are
workload-dependent; identical coefficient interpolation is used in both paths.

## Dense displacement-field spline warps (2026-10-24)

`evaluate_spline_plan_displacement_nd` samples a prepared real or complex
spline plan using a **flat displacement field**. The field has shape
`(source_rank, size(destination))`, in Fortran column-major pixel order.
Displacements are signed, real-valued pixel offsets relative to the one-based
identity map (the default, which requires equal input and output ranks).
For rank-changing maps, pass both `matrix` and `offset`; these use SciPy's
zero-based affine convention, and the resulting source position is

`1 + offset + matmul(matrix, destination_index - 1) + displacement(:,pixel)`.

`evaluate_spline_plan_displacement_grid_nd` takes the more natural dense field
shape `(source_rank, size(destination,1), ..., size(destination,n))`. Fortran's
rank limit restricts that convenience interface to output ranks 1-14. The
flat-field interface also supports rank-15 destinations. Both variants
support noncontiguous outputs, all spline orders 0-5, all eight boundary
modes, and complex source and output data. Existing scalar and batched
geometric callbacks are unchanged.

Example for a two-dimensional flow field:

```fortran
use scifort_ndimage, only : spline_interpolation_plan, &
    prepare_spline_plan_nd, evaluate_spline_plan_displacement_grid_nd
! image(nx,ny), output(nx,ny), flow(2,nx,ny) are allocated by caller.
call prepare_spline_plan_nd(image,plan,order=3,mode='reflect',status=code)
if (code/=0) error stop 'spline preparation failed'
call evaluate_spline_plan_displacement_grid_nd(plan,flow,output,code,batch_size=256)
if (code/=0) error stop 'displacement warp failed'
```

The implementation shares the existing bounded coordinate/interpolation
batches (`batch_size` is 1-8192, default 256). Output arrays remain limited
to 8,000,000 elements. The dense-grid convenience path temporarily flattens
the supplied displacement field, so its memory cost is O(source_rank *
output_pixels); the direct flat-field path avoids that copy. In either case,
a packed output is required for arbitrary-rank/noncontiguous destinations.
Mismatched field shapes or numerical types, incomplete affine parameters,
nonfinite offsets/displacements, invalid batch size, and an identity map
with unequal ranks return invalid-input status (`1`). Destinations are
unspecified on failure.

Native validation is in `test/test_signal_spline_plan_displacements.f90`.
Independent SciPy `ndimage.map_coordinates` comparisons and a rank-changing
affine+flow reference are in
`tools/signal_reference/compare_spline_plan_displacements_scipy.py`.
These are SciFort extensions; SciPy provides numerical references, not source.

## Performance and error reporting

High-dimensional tensor-product interpolation requires up to `(order+1)^rank`
terms per point, before combining repeated mirrored indices. Ranks 1–15 are
supported by the interfaces and filter kernels, but high-order interpolation
at rank 15 is often impractical. To avoid runaway work, coordinate evaluation
returns `spline_nd_status_too_expensive` (value `2`) when the number of distinct
neighbor combinations would exceed 2,000,000 for a point. Singleton or short
axes reduce this count, making many rank-15 datasets practical. A nonzero
status indicates that output values must not be used.

A shape/rank mismatch, empty input dimension, invalid axis/order/mode, or
nonfinite samples/coordinates returns `spline_status_invalid_input` (`1`).
Functions with no numeric failure return `spline_status_success` (`0`).

## Verification and provenance

The library retains its MIT license. The spline kernels implement cardinal
B-spline collocation independently; the corresponding SciPy `ndimage`
functions are used only as numerical references, with acknowledgement in
`CODE_PROVENANCE.md`. No SciPy source is redistributed or imported.

`test/test_signal_spline_nd.f90` checks the public interfaces, integer-knot
reconstruction, complex values, strided sections, rank-15 arrays, and error
reporting. `tools/signal_reference/compare_signal_spline_nd_scipy.py` compares
288 deterministic cases (5,760 interpolation points) with SciPy 1.17.0,
including orders 0–5, ranks 1–5 and 15, complex inputs and boundary modes.
The high-rank reference collapses singleton axes to avoid SciPy's exponential
interpolation workspace requirements. The maximum observed absolute
coefficient and interpolation differences were 1.28e-11 and 3.91e-13,
respectively. Tested with GNU Fortran 14.2.0.

## Callback geometric transforms

`geometric_transform_nd` is also exported by `scifort_ndimage`. Its mapping
callback accepts one-based integer coordinates of an output pixel and fills a
real(dp) vector of one-based sampling coordinates in the input. Mapping is
called in C-order. The caller chooses the rank and extents of the output;
its rank need not equal the input rank. Cubic interpolation with spline
prefiltering, constant boundaries, and zero fill is the default, following
SciPy's `ndimage.geometric_transform` defaults apart from the one-based
coordinate convention. See `docs/ndimage.md` for callback and resource-limit
details.
