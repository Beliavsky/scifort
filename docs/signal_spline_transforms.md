# Rank-15 spline boundary modes and geometric transforms

SciFort provides real and complex Fortran 2018 assumed-rank versions of
`map_coordinates_nd`, `affine_transform_nd`, `shift_nd`, `zoom_nd`, and `rotate_nd`. All accept input
arrays of ranks **1 through 15**, orders **0 through 5**, and one-based array
coordinates. The caller allocates each output; no separate implementation is
needed for each rank. The functions are also exported from `scifort_signal`.

## Zoom and rotate: ranks 1-15 and 2-15

`zoom_nd` and `rotate_nd` extend the same rank-independent sampler with
SciPy-style geometric transformations. They accept `real(dp)` and
`complex(dp)` inputs, caller-allocated outputs, noncontiguous sections,
interpolation orders 0-5, and all eight boundary modes. No new Fortran
package dependency is needed.

| Procedure | Input rank | Notes |
|---|---|---|
| `zoom_shape_nd(shape(x), factors, status)` | 1-15 | Predict output extents using ties-to-even rounding, as NumPy does |
| `zoom_nd(x, factors, y, order, status, mode, prefilter, cval, grid_mode)` | 1-15 | Scalar factor or per-axis array; `grid_mode` defaults to false |
| `rotate_shape_nd(shape(x), angle, axes, reshape, status)` | 2-15 | Predict rotated extents; `reshape` defaults to true |
| `rotate_nd(x, angle, y, axes, reshape, order, status, mode, prefilter, cval)` | 2-15 | Rotate by degrees in the plane of two one-based axes, default `[1,2]` |

For example:

```fortran
use scifort_signal, only : dp, zoom_nd, rotate_nd, &
    zoom_shape_nd, rotate_shape_nd, spline_status_success
real(dp) :: x(8,12)
real(dp), allocatable :: enlarged(:,:), turned(:,:)
integer, allocatable :: sz(:)
integer :: status

! Initialize x before use.
sz = zoom_shape_nd(shape(x), [1.5_dp,0.75_dp], status)
if (status /= spline_status_success) error stop 'invalid zoom size'
allocate(enlarged(sz(1),sz(2)))
call zoom_nd(x, [1.5_dp,0.75_dp], enlarged, order=3, &
    prefilter=.true., mode='mirror', status=status)
if (status /= spline_status_success) error stop 'zoom failed'

sz = rotate_shape_nd(shape(x), 30.0_dp, reshape=.true., status=status)
if (status /= spline_status_success) error stop 'invalid rotation size'
allocate(turned(sz(1),sz(2)))
call rotate_nd(x, 30.0_dp, turned, order=3, prefilter=.true., &
    mode='constant', cval=0.0_dp, status=status)
if (status /= spline_status_success) error stop 'rotation failed'
```

`zoom_nd` calculates its output shape from the input extents times the
requested factors. As in SciPy, the factors determine *only the output
shape*; the actual coordinate scaling is calculated from the input and output
extents. With `grid_mode=.false.`, the first and last sample centers align.
With `grid_mode=.true.`, pixel edges align instead. Each zoom factor must
be finite and positive, and any rounded output extent of zero is rejected,
because the current spline sampling core requires nonempty dimensions.

`rotate_nd` uses the same positive-angle orientation as
`scipy.ndimage.rotate`. When `reshape=.true.`, the two rotated extents
are computed from the outer pixel-boundary corners, as in SciPy. Other axes
remain unchanged. Both the default and explicitly reversed axis pairs are
supported, provided they contain two distinct valid axes.

**SciFort compatibility conventions:** Just like `affine_transform_nd`, the
new transforms default to `mode='mirror'` and `prefilter=.false.`. SciPy
defaults to constant boundaries and prefiltering. Specify the mode and
prefilter arguments to compare directly. High-order `prefilter=.true.`
with `nearest` or `grid-constant` now uses SciPy-style 12-sample prepadding;
short `reflect` arrays retain previously documented numeric differences. A
nonzero status means the output must not be used.

The reproducible reference `tools/signal_reference/compare_signal_zoom_rotate_scipy.py`
checks 1,156 real/complex zoom and rotation cases against SciPy 1.17, with
rank-2 and rank-3 inputs, selected axes, both reshape/grid modes and
interpolation orders 0, 1, 3, and 5. It found 1,116 matches within
`2.47e-15` and 40 short-array `reflect`-prefilter comparisons within
`1.76e-4`; those use a documented `2.5e-4` tolerance. In SciPy 1.17,
complex `zoom` unexpectedly omits `grid_mode` internally, so the comparison
uses separate real and imaginary SciPy zooms when `grid_mode=.true.`.
The permanent Fortran regression also exercises rank-15 arrays,
scalar zoom factors, shape validation, and strided outputs.

## Boundary conditions

The interpolation routines support all eight named `scipy.ndimage` modes:

| Mode | Outside the sampled domain |
|---|---|
| `mirror` | Whole-sample symmetry, reflecting around the first/last sample centers |
| `reflect` | Half-sample symmetry, reflecting around the edges of the end samples |
| `grid-mirror` | Alias of `reflect` |
| `nearest` | Repeat the closest endpoint sample/coefficient |
| `wrap` | Repeat a length-`n-1` interval, with overlapping endpoints |
| `grid-wrap` | Periodic interval of length `n`, without overlapping endpoints |
| `constant` | Return `cval` when the query point lies outside the domain |
| `grid-constant` | Extend by `cval`, **but interpolate across the boundary** |

`constant` and `grid-constant` are deliberately different for orders > 0.
The default remains `mirror`, preserving existing SciFort behavior. This
continues to differ from SciPy's default `mode='constant'`.

For `map_coordinates_nd`, input coordinates have shape `(rank(input), npoints)`
and refer to **one-based** input positions. In contrast, SciPy uses zero-based
coordinates. The output is the caller-supplied rank-one vector of length
`npoints`; all boundary modes accept prefiltered coefficient arrays by default.

### Prefiltering rules and current limits

`prefilter=.false.` remains the default: the input is treated as B-spline
coefficients, with **no** filtering. This path supports **all eight modes** and
all interpolation orders 0–5. The same behavior applies to `affine_transform_nd`
and `shift_nd`.

For `prefilter=.true.` and spline orders 2–5:

- `mirror`, `constant`, `wrap` use the original whole-sample collocation filter.
- `reflect` and `grid-mirror` use exact half-sample symmetric collocation.
- `grid-wrap` solves a periodic, circulant spline system along each axis.
- `nearest` and `grid-constant` prepad each nonsingleton axis by 12 samples
  on each side (edge replication or constant fill, respectively) before
  filtering. The implementation retains the padded spline coefficients for
  interpolation at locations outside the original array. Singleton dimensions
  are handled without allocating 25 elements per singleton axis: `nearest`
  contributes a constant factor, and `grid-constant` uses a one-dimensional
  zero-padded cardinal-spline impulse response. Thus rank-15 arrays with many
  singleton axes remain practical. Orders 0 and 1 bypass prefiltering.
  A prepad requiring more than 8,000,000 packed complex coefficients returns
  `spline_nd_status_too_expensive` (status 2) rather than exhausting memory.

For very short arrays, exact `reflect` collocation differs slightly from SciPy's
finite-precision recursive filter initialization. This was up to `1.40e-3`
for the tested 4-by-5 samples at spline order 5; the difference is not a
boundary-indexing error. In fact SciPy 1.17 can fail to exactly reconstruct
even a two-sample array at its original knots with `reflect` prefiltering,
whereas SciFort retains exact finite-line interpolation. This numerical
difference is deliberately preserved instead of degrading the exact solve. The periodic `grid-wrap` prefilter now uses causal/anticausal pole
recurrences requiring O(n) work per line instead of the former O(n squared)
direct DFT. For repeated transformations, supplying precomputed spline
coefficients still avoids repeated filtering.

## Affine transformations

```fortran
use scifort_signal, only : dp, affine_transform_nd, shift_nd, &
    spline_status_success
real(dp) :: input(64,48), output(80,50), shifted(64,48)
real(dp) :: matrix(2,2)
integer :: status

! Populate input first.
matrix(1,:) = [0.9_dp, -0.15_dp]
matrix(2,:) = [0.1_dp,  0.95_dp]
call affine_transform_nd(input, matrix, [2.0_dp,-1.0_dp], output, &
    order=3, mode='mirror', prefilter=.true., status=status)
if (status /= spline_status_success) error stop 'affine failed'

call shift_nd(input, [1.25_dp,-2.5_dp], shifted, &
    order=1, mode='reflect', status=status)
if (status /= spline_status_success) error stop 'shift failed'
```

The `matrix` has shape `(rank(x), rank(y))`; the input and output must have the
**same rank**, although their shapes can differ. `offset` is a vector of length
`rank(x)`. The input coordinate for an output index vector `o` (one-based) is

`input_coordinate = 1 + offset + matmul(matrix, o - 1)`.

Thus, `offset` measures displacement in sample units in the zero-based
coordinate system used by `scipy.ndimage.affine_transform`. This convention
avoids a hidden extra offset when comparing Fortran and Python examples.
`shift_nd(x, displacement, y)` is equivalent to an identity affine matrix
and `offset=-displacement`; positive shifts move the signal toward higher
array indices.

Both routines support real and complex input/output arrays and noncontiguous
sections. Their output is defined **only when** the status is success. Invalid
matrix/offset shapes, orders, modes, or nonfinite numerical parameters produce
`1` (`spline_status_invalid_input`). The same high-dimensional operation cap
used by `map_coordinates_nd` also applies to transforms.

For memory efficiency, affine transformations generate input coordinates in
batches of 128 destination points, rather than allocating one coordinate
matrix for the entire output. The original input is packed only once per call.

## Verification

- `test/test_signal_spline_modes_transform.f90`: public API, 8 boundary modes,
  real/complex data, spline orders 0–5, identity reconstruction, rank 15,
  noncontiguous sections, shifts, affine mappings, and rejected inputs.
- `tools/signal_reference/compare_signal_spline_modes_scipy.py`: 576
  deterministic reference comparisons against SciPy 1.17.0 (real and complex
  interpolation, affine transformations, and shifts). **528 comparisons**
  agree to `2.89e-15` maximum absolute error; the **48** comparisons with
  short-array `reflect`/`grid-mirror` prefiltering agree within `1.40e-3`.
- `tools/signal_reference/compare_signal_spline_prepad_scipy.py`: 16 additional
  high-order comparisons with a singleton axis, nonzero fill, and complex data
  (maximum difference `9.37e-15`).
- Original spline reference suite still passes 288 cases, including 5,760
  coordinate evaluations and real/complex rank-15 arrays.
- The GNU Fortran 14.2.0 checked signal suite passes 107/107 tests. The
  optimized targeted transform regression passes at `-O3` with bounds checks
  and floating-point traps. The full optimized signal build exceeded the
  execution limit; `fpm test` was not available.

This source is an independent MIT-licensed spline implementation. SciPy's
public API and numerical results are reference oracles, not redistributed
source. For general interpolation formulas see the cardinal B-spline
collocation literature. Relevant SciPy documentation:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.map_coordinates.html
and https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.affine_transform.html.
