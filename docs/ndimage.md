# Multidimensional image filtering and morphology

The `scifort_ndimage` module extends the spline and geometric-transform APIs
with **rank-independent** Gaussian, full-kernel, neighborhood, distance-transform, and binary image-processing
routines. Filtering operations generally retain the input rank and shape;
`geometric_transform_nd` also supports a separately shaped output of a
different rank. Arrays support ranks 1 through 15. All axis parameters are **one-based** and data are stored in
normal Fortran column-major array order. Implementation uses finite separable
kernels rather than copying SciPy source.

## Public procedures

| Procedure | Purpose | Input type |
|---|---|---|
| `geometric_transform_nd(x, mapping, y, order, status, mode, prefilter, cval)` | Output-to-input callback interpolation; output rank can differ | real/complex |
| `prepare_spline_plan_nd(x, plan, order, status, mode, prefilter, cval)` | Prepare reusable spline coefficients for multiple queries | real/complex |
| `evaluate_spline_plan_nd(plan, coords, values, status)` | Query a prepared spline plan | real/complex |
| `evaluate_spline_plan_grid_nd(plan, coords, y, status)` | Query a prepared plan into a shaped output | real/complex |
| `evaluate_spline_plan_geometric_batch_nd(plan, mapping, y, status, batch_size)` | Cached spline warp with one coordinate callback per block (**SciFort extension**) | real/complex |
| `evaluate_spline_plan_displacement_nd(plan, displacements, y, status, matrix, offset, batch_size)` | Cached spline warp with flat per-pixel displacement vectors (ranks 1-15) | real/complex |
| `evaluate_spline_plan_displacement_grid_nd(plan, field, y, status, matrix, offset, batch_size)` | Cached spline warp with leading-component displacement grid (output ranks 1-14) | real/complex |
| `gaussian_filter_nd(x, sigma, y, order, status, mode, cval, truncate, radius, modes, axes)` | Separable Gaussian smoothing or partial derivatives | real/complex |
| `fourier_gaussian_nd(x, sigma, y, status, n, axis)` | Gaussian smoothing applied to complex Fourier coefficients | complex |
| `fourier_ellipsoid_nd(x, sizes, y, status, n, axis)` | Fourier transform of a uniform N-dimensional ellipsoidal kernel | complex |
| `fourier_uniform_nd(x, sizes, y, status, n, axis)` | Fourier multiplier for uniform smoothing | complex |
| `fourier_shift_nd(x, shift, y, status, n, axis)` | Fractional translation in the Fourier domain | complex |
| `gaussian_filter1d_nd(x, sigma, axis, y, order, status, mode, cval, truncate, radius)` | Gaussian kernel along one axis | real/complex |
| `gaussian_gradient_magnitude_nd(x, sigma, y, status, mode, cval, truncate)` | Euclidean norm of Gaussian first derivatives | real |
| `gaussian_laplace_nd(x, sigma, y, status, mode, cval, truncate)` | Sum of Gaussian second partial derivatives | real |
| `generic_gradient_magnitude_nd(x, derivative, y, status, mode, cval, modes, axes)` | Selected-axis gradient using a caller-supplied derivative; complex result uses `sqrt(sum(d**2))` | real/complex |
| `generic_laplace_nd(x, derivative2, y, status, mode, cval, modes, axes)` | Sum of selected-axis second derivatives | real/complex |
| `uniform_filter_nd(x, sizes, y, status, mode, cval, modes, axes)` | Rectangular moving mean | real/complex |
| `uniform_filter1d_nd(x, filter_size, axis, y, status, mode, cval, origin)` | Axis-specific moving mean with origin shift | real/complex |
| `minimum_filter1d_nd(x, filter_size, axis, y, status, mode, cval, origin)` | Axis-specific sliding minimum | real |
| `maximum_filter1d_nd(x, filter_size, axis, y, status, mode, cval, origin)` | Axis-specific sliding maximum | real |
| `minimum_filter_nd(x, sizes, y, status, mode, cval, footprint, origin, modes, axes)` | Rectangular or footprint-defined minimum | real |
| `maximum_filter_nd(x, sizes, y, status, mode, cval, footprint, origin, modes, axes)` | Rectangular or footprint-defined maximum | real |
| `sobel_nd(x, axis, y, status, mode, cval)` | Sobel derivative on a selected axis | real |
| `prewitt_nd(x, axis, y, status, mode, cval)` | Prewitt derivative on a selected axis | real |
| `laplace_nd(x, y, status, mode, cval)` | Sum of second differences along each axis | real |
| `correlate_nd(x, weights, y, status, mode, cval, modes)` | Full N-D correlation with arbitrary weights | real/complex |
| `convolve_nd(x, weights, y, status, mode, cval, modes)` | Full N-D convolution with arbitrary weights | real/complex |
| `correlate1d_nd(x, weights, axis, y, status, mode, cval, origin)` | 1-D correlation along a one-based axis | real/complex |
| `convolve1d_nd(x, weights, axis, y, status, mode, cval, origin)` | 1-D convolution along a one-based axis | real/complex |
| `generic_filter_nd(x, func, y, status, sizes, footprint, origin, mode, cval, modes, axes)` | Caller-supplied reduction over each window | real/complex |
| `vectorized_filter_nd(x, func, y, status, sizes, footprint, origin, mode, cval, modes, axes, batch_size)` | Batched scalar results, including mixed real/complex types; `valid` may shrink | real/complex |
| `vectorized_filter_multi_nd(x, func, y, status, sizes, footprint, origin, mode, cval, modes, axes, batch_size)` | Batched multi-channel results including mixed numeric types; trailing or flat channels | real/complex |
| `generic_filter1d_nd(x, func, filter_size, axis, y, status, mode, cval, origin)` | Call a Fortran procedure once for each padded line along a selected axis | real/complex |
| `median_filter_nd(x, sizes, y, status, mode, cval, footprint, origin, modes, axes)` | N-D median neighborhood | real |
| `rank_filter_nd(x, rank_value, sizes, y, status, mode, cval, footprint, origin, modes, axes)` | Zero-based ordered neighborhood element | real |
| `percentile_filter_nd(x, percentile, sizes, y, status, mode, cval, footprint, origin, modes, axes)` | Nearest-rank neighborhood percentile | real |
| `distance_transform_edt_nd(x, y, status, sampling, feature_indices)` | Exact Euclidean distance and optional nearest-background coordinates | logical -> real/integer |
| `distance_transform_bf_nd(x, y, status, metric, sampling, feature_indices)` | Direct-search Euclidean, taxicab, chessboard distances and nearest-background sites | logical -> real/integer |
| `distance_transform_cdt_nd(x, y, status, metric, feature_indices, axis_weights)` | Taxicab/chessboard, optional per-axis weighted L1 distances | logical -> integer |
| `distance_transform_weighted_cdt_nd(x, weights, y, status, feature_indices)` | Dijkstra shortest path with custom integer neighbor costs | logical + integer -> integer |
| `binary_dilation_nd(x, y, status, connectivity, iterations, border_value, structure, origin, mask, axes)` | Grow a logical mask | logical |
| `binary_erosion_nd(x, y, status, connectivity, iterations, border_value, structure, origin, mask, axes)` | Shrink a logical mask | logical |
| `binary_opening_nd(x, y, status, connectivity, iterations, border_value, structure, origin, mask, axes)` | Erode then dilate | logical |
| `binary_closing_nd(x, y, status, connectivity, iterations, border_value, structure, origin, mask, axes)` | Dilate then erode | logical |
| `binary_fill_holes_nd(x, y, status, connectivity, structure, axes)` | Fill enclosed background along specified axes | logical |
| `binary_propagation_nd(x, mask, y, status, connectivity, structure, axes)` | Geodesic flood-fill along specified axes | logical |
| `binary_hit_or_miss_nd(x, y, status, structure1, structure2, origin1, origin2)` | Locate prescribed foreground/background patterns | logical |
| `grey_dilation_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Grayscale neighborhood maximum | real |
| `grey_erosion_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Grayscale neighborhood minimum | real |
| `grey_opening_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Grayscale erosion followed by dilation | real |
| `grey_closing_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Grayscale dilation followed by erosion | real |
| `morphological_gradient_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Difference of grayscale dilation and erosion | real |
| `morphological_laplace_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Grayscale dilation + erosion - 2 * image | real |
| `white_tophat_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Image minus grayscale opening | real |
| `black_tophat_nd(x, y, status, sizes, footprint, mode, cval, structure, origin, axes, modes)` | Grayscale closing minus image | real |
| `label_nd(x, y, num_features, status, connectivity, structure, axes)` | Label nonzero logical, integer or real foreground | logical/integer/real -> integer |
| `generate_binary_structure_nd(connectivity, structure, status)` | Fill a caller-sized `(3,)*rank` connectivity stencil | logical |
| `iterate_structure_nd(structure, iterations, expanded, status, origin, new_origin)` | Iterated binary structuring element (Minkowski sums) | logical |
| `sum_labels_nd(x, labels, index, totals, status)` | Sum pixels with selected integer region labels | real + integer -> real |
| `mean_labels_nd(x, labels, index, values, status)` | Selected region means | real + integer -> real |
| `variance_labels_nd(x, labels, index, values, status)` | Selected region population variances | real + integer -> real |
| `standard_deviation_labels_nd(x, labels, index, values, status)` | Selected region population standard deviations | real + integer -> real |
| `count_labels_nd(labels, index, counts, status)` | Pixel counts per selected region | integer -> integer |
| `histogram_labels_nd(x, labels, index, low, high, histograms, status)` | Histogram counts for selected regions | real + integer -> integer |
| `labeled_comprehension_nd(x, labels, index, func, default, values, status)` | Apply a user callback to each labeled region | real + integer -> real |
| `center_of_mass_nd(x, labels, index, centers, status)` | Mass-weighted region centers (one-based) | real + integer -> real |
| `extrema_labels_nd(x, labels, index, minimums, maximums, min_positions, max_positions, status)` | Minimum/maximum values and coordinates of labeled regions | real + integer -> real/integer |
| `minimum_labels_nd(x, labels, index, values, status)` | Minimum value for each label | real + integer -> real |
| `maximum_labels_nd(x, labels, index, values, status)` | Maximum value for each label | real + integer -> real |
| `minimum_position_nd(x, labels, index, positions, status)` | One-based position of each label minimum | real + integer -> integer |
| `maximum_position_nd(x, labels, index, positions, status)` | One-based position of each label maximum | real + integer -> integer |
| `bounding_boxes_nd(labels, index, lower, upper, status)` | Inclusive one-based bounds for selected labels | integer -> integer |
| `find_objects_nd(labels, lower, upper, status)` | Inclusive bounds of labels 1..N | integer -> integer |
| `value_indices_nd(x, values, starts, coordinates, status, ignore_value)` | Group integer pixel coordinates by distinct value (one-based CSR representation) | integer -> integer |
| `watershed_ift_nd(x, markers, y, status, structure, legacy_linear_bounds)` | Marker-based IFT segmentation | integer -> integer |

`sigma` accepts either a scalar or a per-axis vector for the all-axis
Gaussian functions, and `sizes` accepts either a scalar or a per-axis integer
vector for the neighborhood filters. `order` and `radius` are scalars in the
scalar-sigma overloads and integer vectors in the vector-sigma overloads.
Gaussian derivative orders 0 through 3 are implemented, with `order=0` by
default. A zero `sigma` skips that axis only when its derivative order is zero;
negative sigma and positive derivative orders with zero sigma are rejected.
The 1-D function requires an explicit one-based `axis`. Gaussian-gradient
magnitude and Gaussian Laplacian require positive sigma in every dimension.

`mode` applies to all axes and defaults to `reflect` (half-sample symmetry),
just as SciPy's filtering routines do. Available modes are `reflect`,
`grid-mirror`, `mirror`, `nearest`, `wrap`, `grid-wrap`, `constant`, and
`grid-constant`. For **filtering**, `wrap` and `grid-wrap` both use period
`size(x,axis)`; `constant` and `grid-constant` are aliases. `cval` defaults to
zero. Unlike geometric spline interpolation, filtering does not interpret
`wrap` with overlapping endpoints. `truncate=4.0_dp` determines Gaussian
radius by default; supplying `radius` overrides it.

`gaussian_filter_nd`, `uniform_filter_nd`, `minimum_filter_nd`,
`maximum_filter_nd`, `correlate_nd`, `convolve_nd`, `generic_filter_nd`,
`median_filter_nd`, `rank_filter_nd`, and `percentile_filter_nd`
accept an optional `modes(:)` argument. Without `axes(:)`, it contains one
character string per Fortran axis; with `axes(:)`, it contains one string
per selected axis in the caller-supplied order. It selects a different
boundary extension for each filtered dimension;
`mode` and `modes` cannot both be provided. This extends the existing scalar
`mode=` API without breaking existing calls. For full N-D kernels and
nonseparable footprints, SciFort supports per-axis modes independently even
where SciPy rejects sequences of modes. `generic_filter_nd` also accepts
complex input/output and a complex-valued callback; its callback receives
active neighbors in C-order (last index varies fastest), matching the real
callback's ordering. Complex generic callbacks are a SciFort extension, not
part of SciPy's standard generic-filter API.


### Axis-selective Gaussian and rectangular filtering

`gaussian_filter_nd`, `uniform_filter_nd`, `minimum_filter_nd` and
`maximum_filter_nd` now accept an optional **one-based** `axes(:)` argument.
Axes must be distinct and lie in `[1,rank(x)]`; the empty vector selects no
axes and returns the input unchanged. As in SciPy, when `axes` is present,
**vector** `sigma`, `sizes`, `order`, `radius`, and `modes` parameters are
supplied for *only those axes*, in the exact order listed. Scalar `sigma`,
`sizes`, `order` and `radius` apply to all selected axes; the other axes are
left untouched. Operations are composed in the order specified in `axes`.
When `axes` is omitted, the existing one-value-per-dimension convention is
unchanged. A common `mode` applies to every selected axis; `modes(:)` instead
specifies one boundary mode for each selected axis. Passing both is invalid.

```fortran
real(dp) :: image(8,9,10), blurred(8,9,10), low(8,9,10)
integer :: status
! Only filter the third and then first axes (the middle axis is unchanged).
call gaussian_filter_nd(image,[0.8_dp,1.4_dp],blurred, &
    axes=[3,1],order=[0,1],modes=[character(len=16) :: 'nearest','reflect'], &
    status=status)
call minimum_filter_nd(image,[4,2],low,axes=[3,1],mode='constant', &
    cval=-5.0_dp,status=status)
```

For `minimum_filter_nd` and `maximum_filter_nd`, `axes` currently applies to
rectangular windows only. Combining `axes` with `footprint` or `origin` is
rejected with status 1; the footprint-based minimum/maximum API remains
available for full-rank windows. Median, rank and percentile filters now
accept `axes` with compact footprints and origins (details below).

Rectangular minimum and maximum filters (including `minimum_filter1d_nd` and
`maximum_filter1d_nd`) now use **monotone-deque sliding windows**. A line of
length `n` with width `w` requires O(`n+w`) time and memory, instead of
O(`n*w`) comparisons. Extensions are generated on demand using the same
boundary modes and shifted-origin conventions as the prior implementation.
The footprint-defined extrema route continues to enumerate active stencil
sites and is not accelerated by the deque. Uniform filtering, including the
origin-shifted `uniform_filter1d_nd`, now uses a similar **running sum**:
a line of length `n` and a window of width `w` takes O(`n+w`) time, not
O(`n*w`). It maintains a separate O(`n`) output line, generates boundary
samples on demand and supports complex values. Successive additions and
subtractions can change last-bit floating-point rounding relative to the
older direct convolution (as in SciPy's running-sum implementation).
The 50,000-site 1-D limit remains
unchanged; the rectangular N-D limit remains `2*max_kernel_radius+1`.

The focused reference tests `compare_selected_axes_scipy.py` and
`compare_gaussian_axes_scipy.py` include mixed boundary modes, even-width
windows, a width of 50,000, Gaussian derivatives, complex filtering, and
validation of malformed axis specifications. The native test
`test_ndimage_selected_axes_deque.f90` additionally checks rank-15 inputs and
strided output arrays.

`compare_order_axes_running_scipy.py` checks another 10,740 real/complex
SciPy values in axis-selected median/rank/percentile filters and running
means, including reordered axes, shifted even-size footprints, all eight
boundary extensions and windows as wide as 50,000. The native
`test_ndimage_order_axes_running.f90` checks invalid arguments, direct
convolution consistency, rank-15 data and noncontiguous output sections.

The Fourier-domain filters operate **directly on precomputed complex Fourier
coefficients**; they do not compute the FFT or inverse FFT. For `n=-1`
(default), every axis uses the full complex FFT ordering. For `n>0`, the
`axis` selects the packed real-FFT dimension of original length `n`, and the
corresponding input/output dimension must have extent `n/2+1` (integer
quotient); all other axes use full FFT frequency ordering. `axis` defaults
to the last axis and is **one-based**, although negative indexing such as
`axis=-1` is also accepted. Fourier Gaussian, uniform, ellipsoid, and shift parameters
accept either a scalar or one value per dimension. Outputs retain the
original complex spectrum's dimensions, ready for an inverse transform.
`fourier_ellipsoid_nd` multiplies each spectrum coefficient by the normalized
Fourier transform of a uniform ellipsoid with specified axis diameters. SciPy
implements this operation only in ranks 1--3; SciFort extends the same
Bessel-function formula to ranks 4--15. Both positive and negative sizes are
accepted because their geometric Fourier response depends on squared radii.
The routine rejects nonfinite parameters and sizes exceeding 1e6 in absolute
value, and follows the packed real-FFT `n`/`axis` conventions above.

`correlate_nd` and `convolve_nd` accept arbitrary real/complex kernels with the
same rank as `x`, **including asymmetric and even-sized kernels**. Complex
correlation conjugates the complex coefficients, following SciPy. The current
overloads require real kernels for real images and complex kernels for complex
images; a real/complex mixed expression can first be converted with `cmplx`.
`median_filter_nd`, `rank_filter_nd`, and `percentile_filter_nd` accept a common
scalar window width or a vector of widths. A negative `rank_value` counts
backwards from the largest neighborhood element (`-1` selects its maximum).
Percentiles are limited to `[0,100]` and follow SciPy's discrete nearest-rank
convention rather than interpolating values. These three order-statistic APIs
also accept `footprint(..)` and `origin(:)`, including asymmetric or even-sized
stencils; the active footprint cells alone determine the sample count and
requested rank. `footprint` overrides `sizes`. A vector `sizes(:)` may be
omitted when supplying a footprint; then pass the required output by keyword,
for example `median_filter_nd(image, y=filtered, footprint=stencil)`.
Origins have valid range `[-width/2,(width-1)/2]` on each axis. The optional
`axes(:)` selects a subset of one-based axes, including an empty subset
(which returns the input unchanged). With `axes` present, `sizes(:)` and the
rank/shape of `footprint(..)` cover only the selected dimensions, not the
full image. To match SciPy's *nonseparable rank-filter* convention, these
compact stencil dimensions are assigned to the selected axes **in ascending
axis order**, regardless of the order in `axes(:)`. By contrast,
`origin(:)` entries refer to the selected axes **in the order supplied**;
this distinction matters for unsorted axes and even-sized footprints.
`modes(:)` similarly refers to the selected axes in the supplied order, as
in the other SciFort axis-selective filters. A scalar `mode` applies to all
axes. Providing both `mode` and `modes` is invalid. SciPy may reject mode sequences
for nonseparable footprints; SciFort applies each axis extension independently.
A footprint without active entries or a rank outside the selected entries
is rejected. Selected footprints are created in compact rank before their
offsets are embedded in the original rank-1 through rank-15 image.

```fortran
real(dp) :: image(8,9,10), result(8,9,10)
integer :: status
! SciPy's rank-filter convention: widths 4,3 apply to axes 1,3
! in ascending order, even though the axes are supplied as [3,1].
call median_filter_nd(image,[4,3],result,axes=[3,1], &
    mode='nearest',status=status)
```


`generic_filter1d_nd` accepts a one-based `axis` and a Fortran callback
with interface `subroutine func(input_line, output_line)`, both arrays being
`real(dp)` for real images or both being `complex(dp)` for complex images. Each `input_line` is the selected axis extended by `filter_size-1`
samples; `output_line` has the original axis extent. The callback is invoked
once per line rather than once per pixel, and is responsible for filling every
output element. Its padded argument matches SciPy's `generic_filter1d`
convention, including even-sized filters, signed `origin`, all eight scalar
boundary modes, and constant `cval`. The filter size is limited to 100,000
samples; invalid arguments and resource-limit failures return status 1 and 2.

`minimum_filter_nd` and `maximum_filter_nd` now support optional Boolean
`footprint(..)` and per-axis `origin(:)` with the same origin convention as
`median_filter_nd`. The vector overload accepts omitted `sizes` when a
footprint is given (write `y=` to identify the required output); a default
3-by-3-by-... neighborhood is used if both are omitted. Rectangular size-only filtering uses the separable monotone-deque implementation.
With a footprint or origin, the code enumerates the selected stencil sites
in Fortran column-major order. Unlike `grey_dilation_nd`, `maximum_filter_nd`
uses the **unreversed** stencil offsets, matching SciPy's `maximum_filter`
for asymmetric/even-sized footprints. Empty footprints, invalid origins, and
footprints above the 50,000-site limit are rejected.

`distance_transform_cdt_nd` provides SciPy-style *discrete* distances for
`metric='taxicab'` (the default, also `'cityblock'` or `'manhattan'`) and
`metric='chessboard'`. It uses separable integer distance transforms; for the
chessboard metric each 1-D line is processed with a range-min segment tree.
The input foreground mask is logical; the caller provides a same-shaped
**integer** distance output. A zero distance is returned at background pixels.
If no background exists, the returned distances are `-1`, following SciPy.
The optional `feature_indices(rank(x),size(x))` contains **one-based**
nearest-background coordinates in column-major flattened order, as for EDT.
When multiple background pixels tie for minimum distance, the chosen site is
valid but its index need not match SciPy's tie-breaking order. For an
all-foreground array, feature indices remain the pixel's own coordinates.
Weighted axis costs and general integer chamfer stencils are available separately, as described below.

`count_labels_nd` returns counts for arbitrary integer labels, including zero,
negative, repeated, and absent labels. `histogram_labels_nd` writes a
`(number_of_bins,number_of_selected_labels)` integer array; its equally
spaced bins cover `[low,high]` and include `high` in the last bin, like
NumPy/SciPy histogram semantics. Missing labels produce zero histograms
(instead of SciPy's `None`), and out-of-range pixel values are ignored.
`standard_deviation_labels_nd` takes the square root of the population
variance (divisor `N`), returning NaN for absent labels. All these procedures
accept images or labels of ranks 1 through 15.

For example:

```fortran
use scifort_kinds, only : dp
use scifort_ndimage, only : minimum_filter_nd, distance_transform_cdt_nd, &
    histogram_labels_nd
real(dp) :: image(6,7), lowpass(6,7)
logical :: footprint(3,3), foreground(6,7)
integer :: distances(6,7), nearest(2,42), labels(6,7), hist(8,2), status
! Initialize image, footprint, foreground, and labels first.
call minimum_filter_nd(image,y=lowpass,status=status, &
    footprint=footprint,origin=[-1,0])
call distance_transform_cdt_nd(foreground,distances,status, &
    metric='chessboard',feature_indices=nearest)
call histogram_labels_nd(image,labels,[1,2],0.0_dp,8.0_dp,hist,status)
```

`distance_transform_bf_nd` implements SciPy's independent brute-force
reference algorithm: for each foreground cell it searches all background cells.
The `metric` keyword is `euclidean` (default), `taxicab`, `chessboard`,
`cityblock`, or `manhattan`; `sampling(:)` is a positive per-axis physical
spacing used for Euclidean distances only. Unlike the faster integer CDT
routine, the caller allocates a **real(dp)** output for every metric (including
integer-valued taxicab/chessboard distances). Optional
`feature_indices(rank(x),size(x))` uses the existing one-based Fortran
column-major convention. Equidistant sites may be selected differently than
in SciPy, but are valid closest-background witnesses. For a mask without
background, SciPy's brute-force sentinel is returned: `sqrt(huge(1.0_dp))`
for Euclidean distance or `4294967295.0_dp` otherwise, with all feature
indices zero. This differs from the virtual-site convention of the EDT API.
The brute-force implementation intentionally imposes a 30,000,000
`uniform_filter1d_nd`, `minimum_filter1d_nd` and
`maximum_filter1d_nd` process only the selected **one-based** axis, leaving
all other dimensions unchanged. Their scalar `filter_size` is positive, and
`origin` follows SciPy's valid range `[-filter_size/2, (filter_size-1)/2]`
(integer division). All eight extension modes are supported. Window widths
above 50,000 return resource status 2. Uniform filtering also accepts
`complex(dp)` input and `cval`. The new running-sum implementation is
linear in line length plus window width; memory is linear in line length.
These operations retain ranks 1 through 15.

`size(x) * number_of_background_pixels` comparison limit; it returns status
`ndimage_status_too_expensive` (2) when exceeded, rather than launching
unbounded quadratic work. Use `distance_transform_edt_nd` or
`distance_transform_cdt_nd` on large masks.

`generic_gradient_magnitude_nd` and `generic_laplace_nd` accept a derivative
subroutine for each axis. Callbacks have the explicit interface
`ndimage_derivative_func`:

```fortran
subroutine derivative(x,axis,y,mode,cval,status)
    real(dp), intent(in) :: x(..)
    integer, intent(in) :: axis
    real(dp), intent(out) :: y(..)
    character(len=*), intent(in) :: mode
    real(dp), intent(in) :: cval
    integer, intent(out) :: status
end subroutine derivative
```

The optional `axes(:)` argument restricts the derivative callbacks to a
nonempty set of unique **one-based** axes. The previous default uses all
axes. `modes(:)` may contain one string per selected axis in the same order
when a strict subset is selected, or one string per full image axis (the
existing SciFort convention). Invalid/duplicate axes and invalid modes
return `ndimage_status_invalid_input`. For example, `axes=[1,3]` on a
3-D image differentiates only the first and third axes, while
`modes=['nearest','constant']` applies the respective extensions to them.
The empty-axis case returns status 1, unlike SciPy's empty-axis shortcut.

The callback receives the *original* image, a **one-based** axis and a
same-shaped scratch output, and sets `status=0` on success. The gradient
operator accumulates the sum of squared first derivatives and takes its
square root; the generic Laplacian sums the second derivatives. The optional
`modes(:)` argument selects a boundary mode per axis; `mode=` and `modes=`
are mutually exclusive, and the selected mode is passed to each callback.
Any nonzero callback status is propagated. For example, a caller can
implement the callback by passing `x`, `axis`, `y`, `mode` and `cval` to
`correlate1d_nd` with weights `[-0.5_dp,0.0_dp,0.5_dp]` for the first
derivative or `[1.0_dp,-2.0_dp,1.0_dp]` for the second.

`distance_transform_edt_nd` returns the exact Euclidean distance from every
foreground (`.true.`) cell to its nearest background (`.false.`) cell, with
optional positive per-axis physical spacing `sampling(:)`. It computes
squared-distance lower envelopes along each axis and takes square roots at
the end; work scales linearly with pixel count times rank. If the input is
entirely foreground, it follows SciPy's special convention equivalent to a
single virtual background site at zero-based `(-1,0,...,0)`. The optional `feature_indices(rank(x),size(x))` returns the one-based coordinates
of the nearest background pixel for every flattened Fortran-order input
position: `feature_indices(:,i)` describes `x` at linear position `i`.
This two-dimensional layout supports input ranks up to 15 without requiring
a rank-16 integer array. For all-foreground images, the virtual site is
represented as `(0,1,1,...)` in one-based coordinates. In tie cases, more
than one site can minimize the distance; the implementation chooses a site
using its separable envelope order.

```fortran
use scifort_kinds, only : dp
use scifort_ndimage, only : gaussian_filter_nd, gaussian_gradient_magnitude_nd, &
    maximum_filter_nd, binary_dilation_nd, label_nd
real(dp) :: image(64,48), blurred(64,48), edges(64,48), peaks(64,48)
logical :: mask(64,48), expanded(64,48)
integer :: components(64,48), ncomponents, status

! Populate image before these calls.
call gaussian_filter_nd(image,[1.3_dp,0.8_dp],blurred,status=status)
call gaussian_gradient_magnitude_nd(image,1.0_dp,edges,status=status)
call maximum_filter_nd(edges,[3,3],peaks,status=status)
mask=edges>0.25_dp
call binary_dilation_nd(mask,expanded,status=status,iterations=2)
call label_nd(expanded,components,ncomponents,status=status)
```

## Binary connectivity

`generate_binary_structure_nd(connectivity, structure, status)` generates a
Boolean `(3,3,...,3)` stencil using the rank of the caller-provided output,
with face, edge, or full connectivity. `connectivity` is clamped to the
interval `[1,rank]`, as with SciPy. `iterate_structure_nd(structure,
iterations, expanded, status, origin, new_origin)` applies repeated Minkowski
sums of an arbitrary (including even-sized/asymmetric) binary footprint.
For `iterations>=2`, every output dimension must be
`(shape(structure)-1)*iterations+1`; for lower counts, the original shape
is copied. `origin`/`new_origin` are optional vectors for propagating the
per-axis footprint origins: `new_origin = origin*max(1,iterations)`. The
iterative path uses bounded work and returns status 2 if its operation or
memory limit is exceeded. Both procedures accept ranks 1 through 15.

`label_nd` now also accepts default-kind integer and `real(dp)` source
arrays. Any nonzero element is foreground and all disconnected foreground
pixels are segmented into positive labels; the original logical input
remains supported without any signature change. The labels' integer IDs
are assigned in Fortran storage order, not SciPy's row-major traversal
order. The partition and number of components match SciPy.

The logical morphology functions and `label_nd` use a centered radius-one
Boolean neighborhood. Default `connectivity=1` includes face-adjacent pixels
(along one axis at a time). At rank `d`, `connectivity=d` includes all
`3**d` immediate neighbors, including diagonals. Intermediate connectivity
values include offsets with at most that many nonzero components, matching
`scipy.ndimage.generate_binary_structure(d, connectivity)`. Background pixels
are labeled zero; distinct components receive positive integer labels.
Component **numbers** follow Fortran storage order and need not have SciPy's
numerical identifiers, but the connectivity partition agrees. `iterations`
for dilation and erosion defaults to 1; values of zero or less run until the
result stabilizes, with a finite-step resource guard returning status 2 for
nonconvergent custom stencils. `border_value` defaults to `.false.`.

Opening applies the specified number of erosions and then dilations;
closing applies dilations and then erosions. Hole filling floods the
background connected to any outer face of the array, then fills the remaining
background. `binary_propagation_nd` floods from all initial true seeds through
`mask`, retaining seeds even if they are outside the allowed mask (as SciPy
does). Propagation uses a queue and converges without an iteration limit.

The binary operations also accept an optional Boolean `structure(..)` with the
same rank as the image; it overrides `connectivity` when supplied. Arbitrary
odd- and even-width footprints are allowed. Each axis is centered at
`shape(structure)/2` (integer division, in zero-based coordinates). Dilation
uses the reflected element, while erosion uses the supplied element. The
same custom footprints work with opening, closing, propagation, hole filling,
and connected-component labeling. Labeling requires a centrosymmetric
footprint and returns an invalid-input status otherwise. The connected-component
identifiers still follow Fortran storage order.

`binary_hit_or_miss_nd` matches the foreground pattern `structure1` and
rejects foreground at every active position of `structure2`. By default,
`structure1` has face connectivity and `structure2` is its complement in the
same centered box. If only `structure1` is given, its complement is used for
`structure2`; both structures may be specified separately.

Binary dilation, erosion, opening, and closing now accept a per-axis integer
`origin(:)` argument (vector length equals input rank); `origin=0` is the
centered default. The valid interval for an axis with stencil width `w` is
`[-w/2, (w-1)/2]` with integer division. For binary hit-or-miss, `origin1(:)`
and `origin2(:)` shift foreground and background stencils separately;
`origin2` defaults to `origin1`. These optional arguments follow existing
positional parameters, preserving older calls. Origins do not apply to
`label_nd` or geodesic propagation.

Dilation, erosion, opening, and closing accept an optional logical `mask(..)`
with the same rank and shape as the image. Pixels where this mask is false keep
their current value **on every iteration**, including both stages of opening
and closing, matching SciPy's editable-pixel semantics. When `iterations<=0`,
each stage converges independently. This `mask` argument is unrelated to the
second positional `mask` argument of `binary_propagation_nd`, which instead
restricts reachable foreground. Labeling currently accepts logical masks
(not general numeric arrays). For high-rank hit-or-miss calls
using the default complement, the size of the `3**rank` stencil may exceed
the 50,000-element safety limit.

## Grayscale morphology

`grey_dilation_nd`, `grey_erosion_nd`, `grey_opening_nd`, `grey_closing_nd`,
`morphological_gradient_nd`, `morphological_laplace_nd`, `white_tophat_nd`, and
`black_tophat_nd` process real(dp) images of ranks 1 through 15. All take
caller-allocated same-shape output arrays. Choose `sizes(:)` for a filled
rectangular neighborhood or `footprint(..)` for a Boolean custom stencil;
when none is supplied, a full width-three neighborhood is used. A real(dp)
`structure(..)` provides non-flat additive weights: dilation takes a maximum
of `input + structure`, erosion a minimum of `input - structure`. An optional
Boolean `footprint(..)` selects active entries; when both are present their
shapes must agree. `structure` or `footprint` takes precedence over `sizes`.
The optional `origin(:)` provides the SciPy-style per-axis shift; its allowed
range is `[-w/2, (w-1)/2]` for a width-`w` stencil. A footprint with no
active elements is invalid. Even-width and asymmetric footprints are supported.
All eight numeric filtering boundary modes are supported, with `reflect` as
the default and `cval=0.0_dp` for constant extension. Optional `axes(:)`
restricts processing to selected one-based dimensions. In that case, `sizes`,
`footprint`, `structure`, and `origin` describe the compact selected-axis
stencil, without singleton dimensions; unselected image axes are left intact.
Reordered axes follow SciPy's separable-window convention for rectangular
sizes and SciPy's ascending-axis expansion for nonseparable footprints and
non-flat structures, including the special even-width dilation origins.
The new `modes(:)` argument accepts one boundary string per selected axis
(in the supplied `axes` order). For nonseparable stencils this mixed-mode
behavior extends SciPy, which accepts only one mode in that case. Supplying
both `mode` and `modes` is invalid; existing scalar `mode` calls are unchanged.
`axes=[integer::]` with `sizes=[integer::]` acts as the identity filter.
Grayscale morphology
requires ordered real arrays (not complex); NaN or infinite stencil weights
are rejected with `status=1`.

```fortran
use scifort_kinds, only : dp
use scifort_ndimage, only : grey_opening_nd, white_tophat_nd, &
    distance_transform_edt_nd
real(dp) :: image(5,5),opened(5,5),bright_peaks(5,5),distance(5,5)
logical :: footprint(3,3),mask(5,5)
integer :: feature_indices(2,25),status

footprint=.true.
footprint(1,1)=.false.
call grey_opening_nd(image,opened,status,footprint=footprint)
call white_tophat_nd(image,bright_peaks,status,footprint=footprint)
! On a 3D image, compact 2D footprint dimensions can target axes 3 and 1:
! call grey_dilation_nd(image3,filtered3,status,sizes=[2,3],axes=[3,1])
mask=image>0.0_dp
call distance_transform_edt_nd(mask,distance,status,feature_indices=feature_indices)
```

## Generic callbacks, one-axis kernels, and weighted chamfer metrics

`correlate1d_nd` and `convolve1d_nd` apply real or complex coefficients along
one **one-based** axis. Odd/even kernel centers and the optional scalar `origin`
match `scipy.ndimage.correlate1d`/`convolve1d`; valid origins are from
`-size(weights)/2` through `(size(weights)-1)/2` (integer division). Complex
correlation conjugates the coefficients; convolution does not. Both accept
all eight filtering modes and a constant `cval` where relevant.

`generic_filter_nd` calls a Fortran function for every footprint window,
with the callback interface `real(dp) function func(values)` or a complex(dp)
function for complex images, where `values(:)` contains the selected samples
in **C order**, as in SciPy, even though the input and output arrays are
stored in Fortran order. Supply `sizes(:)` or `footprint(..)` (optional
per-axis `origin(:)`); a footprint overrides sizes. Optional `axes(:)`
selects one-based dimensions and allows a compact footprint whose rank is the
number of selected axes. As in SciPy's generic filter, compact footprint
shapes correspond to the **ascending selected axes**, whereas `origin(:)`
and `modes(:)` correspond to the **supplied axes order**. This distinction
matters when axes are reordered. Empty `axes` and empty `sizes` reduce each
one-element neighborhood to the original pixel. The complex callback takes
and returns `complex(dp)`; complex filtering with arbitrary callbacks and
mixed boundary modes is an additional SciFort capability.
For example:

```fortran
module user_callbacks
    use scifort_kinds, only : dp
    implicit none
contains
    real(dp) function range_of_values(values) result(value)
        real(dp), intent(in) :: values(:)
        value=maxval(values)-minval(values)
    end function range_of_values
end module user_callbacks

! Elsewhere in a program, with real(dp) arrays image and filtered:
! use user_callbacks, only : range_of_values
! use scifort_ndimage, only : generic_filter_nd
! call generic_filter_nd(image,range_of_values,filtered,status,sizes=[3,3])
```

`labeled_comprehension_nd` accepts a callback with the same signature and
runs it once for each **nonempty** ID in `index(:)`. Missing IDs receive
`default` without invoking the callback. Repeated IDs yield repeated results.
It preserves the C-order selection of pixels. SciPy 1.17 sorts grouped labels
with an **unstable** sort, so order-dependent callbacks can produce different
results, even though symmetric reductions (sum, extrema, variance, etc.) match.
For zero selected IDs, the output is empty and no callback is invoked. Only
real images/results are currently supported.

`distance_transform_cdt_nd` accepts optional positive `axis_weights(:)` for
weighted **taxicab** distances. This preserves the fast separable algorithm
and ranks 1–15; a chessboard metric cannot be combined with axis weights.
For general chamfer costs, use `distance_transform_weighted_cdt_nd` with an
integer stencil whose shape is `(3,3,...)` and rank equals that of the mask.
A positive weight at offset `(+1,0,...)` means that moving from a pixel to
its positive-axis neighbor costs that many units; zero forbids the move.
The center weight must be zero; **asymmetric/directed costs** are supported.
The implementation uses a multi-source binary-heap Dijkstra algorithm,
choosing the first Fortran-linear background site on equal-distance ties.
Unreachable foreground cells receive distance `-1` and all-zero witness
coordinates, while background sites receive distance zero and their own
one-based coordinates. `feature_indices` has shape `(rank(x),size(x))`.
Unlike SciPy's custom `distance_transform_cdt(metric=...)` array, which
encodes permitted neighbors, SciFort's custom array encodes **actual positive
integer costs**. The full-stencil algorithm caps the metric at 50,000 cells
and work at 50 million pixel/edge visits; weighted taxicab remains available
for every supported rank. These are resource guards, not distance definitions.

## Region measurements and marker-controlled watershed

The scalar-valued region measurements and extrema use a real(dp) image `x(..)` and an integer
`labels(..)` array of the same shape. The requested label IDs are given in
`index(:)`, including zero or negative IDs if desired; missing labels yield
sum 0 or NaN for mean, variance, and center of mass. Outputs for sums, means,
and variances are vectors with `size(index)` entries. `center_of_mass_nd`
returns a `(rank(x),size(index))` matrix, using **one-based Fortran pixel
coordinates** (add one to SciPy's zero-based result). Input mass may be
negative. Zero-total-mass regions produce NaN or signed infinity on each
coordinate, rather than raising divide-by-zero traps. Variances use population
normalization (dividing by the number of pixels) with stable online
accumulation. These interfaces support ranks 1 through 15; scanning currently
scales with the number of selected indices.

`extrema_labels_nd` combines SciPy-style minimum and maximum label statistics
with their positions, with individual wrappers for minimum/maximum values and
positions. For missing labels, values and coordinate vectors are filled with
zero. Each position array has shape `(rank(x),size(index))`; positions are
**one-based**. When equal extrema occur at several pixels, the first location
in C-order (axis 1 checked before axis 2, etc.) is returned, following SciPy.

`bounding_boxes_nd` returns an inclusive one-based `lower(rank,size(index))`
and `upper(rank,size(index))` pair for each requested label, with `(0,0,...,0)`
for a missing label. `find_objects_nd` returns the same bounds for all positive
labels `1..size(lower,2)` and ignores background and negative labels. Unlike
SciPy's zero-based Python slices, `upper` is **inclusive**, so a region can be
sliced in Fortran as `image(lower(1,j):upper(1,j),...)`. All these procedures
support ranks 1 through 15 and preserve repeated label IDs in `index`.

```fortran
use scifort_kinds, only : dp
use scifort_ndimage, only : sum_labels_nd, center_of_mass_nd, watershed_ift_nd
real(dp) :: intensities(6,7), masses(6,7), sums(2), centers(2,2)
integer :: markers(6,7), segmented(6,7), status
! Populate intensities, masses, and markers first. Then:
call watershed_ift_nd(nint(intensities),markers,segmented,status)
call sum_labels_nd(masses,segmented,[1,2],sums,status)
call center_of_mass_nd(masses,segmented,[1,2],centers,status)
```

Example of one-based selected-region bounds and masked neighborhood operations:

```fortran
use scifort_kinds, only : dp
use scifort_ndimage, only : bounding_boxes_nd, median_filter_nd, binary_dilation_nd
real(dp) :: image(6,7),filtered(6,7)
logical :: stencil(2,3),binary(6,7),editable(6,7),expanded(6,7)
integer :: labels(6,7),lo(2,2),hi(2,2),status

! Populate the image, footprint, masks, and integer region IDs first.
call median_filter_nd(image,y=filtered,status=status, &
    footprint=stencil,origin=[0,1])
call binary_dilation_nd(binary,expanded,status,iterations=0,mask=editable)
call bounding_boxes_nd(labels,[1,2],lo,hi,status)
! Region j is included in image(lo(1,j):hi(1,j),lo(2,j):hi(2,j))
! when lo(:,j) is nonzero.
```

`watershed_ift_nd` accepts **integer** intensities in the interval 0..65535,
corresponding to SciPy's `uint8` and `uint16` input types, and arbitrary
signed integer markers (zero means unlabeled). The image foresting transform
propagates each marker using the minimax of absolute intensity differences
between neighboring pixels. A supplied Boolean stencil must have extent
three on every input axis; otherwise the face-connected stencil is used.
Positive markers take priority over negative markers at equal costs. The
implementation supports ranks 1 through 15 and uses bucket queues rather
than recursively traversing the image. It enforces each axis' spatial bounds
by default. `legacy_linear_bounds=.true.` instead emulates the linearized
neighbor-boundary check in SciPy 1.17's C code; this can permit apparent
cross-edge connections, and tie labels may still differ from SciPy. Reference
checks verify the bounded implementation separately from the legacy option.

## Axis-selective binary morphology, complex line callbacks, and value indices (2026-10-23)

`binary_dilation_nd`, `binary_erosion_nd`, `binary_opening_nd`, and
`binary_closing_nd` accept optional one-based `axes(:)` alongside their
existing mask, iteration, border, and footprint arguments. When `axes` is
present, the Boolean `structure(..)` has **rank `size(axes)`** instead of the
full image rank. Its footprint dimensions correspond to selected image axes
in **ascending order**, regardless of `axes(:)` ordering; `origin(:)` follows
the order in which axes were supplied. With no structure, connectivity applies
to the selected axes. An empty `axes` array leaves the mask unchanged, and
repeated or out-of-range axes return `status=1`. The implementation embeds
small selected-axis offsets directly into full-rank coordinates, without
allocating an expanded rank-15 structuring element. As with earlier calls,
nonselected dimensions are preserved and each axis is one-based.

`generic_filter1d_nd` now overloads its callback for `complex(dp)` images,
with a corresponding public `ndimage_complex_line_func` abstract interface.
Its callback receives a padded complex line and fills a complex output line.
`cval` is complex for the new overload. Input/output shapes, boundary modes,
origin validation, line ordering, and the 100,000-sample safety limit are
identical to the existing real interface. SciPy 1.17 rejects complex input to
`generic_filter1d`, so the tests use independently decomposed real SciPy
callbacks as the numerical reference for the complex extension.

`value_indices_nd` returns the same region-coordinate information as
`scipy.ndimage.value_indices` but in a Fortran-friendly compressed form:

```fortran
integer :: labels(4,5),status,j
integer, allocatable :: values(:),starts(:),coordinates(:,:)
! Populate labels before calling:
call value_indices_nd(labels,values,starts,coordinates,status,ignore_value=0)
if (status==0) then
    do j=1,size(values)
        ! All one-based coordinates of label values(j):
        print *, values(j),coordinates(:,starts(j):starts(j+1)-1)
    end do
end if
```

`values(:)` lists distinct integer IDs in ascending order. `starts(:)` has
length `size(values)+1`, and `coordinates(:,:)` has shape `(rank(x),N)`;
`starts(j):starts(j+1)-1` selects all coordinates of key `values(j)`.
Coordinates are one-based and appear in C-order within each group to match
SciPy's occurrence order. Negative and zero label values are supported, and
`ignore_value` removes a specified key. If every pixel is ignored,
`values(:)` and `coordinates(:,:)` are empty and `starts(:)=[1]`.
The implementation uses an iterative heapsort with O(N log N) time and
O(N) additional storage, not a potentially enormous `max(label)` table.

`compare_complex_line_binary_axes_scipy.py` and
`compare_value_indices_scipy.py` validate the additions; the paired native
regressions include invalid arguments and high-rank inputs.

## Validation, limits, and provenance

All functions use `status=0` for success, `1` for invalid inputs, and `2`
for a kernel or neighborhood that exceeds an explicit safety limit. The
numerical filters cap kernel radius at 100,000 samples, full kernels and
order-statistic windows at 50,000 elements, and the morphology/grayscale stencil at
50,000 offsets. These are resource guards, not SciPy limitations.
Invalid calls leave the caller's output undefined; inspect status first.
Values containing NaNs or infinities are rejected by the numeric filters.

`test/test_ndimage_filters.f90`, `test/test_ndimage_morphology.f90`, and
`test/test_ndimage_extended.f90`, `test/test_ndimage_footprints.f90`,
`test/test_ndimage_advanced.f90`, `test/test_ndimage_regions_masks.f90`, and
`test/test_ndimage_chamfer_regions.f90`, and
`test/test_ndimage_callbacks_weighted.f90`
cover scalar/vector options, separability,
complex coefficients, status and resource checks, binary connectivity,
strided targets, and rank-15 inputs. Reproduce the SciPy 1.17 cross-checks:

```text
python tools/ndimage_reference/compare_ndimage_scipy.py
python tools/ndimage_reference/compare_morphology_scipy.py
python tools/ndimage_reference/compare_extended_filters_scipy.py
python tools/ndimage_reference/compare_extended_morphology_scipy.py
python tools/ndimage_reference/compare_distance_scipy.py
python tools/ndimage_reference/compare_distance_random_scipy.py
python tools/ndimage_reference/compare_footprint_feature_scipy.py
python tools/ndimage_reference/compare_feature_random_scipy.py
python tools/ndimage_reference/compare_advanced_gray_scipy.py
python tools/ndimage_reference/compare_binary_origin_scipy.py
python tools/ndimage_reference/compare_measurements_scipy.py
python tools/ndimage_reference/compare_watershed_scipy.py
python tools/ndimage_reference/compare_regions_masks_rank_scipy.py
python tools/ndimage_reference/compare_extrema_histogram_scipy.py
python tools/ndimage_reference/compare_chamfer_random_scipy.py
python tools/ndimage_reference/compare_callbacks_weighted_scipy.py
python tools/ndimage_reference/compare_weighted_stencil_random.py
python tools/ndimage_reference/compare_structure_axes_scipy.py
python tools/ndimage_reference/compare_axis_extrema_scipy.py
python tools/ndimage_reference/compare_complex_line_binary_axes_scipy.py
python tools/ndimage_reference/compare_value_indices_scipy.py
```

Run every reference check and the focused Fortran regression programs
(checked debug and optimized builds) with one command when GNU Fortran,
NumPy, and SciPy are installed:

```text
python tools/ndimage_reference/run_all_checks.py
```

The reference suite validates **818 distinct numerical/morphological cases**:
128 original filtering, 45 original binary morphology/labeling, 152 full-kernel
and order-statistic filtering (including Prewitt), 50 morphology composites
and flood-fills, and 147 exact distance transforms (including 135 randomized
shapes, masks, and anisotropic spacings), plus 168 custom-footprint,
grayscale, hit-or-miss, and feature-index cases, and 128 randomized exact
nearest-background feature maps. Eight filtering boundary
modes are exercised. The largest measured numerical-filter absolute error
was `7.11e-15`; the new EDT comparisons agreed within `2.22e-16`.
An additional 256 non-flat grayscale, 20 shifted binary, and 77 region
measurement SciPy comparisons pass. The watershed reference suite checks 36
bounded and 36 legacy-graph cases; 34 of the legacy results match SciPy 1.17
exactly, with two tie-label differences. The region/mask/footprint extension
adds another 8,526 successful 2D/3D SciPy comparisons, including negative and
zero iteration counts, C-order extrema ties, every filtering boundary mode,
asymmetric/even-sized footprints, and missing labels. These implementations are MIT-licensed
SciFort code; additional provenance is recorded in `CODE_PROVENANCE.md`.

The newest extension adds another **4,620 SciPy scalar comparisons** of
footprint/origin minimum and maximum filters plus labeled counts, histograms,
and population standard deviations, and **17,280 randomized SciPy discrete
chamfer distance values** (with nearest-background witness validation).
The seven focused Fortran regression executables pass with both checked
`-O0` and `-O3` GNU Fortran 14.2.0 builds. The 15 reference scripts pass
when run individually; the combined runner can exceed execution limits in
restricted environments. These additions require no new FPM dependencies.

The callback/weighted-metric extension additionally verifies **21,624**
SciPy/independent-reference elements (including real/complex 1-D kernels,
C-order generic neighborhoods, labeled callbacks, weighted axis distances,
and general weighted Dijkstra graphs) and **12,672** randomized directed-chamfer
distances/feature indices. The new native regression includes invalid-input
checks and rank-15 callbacks/anisotropic distances. No runtime dependencies
have been added.

## Mixed boundary modes and Fourier reference tests (2026-10-22)

`compare_mixed_complex_scipy.py` checks 19,584 numerical entries covering eight
boundary modes and their 2-D combinations for real/complex Gaussian,
uniform, extrema, full-kernel and generic filters, plus 3-D Gaussian and
uniform combinations. Mixed-mode multidimensional correlation/convolution
and complex generic reductions are additionally checked against independent
coordinate-space references because SciPy's public routines do not accept
those exact combinations.

`compare_fourier_scipy.py` checks 402 complex coefficients against SciPy
1.17.0 for real/complex FFT frequency conventions, including a nonterminal
packed real-FFT axis, odd real-FFT lengths, Nyquist frequencies, and
fractional shifts. `test_ndimage_mixed_fourier.f90` checks status validation,
composition identities, even-length Fourier frequencies, and rank-15 arrays.
The new implementations have no additional FPM dependencies.

## Fourier ellipsoid, mixed order modes, and generic line callbacks (2026-10-23)

`compare_ellipsoid_order_line_scipy.py` checks 9,134 numerical output
components: SciPy 1.17.0 Fourier ellipsoid cases in one, two, and three
dimensions (including packed real-FFT axes); independent Bessel-function
references for ranks 4 and 15; mixed-mode median, rank and percentile
filters in two and three dimensions; and `generic_filter1d` callbacks across
all eight boundary modes. Matching scalar-mode cases are also checked against
SciPy to validate the independent mixed-mode neighborhood reference.
`test_ndimage_ellipsoid_order_line.f90` adds invalid-input, asymmetric-origin,
rank-15, and exact Fourier-response tests. No new dependencies are required.

## Complex generic derivatives, axis-selective segmentation, and indexed region reductions (2026-10-23)

`generic_gradient_magnitude_nd` and `generic_laplace_nd` now overload on
`complex(dp)` images and user callbacks matching the public
`ndimage_complex_derivative_func` interface. A callback receives the complex
source, a **one-based** axis, same-shape complex destination, a boundary mode,
complex `cval`, and an integer status. Complex gradients follow SciPy's
**complex sum-of-squares** convention, `sqrt(sum(derivative**2))`; they are
*not* the real, Hermitian norm `sqrt(sum(abs(derivative)**2))`. The Laplacian
sums complex second derivatives. Optional `axes` and `modes` keep the same
conventions as real generic derivatives; a selected-axis mode vector is in
the order of `axes`.

`binary_propagation_nd` and `binary_fill_holes_nd` accept optional one-based
`axes(:)` and compact structures over selected dimensions. This extends the
SciFort API beyond SciPy's full-rank `binary_propagation` and
`binary_fill_holes` calls. A fill-holes boundary is evaluated **only along
selected axes**, so unselected dimensions define independent slices. An
empty `axes` list is the identity for both functions. Existing calls without
`axes` are unchanged.

The region reductions `sum_labels_nd`, `mean_labels_nd`, `variance_labels_nd`,
`standard_deviation_labels_nd`, `center_of_mass_nd`, `count_labels_nd`,
`histogram_labels_nd`, `extrema_labels_nd`, and `bounding_boxes_nd` now sort
requested label IDs and use binary search per image pixel. This changes
request lookup from `O(N*K)` to `O(N*log(K) + K*log(K))` for `N` pixels and
`K` requested labels, with only `O(K)` additional lookup storage. Repeated
indices retain their original result positions and values, including missing
labels and negative or large label IDs. `labeled_comprehension_nd` is now accelerated in the next release section.

The reference programs
`compare_complex_derivatives_labels_scipy.py` and
`compare_propagation_axes_scipy.py` compare complex derivatives and region
statistics directly to SciPy and selected-axis binary segmentation to
independent SciPy slices. The native regression tests include rank-15
arrays, noncontiguous output sections, repeated/missing labels, and invalid
arguments. No additional dependencies are required.

## Batched, vectorized image filtering (2026-10-23)

`vectorized_filter_nd` implements the scalar-output subset of
`scipy.ndimage.vectorized_filter` using **bounded-memory window batches**.
For example, with `real(dp) :: image(128,64), filtered(128,64)`:

```fortran
call vectorized_filter_nd(image, window_mean, filtered, status, &
    sizes=[3,5], mode='reflect', batch_size=1024)

! Place the callback inside a module or the caller's CONTAINS section.
subroutine window_mean(windows, results)
    real(dp), intent(in) :: windows(:,:) !! window value in row, output pixel in column
    real(dp), intent(out) :: results(:) !! one result per column
    results=sum(windows, dim=1)/real(size(windows,1),dp)
end subroutine window_mean
```

Two generic overloads accept real and complex input, callbacks, and output.
The callback receives `windows(neighbor, pixel)`, with active footprint
entries in NumPy/SciPy **C-order** within each window. Batch columns follow
**Fortran array storage order**, rather than the callback invocation order of
SciPy. Optional `footprint`, `origin`, `axes`, `mode`, `cval`, and `modes(:)`
follow the existing SciFort neighborhood-filter convention. A compact
footprint is mapped to `axes` **in the supplied order**. As elsewhere in
SciFort, axes are **one-based**.

`mode='valid'` returns only complete windows: for a width `w`, an input axis
of length `n` produces `n-w+1` output pixels, and the caller must allocate
`y` with the corresponding shape. Only selected axes shrink. Shifted origins
and mixed `modes(:)` are incompatible with `valid`. Unlike
`generic_filter_nd`, this operation permits NaNs in the source and constant
fill so callers can implement NaN-aware callbacks, as in SciPy's vectorized
interface. Scalar and multi-valued
callbacks are supported, including mixed real/complex output types as
described below.

`batch_size` is a positive **maximum number of output pixels** passed per
callback invocation; the implementation reduces it as needed so each
materialized window matrix contains at most two million elements. It does
not build a padded full-image array. The callback may be invoked multiple
times and should not rely on NumPy's C-order batch traversal or side effects.
The implementation is a general-purpose interface, not an automatic
acceleration for inexpensive reductions; specialized separable filters are
preferable for means, extrema, and Gaussian operations.

The new native test `test_ndimage_vectorized.f90` covers rank-15 arrays,
noncontiguous outputs, empty axes, mixed boundaries, NaN callbacks, batching
and input guards. `compare_vectorized_filter_scipy.py` checks real and
complex arrays, shifted origins, sparse reordered footprints, and `valid`
outputs against SciPy 1.17.0. `map_coordinates_grid_nd` is also now
re-exported directly from `scifort_ndimage`.

## Multi-valued vectorized filtering and interior-window acceleration (2026-10-23)

`vectorized_filter_multi_nd` adds a batched Fortran callback that returns
**several results from each neighborhood**, with real or complex input and
output. The callback takes `windows(neighbor,batch)` and fills
`results(channel,batch)`; channel count comes from the last dimension of
`y`. For example, with `real(dp) :: x(128,64), result(128,64,2)`:

```fortran
call vectorized_filter_multi_nd(x, mean_and_square, result, status, &
    sizes=[3,5], mode='reflect', batch_size=256)

! A callback in a CONTAINS section or a module.
subroutine mean_and_square(windows, results)
    real(dp), intent(in) :: windows(:,:) !! neighbors by pixel batch
    real(dp), intent(out) :: results(:,:) !! channels by pixel batch
    results(1,:)=sum(windows,dim=1)/real(size(windows,1),dp)
    results(2,:)=sum(windows**2,dim=1)/real(size(windows,1),dp)
end subroutine mean_and_square
```

For input shape `(n1,...,nr)` and `c` returned channels, the standard
output shape is `(m1,...,mr,c)`. The first `r` axes retain the input
axis order and reduce under `mode='valid'` as in the scalar callback API.
For **rank-15 input**, where an additional axis would exceed Fortran's
rank limit, use the flat output shape `(product([m1,...,mr]), c)` instead.
This flat shape is also accepted at lower ranks, with pixels in Fortran
array-element order. Results and input are permitted to be noncontiguous.
The procedure has the same `sizes`, `footprint`, `origin`, `mode`, `modes`,
`axes`, and `cval` options as `vectorized_filter_nd`.

The callback is invoked in bounded batches and sees footprint entries in
SciPy's C-order; batch columns follow Fortran pixel order. Its output
has shape `(channels, current_batch_size)` on every call. Mixed real-to-complex and complex-to-real callbacks are now supported for
both this API and the single-result API. More than one new non-image result
axis is not supported. The window-plus-channel buffers use
at most two million elements; the complete output is allocated before
scattering, so total result storage is still proportional to output size.

The original single-result `vectorized_filter_nd` and the new
multi-result filter now share a neighborhood gather routine. A fast path
precomputes linear stencil displacements for pixels whose full footprint
is in bounds; only border pixels require per-axis boundary indexing.
The operation is independent of the window statistic and preserves all
supported modes. A small GNU Fortran 14.2.0 `-O1` benchmark (`160 x 128`
pixels, `7 x 7` windows, four repeats) measured 0.0546 CPU seconds for two
single-result passes versus 0.0279 seconds for one two-channel pass,
approximately 1.96x faster. This comparison is workload-specific, not a
guaranteed universal speedup. Reproduce it with
`benchmarks/bench_ndimage_vectorized_multi.f90`.

The SciPy 1.17.0 reference exhibits an unusual output-axis permutation
with **multi-result** callbacks when `axes` is reordered. SciFort
intentionally retains `(...,channels)` rather than reproducing this
surprising layout; comparisons in that case use two independent SciPy
single-result calls followed by `np.stack`.

`test_ndimage_vectorized_multi.f90` exercises channel layout,
noncontiguous destinations, valid shapes, complex arrays, invalid
arguments, batching, and rank-15 flat output.
`compare_vectorized_multi_scipy.py` validates numerical parity with
SciPy, while the earlier scalar vectorized-filter comparisons protect the
shared neighborhood-extraction path.

## Mixed-type batched neighborhood callbacks (2026-10-23)

Both `vectorized_filter_nd` and `vectorized_filter_multi_nd` now support
**real input with complex output** and **complex input with real output**,
in addition to their original same-type variants. The new callback interfaces
are `ndimage_batch_real_to_complex_func`,
`ndimage_batch_complex_to_real_func`, and their `_multi_func` counterparts.
For example, computing a power estimate from a complex image:

```fortran
complex(dp) :: image(40,30)
real(dp) :: power(40,30)
integer :: status
call vectorized_filter_nd(image, average_power, power, status, &
    sizes=[3,5], mode='reflect')

! Define average_power in a CONTAINS section or a module.
subroutine average_power(windows, results)
    complex(dp), intent(in) :: windows(:,:) !! complex windows by output batch
    real(dp), intent(out) :: results(:) !! average squared magnitude
    results=sum(abs(windows)**2,dim=1)/real(size(windows,1),dp)
end subroutine average_power
```

The output type selects the overload and is not inferred at runtime from the
callback. The boundary `cval` uses the **input** type; `valid` mode and all
existing footprint, origin, selected-axis, and batching options are retained.
Multi-valued callbacks still fill `results(channel,batch)` and use a trailing
channel dimension or the rank-15 flat-output representation. For mixed-type
callbacks this is an intentional **SciFort extension** of the vectorized
Fortran interface, with semantics numerically verified against SciPy 1.17's
vectorized filtering and custom NumPy reductions. The reference script
`tools/ndimage_reference/compare_vectorized_mixed_scipy.py` verifies 2,016
values over 26 configurations, and the native regression
`test/test_ndimage_vectorized_mixed.f90` covers noncontiguous outputs, invalid
arguments, `valid` windows, and rank-15 arrays.

## Callback geometric transforms, batched labels, and axis-selective labeling (2026-10-23)

`geometric_transform_nd(x, mapping, y, order, status, mode, prefilter, cval)`
interpolates at source coordinates returned by a user-supplied mapping:

```fortran
subroutine mapping(output_coords, input_coords)
    integer, intent(in) :: output_coords(:) !! one-based coordinates in y
    real(dp), intent(out) :: input_coords(:) !! one-based coordinates in x
end subroutine mapping
```

The callback is invoked once per output pixel, in **C coordinate order** (last
output axis varying fastest). The input and output may have **different ranks**
within 1..15. Coordinates in both arrays are one-based, unlike SciPy's
zero-based callback contract. The default interpolation is cubic (`order=3`),
with `mode='constant'`, `cval=0`, and `prefilter=.true.`; users may request
orders 0 through 5 and all eight existing interpolation boundary modes.
Real and complex overloads are supported. A guard caps the allocated
coordinate table at four million real entries: exceeding this limit returns
status 2 before invoking the callback. Nonfinite coordinates are rejected.
Higher-order `reflect` and `grid-mirror` interpolation on short input axes
retains the previously documented spline prefilter discrepancy; the new
reference suite bounds this separately rather than claiming exact agreement.

`labeled_comprehension_nd` now builds an indexed bucket of selected pixels
in one C-order traversal instead of rescanning the entire input per requested
label. Sorting selected IDs and binary-searching the sorted requests gives
**O(N log K + K log K)** grouping time plus callback work, using O(N+K)
storage, for N pixels and K requested labels. A callback still runs **once
per nonempty request**, including repeated IDs (side effects are preserved).
Missing labels return the supplied default without invoking the callback.
Callback pixel order remains the previous deterministic C coordinate order,
which may differ from SciPy's internal grouping order for order-sensitive
callbacks. Negative, large, missing and duplicated label IDs are accepted.

`label_nd` accepts `axes(:)` for logical, integer and real foreground masks.
This is a **SciFort extension**: SciPy's `ndimage.label` has no `axes`
parameter. Connected components are restricted to the specified one-based
axes; all other axes act as independent slices, with **globally unique**
component IDs across slices. When supplied, the symmetric compact
`structure` footprint has the selected rank and its dimensions correspond
to the selected axes in ascending numerical order. Duplicate or out-of-range
axes are rejected, while an empty axes list treats every foreground pixel as
a separate component. Calls without `axes` retain their original behavior.

An optional performance experiment,
`benchmarks/bench_ndimage_labeled_comprehension.f90`, compares bucketed
reductions against an independent repeated-scan implementation. On this test
machine, 80,000 pixels and 1,200 requested labels produced a 13.1x speedup
at `-O2` (0.00468 vs 0.0613 seconds per iteration); this is a single workload,
not a general performance guarantee.

The reproducible checks are `compare_geometric_comprehension_scipy.py`,
`compare_label_axes_scipy.py`, `test_ndimage_geometric_comprehension.f90`, and
`test_ndimage_label_selected_axes.f90`. The independent SciPy slice-partition
checks do not assume arbitrary component-ID numbering. No new dependencies.
