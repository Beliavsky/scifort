# Multidimensional relative extrema (SciFort extension)

SciFort provides SciPy-compatible comparisons along a selected dimension of a
real signal through rank 15. Unlike `scipy.signal.argrelextrema`, which returns
a tuple of coordinate arrays, these **SciFort extension** procedures write a
logical mask of the same shape as the input:

```fortran
use scifort_signal, only : dp, argrelmax_mask_into, &
    argrelmin_mask_into, argrelextrema_mask_into, extrema_status_success
real(dp) :: samples(3, 128, 4)
logical :: peaks(3, 128, 4)
integer :: status

! Detect maxima separately in each time series along dimension two.
call argrelmax_mask_into(samples, peaks, axis=2, order=3, &
    mode='clip', status=status)
if (status /= extrema_status_success) error stop 'invalid extrema options'
```

`argrelmax_mask_into` and `argrelmin_mask_into` require strict greater-than
and less-than comparisons. `argrelextrema_mask_into(x, mask, comparator, ...)`
accepts a user-supplied pure scalar logical comparator with arguments of type
`real(dp)`. Custom comparators need not be strict, so they can retain tied
values; for example, a greater-than-or-equal comparison returns a true mask
for a constant singleton axis.

The optional arguments are `axis` (one-based, default 1; negative axes count
from the end), `order` (number of neighbors to test on *both* sides, default 1,
positive), `mode` (`'clip'` by default or `'wrap'`) and `status` (integer zero
on success). For the same physical axis as SciPy, use `axis = scipy_axis + 1`
for nonnegative indices. If `mode='clip'`, out-of-range neighbors repeat the
endpoint; if `'wrap'`, they wrap periodically. The output shape and rank must
exactly match the input. The input may be strided, and so may the output mask.

No assumptions are made about finiteness. NaNs compare according to the
chosen comparator (strict comparisons with a NaN are false). Empty input
arrays succeed and yield empty output masks. Invalid rank, shape, axis,
`order` or `mode` yields `extrema_status_invalid_input`, leaving the entire
caller-provided mask unchanged. Each call stages results before updating the
mask; it does not modify the input signal.

The implementation currently copies inputs to contiguous Fortran-order
storage and allocates a logical mask of the same element count, using O(N)
additional memory and O(N * min(order, axis_length)) comparisons. The legacy
`argrelmax`, `argrelmin` and `argrelextrema` remain available for rank-one
signals and return one-based integer index arrays.

## Verification

`test/test_signal_extrema_nd.f90` covers errors, custom comparators,
rank 15, empty and noncontiguous arrays. The independent reference command:

```sh
python tools/signal_reference/verify_extrema_nd_scipy.py
python tools/signal_reference/verify_extrema_nd_scipy.py --optimized
```

compares 720 deterministic seeded cases against SciPy 1.17 per build. The
standalone public-module test can be selected with:

```sh
python tools/run_signal_standalone.py checked ./build_signal_checked --match extrema
```
