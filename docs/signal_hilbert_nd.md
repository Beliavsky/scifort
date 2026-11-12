# Rank-independent Hilbert analytic signals

SciFort's public `scifort_signal` module exports two caller-output routines:

- `hilbert_into(x, y, axis, n, status)`: the 1-D analytic signal along a
  selected axis of a real `real(dp)` array, rank 1 through 15.
- `hilbert2_into(x, y, axes, n, status)`: the 2-D (single-orthant) analytic
  signal along **two** selected axes of a real array, rank 2 through 15.

The analytic signal is `complex(dp)`. Both routines use the same FFT frequency
masks and normalization as SciFort's existing vector `hilbert` and matrix
`hilbert2` functions. In particular, the two-dimensional mask follows SciPy
1.17's convention of **suppressing an even-length Nyquist bin** along either
axis; the one-dimensional mask retains its Nyquist bin. The two-axis transform
is not generally equal to independently calling the one-dimensional Hilbert
transform and summing the results.

```fortran
use scifort_signal, only : dp, hilbert_into, hilbert2_into
real(dp) :: samples(4,128,3)
complex(dp) :: analytic(4,256,3), analytic2(4,256,96)
integer :: st

call hilbert_into(samples, analytic, axis=2, n=256, status=st)
call hilbert2_into(samples, analytic2, axes=[2,3], n=[256,96], status=st)
```

## Argument conventions

The `axis` argument defaults to the last dimension and `axes` defaults to the
last two dimensions in increasing order. **Positive axes are one-based**, as
in other SciFort interfaces; negative axes count from the end. The two axes
of `hilbert2_into` must differ. Explicit transform length `n` is a scalar for
`hilbert_into` and a two-element integer array for `hilbert2_into`, in the
same order as `axes`. If omitted, lengths equal the input extents.

A transform length may exceed or be less than the input axis extent, padding
with trailing zeros or truncating that axis, respectively. The output must
have the **same rank** as the input, with each transformed extent equal to the
requested length and all other extents unchanged. Empty channel/batch extents
are accepted when the transformed extents are positive.

Inputs must be real and finite. The output is allocated by the caller and may
be noncontiguous. If an axis, transform length, shape, or input value is
invalid, `status` is `hilbert_status_invalid_input` and **none of the output
is changed**. Success reports `hilbert_status_success=0`. The routines make
one flattened staging copy before committing results; peak temporary memory
scales with the total requested result size. Reusable FFT plans are local to
one transform call and are not shared across threads.

The 2-D axis-selection capability, assumed-rank caller-output API, and
transactional failure behavior are **SciFort extensions**. SciPy's
`scipy.signal.hilbert2` transforms the last two axes; SciPy's `hilbert`
accepts a single zero-based axis. The older SciFort `hilbert` and `hilbert2`
array-returning entry points are unchanged.

## Validation

`test/test_signal_hilbert_nd.f90` compares the new operations against the
existing rank-one/rank-two kernels and tests rank 15, negative axes,
noncontiguous array sections, even-length Nyquist bins, and invalid inputs.
`test/test_signal_hilbert_nd_public.f90` exercises both new interfaces from
`scifort_signal`.

The independent reference pair
`validation/reference_drivers/signal_hilbert_nd_reference_driver.f90` and
`test/compare_signal_hilbert_nd_scipy.py` compares 500 seeded SciPy 1.17
cases and 191,284 complex output values across ranks 1, 2, 3, 5, and 15,
including arbitrary nonadjacent axes and zero padding/truncation:

```sh
python test/compare_signal_hilbert_nd_scipy.py /path/to/signal_hilbert_nd_reference_driver
```

The comparison uses SciPy's `hilbert2` on the last two dimensions after moving
selected axes there, then restores their original ordering.
