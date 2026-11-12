# Decimation with caller-supplied FIR or IIR filters

The `decimate` interface now allows a user-designed antialias filter, rather
than requiring one of SciFort's standard Hamming FIR or Chebyshev-I IIR
prototypes. These are *trailing optional keywords*, so all existing positional
calls remain compatible:

| Keyword | Type and shape | Meaning |
| --- | --- | --- |
| `taps` | `real(dp) :: taps(:)` | Finite FIR coefficients, in causal order; at least one coefficient |
| `sos` | `real(dp) :: sos(:,:)` | Real IIR sections of shape `(nsections,6)`, with each `sos(:,4) == 1.0_dp` |

Exactly one of `taps=` and `sos=` can be supplied. When supplying either,
omit `n=` and `ftype=`; `zero_phase=` remains available. The new arguments
are supported with real and complex signals in the following APIs:

- `decimate(x, q, ..., taps=..., sos=...)` for vectors, and with an `axis=`
  keyword for rank-two and rank-three arrays;
- `call decimate_into(x, y, q, ..., taps=..., sos=...)` for caller-allocated
  rank-one through rank-15 arrays with an `axis=` keyword.

For example:

```fortran
use scifort_signal, only : dp, decimate, decimate_into, butter, butter_result
implicit none
real(dp) :: x(100), h(5)
type(butter_result) :: filter
real(dp) :: multichannel(3,100), downsampled(3,25)
real(dp), allocatable :: y(:)
integer :: status, i

x = [(sin(0.1_dp*real(i,dp)), i=1,100)]
multichannel(1,:) = x
multichannel(2,:) = 2.0_dp*x
multichannel(3,:) = -x
h = [0.05_dp, 0.2_dp, 0.5_dp, 0.2_dp, 0.05_dp]
y = decimate(x, 4, taps=h, status=status)
y = decimate(x, 4, taps=h, zero_phase=.false., status=status)

! A caller can also use a previously designed stable SOS filter:
filter = butter(4,0.24_dp)
! filter%sos(n,1:6) holds (b0,b1,b2,a0,a1,a2), with a0=1.
call decimate_into(multichannel, downsampled, 4, sos=filter%sos, &
                   axis=2, status=status)
```

**Filtering conventions:** For FIR designs, `zero_phase=.true.` (the default)
uses delay-compensated polyphase filtering, equivalent to SciPy's
`signal.decimate(..., ftype=signal.dlti(taps, [1.0]))` with zero phase.
`zero_phase=.false.` uses causal `upfirdn` followed by downsampling. For
SOS designs, the zero-phase default uses `sosfiltfilt` and causal mode uses
`sosfilt`, with downsampling by `q` after filtering. Output length is
`ceiling(size(x,axis)/q)` in both modes. A real SOS can be obtained from
SciFort's existing IIR filter-design routines (for example, `butter` and
`cheby1`) or supplied directly.

**Validation:** Both FIR and SOS arrays must contain finite coefficients.
Invalid shapes, non-normalized `a0` coefficients, ambiguous filter arguments,
nonpositive downsampling factors, and signals too short for zero-phase IIR
padding return `resample_status_invalid_input`. Vector and rank-2/3
array-returning calls return empty arrays; `decimate_into` preserves its
caller-provided output on failure. As with normal IIR filtering, users are
responsible for choosing a stable IIR filter and suitable antialias cutoff.

**Compatibility limits:** SciPy accepts an arbitrary `dlti` object in `ftype`;
Fortran instead accepts its FIR coefficients or a *real*, normalized SOS
representation explicitly. Complex SOS coefficients are now supported by
additional `sos=` overloads with complex output, as described in
[`signal_decimate_complex_sos.md`](signal_decimate_complex_sos.md). Numeric
transfer-function `(b,a)` decimators, including complex IIR coefficients,
are also available via `b=` and `a=` as described in
[`signal_decimate_transfer_function.md`](signal_decimate_transfer_function.md).
Python DLTI objects themselves are not accepted. SciFort's zero-phase SOS algorithm also uses its existing padding
rules, which can differ from a SciPy `dlti` object that is converted to a
different sequence of SOS sections.

## Validation

From the project root, with GNU Fortran, NumPy, and SciPy installed:

```sh
python tools/run_signal_decimate_custom_filters.py
```

The independent SciPy 1.17 comparison includes 336 real and complex
custom-filter cases / 10,500 decimated samples. The maximum absolute
difference observed was `1.288e-14`. The native suite checks 1D, rank-2,
rank-3, rank-5 and rank-15 arrays, noncontiguous output slices, both phase
options, inconsistent parameters,
invalid coefficient shapes, negative-axis handling and output preservation.
The standalone checked build and targeted `-O3` build passed. The existing
53-case multidimensional and 184-case complex multirate SciPy comparison
suites also passed in both builds. Unmodified large dependencies were compiled
in checked mode for the optimized build. The complete `fpm test` suite has not
been run because FPM is unavailable in the validation environment.

This code reuses SciFort's independent FIR and SOS filtering engines. SciPy
is a test-only reference dependency; no new library or runtime dependency is
introduced.
