# Complex second-order-section filters for decimation

`decimate` and `decimate_into` accept **complex-valued SOS coefficients** using
`sos=`. The implementation reuses SciFort's existing complex `sosfilt` and
`sosfiltfilt` kernels and supports real or complex signals with complex output.
The previously supported real `sos=` argument, FIR `taps=`, numeric `b=`/`a=`
transfer functions, and built-in filter designs remain unchanged.

The SOS layout is `(nsections,6)` with one row per section:
`[b0,b1,b2,a0,a1,a2]`. Each `a0` must be exactly `(1.0_dp,0.0_dp)`;
all coefficient components must be finite and at least one section is
required. Convert non-unit denominators to normalized form before calling.
As in SciPy's SOS filters, complex coefficients are accepted in causal and
forward-backward filtering, but users remain responsible for filter stability
and adequate antialias attenuation.

```fortran
use scifort_signal, only : dp, decimate, decimate_into
implicit none
real(dp) :: x(100), channels(2,100)
complex(dp) :: sections(1,6), filtered(2,34)
complex(dp), allocatable :: y(:)
integer :: status

sections(1,:) = [cmplx(0.25_dp,0.05_dp,dp), &
                 cmplx(0.10_dp,-0.02_dp,dp), &
                 cmplx(0.04_dp,0.0_dp,dp), &
                 (1.0_dp,0.0_dp), &
                 cmplx(-0.45_dp,0.04_dp,dp), &
                 cmplx(0.12_dp,-0.01_dp,dp)]
y = decimate(x,3,sos=sections,status=status)
call decimate_into(channels,filtered,3,sos=sections, &
                   axis=2,zero_phase=.false.,status=status)
```

Both forms return `ceiling(size(x,axis)/q)` samples along the selected axis.
`q` must be positive. The default `zero_phase=.true.` applies
`SOSFILTERFORWARD -> reverse -> SOSFILTERBACKWARD`, with SciFort's existing
odd-padding/steady-state initial conditions, before taking every `q`th
sample. `zero_phase=.false.` applies a causal SOS cascade, then downsamples.
Very short signals may not permit the default zero-phase padding.

Supported APIs:

- `decimate(x,q,sos=...)` for rank-one real or complex `x`, with a complex
  vector result. Rank-two and rank-three functions also accept `axis=`, using
  one-based positive axes or negative indices (default last axis).
- `decimate_into(x,y,q,sos=...)` for rank-one through rank-15 assumed-rank
  real or complex input `x` with a caller-provided **complex** output `y`.
  The output must have the same rank, unchanged other extents, and selected
  dimension of length `ceiling(size(x,axis)/q)`. Strided and noncontiguous
  slices are supported.

The new complex-SOS overloads require `sos=` and support `zero_phase=`,
`axis=` where appropriate, and `status=`. To specify additional filter-design
arguments (`n=` or `ftype=`), use the corresponding existing designed-filter
call instead; Fortran's generic interface rejects these combinations with
complex-SOS overloads at compile time. The existing same-type overloads
continue to reject ambiguous `taps=` / real `sos=` / `b=` / `a=` combinations
at runtime.

Invalid shapes, nonfinite coefficients or input, nonunit `a0`, invalid `q`
or `axis`, or insufficient samples for zero-phase filtering return
`resample_status_invalid_input` and an empty array for the array-returning
functions. The assumed-rank `decimate_into` procedure preserves all elements
of the caller's output on failure.

## Validation

Run from the project root with GNU Fortran, SciPy 1.17 and NumPy installed:

```sh
python tools/run_signal_decimate_complex_sos.py
```

The independent reference compares SciFort against public SciPy
`signal.sosfilt` and `signal.sosfiltfilt` results followed by `[::q]`.
It covers 288 seeded configurations (7,434 complex-valued samples), real or
complex signals, real or complex SOS values, one to three cascaded sections,
multiple downsampling factors, two phase choices, and first-order sections.
All comparisons passed; the maximum absolute difference observed was zero
at the printed binary64 precision for this test matrix. Native checks cover
rank-two/three and rank-15 arrays, negative/positive axes, noncontiguous
slices, invalid inputs, and transactional output preservation.

The checked build and targeted `-O3` build of both modified modules passed,
as did the earlier 576-case transfer-function, 336-case custom-filter,
53-case multidimensional and 184-case complex-resampling comparisons.
The full FPM suite was not run in the validation environment.

SciPy is a test-only dependency, not an additional Fortran runtime dependency.
