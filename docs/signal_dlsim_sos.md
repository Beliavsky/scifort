# Arbitrary-input simulation using digital SOS cascades

`dlsim_sos` is a **SciFort extension** to the `scipy.signal.dlsim`-style
API. SciPy 1.17 accepts transfer-function, ZPK, or state-space LTI models
in `dlsim`, but has no public SOS-native discrete-time simulation entry
point. SciFort now simulates a sampled SISO input directly through a cascade
of digital second-order sections (SOS), without expanding the cascade into
high-order transfer-function polynomial coefficients.

## Interfaces

The public `scifort_signal` and `scifort_signal_discrete_simulation` modules
export:

```fortran
simulation = dlsim_sos(sos, u, dt, t, zi)
simulation = dlsim(sos, u, dt, t, zi)
```

The second call is an overload of the existing `dlsim` generic.

- `sos(nsections,6)` contains `[b0,b1,b2,1,a1,a2]` per row, in the same
  order and normalization as `scipy.signal.sosfilt`. `nsections >= 1`.
- `u(:)` is the input sampled at intervals `dt` (if `t` is absent) or at
  the timestamps given by `t`. It must not be empty.
- `dt` is the strictly positive finite *output* sampling interval.
- `t(:)` is optional, must match `size(u)` and be finite and strictly
  increasing. For this option at least two samples are required.
- `zi(nsections,2)` is an optional initial transposed direct-form-II
  state matrix. If omitted, initial filter states are zero. Real filters
  and real inputs use a real `zi`; when either is complex, a complex `zi`
  is accepted. State pairs are in the same order as in `sosfilt`.

All four combinations of real and complex filters/signals are supported.
Output values and states use the complex-valued
`type(discrete_simulation_result)`, even for real calculations.

```fortran
use scifort_signal, only : dp, butter, butter_result, dlsim_sos, &
    discrete_simulation_result, discrete_sim_success

real(dp) :: samples(256)
type(butter_result) :: filter
type(discrete_simulation_result) :: result

! Populate samples ...
filter = butter(8,0.2_dp)
result = dlsim_sos(filter%sos,samples,dt=0.01_dp)
if (result%status /= discrete_sim_success) error stop 'SOS simulation failed'
! result%y(:,1): filtered signal
! result%x(:,1:): pre-sample filter states
! result%zf(:,:): SOS states after the last sample
```

## State layout and streaming

`%t(nt)` contains output times; `%y(nt,1)` contains SISO outputs;
`%x(nt,2*nsections)` contains states *before* each sample is applied,
interleaved by section:
`[section1_delay1, section1_delay2, section2_delay1, section2_delay2,...]`.
`%zf(nsections,2)` contains states **after** the final sample, in `sosfilt`
shape. Feed `%zf` to a subsequent call as `zi` to resume filtering without
resetting the filter. For real SOS+input, use `real(result%zf,dp)` when calling
the real-input overload; for a mixed or complex case `%zf` can be supplied
directly. `%zf` is unallocated for the existing transfer-function, ZPK and
state-space variants of `dlsim`.

The state recurrence for each section (with input `v` and delays `z1,z2`) is

\[
  y=b_0 v+z_1, \quad z_1'=b_1 v-a_1 y+z_2, \quad z_2'=b_2 v-a_2 y.
\]

This is the transposed direct-form-II filter used in SciFort's independent
`sosfilt` implementation. High-order SOS cascades are evaluated section by
section, avoiding high-order polynomial conversion.

## Irregular input timestamps

When `t` is supplied, output times are

`linspace(0,t(end),floor(t(end)/dt)+1)`.

The signal is linearly interpolated (or linearly extrapolated at the edges)
from the input timestamp grid onto that output grid, matching SciFort's
existing state-space `dlsim` convention. With no `t`, the output has exactly
`size(u)` samples at `(k-1)*dt`. Time spacing changes only the time labels
when `t` is omitted, not the discrete filter coefficients or response.

Invalid SOS dimensions, nonunit `a0`, nonfinite coefficients or inputs,
invalid `dt`/`t`, and malformed/nonfinite `zi` return
`discrete_sim_invalid_input`, leaving the result arrays unallocated.
Nonfinite computed trajectories return `discrete_sim_numerical_failure`.

## Reproducible validation

```bash
python tools/run_signal_dlsim_sos.py
```

The runner independently compares 240 real/complex and high-section-count
simulation cases against SciPy 1.17 `signal.sosfilt` (9,120 output samples
and 310,080 state-history values); it also reruns the previous 288-case SOS
step/impulse reference suite. These 240 cases include optional initial states,
input timestamp interpolation, 1--64 SOS sections, and four real/complex
combinations. Native Fortran tests cover continuation, noncontiguous inputs,
generic overloading, invalid arrays, and numerical overflow reporting.

The `-O3` stage compiles the modified discrete-simulation module and test
programs at `-O3`, with unchanged dependencies built with runtime checks.
The standalone runner does not replace `fpm test`.

The implementation and tests are original MIT-licensed SciFort code. SciPy
1.17 and NumPy are numerical references only, not runtime dependencies.
