# Parallel discrete-time SOS simulation (`dlsim_sos_bank`)

`dlsim_sos_bank` extends the scalar SOS simulator to **independent, parallel
channels** while keeping the filters in second-order-section representation.
This is a **SciFort extension**, not a separate public SciPy function. It uses
the same transposed-direct-form-II recurrence as the one-dimensional
`dlsim_sos` implementation and avoids forming high-order polynomials.

## Public Fortran API

```fortran
use scifort_signal, only : dp, dlsim_sos_bank, discrete_sos_bank_result

type(discrete_sos_bank_result) :: simulation

! Shared SOS filter: all channels use the same cascade.
simulation = dlsim_sos_bank(sos, inputs, dt=0.01_dp)

! Independent SOS filters, one cascade per channel.
simulation = dlsim_sos_bank(sos_per_channel, inputs, dt=0.01_dp)
```

`dlsim_sos(sos,inputs,dt,t,zi)` and the generic
`dlsim(sos,inputs,dt,t,zi)` accept these same new rank-two input overloads,
without changing the existing SISO overloads. All eight combinations of
shared/per-channel SOS filters and real/complex SOS/input arrays are supported.

- `inputs(nsamples,nchannels)` stores one time series in each column.
- **Shared filter:** `sos(nsections,6)` applies the same filter to every
  channel. **Per-channel filter bank:** `sos(nsections,6,nchannels)` provides a
  different cascade for each channel. The section count must be the same
  across all channels; different section counts can be handled by separate
  calls or padding with identity SOS sections.
- Each section is normalized as `[b0,b1,b2,1,a1,a2]`, matching SciFort
  `sosfilt` and SciPy's SOS order. All dimensions must be nonzero.
- `dt` is a positive finite output sampling interval. With no `t` argument,
  each input row is sampled at `(i-1)*dt` and the result has exactly
  `nsamples` rows.
- Optional `t(nsamples)` contains finite, strictly increasing sample times.
  As in the existing SISO `dlsim_sos`, the output time grid is
  `linspace(0,t(end),floor(t(end)/dt)+1)` and each channel is independently
  linearly interpolated to that grid. An explicitly provided time vector
  requires at least two input rows.
- Optional `zi(nsections,2,nchannels)` holds per-channel initial SOS delay
  states. If not provided, they are zero. A real `zi` is used for real
  filters **and** real signals; otherwise `zi` must be complex.

The result has type `discrete_sos_bank_result`:

| Component | Shape | Meaning |
|---|---|---|
| `%t` | `(nt)` | Output time grid |
| `%y` | `(nt,nchannels)` | Complex-valued filtered signal, also for real inputs |
| `%x` | `(nt,2*nsections,nchannels)` | Filter states **before** each input sample |
| `%zf` | `(nsections,2,nchannels)` | Final filter states, suitable for the next call's `zi=` |
| `%status` | scalar integer | `discrete_sim_success`, `discrete_sim_invalid_input`, or `discrete_sim_numerical_failure` |

In `%x`, state coordinates are interleaved by section, i.e.
`[sec1_delay1,sec1_delay2,sec2_delay1,sec2_delay2,...]`.

## Streaming without resetting filter history

```fortran
use scifort_signal, only : dp, dlsim_sos_bank, &
    discrete_sos_bank_result, discrete_sim_success

real(dp) :: sos(2,6), input(500,3)
type(discrete_sos_bank_result) :: first, second

! Initialize sos and input; each of the three columns is one channel.
first = dlsim_sos_bank(sos, input(:250,:), dt=0.02_dp)
if (first%status /= discrete_sim_success) error stop 'first block'

second = dlsim_sos_bank(sos, input(251:,:), dt=0.02_dp, &
                        zi=real(first%zf,dp))
if (second%status /= discrete_sim_success) error stop 'second block'
```

`%zf` preserves the filter state across blocks. The example's local time
arrays both start at zero; assemble an absolute time axis separately for
streamed blocks. The example uses a real initial state array because both
filters and samples are real. For complex signals or filters, pass `%zf`
without the `real(...)` conversion.

## Validation and error behavior

The routine checks SOS dimensions/normalization, input channel count, all
coefficients, samples, optional initial states, sample interval and timestamp
array. An invalid request returns `discrete_sim_invalid_input` with `%y`
and `%zf` unallocated. Nonfinite numerical trajectories are reported as
`discrete_sim_numerical_failure` rather than success. The input arrays are
not overwritten. Strided/noncontiguous array sections are accepted.

Run independent numerical comparisons and existing regressions with:

```bash
python tools/run_signal_dlsim_sos_bank.py
```

The reference driver generates 256 parameter combinations covering shared
and per-channel SOS, one through 24 sections, one and three channels,
real/complex combinations, initial states and irregular timestamps. It
compares 9,808 output samples and 171,680 pre-sample state components against
**independent per-channel calls to SciPy 1.17 `signal.sosfilt`**, including
final SOS states. Maximum measured absolute error was `7.22e-16` in both
checked and targeted `-O3` builds. The runner also repeats the previous
240-case arbitrary-input SOS and 288-case SOS step/impulse reference suites.

This code is original MIT-licensed SciFort work. NumPy and SciPy are used
only by the optional comparison scripts, not at runtime.
