# SOS-native discrete impulse and step responses

SciFort extends its `scipy.signal.dstep` / `scipy.signal.dimpulse`-style
API with **direct digital second-order-section (SOS) responses**.  This
is a **SciFort extension**: SciPy 1.17 does not accept an SOS matrix directly
in its public `dstep`/`dimpulse` functions. The response is computed using
SciFort's existing transposed-direct-form-II `sosfilt` kernel, section by
section, rather than expanding an entire SOS cascade to one high-degree
numerator/denominator pair. This avoids a common source of catastrophic
coefficient cancellation for high-order digital filters.

## Fortran API

The following procedures accept real or complex `sos(:,:)` arrays and return
`type(discrete_response_result)`:

- `dstep_sos(sos, dt, n, zi)` and `dimpulse_sos(sos, dt, n, zi)`
- `dstep(sos, dt, n, zi)` and `dimpulse(sos, dt, n, zi)` are **overloads**
  of the existing public generic names, not replacements.

Every SOS row is `[b0,b1,b2,a0,a1,a2]`, with **`a0=1`**, as in
`scipy.signal.sosfilt`. Numerator and denominator coefficients may be complex;
real SOS matrices yield real-valued responses stored in the existing complex
result container. `dt` is a positive finite sampling interval. The optional
`n` gives the sample count (default 100). The optional `zi` is a real/complex
array of shape `(nsections,2)` containing **per-section transposed-DF-II
initial delay states**, rather than the controller-canonical `x0` vector
used by `dstep` for transfer functions and state-space models.

```fortran
use scifort_signal, only : dp, butter, butter_result, dstep_sos, dimpulse, &
    discrete_response_result, discrete_sim_success

type(butter_result) :: filt
type(discrete_response_result) :: step_response, impulse_response

filt = butter(10, 0.25_dp)
step_response = dstep_sos(filt%sos, dt=0.01_dp, n=200)
impulse_response = dimpulse(filt%sos, dt=0.01_dp, n=200)

if (step_response%status /= discrete_sim_success) error stop 'bad SOS response'
! Time:       step_response%t(k)        = (k-1)*0.01_dp
! Output:     step_response%y(k,1,1)   = step response at k
```

The result uses the same existing `discrete_response_result` as other
`dstep`/`dimpulse` overloads: real `%t(n)` and complex `%y(n,1,1)`.
The discrete unit impulse is **1 at sample zero**, not `1/dt`;
therefore the sample values depend only on SOS coefficients and the input
sequence, not on the time-grid spacing. An initial-state-free step response
is the cumulative sum of the impulse response, subject to floating-point
roundoff.

Nonfinite `dt`, `dt <= 0`, `n < 1`, nonfinite SOS entries, incorrect matrix
shape, nonunit `a0`, and wrong/invalid initial state arrays yield
`discrete_sim_invalid_input` without allocating results. Nonfinite filter
outputs yield `discrete_sim_numerical_failure`. No `t` or `x0` keywords are
accepted in SOS overloads; those retain their existing state-space meanings.

## Validation

Reproduce with:

```bash
python tools/run_signal_sos_responses.py
```

The clean runner builds with `-fcheck=all` and with `-O3` using GNU Fortran,
executes a native regression program, and compares 288 independently generated
real/complex SOS cases against SciPy 1.17 `sosfilt`, including random stable
biquads, Butterworth filters up to **64 sections**, unit impulses, unit
steps, optional per-section initial states, noncontiguous filter arrays,
and several sample spacings. Low-order, zero-state cases are additionally
checked against SciPy's public discrete-time `dstep`/`dimpulse` after
conversion to a transfer function (only where that conversion is stable).
Both builds pass 15,153 compared response samples, with zero observed
absolute difference at the comparison's reported precision. The native
regression tests also check invalid inputs and the public overloads.

The full FPM suite has not been run by this standalone runner.
