# SOS group delay (2026-10-28)

SciFort now provides **`group_delay_sos(sos, worn, whole, fs, w)`**, a
numerically stable method for calculating group delay of a digital filter
represented as second-order sections (SOS). The `group_delay` generic also
accepts a two-dimensional SOS argument, without changing existing calls
with one-dimensional transfer-function coefficients.

This is a **SciFort extension**: SciPy 1.17's `scipy.signal.group_delay`
accepts `(b, a)` coefficient pairs rather than a matrix of SOS sections.
The computed quantity and frequency-grid conventions otherwise match SciPy.

```fortran
use scifort_signal, only : dp, group_delay_sos, group_delay_result

real(dp) :: sos(2,6)
type(group_delay_result) :: delay

! A two-section example, arranged as (sections, [b0,b1,b2,a0,a1,a2]).
sos(1,:) = [0.3_dp,0.2_dp,0.1_dp,1.0_dp,-0.4_dp,0.1_dp]
sos(2,:) = [1.0_dp,0.5_dp,0.0_dp,1.0_dp,-0.3_dp,0.0_dp]
delay = group_delay_sos(sos, worn=256, fs=2000.0_dp)
if (delay%status == 0) print *, delay%frequencies, delay%delay
```

## Arguments and output

- `sos`: real(dp) or complex(dp) rank-two matrix `(nsections,6)` with
  `nsections >= 1`, finite coefficients, `a0 == 1` in column 4, and at
  least one nonzero numerator coefficient per section. All six coefficients
  are expressed in **ascending powers of unit delay** (`z^-1`).
- `worn`: optional positive number of frequencies (default 512).
- `whole`: include the entire frequency circle when generating a grid
  (default `.false.`). The endpoint is excluded.
- `fs`: finite positive sampling frequency (default `2*pi`).
- `w`: optional explicit array of frequencies in `fs` units; overrides
  `worn` and `whole`.
- Returns a `group_delay_result` with allocatable `frequencies` and `delay`
  arrays and integer `status`.

Both real and complex section coefficients work. Group delay is always real
and is expressed **in samples** even when `fs` is given in Hertz.
A missing or malformed SOS matrix, bad frequency grid, or invalid sampling
frequency sets `design_status_invalid_input=1` without allocating a delay
array. If a zero or pole lies on the evaluated unit circle (or the result
cannot be represented), the affected delay is IEEE NaN and the overall
status is `design_status_singular=2`; other bins remain usable. Unlike
SciPy's `group_delay`, SciFort deliberately does not replace undefined
values with zeros.

## Numerical method

For each SOS, write `B(q) = b0+b1*q+b2*q^2` and
`A(q) = a0+a1*q+a2*q^2`, with `q = exp(-i*omega)`. Group delay is

`Re[ q B'(q)/B(q) - q A'(q)/A(q) ]`.

The group delay of a cascade is the sum of the group delays of its sections.
SciFort evaluates and normalizes each section separately. **It does not
multiply the SOS sections into a global polynomial**; this avoids severe
coefficient cancellation for high-order designs. Gains of individual
numerators are independently scaled to handle extreme values without
affecting their group delay. The method still returns undefined delay where
a section has a zero or pole on the evaluation circle.

## Verification

From the repository root:

```bash
python tools/run_signal_group_delay_sos.py
```

The runner builds the independent SOS implementation in two GNU Fortran
configurations (`-fcheck=all` and `-O3`) and checks native input-validation
cases, existing transfer-function group-delay tests, and seeded numerical
comparisons against SciPy 1.17. For short cascades, the SciPy reference uses
`sos2tf` followed by `group_delay`; for long cascades, it adds independent
SciPy group delays of individual sections. The Python tests also compute
an analytic derivative independently, including stable and complex SOS
filters through 128 sections. No FPM or additional Fortran dependencies
are required for this validation runner.

## Reference

- Julius O. Smith III, *Introduction to Digital Filters with Audio
  Applications*, "Numerical Computation of Group Delay" (as referenced in
  the SciPy 1.17 `signal.group_delay` documentation).
- Richard G. Lyons, *Understanding Digital Signal Processing*, 3rd ed.,
  p. 830 (also cited by SciPy's group-delay documentation).
