# Fast Hankel transforms on logarithmic grids

Import `fht`, `ifht`, and `fhtoffset` from `scifort_fft`.
They implement FFTLog using binary64 real samples, SciFort FFT plans, and
complex log-gamma kernels. The implementation adapts SciPy v1.15.3 under
BSD-3-Clause; see [provenance](../CODE_PROVENANCE.md) and
[the retained license](../LICENSES/scipy-fftlog-BSD-3-Clause.txt).

## Mathematical convention and grids

The forward convention follows SciPy:

```text
A(k) = integral_0^infinity a(r) J_mu(k*r) k dr
a(r) = integral_0^infinity A(k) J_mu(k*r) r dk
```

The forward measure is `k dr`, unlike the commonly used `r dr` Hankel
convention. For a conventional transform of `f(r)` with measure `r dr`,
sample `a(r)=r*f(r)` and divide the forward result by the positive grid
values `k`.

For zero-based mathematical index `j=0,...,n-1`, set
`center=(n-1)/2`. The grids are

```text
r(j) = r_center * exp((j-center)*dln)
k(j) = k_center * exp((j-center)*dln)
offset = log(r_center*k_center)
k(j) = exp(offset)/r(n-1-j)
```

The supplied samples have uniform natural-log spacing `dln`, which can
be positive or negative but must be finite and nonzero. Grid coordinates
are not passed to the routines. `mu` is the real Bessel order, not a
dimension or a transform type.

FFTLog is a discrete transform of a logarithmically periodic sequence.
It approximates an integral on `(0,infinity)` only when the chosen range,
resolution, and bias resolve the function adequately. Ringing and aliasing
from truncation and periodic replication are possible.

## One-dimensional API

```text
b = fht(a, dln, mu [, offset, bias, status])
a = ifht(b, dln, mu [, offset, bias, status])
offset = fhtoffset(dln, mu [, initial, bias, status])
```

Input arrays are finite, real, rank-one sequences of positive length; strided
sections are accepted and preserved. Results are allocatable arrays of the
same length. Rejected transform calls return an allocated zero-length result.
Length one is supported. Complex samples are not accepted, matching SciPy's
real-data interface.

`offset` and `bias` default to zero. Use the same `dln`, `mu`, `offset`,
and `bias` for an inverse pair. No separate FFT normalization is required.

The bias `q=bias` removes `exp(q*(j-center)*dln)` from the input before the
discrete transform and restores the corresponding output power law. It
changes which biased sequence is assumed periodic. It can improve a
continuous-integral approximation, but a large bias span amplifies rounding.

`fhtoffset` returns the nearest low-ringing offset to `initial` (default
zero), within half a logarithmic grid step, subject to binary64 rounding.
It imposes a real endpoint on the continuous Mellin coefficient at
`omega=pi/dln`. Nearest-even ties follow SciPy's rounding convention.
The low-ringing choice improves the periodic join of the coefficients;
it cannot repair unresolved input or eliminate all ringing.

## Batches and selectable axes

```text
call fht_nd(a, b, dln, mu [, offset, bias, axis, status])
call ifht_nd(a, b, dln, mu [, offset, bias, axis, status])
```

These subroutines accept real arrays of rank 1-15 and a preallocated output
of identical shape. Input and output must not overlap. By default, the last
dimension is transformed, as in SciPy's batched `fht`/`ifht`.
The optional `axis` is a **SciFort extension** selecting any positive,
one-based dimension. The transform extent must be positive; an untransformed
zero extent is a valid empty batch.

Fourier plans and Mellin weights are computed once per call and reused across
all lines. There is no global cache. Strided sections use compiler-managed
temporaries. The routines also allocate a full result buffer, coefficients,
and line work arrays; all element counts must fit default integers.
Failed calls leave representable nonempty rank 1-15 outputs filled with NaNs.

## Singularities, range, and accuracy

Optional status values are:

| Value | Meaning |
| --- | --- |
| 0 | Successful call |
| 1 | Invalid shape, axis, empty transform, zero spacing, or nonfinite input/parameter |
| 2 | Caught allocation or FFT setup/execution failure |
| `fft_status_singular = 3` | Singular transform direction or nearly zero even Nyquist projection |
| `fft_status_range_error = 4` | Parameters, multiplier spread, bias span, or output exceed guarded numerical range |

`fhtoffset` returns a quiet NaN on failure. It uses statuses 0, 1, and 4.
An empty batch has no data/coefficients to evaluate; parameter and shape
validation still apply.

The Mellin multiplier is

```text
u(omega) = exp(-i*omega*offset) * 2**(q+i*omega)
           * Gamma((mu+1+q+i*omega)/2) / Gamma((mu+1-q-i*omega)/2)
```

At zero frequency, a numerator gamma pole makes the forward direction
singular; a denominator pole makes the inverse singular. SciPy warns and
omits the divergent constant in such directions. SciFort instead rejects
the call with status 3. The well-defined opposite direction has a zero
DC coefficient and is supported. It loses that component and therefore
does not have an invertible round trip. Simultaneous numerator/denominator
poles have a finite residue ratio and are supported.

For even lengths, FFTLog uses the real part of the last multiplier.
SciFort rejects the call with status 3 when its cosine factor has magnitude
at most `32*epsilon(1.0_dp)`, rather than divide by a numerically vanishing
coefficient. Choosing `fhtoffset` avoids this endpoint degeneracy in ordinary
parameter ranges.

Amplitude normalization and logarithmic restoration avoid unnecessary
intermediate overflow for large/small inputs, including tested amplitudes
`1e-200` and `1e200`. They do not create additional precision or prevent
all underflow. The following conservative guards return status 4:

- `abs(mu)` or `abs(bias)` exceeds `real(huge(0),dp)/4`, protecting integer
  recurrence-index conversions in the existing log-gamma kernel.
- `abs(dln)`, `abs(mu)`, `abs(bias)`, `abs(offset)`, or `abs(initial)` exceeds
  `sqrt(huge(1.0_dp))/64`, or `abs(dln)` is below its reciprocal.
- For a transform, `abs(bias)*((n-1)*abs(dln)/2+abs(offset))` exceeds
  `log(huge(1.0_dp))`.
- The nonzero log-magnitude spread of the multipliers exceeds that same
  logarithmic limit, or restoring a result would overflow. A small rounding
  margin is retained at the upper representable output limit.

Very large orders and frequencies also reduce phase and gamma-ratio accuracy
before a range guard is reached. Wide bias spans can amplify tiny coefficient
and FFT errors through the two power-law scalings. In a tested sensitive
case (`n=31, mu=100, dln=0.7, abs(bias)=0.8`), raw scaled errors against
90-digit direct sums were about `1.7e-9`; SciPy also lost accuracy in this
regime. This is not a uniform machine-precision guarantee. The validation
record reports raw errors alongside bounds on rounding amplification.
Prefer modest bias spans, check resolution/range convergence, and assess
continuous-integral accuracy independently of a discrete inverse round trip.

## Example

```fortran
program hankel_example
    use scifort_fft, only : dp, fht, ifht, fhtoffset
    implicit none
    integer, parameter :: n = 256
    real(dp) :: r(n), k(n), a(n), dln, offset
    real(dp), allocatable :: transformed(:), recovered(:)
    integer :: j, status

    dln = 8.0_dp*log(10.0_dp)/real(n-1,dp)
    do j = 1, n
        r(j) = exp(-6.0_dp*log(10.0_dp)+real(j-1,dp)*dln)
    end do
    offset = fhtoffset(dln,0.0_dp,initial=-4.0_dp*log(10.0_dp),status=status)
    if (status /= 0) error stop 'offset failed'
    k = exp(offset)/r(n:1:-1)
    a = r*exp(-0.5_dp*r*r)
    transformed = fht(a,dln,0.0_dp,offset,status=status)
    if (status /= 0) error stop 'fht failed'
    recovered = ifht(transformed,dln,0.0_dp,offset,status=status)
    if (status /= 0) error stop 'ifht failed'
end program hankel_example
```

## Verification and sources

See the [validation record](../validation/fftlog_2026-11-12.md) for exact
commands, compiler, reference comparisons, and error metrics.
`test/test_fftlog.f90` covers inverse pairs, norm preservation, power laws,
the Gaussian integral identity, negative spacing/orders, pole limits,
batched/scalar agreement, strided inputs/outputs, rank 15, and failure cases.

The defining specifications are
[SciPy fht](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fht.html),
[ifht](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.ifht.html),
[fhtoffset](https://docs.scipy.org/doc/scipy/reference/generated/scipy.fft.fhtoffset.html),
and [Hamilton's FFTLog derivation](https://jila.colorado.edu/~ajsh/FFTLog/).
The established original Fortran import remains deferred because its
adaptation permissions are unclear; the precise issue and inspected sources
are recorded in [provenance](../CODE_PROVENANCE.md).
