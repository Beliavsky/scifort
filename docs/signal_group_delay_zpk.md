# Group delay from zeros, poles, and gain (2026-10-28)

SciFort provides **`group_delay_zpk(z, p, k, worn, whole, fs, w)`**, an
analytic group-delay calculation for a digital filter in zero-pole-gain
(ZPK) form. The `group_delay` generic accepts this three-argument form as
well. It complements `group_delay(b,a,...)` and `group_delay_sos(sos,...)`
without converting roots into potentially ill-conditioned transfer-function
polynomials.

This is a **SciFort extension**: SciPy 1.17's `scipy.signal.group_delay`
accepts numerator/denominator coefficients, not ZPK data directly.
SciFort matches the resulting group delay for equivalently represented
filters, but computes it root by root for improved accuracy at high order.

```fortran
use scifort_signal, only : dp, group_delay_zpk, group_delay_result

real(dp) :: zeros(2), poles(2)
type(group_delay_result) :: delay

zeros = [0.25_dp, -0.40_dp]
poles = [0.70_dp, -0.60_dp]
delay = group_delay_zpk(zeros, poles, 1.0_dp, worn=512, fs=2000.0_dp)
if (delay%status == 0) print *, delay%frequencies, delay%delay
```

## Interface and conventions

- `z`, `p`: vectors of zeros and poles in the **z-plane**; each may be empty.
  Real and complex vectors are accepted in either combination.
- `k`: a required, finite, nonzero, real or complex scalar gain. Group delay
  is independent of the gain's magnitude and constant phase.
- `worn`: optional grid count (positive, default 512).
- `whole`: use frequencies around the full circle; default `.false.`.
  The endpoint is excluded in both modes.
- `fs`: finite positive sample frequency (default `2*pi`).
- `w`: optional explicit real frequency grid, in `fs` units; overrides the
  generated `worn`/`whole` grid.

`group_delay_result%frequencies` contains the requested or generated grid;
`%delay` is **always real and measured in samples**, not seconds, even when
`fs` is in Hz. `%status` is `design_status_success=0` on successful
calculation, `design_status_invalid_input=1` for malformed roots, gain, or
grid, and `design_status_singular=2` when at least one grid point is
undefined. Undefined bins contain IEEE quiet NaN; other bins remain valid.
An invalid request returns no allocated delay array.

Both explicit `group_delay_zpk` and generic `group_delay(z,p,k,...)`
support the eight combinations of real/complex roots and gain.

## Formula and numerical stability

For the transfer function

\[
H(z)=k\frac{\prod_{j=1}^{N_z}(z-z_j)}{\prod_{j=1}^{N_p}(z-p_j)},\qquad z=e^{i\omega},
\]

the group delay is

\[
\tau(\omega) = N_p-N_z
 - \operatorname{Re}\sum_{j=1}^{N_z}\frac{z_j}{e^{i\omega}-z_j}
 + \operatorname{Re}\sum_{j=1}^{N_p}\frac{p_j}{e^{i\omega}-p_j}.
\]

The integer `Np-Nz` term matters when the counts differ: for example,
`z=[]`, `p=[0]` represents one unit of sample delay and has `tau=1`.
For equal-degree systems, the expression reduces to the usual derivative
of the transfer-function phase. SciFort uses compensated summation and
checks for roots on the evaluation unit circle. A zero and pole that cancel
algebraically on the circle are still treated as an undefined point in the
supplied uncanceled representation.

Because the roots are never multiplied into a polynomial, the algorithm
avoids numerical cancellation in high-order filters. It does not guarantee
that delay is well-conditioned arbitrarily close to unit-circle roots;
the result can become large as a pole approaches a frequency bin.

## Verification

From the project root:

```bash
python tools/run_signal_group_delay_zpk.py
```

The standalone runner uses GNU Fortran with runtime checks and `-O3` and
compares 221 seeded real/complex ZPK designs (6,774 frequency points)
against independently summed SciPy 1.17 one-root transfer-function group
delays, including designs with hundreds of roots and selected Butterworth
and near-unit-circle all-pass filters. Low-order equal-degree designs also
use SciPy's `zpk2tf` followed by `group_delay`. The runner repeats the
existing SOS and transfer-function group-delay native regressions. The
separate public-module test is
`test/test_signal_group_delay_zpk_public.f90`.

Neither NumPy nor SciPy is a Fortran runtime dependency; they are used only
for independent numerical comparisons.
