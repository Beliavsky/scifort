# Complex Riemann zeta

`zeta(s)` in `scifort_special` and `scifort_zeta` accepts real or complex
binary64 arguments. The complex overload is `pure elemental`, so it also
accepts arrays. It computes the analytic continuation of the Riemann zeta
function, with its single pole at one. The existing real overload and its
optional real Hurwitz shift `q` are unchanged. Complex Hurwitz arguments and
complex `zetac` are not provided.

```fortran
program example
    use scifort_kinds, only : dp
    use scifort_special, only : zeta
    implicit none
    complex(dp) :: s, value
    s = cmplx(2.0_dp,3.0_dp,dp)
    value = zeta(s)
    print *, value
    ! Approximately (0.798021985146276, -0.113744308052939)
end program example
```

## Numerical contract and limits

- Exactly real complex inputs delegate to the existing real Riemann routine,
  including the value `-0.5` at zero, exact zeros at negative even integers,
  positive infinity at one, and a value of one at positive infinity. Their
  imaginary component retains the input signed zero. This axis path has the
  existing real routine's domain and range behavior.
- Nonreal inputs require finite components, `abs(real(s)) <= 10000`, and
  `abs(aimag(s)) <= 100000`. Inputs outside this rectangle return NaNs in both
  components. This is an implementation limit, narrower than SciPy's domain.
- Within the rectangle, `abs(s-1) < 16/huge(1.0_dp)` returns complex NaNs to
  prevent overflow in the reciprocal pole. Negative-real-part continuation
  also returns complex NaNs when its computed log magnitude exceeds
  `log(huge(1.0_dp))-1`; the safety margin conservatively excludes some
  representable values near the overflow boundary. Results below the smallest
  subnormal magnitude are returned as complex zero.
- Other nonreal inputs use Euler-Maclaurin evaluation for `real(s) >= -0.25`
  and the logarithmic functional equation to the left. Conjugate inputs have
  conjugate values to rounding accuracy. No branch cut is introduced for zeta.
- Work grows linearly with the exponent magnitude, up to roughly 100,531
  direct terms at the rectangle's corners; temporary storage is constant.
  There is no Riemann-Siegel or other high-frequency acceleration.

The [validation record](../validation/complex_zeta_2026-11-13.md) measures
`abs(result-reference)/(1+abs(reference))`, providing absolute control near
zeros and relative control for large results. In the 248-point reference set,
the maximum was about `2.45e-13` for `abs(aimag(s)) <= 200` and `2.34e-10`
overall, at `s=-10+100000i`. Acceptance bounds are `5e-13` for the former and
`3e-9` for the higher-imaginary cases. These are sampled validation bounds,
not uniform error guarantees throughout the rectangle. Relative accuracy
necessarily deteriorates near zeros, and phase error grows at large imaginary
arguments. SciPy 1.15.3 itself loses accuracy at several tested near-zero or
near-pole inputs; mpmath remains the independent oracle there.

## Sources

The Euler-Maclaurin kernel and 25 Bernoulli coefficients are adapted from
Alexey Kuznetsov's BSD-3-Clause
[Riemann_zeta-and-Gamma](https://github.com/Alexey-Kuznetsov-math/Riemann_zeta-and-Gamma)
Fortran implementation. SciFort uses binary64, compensated direct summation,
and its existing complex log-gamma, and retains the
[full license](../LICENSES/kuznetsov-zeta-BSD-3-Clause.txt).
See [provenance](../CODE_PROVENANCE.md) for the exact source commit and changes.
The formulas are [DLMF 25.2.10](https://dlmf.nist.gov/25.2.E10) and
[DLMF 25.4.2](https://dlmf.nist.gov/25.4.E2).
