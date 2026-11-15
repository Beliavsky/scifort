# Modern smoothing near the least-squares limit

`make_splrep` and `make_splprep` minimize derivative jumps while meeting a
shared weighted squared-residual target. Supplied knots request smoothing,
not just unpenalized least squares. Automatic constructors use the same solve
after selecting knots within the modern generator limits.

For simple clamped ordinary knots and canonical periodic knots, feasible
positive targets just above the least-squares residual now enter the finite-
penalty rational search even when that residual already lies within the
`0.001*s` acceptance band. Previously these fits could return the unpenalized
coefficients immediately. The extra search follows the inspected modern
reference behavior and can change coefficients and higher derivatives while
both old and new residuals satisfy the tolerance.

The change applies to scalar and shared-channel fits, degrees 1-5. Public
arguments, knot retention, weights, parameter generation, periodic closing-
sample rules and the residual tolerance are unchanged. Legacy tuple fitting
and continuation retain their previous stopping behavior. Polynomial/constant
limits and existing acceptance allowances at or immediately below the least-
squares residual also remain unchanged.

This is a targeted compatibility correction. General nonclamped, repeated-
knot and rank-deficient fallback solves still use their documented rank-aware
search. Canonical periodic selection still intentionally differs from stale
exterior-knot reference paths. Exact root histories for all inputs, general
conditioning guarantees and matching invalid reference results are not
promised. C ABI and wrappers remain deferred.

The [validation record](../validation/modern_root_2026-11-14.md) contains
reference versions, source hashes, tests and numerical comparisons. See also
[modern selection](interpolate_modern_selection.md),
[general fixed knots](interpolate_fixed_general.md), and
[provenance](../CODE_PROVENANCE.md).
