# Broader supplied knots and rank-deficient smoothing

Native `make_splrep` and `make_splprep` accept nonclamped full knot vectors
for positive `s`, with ordinary or periodic boundaries, degrees 1-5.
They also support rank-deficient observation matrices when the shared
derivative-jump smoothing problem can be solved numerically.

## Accepted knot vectors

The full vector `t` must be finite and nondecreasing, with
`t(k+1)` equal to the left bound and `t(size(t)-k)` equal to the right bound.
The base interval, from the first bound through the last, must have strictly
increasing knots: repeated interior knots remain invalid. Exterior knots may
be distinct, partially repeated, or clamped. Periodic exterior knots need not
be cyclic extensions. The supplied vector is returned exactly.

Lengths remain `2*(k+1)` through `m+k+1` for ordinary fits and through
`m+2*k` for periodic fits, where `m` includes the closing periodic sample.
`t` excludes `nest` and requires positive `s`. Ordinary enclosing bounds and
periodic closed-input, final-weight, residual, and extrapolation rules remain
as documented in the constructor guides.

```fortran
! Cubic full knots with nonclamped exterior knots and simple interiors.
fit = make_splrep(u, y, k=3, s=0.1_dp, &
    t=[-0.3_dp,-0.2_dp,-0.1_dp,0.0_dp,0.25_dp,0.5_dp,0.75_dp, &
       1.0_dp,1.1_dp,1.2_dp,1.3_dp], bc_type='periodic', status=status)
! y must be closed; check status == 0 and fit%valid() before evaluation.
```

For periodic fits, the solver enforces equality of endpoint values and
derivatives through order `k-1` on the supplied basis. Cyclic coefficient
copying alone would not enforce those conditions for arbitrary exterior knots.

## Solver and numerical contract

Conventional clamped or cyclic knots first use the existing attributed banded
FITPACK solver. Other exterior forms and numerical failures of that path use
a dense fallback built from existing native B-spline, QR, and SVD procedures.
No new dependency is required. Automatic knot selection is unchanged.

The fallback projects onto the null space of the periodic endpoint constraints
when applicable. It minimizes squared kth-derivative jumps at the interior
knots, subject to the existing weighted squared-residual bound, with one shared
multiplier across all coordinates. Augmented least-squares systems avoid
forming normal equations. Numerical rank uses the existing QR/SVD relative
cutoff `max(rows,columns)*epsilon`.

A rank-deficient observation matrix no longer causes automatic rejection.
A feasible regularized problem can succeed because the penalty determines
otherwise unobserved coefficients. This rank-aware fallback is a **SciFort
extension** to the shared smoothing interface; exact SciPy failure behavior
on rank-deficient inputs is not reproduced. SVD least-squares steps use
minimum coefficient norm to resolve numerical null directions.

The large-target limit is a degree-k polynomial for ordinary fits and a
weighted constant for periodic fits. Polynomial B-spline coefficients come
from Marsden's identity, rather than an ill-conditioned numerical penalty
null space. Supplied knots are retained even at those limits. Infeasible
residual targets and failed numerical checks return status 3.

The dense fallback is limited to **256 coefficients**, **10,000 sites**, and
**256 coordinate channels**. These limits apply only when the fallback is
needed. Its normalized exterior knots must lie within magnitude `1e4` of
the normalized left bound. Existing finite-range, normalized-gap, relative
weight, physical-gap, and physical-residual checks also apply. A size or
range failure returns status 3.

Rank handling does not remove ill-conditioning. Clustered degree-five knots
can produce augmented systems with condition numbers above `1e12`, reducing
coefficient accuracy in binary64 even when backward errors and residuals are
small. Derivatives at clustered knots or outside the base interval can
amplify rounding further. No uniform forward-accuracy guarantee is claimed.

The fixed nonperiodic banded kernels also use a general interval search when
checking smoothing residuals. This handles multiple supplied knots between
sample sites; their automatic-knot paths retain their existing behavior.

## Validation and references

See the [validation record](../validation/fixed_general_2026-11-14.md),
[scalar constructor](interpolate_make_splrep.md),
[parametric constructor](interpolate_make_splprep.md), and
[periodic constructor](interpolate_periodic_constructors.md).

The mathematical objective follows
[SciPy's constructor specification](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.make_splrep.html)
and P. Dierckx, *Curve and Surface Fitting with Splines*, Oxford University
Press, 1993. Polynomial reproduction follows Carl de Boor,
[*B(asic)-Spline Basics*, Section 4](https://ftp.cs.wisc.edu/Approx/bsplbasic.pdf),
Marsden's identity and polynomial coefficient functionals. These mathematical
references are not imported implementation code; see
[provenance](../CODE_PROVENANCE.md).
