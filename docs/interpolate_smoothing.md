# Natural cubic smoothing splines

SciFort's `scifort_interpolate` public module now exports
`make_smoothing_spline`, which implements the mathematical problem solved by
`scipy.interpolate.make_smoothing_spline`:

\[
\min_f \left\{ \sum_{i=1}^{n}w_i (y_i-f(x_i))^2 +
\lambda \int_{x_1}^{x_n}[f''(t)]^2\,dt \right\}.
\]

The solution is a **natural cubic smoothing spline**, not an interpolating
spline unless `lam=0`. Like SciPy, `n >= 5`, `x` must be strictly increasing,
`w` must be positive (and defaults to 1), and **weights are not squared** in
the least-squares objective.

## API

```fortran
use scifort_kinds, only : dp
use scifort_interpolate, only : make_smoothing_spline, ppoly_type

type(ppoly_type) :: s
real(dp) :: x(6), y(6), w(6), selected_lambda
integer :: status

x = [0.0_dp, 0.4_dp, 1.0_dp, 1.3_dp, 2.4_dp, 3.0_dp]
y = sin(x) + [0.1_dp, -0.1_dp, 0.05_dp, 0.1_dp, -0.05_dp, 0.02_dp]
w = 1.0_dp

! Choose lambda by generalized cross-validation (GCV)
s = make_smoothing_spline(x, y, w=w, status=status, &
    selected_lam=selected_lambda)
if (status /= 0) error stop 'smoothing spline failed'
print *, s%value(1.5_dp)
print *, s%value(1.5_dp, nu=1)
print *, s%integral(x(1), x(6))

! Fixed lambda, for reproducible regularization strength
s = make_smoothing_spline(x, y, w=w, lam=0.05_dp, status=status)
```

For multiple outputs, pass a `y(n,nchannels)` matrix. A scalar `lam` applies
to all channels; if absent, GCV selects a separate penalty for each channel.
`selected_lam` is then an output array with length `nchannels`.

**Intentional interface difference:** SciPy returns a `BSpline`. SciFort
returns `ppoly_type` with cubic coefficients in the local *power* basis,
reusing the existing `%evaluate`, `%value`, `%derivative`, `%integrate`,
`%integral`, `%coefficients`, and `%knots` methods. Spline values and
calculus operations are mathematically equivalent.

## Algorithm

Write `h_i=x_(i+1)-x_i`, form the three-point divided-difference matrix
`Q` and tridiagonal cubic roughness matrix `R`. The interior second
derivatives solve

\[
(R + \lambda Q^T W^{-1}Q)M = Q^T y,
\quad f(x_i) = y_i - \lambda[W^{-1}QM]_i.
\]

The system is symmetric positive definite and has **half-bandwidth two**.
SciFort uses banded Cholesky factorization and substitution, with
`O(n)` factorization and storage. It then constructs the exact natural cubic
piecewise-polynomial coefficients directly from fitted values and second
derivatives, avoiding a second spline solve.

When `lam` is omitted, a bounded golden-section minimization of GCV over
`[0,n]` selects the penalty. The GCV score is evaluated in `O(n)` using a
selected-inverse recursion for the diagonal and two upper bands of the
inverse of the smoothing system. For numerical stability, the explicit
`lambda^2` factors cancel from the GCV numerator and denominator.

**Optimizer tolerance:** SciPy's bounded scalar minimizer uses a default
absolute tolerance of about `1e-5`; SciFort uses a tighter interval
stopping criterion. Thus GCV-selected values can differ slightly when the
optimal lambda lies very close to zero or to `n`, despite the *fixed-lambda*
solutions agreeing at nearly double precision. The reference suite compares
both outputs at the same penalty and separately reports deviations from
SciPy's default automatic-penalty selection.

The method has no FITPACK, BLAS, LAPACK, Qhull or other new dependency.

## Supported behavior and limitations

- Strictly increasing, finite, one-dimensional `x` with at least five points.
- Finite scalar or matrix-valued real observations with the data axis first.
- Positive finite point weights, nonnegative finite `lam`, or automatic GCV.
- Natural endpoint second derivatives; values, derivatives, integration,
  and cubic coefficients through the existing `ppoly_type` implementation.
- No arbitrary-rank ordinate arrays yet (`y` may be rank 1 or 2).
- Unlike `UnivariateSpline`/`splrep`, the method uses a **roughness penalty**
  controlled by lambda, not a target weighted sum-of-squares `s`.
- Very large smoothing penalties or extreme coordinate/weight ratios may
  require input rescaling to avoid floating-point overflow or loss of
  significance; such errors return `smoothing_status_singular` where detected.

## Reproduction

From the repository root, with GNU Fortran, NumPy, and SciPy installed:

```bash
python tools/interpolate_reference/compare_smoothing_scipy.py
python tools/interpolate_reference/compare_smoothing_scipy.py --seed 8675309 --cases 95 --mode optimized
python tools/interpolate_reference/run_all_native_interpolation.py
```

These commands do not need FPM. The native test covers exact linear
reproduction, zero-penalty natural interpolation, natural endpoint
curvature, affine weighted residual orthogonality, multiple channels,
invalid inputs, GCV, and a 1,500-point banded fit. The Python suite
compares to the independently installed SciPy 1.17 implementation.

## Provenance

The mathematical roughness-penalty formulation and `Q`/`R` banded identity
are standard cubic-smoothing-spline mathematics. The implementation,
banded Cholesky and selected-inverse code are original SciFort code under
the MIT license, not a translation of SciPy's B-spline-basis GCV
implementation. SciPy is used only for independent numerical comparisons.
