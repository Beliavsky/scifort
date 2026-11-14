# Repeated interior knots in smoothing constructors

Repeated interior knots in positive-target `make_splrep` and `make_splprep`
are a **SciFort extension**, for ordinary and periodic boundaries. SciPy
v1.18.0's `_make_splrep_impl` calls `fpcheck`, which requires distinct interior
knots; its derivative-jump kernel also assumes simple interiors. These
repeated-knot cases were previously misclassified as missing shared SciPy
functionality. SciPy's separate `make_lsq_spline` accepts suitable repeated
knots and minimizes the residual without this smoothing penalty.

## Knot and continuity contract

Degrees 1-5 are supported. Interior knots must be strictly inside the base
interval, with multiplicity `r` from 1 through `k+1`. The full knot vector is
nondecreasing, finite and retained exactly. Base endpoints and length rules
remain those of the [supplied-knot interface](interpolate_fixed_general.md).
Extra copies of a base endpoint inside the interior slice and multiplicity
above `k+1` return invalid-input status 1.

A degree-k basis with multiplicity r guarantees continuity through derivative
order `k-r`; multiplicity `k+1` permits a discontinuity. The fitted function
may be smoother than that minimum. At an interior knot, native evaluation
uses the right-hand piece; the right base endpoint uses the final left-hand
piece. Pointwise derivative queries use those piecewise derivatives even
where a classical two-sided derivative does not exist. `derivative_spline`
can reject derivative objects when knot multiplicities prevent global
classical differentiability; check its status. Integration sums the ordinary
piecewise-polynomial integrals; jumps have no extra point-mass contribution.

## Objective and zero-penalty limit

The objective remains the sum over channels of squared jumps of the **kth
piecewise derivative**, counted **once per distinct interior knot**, under
the shared weighted squared-residual bound. Lower-order jumps are not added
to the penalty. It is not a distributional derivative penalty. Thus repeated
knots can preserve corners and discontinuities even at zero penalty.

For example, a cubic with a double interior knot can represent a piecewise
quadratic with a curvature change. Its third derivative can vanish on both sides, so its
third-derivative-jump penalty is zero. Assuming the zero-penalty space is
only global cubic polynomials would solve a different problem.

The solver computes the numerical null space of the jump operator after
applying any periodic endpoint constraints. If the least-squares fit within
that space meets `s`, it returns that fit without forcing residual equality.
Otherwise it uses the existing augmented least-squares penalty search with
the same `0.001` relative residual tolerance. Minimum coefficient norm
resolves numerical null directions. Periodic fits still enforce endpoint
values and derivatives through `k-1`; interior discontinuities do not relax
those boundary constraints. Final periodic samples and weights remain ignored
in the residual, following the constructor contract.

## Limits and validation

Repeated interiors always use the existing rank-aware dense solver, limited
to **256 coefficients, 10,000 sites and 256 channels**, with its existing
normalization, range and physical-residual checks. Size/range, unresolved
numerical failures and infeasible residuals return status 3 and invalid
objects; failed parametric calls leave `parameters` unallocated. Clustered
high-multiplicity knots can be severely ill-conditioned: small residuals
and backward errors do not guarantee accurate coefficients or derivatives.

Simple-interior banded paths, automatic generation, interpolation, public
signatures, legacy tuple knot rules and C ABI are unchanged. See the
[validation record](../validation/repeated_knots_2026-11-14.md).

References: SciPy v1.18.0
[`fpcheck`](https://github.com/scipy/scipy/blob/v1.18.0/scipy/interpolate/_bsplines.py)
and [`_make_splrep_impl` / `disc`](https://github.com/scipy/scipy/blob/v1.18.0/scipy/interpolate/_fitpack_repro.py);
Carl de Boor, [*B(asic)-Spline Basics*](https://ftp.cs.wisc.edu/Approx/bsplbasic.pdf)
for multiplicity, continuity and piecewise differentiation. The extension
uses SciFort's existing derivative-jump objective and native SVD/QR solvers;
it does not translate a SciPy repeated-knot smoothing implementation.
