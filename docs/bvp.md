# Boundary-value ordinary differential equations

`scifort_integrate` exports the real/complex generic `solve_bvp`. It solves
first-order two-point systems on a finite interval:

```text
y' = fun(x,y,p) + S*y/(x-a)
bc(y(a),y(b),p) = 0
```

The singular term is optional. There are `n` state variables, `k` unknown
parameters, and exactly `n+k` boundary equations, including nonseparated
conditions that involve both endpoints. Convert higher-order equations to
first-order systems before calling the solver. All real values use binary64;
complex systems require holomorphic equations and boundary conditions. Split
nonholomorphic systems into real and imaginary components yourself.

## Calling convention

```fortran
call solve_bvp(fun, bc, x, y, res, p=p, s=s, fun_jac=fun_jac, &
    bc_jac=bc_jac, tol=1e-6_dp, max_nodes=1000, bc_tol=1e-8_dp)
```

Only the first five arguments are required. `x(:)` must contain at least two
finite, strictly increasing coordinates with finite positive interval widths.
`y(:,:)` is a finite initial guess with shape `(n,size(x))`, `n>=1`.
Input arrays are not modified. Use `solve_bvp_result` for real `y` and
`solve_bvp_complex_result` for complex `y`; `x` remains real in both cases.

`p(:)` supplies the finite initial parameter guess; omit it when `k=0`.
`s(:,:)` is a finite `(n,n)` real or complex matrix matching the state type.
The solver projects the left endpoint onto `S*y(a)=0`, applies the corresponding
pseudoinverse correction to its derivative, and enforces that constraint during
Newton iteration. Boundary conditions must be consistent with this constraint.

`tol` defaults to `1e-3`; `bc_tol` defaults to the effective `tol`.
Both must be finite and positive. The effective `tol` is at least
`100*epsilon(1.0_dp)` and is returned in `res%tol`. `max_nodes` defaults to 1000
and must be at least the initial mesh size. Size combinations that exceed
default-integer indexing are rejected before allocation.

Callbacks are subroutines with assumed-shape arrays and explicit interfaces:

```fortran
subroutine fun(x,y,p,f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:,:), p(:)
    real(dp), intent(out) :: f(:,:)
end subroutine

subroutine bc(ya,yb,p,residual)
    real(dp), intent(in) :: ya(:), yb(:), p(:)
    real(dp), intent(out) :: residual(:)
end subroutine
```

For complex problems replace state, parameter, and output declarations by
`complex(dp)`, retaining real `x`. Callbacks always receive `p`, including a
zero-length array when parameters are absent. `fun` must fill its entire
`(n,m)` output and evaluate independent columns, including batches with
different sizes during refinement. `bc` fills `n+k` residuals. Nonfinite
outputs fail the solve or reject a trial step; callbacks must not perform
invalid arithmetic internally.

Optional analytic callbacks have signatures
`fun_jac(x,y,p,df_dy,df_dp)` and
`bc_jac(ya,yb,p,dya,dyb,dp_bc)`. Output shapes are respectively
`(n,n,m)`, `(n,k,m)`, `(n+k,n)`, `(n+k,n)`, `(n+k,k)`.
Fill the complete arrays, including zero entries; zero-length parameter
dimensions are allowed. Supply derivatives of the regular `fun` only: the
solver adds derivatives of the singular term. Missing Jacobians are estimated
by forward differences with realized perturbations of size approximately
`sqrt(epsilon)*(1+abs(value))`. Complex differences perturb the real direction.

## Results and dense evaluation

Results contain `x`, `y`, `yp`, `p`, `rms_residuals`, `bc_residuals`, `sol`,
`status`, `success`, `message`, `tol`, `niter`, `nfev`, and `njev`.
`yp` contains the modified RHS, including any singular term. `niter` counts
mesh iterations, `nfev` counts vectorized RHS calls, and `njev` counts ODE
Jacobian batches, including estimated batches; these counts need not match
SciPy's implementation. Boundary calls are not included in `nfev`.

The retained C1 cubic solution owns its mesh and data independently of the
public result arrays. Evaluate a scalar coordinate into a vector or a vector
of coordinates into an `(n,nquery)` array:

```fortran
call res%sol%evaluate(0.5_dp, value, status)
call res%sol%evaluate(points, values, status, nu=1, extrapolate=.false.)
```

Queries can be unsorted. `nu` is a nonnegative derivative order, default zero;
orders above three return zero. Second and third derivatives need not be
continuous; at an interior knot the right interval is selected, and at the
last knot the left interval is selected. Extrapolation extends endpoint cubics
by default. Invalid shape, nonfinite coordinate, negative order, disabled
extrapolation outside the interval, or an uninitialized solution returns NaNs
and status -1. Valid evaluation returns status zero. This evaluation status
describes argument validation, not an accuracy or floating-point range bound
for arbitrarily distant extrapolation.

| Solver status | Meaning |
| --- | --- |
| `bvp_status_success` (0) | Interval RMS defects and absolute boundary residuals meet their tolerances. |
| `bvp_status_max_nodes` (1) | Required refinement exceeds the mesh budget. |
| `bvp_status_singular` (2) | Numerically singular collocation Jacobian. |
| `bvp_status_bc_tolerance` (3) | Boundary tolerance remains unmet after ten mesh iterations without further required refinement. |
| `bvp_status_invalid_argument` (-1) | Invalid input or unsupported indexing size. |
| `bvp_status_numerical_failure` (-2) | Nonfinite evaluation, unusable Newton step, or mesh stagnation. |
| `bvp_status_setup_failure` (-3) | Caught work-array allocation failure or failed singular-term pseudoinverse. |

Check `success` before relying on a result. Convergence failures retain the
last iterate, with dense output when its mesh derivatives were evaluated
successfully. RMS residuals are returned on success and the mesh-budget and
boundary-tolerance exits. Invalid initial arguments leave solution arrays
unallocated. Setup errors may return only partial diagnostics. Allocations
inside reused linear-algebra routines are subject to their existing behavior.

## Method and limits

The solver adapts SciPy 1.15.3's fourth-order Lobatto IIIA collocation,
damped Newton correction, five-point Lobatto defect estimate, and adaptive
two-/three-way interval subdivision. The residual is normalized componentwise
by `1+abs(fun)` and integrated as an interval RMS, as in SciPy. This is a
local defect contract; it does not bound global solution error in an
ill-conditioned BVP. Different initial guesses can select different nonlinear
branches or fail to converge.

SciFort uses its native row-pivoted sparse LU with row-scaled collocation
equations and no external LAPACK dependency. It refreshes the Jacobian on
each of at most eight Newton corrections per mesh iteration; SciPy also
reuses frozen Jacobians. Complex line search uses a real squared magnitude
norm. Consequently iteration counts, meshes, convergence paths, and failure
classifications can differ. Sparse pivot searches can become expensive on
large systems; highly ill-conditioned matrices can be classified singular
at the backend's binary64 pivot threshold. Arithmetic must remain in the
binary64 range. The solver is intended for finite intervals and moderate
systems, not as an unrestricted replacement for specialized large BVP codes.

There is no verbose output or interactive input, callback event facility,
periodic/multipoint interface, continuation driver, or C ABI addition.
The underlying two-point boundary callback can express periodic conditions
when they supply a well-posed system. Wrapper and C ABI expansion are deferred.

See the [validation record](../validation/bvp_2026-11-12.md),
[provenance](../CODE_PROVENANCE.md), and
[SciPy API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html).
