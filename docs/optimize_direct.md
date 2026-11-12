# `scipy.optimize.direct`-style deterministic global optimization

SciFort adds `scifort_optimize_direct` and exports `direct` through
`scifort_optimize`.

This is an independent MIT-licensed implementation of the DIRECT
(DIviding RECTangles) algorithm described by Jones, Perttunen, and Stuckman
(1993), *Lipschitzian Optimization Without the Lipschitz Constant*, JOTA
79:157-181, DOI [10.1007/BF00941892](https://doi.org/10.1007/BF00941892).
The normalization, splitting procedure, and potentially optimal rectangle
criteria are defined by that paper. The SciPy public API documentation was
consulted for names and option meanings, not copied as implementation code.

## Basic use

```fortran
use scifort_kinds, only : dp
use scifort_optimize_direct, only : direct, direct_result
! ...
real(dp) :: bounds(2,2)
type(direct_result) :: res
bounds(:,1) = [-5.0_dp,-5.0_dp]
bounds(:,2) = [ 5.0_dp, 5.0_dp]
res = direct(objective, bounds, maxfun=3000, maxiter=500, &
             locally_biased=.false.)
! objective(x) is a function accepting real(dp) x(:) and returning real(dp)
```

`direct` accepts a scalar objective callback and a shape `(n,2)` bounds
array with **finite, strictly increasing** limits. It returns best solution
`x`, best value `fun`, number of function calls `nfev`, iteration count `nit`,
number of rectangles `nrect`, a status code, `success`, and a message.

| Option | Default | Meaning |
|---|---|---|
| `eps` | `1e-4` | Potential-optimality improvement requirement |
| `maxfun` | `min(1000*n,1000000)` | **Strict** objective evaluation limit |
| `maxiter` | 1000 | Rectangle-selection iterations |
| `locally_biased` | `.true.` | DIRECT-L-style longest-edge size metric; otherwise half-diagonal |
| `f_min` | absent | Known global objective value for stopping |
| `f_min_rtol` | `1e-4` | Relative target tolerance (unit-scaled at zero) |
| `vol_tol` | `1e-16` | Normalized volume stopping criterion |
| `len_tol` | `1e-6` | Best-rectangle half-length stopping criterion |
| `callback` | absent | `subroutine callback(x, value, stop)` after a completed iteration |

Stop statuses: `direct_maxfun=1`, `direct_maxiter=2`, `direct_vol_tol=3`,
`direct_len_tol=4`, `direct_target=5`, `direct_callback_stop=6`,
`direct_invalid=-1`, `direct_numerical=-2`.
Only tolerance/target stops set `success=.true.`. Budget and callback exits
return a usable best point but are **not** reported as successful convergence.

## Algorithm and compatibility

The hyperrectangle centers and edge lengths are normalized to the unit
hypercube. Each iteration chooses low-valued rectangles on a slope-based
lower hull of size versus center objective, using the improvement condition
from the original algorithm. For each selected rectangle, it evaluates the
two off-center points along every longest edge (or as many complete pairs
as fit in the remaining evaluation budget), orders axes by their better
objective value, and trisects along those axes.

A single minimum-valued rectangle is chosen for each numerical size class.
This is an original standalone implementation, **not** a byte-for-byte port
of SciPy's C/Fortran backend, so subdivision ties and trajectories may differ.
This release uses the standard maximum-length/half-diagonal size choices;
its locally biased variant is DIRECT-L-*style*, not an assertion of exact
SciPy DIRECT_L sample ordering or stopping codes.

Other deliberate differences: evaluation budgets are enforced strictly,
while SciPy treats `maxfun` as approximate; nonfinite objective values are
reported as a numerical failure rather than being permitted as penalties;
and there are no integer or general nonlinear constraints. The solver is
deterministic and does not accept a seed. Each call has its own working arrays.

## Validation

Run `bash scripts/validate_optimize_direct.sh` to compile and test both checked
and optimized GNU Fortran builds. The test runner executes native tests,
compares ten benchmark functions at three evaluation budgets in two selection
modes with `scipy.optimize.direct`, and runs the example.

The standalone validation does not replace the complete FPM test suite.
No external code or dependencies were added for this module.
