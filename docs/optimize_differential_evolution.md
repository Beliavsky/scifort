# `scipy.optimize.differential_evolution` in SciFort

The new `scifort_optimize_differential_evolution` module provides bounded,
derivative-free population optimization. It is exported through the existing
`scifort_optimize` facade (as well as directly from the new module).

## Example

```fortran
use scifort_kinds, only : dp
use scifort_optimize, only : differential_evolution, differential_evolution_result
type(differential_evolution_result) :: res
real(dp) :: bounds(2,2)
bounds(:,1) = [-5.12_dp, -5.12_dp]
bounds(:,2) = [ 5.12_dp,  5.12_dp]
res = differential_evolution(rastrigin, bounds, strategy='best1bin', &
                             seed=2026, popsize=15, maxiter=250)
print *,res%x,res%fun,res%success,res%nfev
```

`rastrigin(x)` is an ordinary real(dp) scalar function accepting `real(dp) x(:)`.
A complete example is in `example/optimize_differential_evolution_example.f90`.

## Inputs

- `bounds(n,2)` is mandatory, **lower bounds in column 1, upper in column 2**.
  Both must be finite, lower may equal upper to fix a variable, and the width
  must be representable in real(dp).
- `strategy='best1bin'` by default. The six donor strategies `best1`, `rand1`,
  `best2`, `rand2`, `randtobest1`, `currenttobest1` each accept `bin` or `exp`
  crossover, corresponding to SciPy's 12 standard string strategies.
- `maxiter=1000`, `popsize=15` (multiplier, not total population). The total
  population is `max(6,popsize*n_free)` for unfixed variables; a fully fixed
  problem requires just one objective evaluation.
- `tol=0.01_dp`, `atol=0.0_dp`: population-energy convergence condition
  `std(energies) <= atol + tol*abs(mean(energies))`. Return status 0 on
  convergence. A maximum-generation stop is not reported as successful.
- `mutation`: constant real factor in `[0,2)`; **or** `mutation_range(2)` for
  per-generation dithering. If neither is specified, the default interval
  is `[0.5_dp,1.0_dp]` as in SciPy. Do not supply both options.
- `recombination=0.7_dp`: binomial/exponential crossover probability.
- `init='latinhypercube'` (default) or `init='random'`; `x0(:)` can replace
  the first initialized population member.
- `updating='immediate'` (default) lets later candidates within a generation
  use improvements from earlier candidates; `'deferred'` freezes the population
  during a generation and performs simultaneous replacement at its end.
- `seed` is an optional default-kind integer. A **local Park–Miller RNG**
  provides repeatable results without changing Fortran's global random state.
  The same numeric seed does **not** reproduce NumPy/SciPy's random draws.
- `maxfun`: optional total objective-evaluation budget; in this first version
  it must be at least the initial population size. `callback(x,fun,nit,stop)`
  may request termination after any completed generation.

### Integer coordinates, constraints, and polishing (2026-11-10)

- `integrality(n)` is optional `logical` and marks integer-valued coordinates.
  On those axes, all population and trial values are rounded and clipped to
  integer values in `[ceiling(lb), floor(ub)]`. A bound interval containing no
  integer is rejected; integer bounds above 2^53-1 in magnitude are rejected
  to avoid silently losing integral precision in `real(dp)`. The user-provided
  `x0` is projected to the integer lattice as well.
- Linear constraints use `linear_a(m,n)`, `linear_lb(m)`, `linear_ub(m)` with
  `linear_lb <= matmul(linear_a,x) <= linear_ub`. Supply all three arrays.
- Nonlinear constraints use `nonlinear(x,c)` where `real(dp) x(:)` is the
  candidate, `real(dp) c(:)` is the output of length `m`, and bounds are in
  `nonlinear_lb(m)` and `nonlinear_ub(m)`. Supply both arrays together with
  the callback. A row may be equality-constrained by equal lower/upper bounds.
  Infinite constraint bounds are supported; NaN bounds are rejected. A
  nonfinite callback output is treated as an infeasible candidate.
- Constraints are enforced by **feasibility-first selection**: feasible trials
  outrank infeasible ones; among feasible trials the smaller objective wins;
  among infeasible trials the smaller total constraint violation wins. The
  implementation is similar to but not identical to SciPy's Lampinen rule.
  Infeasible candidates do not cause objective evaluations.
- `constraint_tol` (default `1.0e-10_dp`) is an absolute tolerance on the
  **maximum individual constraint violation**. It is not a scaled relative
  tolerance; demanding exact equalities may require more generations.
- `polisher` is an **optional procedure**. The module-level
  `scifort_optimize_differential_evolution_polish` provides
  `de_lbfgsb_polisher` to reuse SciFort's existing L-BFGS-B solver. Pass
  `polisher=de_lbfgsb_polisher` (also exported by `scifort_optimize`). The
  polisher is called only for unconstrained, feasible runs that contain at
  least one continuous variable, with integer variables fixed to the selected
  integers. A polish result is accepted only when feasible with respect to
  those bounds and strictly improves the objective. Constraint polishing
  is **not yet implemented**: SciPy instead uses `trust-constr` when nonlinear
  or linear constraints are supplied. Polishing is disabled by default here,
  unlike SciPy. It can make extra objective calls beyond `maxfun`.

All objective functions should return finite values on attainable points.
Nonfinite objective values are penalized and never counted as convergence.

## Result

`differential_evolution_result` has `x(:)`, `fun`, `nit`, `nfev`, `status`,
`success`, `message`, `population(n,np)` and `population_energies(np)`.
For constraints it additionally provides `feasible`, `maxcv` (maximum
component violation at `x`), `ncev` (candidate constraint checks), and
`population_constraint_violation(np)`. `nfev` counts **only actual objective
calls**; infeasible candidates are not evaluated. An infeasible result has the
least-violating candidate as `x`, a large placeholder `fun`, and `success` false.
Status codes: 0 converged, 1 max generations, 2 callback termination,
3 maximum evaluations, 4 no feasible solution found, -1 invalid inputs.
There is no exception for ordinary input errors; the result is left
unallocated and contains a readable message.

## Scope and differences

The implementation covers core DE strategies, integer coordinates, linear
and nonlinear constraints, and optional bounded L-BFGS-B polishing.
It does **not** yet offer SciPy's Sobol/Halton initialization, user-defined
population arrays, arbitrary collections of nonlinear constraint objects,
constrained `trust-constr` polishing, vectorized objectives, or parallel workers.
Random trajectories and iteration counts are not identical to SciPy.

## Testing

With GNU Fortran, Python, NumPy, and SciPy installed:

```text
python tools/run_differential_evolution_tests.py
```

The script compiles the module and both native test programs in checked and
optimized modes. It runs 72 original unconstrained benchmarks **and 80 new
integer/constrained benchmark cases per build** against SciPy. The constrained
benchmarks verify independent objective evaluations and true mathematical
feasibility rather than expecting identical random trajectories. The standalone
tests exercise the `de_polisher` callback; the L-BFGS-B bridge builds against
the full SciFort FPM dependency graph. The existing pinned `lbfgsb` dependency
is reused; there are **no new runtime FPM dependencies**.

## Attribution

Numerical DE/rand/1/bin with dithering follows the standard Storn–Price
algorithm and the MIT-licensed example described by Ian Morrall:
https://github.com/ian-mmm/differential-evolution_f95 . The SciFort module
is independently implemented with corrected distinct-donor selection,
explicit bounds handling, and broad strategy coverage; no Fortran source
from that project was inserted verbatim. We retain its MIT copyright and
permission notice in `LICENSES/ian-morrall-differential-evolution-MIT.txt`.
The MIT-licensed StochOptim project was also reviewed, but its MPI-oriented
framework was not incorporated:
https://github.com/keurfonluu/StochOptim .
