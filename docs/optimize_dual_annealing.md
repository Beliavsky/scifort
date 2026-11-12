# Dual annealing (`scipy.optimize.dual_annealing`)

SciFort provides generalized simulated annealing in two independent Fortran
modules, compiled with the existing `real(dp)` kind:

- `scifort_optimize_dual_annealing`: the standalone `dual_annealing_core`
  numerical search engine, optionally taking a local-search callback.
- `scifort_optimize_dual_annealing_minimize`: the public `dual_annealing`
  adapter, normally invoking SciFort's existing bounded L-BFGS-B (or Powell)
  solver for local polishing.

Both interfaces are re-exported by `scifort_optimize`. Their implementation
is original MIT-licensed code based on the published generalized-annealing
formulas of Xiang et al. and the public SciPy API. No third-party sources
have been copied and no new dependencies are introduced.

## Usage

The **public** interface with local polishing (requires normal SciFort FPM
build with its existing L-BFGS-B dependency):

```fortran
use scifort_kinds, only : dp
use scifort_optimize, only : dual_annealing, dual_annealing_result
implicit none
real(dp) :: bounds(2,2)
type(dual_annealing_result) :: res
bounds(:,1) = [-4.0_dp,-4.0_dp]
bounds(:,2) = [ 4.0_dp, 4.0_dp]
res = dual_annealing(objective, bounds, seed=49, maxiter=300)
print *, res%fun, res%x, res%nfev
```

For stand-alone (no local minimization), call `dual_annealing_core(objective,
bounds, maxiter=300, seed=49)`, or use the public adapter with
`no_local_search=.true.`. The complete runnable stand-alone example is
[`example/optimize_dual_annealing_example.f90`](../example/optimize_dual_annealing_example.f90).

`objective` is a scalar function accepting one `real(dp), intent(in) :: x(:)`
argument and returning `real(dp)`.

## Control parameters

| Keyword | Default | Description |
|---|---:|---|
| `maxiter` | `1000` | Maximum temperature steps; zero evaluates only the initial point |
| `initial_temp` | `5230.0_dp` | Starting temperature in (0.01, 50000] |
| `restart_temp_ratio` | `2.0e-5_dp` | Restart when temperature is below this fraction of initial temperature |
| `visit` | `2.62_dp` | Heavy-tail visiting parameter (currently 1 < visit < 3) |
| `accept` | `-5.0_dp` | Generalized acceptance parameter in [-10000,-5] |
| `maxfun` | `10000000` | Objective-evaluation cap; soft when local minimization is active |
| `seed` | time-dependent | Private Park–Miller generator; independent of Fortran's global `random_seed` |
| `x0` | random in bounds | Optional feasible starting vector |
| `no_local_search` | `.false.` in public adapter | Disable local polishing |
| `method` | `'l-bfgs-b'` | Local method `'l-bfgs-b'` or `'powell'` (public adapter) |
| `local_maxiter` | `200` | Local-minimizer iteration cap (public adapter) |
| `local_interval` | `20` | Also try polishing after this many outer steps; always polish on improvements |

The standalone core performs annealing only if no `minimizer=` callback is
supplied. If a callback *is* supplied it is used for polishing unless
`no_local_search=.true.`. The public adapter defaults to bounded L-BFGS-B.
A supplied `minimizer=` callback overrides the named local method.

### Local minimizer and callback interfaces

```fortran
subroutine my_local(fun, lower, upper, start, x, value, nfev, success)
    ! fun is the objective procedure; lower,upper,start are input vectors
    ! x(:), value, nfev, success are returned by the local optimizer
end subroutine

subroutine my_callback(x, value, context, stop)
    ! x(:) and value are a newly improved global optimum
    ! context = 0 for annealing, 1 for local search
    ! set stop=.true. to stop the run
end subroutine
```

Local minimizers must return accurate objective evaluation counts, finite
outputs satisfying bounds and a truthful convergence indicator. Failed local
solves are ignored rather than treated as improvements. If an annealing
objective evaluation is nonfinite, that candidate is skipped.

## Results

`dual_annealing_result` provides `x(:)`, `fun`, `current_x(:)`,
`current_fun`, `nit`, `nfev`, `nlocal`, `naccepted`, `nrestarts`,
`nbest_updates`, `status`, `success`, and `message`.

- `da_maxiter = 0` — requested annealing iterations completed (`success=.true.`).
- `da_maxfun = 1` — objective evaluation budget reached.
- `da_callback_stop = 2` — requested stop from callback.
- `da_invalid = -1` — rejected invalid bounds, controls or callback configuration.
- `da_no_initial = -2` — no finite initial objective could be found.

`success` means normal completion, **not** mathematical proof that the global
minimum was found. With a local solver the `maxfun` limit is soft, matching
SciPy's allowance for a local search to finish its ongoing iteration.

## Algorithm and differences from SciPy

The generalized visiting distribution uses the distorted Cauchy–Lorentz
Gaussian-ratio construction; moves change either all coordinates or one
coordinate at a time. The cooling law and acceptance probability use the
Tsallis/Xiang generalized simulated-annealing formulas. When the temperature
falls below the requested restart fraction, the chain restarts from another
feasible point while retaining the best objective value.

The numerical machinery is a **new implementation, not a line-by-line port**.
Its random-number generator and exact strategy for triggering local searches
are different from SciPy's, so equal seeds do not produce equal trajectories.
At the singular visiting limit `visit=3`, SciFort currently rejects the input
rather than attempting SciPy's formally allowed endpoint. Very small visiting
parameters close to one may have different numerical normalization conventions
because the gamma factors approach poles. In particular, the extreme visiting
parameter regimes have not been validated to match SciPy's distribution.
The callback currently produces contexts 0 and 1, not SciPy's optional context
2. Other SciPy options such as `minimizer_kwargs` dictionaries are represented
by a Fortran callback or the explicit `method`/`local_maxiter` options.

## Verification and limits

Run `./scripts/validate_optimize_dual_annealing.sh` with GNU Fortran, Python,
NumPy and SciPy installed. It builds the independent engine with runtime
checks and under `-O3`, checks native error paths, and compares eight benchmark
functions with SciPy `dual_annealing(no_local_search=True)`. Comparisons check
feasibility, actual objective values, function-evaluation counts, and ability
to approach known minima. They deliberately do not compare RNG draws.

The public adapter is compiler-tested with a clearly labeled **test-only mock**
optimizer. This is not a numerical test of the actual pinned L-BFGS-B or Powell
backends. A complete FPM build/test with the project's external dependencies
is still needed to verify the real local-search integration.
