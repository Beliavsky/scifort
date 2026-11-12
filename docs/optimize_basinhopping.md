# Basin hopping (`scipy.optimize.basinhopping`)

SciFort implements basin hopping in two modules:

- `scifort_optimize_basinhopping`: dependency-light global-search engine `basinhopping_core`, accepting any user-supplied local optimizer.
- `scifort_optimize_basinhopping_minimize`: public `basinhopping` wrapper using SciFort's Nelder–Mead, Powell, or BFGS optimizer by default. Both are reexported by `scifort_optimize`.

The algorithm starts by locally minimizing `x0`, then repeats randomized coordinate displacements, local minimization, Metropolis acceptance/rejection, and tracking the best successful minimum. A local-minimization failure rejects the corresponding trial; it is **never treated as an improved minimum**.

## Basic example

```fortran
use scifort_kinds, only : dp
use scifort_optimize, only : basinhopping, basinhopping_result
implicit none
type(basinhopping_result) :: res
res = basinhopping(objective, [2.4_dp], niter=80, seed=47, &
    temperature=0.4_dp, stepsize=1.5_dp, method='nelder-mead')
print *, res%fun, res%x, res%nit, res%naccepted
```

Here `objective` is a procedure taking `real(dp), intent(in) :: x(:)` and returning a scalar `real(dp)`.

## Main arguments

| Argument | Default | Meaning |
|---|---:|---|
| `niter` | `100` | Outer attempts after initial minimization |
| `temperature` | `1.0_dp` | Nonnegative Metropolis temperature (zero: reject all strictly uphill moves) |
| `stepsize` | `0.5_dp` | Uniform random displacement range `[-stepsize,+stepsize]` for each coordinate |
| `interval` | `50` | Step-size adaptation interval |
| `target_accept_rate` | `0.5_dp` | Acceptance-rate target, strictly between 0 and 1 |
| `stepwise_factor` | `0.9_dp` | Multiply or divide displacement amplitude during adaptation |
| `seed` | system clock | Independent local Park–Miller random state; `seed` repeats a trajectory on the same implementation |
| `niter_success` | disabled | Stop after more than this many consecutive steps without improving the best minimum |
| `method` | `nelder-mead` | Built-in local method: `nelder-mead`, `powell`, or `bfgs` |
| `local_maxiter` | local default | Maximum iterations per local solve |

`minimizer=` overrides `method=` with an explicitly supplied local solver. `take_step=` accepts a trial vector and the current `stepsize` and modifies the vector. `accept_test=` may return a decision of `-1` (force reject), `0` (use the Metropolis test), or `1` (force accept). `callback=` receives `(x,value,accepted,iteration,stop)`; it runs at the initial minimum (`iteration=0`) and at each **successful trial local minimization**. Custom step code is responsible for its own random-number generation, if any. The global `random_seed` state is not modified by the default stepper.

`basinhopping_core` has the same controls but takes an explicit third argument `minimizer`, with this interface:

```fortran
subroutine my_minimizer(fun, x0, x, value, nfev, success)
    ! fun: scalar objective procedure; x0: input vector
    ! x: output minimum; value: its objective; nfev: objective calls
    ! success: true only if the returned local minimum is reliable
end subroutine
```

## Result and status codes

`basinhopping_result` exposes the best `x(:)` and `fun`, the final accepted chain state `current_x(:)` and `current_fun`, final adaptive `final_stepsize`, `nit`, `nlocal`, `nfev`, `naccepted`, `nrejected`, `nfailed_local`, `nbest_updates`, `success`, `status`, and `message`. If the initial local search succeeds, a zero-iteration call still returns that minimum with `nlocal=1`.

| Code | Constant | Meaning |
|---|---|---|
| 0 | `bh_completed` | Requested moves attempted |
| 1 | `bh_no_improvement` | User's no-improvement threshold met |
| 2 | `bh_callback_stop` | Callback requested early termination |
| -1 | `bh_invalid` | Invalid control/input or callback decision |
| -2 | `bh_local_failure` | Initial local minimization failed |

The engine counts objective evaluations as *reported by the local minimizer*. A custom minimizer must report its own evaluations accurately. For a failed trial, no objective value from that local solve is used as a minimum. Failed trial minimizations count toward the no-improvement stopping rule.

## Comparison with SciPy

The algorithm follows the documented basin-hopping/Metropolis/adaptive-step paradigm, but it is **not a bit-for-bit port of SciPy**. RNG streams, local-optimizer stopping rules, iteration counting for early-stop conditions, and local factorization behavior may differ. Unlike Python `minimizer_kwargs`, bounds and constraints are not directly forwarded by this wrapper; users can supply a bounded or constrained `minimizer=` callback. The public wrapper accepts exactly three local method names, not all SciPy `minimize` options. Custom step callbacks receive a numeric `stepsize`; their implementation can apply it to adapt their own step magnitude.

Basin hopping does not guarantee the true global minimum has been reached. Independent seeded runs are advisable for nonconvex functions.

## Reproducing validation

From the project root, with GNU Fortran, Python, NumPy and SciPy installed:

```sh
./scripts/validate_optimize_basinhopping.sh
```

This compiles the algorithm in checked and optimized modes, runs its native tests and an independent 48-case SciPy comparison per build. It also checks that the public wrapper compiles against a **test-only substitute** for the large `scifort_optimize_minimize` dependency. That substitute is **not** a test of the real FPM minimizer. To run the actual full SciFort integration, including `test/optimize_basinhopping_public_reference.f90`, use `fpm test` in an environment with all pinned dependencies. The complete FPM suite has not been executed in the standalone GNU Fortran environment.

## Attribution

The mathematical basin-hopping method is associated with Wales and Doye, *J. Phys. Chem. A* **101** (1997), 5111–5116. SciPy's public reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html . All added Fortran implementation, native tests and examples are original MIT-licensed SciFort contributions; no external source code is copied and no runtime dependencies are added.
