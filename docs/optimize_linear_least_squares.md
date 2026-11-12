# SciFort `scipy.optimize.nnls` and `lsq_linear` (2026-11-10)

These are new **dense real** least-squares interfaces in `scifort_optimize`.
They do not modify existing SciFort least-squares or optimization APIs.

## Native functions

```fortran
use scifort_optimize, only: nnls, lsq_linear, linear_lsq_result
real(dp) :: a(m,n), b(m), x(n), rnorm
integer :: status
type(linear_lsq_result) :: result
call nnls(a, b, x, rnorm, status, maxiter=300)
result = lsq_linear(a, b, lower=lb, upper=ub, method='bvls')
```

`nnls` returns `status=0` on success and negative status on input, iteration,
or numerical failure; `rnorm` is `||A*x-b||_2`.
`lsq_linear` returns a `linear_lsq_result` with `x`, `fun`, `grad`,
`active_mask` (-1 lower, 0 free, +1 upper), `cost`, `optimality`, `nit`,
`status`, `success`, and `message`. It reports `status=3` if the
unconstrained SVD solution is feasible, `status=1` if KKT conditions are
reached by active-set iterations, and negative status on failure.

`lower` and `upper` are optional rank-one arrays of length `n`; absent
bounds are unbounded. Infinite bounds are accepted. The interface supports
zero-row or zero-column inputs as SciFort extensions and accepts equal
lower/upper bounds to freeze a variable. `tol` defaults to `1e-10` and
`maxiter` defaults to `max(100,20*n)` active-set releases.

**Current scope**: Only `method='bvls'` is implemented. SciPy's default
trust-region reflective (`method='trf'`), sparse matrices, linear operators,
`lsmr_tol`, and arbitrary batches are **not** implemented here. A request
for method `trf` is rejected, never silently substituted. SciFort's
iteration counts and stopping tolerances are not identical to SciPy's.
The native BVLS method is an independent active-set implementation;
free-variable least-squares systems use the existing SciFort Jacobi SVD
rather than forming ill-conditioned normal equations. This implementation
has no new dependencies or third-party copied code.

## Existing SLSQP backend

The separate `nnls_slsqp` and `lsq_linear_slsqp` functions call **the actual
`bvls` subroutine** in `bvls_module`, supplied by SciFort's existing pinned
`jacobwilliams/slsqp` v1.6.0 FPM dependency. That source is not copied
into this archive. Its historical Lawson–Hanson code and Jacob Williams's
BSD-style provenance and license remain in the upstream dependency.

These bridge functions are exported by `scifort_optimize` and require the
FPM build to link SLSQP. The SciFort bridge returns a common result type.
The underlying BVLS routine overwrites its matrix and right-hand side;
the bridge copies them to preserve the caller's inputs. The dependency
returns no iteration counter, so `nit=-1` in this backend.

**Verification qualification:** The bridge was type-checked against a
source-compatible mock module. A real SLSQP/FPM compilation and runtime
verification were not possible in the local environment; run the included
`test/optimize_linear_least_squares_slsqp_reference.f90` with FPM before
relying on that backend. The native implementation was compiled and
numerically checked independently against SciPy 1.17.0.

## Standalone validation

From the project root, with `gfortran`, NumPy, and SciPy installed:

```sh
python tools/compare_optimize_linear_least_squares_scipy.py --count 480 --opt=-O0
python tools/compare_optimize_linear_least_squares_scipy.py --count 480 --opt=-O3
```

For native tests, compile `src/kinds.f90`, `src/linalg.f90`,
`src/linalg_svd.f90`, `src/optimize/linear_least_squares.f90`, and
`test/optimize_linear_least_squares_reference.f90` into an executable.
The example directly imports `scifort_optimize_linear_least_squares` and can
be compiled with only the four source modules listed above. Users building
with FPM can equivalently import these routines through `scifort_optimize`.
